/* Fetching
 * Elaine Angelino, Eddie Kohler
 * Copyright (c) 2012-2014 President and Fellows of Harvard College
 *
 * Permission is hereby granted, free of charge, to any person obtaining a
 * copy of this software and associated documentation files (the "Software"),
 * to deal in the Software without restriction, subject to the conditions
 * listed in the Fetching LICENSE file. These conditions include: you must
 * preserve this copyright notice, and you cannot mention the copyright
 * holders in advertising related to the Software without their permission.
 * The Software is provided WITHOUT ANY WARRANTY, EXPRESS OR IMPLIED. This
 * notice is a summary of the Fetching LICENSE file; the license in that file
 * is legally binding.
 */
#ifndef MW_MCMC_WORKER_HH
#define MW_MCMC_WORKER_HH
#include "protocol.hh"
#include <vector>
#include <deque>
#include <boost/mpi/communicator.hpp>
#include <sys/resource.h>
#include "time.hh"
#include <unistd.h>


namespace Worker {

template <typename Executor>
void setup_job(Executor& executor, JobInfo& info, ProposalInfo& proposal) {
    size_t pos = 0;
    while (1) {
        proposal.proposal = executor.load(info);
        proposal.log_prior = executor.log_prior();
        proposal.random_position = executor.random_position();
        executor.prepare();
        proposal.next_random_position = executor.random_position();

        if (info.has_proposal() || pos == info.path.length()) {
            proposal.id = info.id;
            break;
        }

        ++pos;
        info.random_position = proposal.next_random_position;
        info.state = JobInfo::s_need_proposal;
        if (info.path[info.path.length() - pos] == 'a')
            info.theta = proposal.proposal;
    }
}

} // namespace Worker

template <typename Executor>
double worker(Executor& executor, boost::mpi::communicator& world,
              size_t data_size, size_t batch_size=1000, size_t update_style=INCREMENTAL,
              useconds_t usleep_time=0) {
    JobInfo info;
    MetricChange change;
    ProposalInfo proposal;
    AbandonInfo abandon;
    std::deque<boost::mpi::request> reqs;
    change.id = 0;
    change.step_time = 0.0;
    double wait_time = 0.0;

    world.send(0, MSG_W2M_REGISTER, data_size);

    while (1) {    
        //if (change.id == 0)
            //std::cerr <<  world.rank() << " abandon " << info.utility << "\n";
        if (change.id == 0 || change.last_position == info.last_position) {
            // request new work
            double tw = timestamp();
            world.send(0, MSG_W2M_WANTWORK);
            //fprintf(stderr, "%f want work on rank: %d\n", timestamp(), world.rank());
            world.recv(0, MSG_M2W_HAVEWORK, info);
            //fprintf(stderr, "%f have work on rank: %d\n", timestamp(), world.rank());
            //std::cerr << world.rank() << " got " << info.id << " " << info.utility << "\n";
            if (info.id == 0) {
                world.send(0, MSG_W2M_QUIT);
                if (world.rank() == 1) {
                    struct rusage r;
                    struct timeval tv;
                    getrusage(RUSAGE_SELF, &r);
                    tv = r.ru_utime;
                    double tu = (double) tv.tv_sec + (double) tv.tv_usec / 1000000.0;
                    tv = r.ru_stime;
                    double ts = (double) tv.tv_sec + (double) tv.tv_usec / 1000000.0;
                    std::stringstream sbuf;
                    sbuf << "user: " << tu << "\n"
                         << "sys: " << ts << "\n"
                         << "ru_maxrss: " << r.ru_maxrss << "\n"
                         << "ru_ixrss: " << r.ru_ixrss << "\n"
                         << "ru_idrss: " << r.ru_idrss << "\n"
                         << "ru_isrss: " << r.ru_isrss << "\n"
                         << "ru_ixrss: " << r.ru_ixrss << "\n"
                         << "ru_minflt: " << r.ru_minflt << "\n"
                         << "ru_majflt: " << r.ru_majflt << "\n"
                         << "ru_nswap: " << r.ru_nswap << "\n"
                         << "ru_inblock: " << r.ru_inblock << "\n"
                         << "ru_oublock: " << r.ru_oublock << "\n"
                         << "ru_msgsnd: " << r.ru_msgsnd << "\n"
                         << "ru_msgrcv: " << r.ru_msgrcv << "\n"
                         << "ru_nsignals: " << r.ru_nsignals << "\n"
                         << "ru_nvcsw: " << r.ru_nvcsw << "\n"
                         << "ru_nivcsw: " << r.ru_nivcsw << "\n";
                    std::cerr << sbuf.str();
                }
                fflush(stderr);
                return wait_time;
            }
            wait_time += timestamp() - tw;

#if 0
            std::stringstream strs;
            strs << world.rank() << " exec " << info.id << " from [" << info.first_position << ", " << info.last_position << ") " << info.path << "\n";
            std::cerr << strs.str();
#endif

            double tp = timestamp();

            Worker::setup_job(executor, info, proposal);

            tp = timestamp() - tp;
            if (info.state != JobInfo::s_has_proposal) {
                proposal.proposal_time = tp;
                reqs.push_back(world.isend(0, MSG_W2M_SETPROPOSAL, proposal));
            }

            change.id = info.id;
            change.first_position = 0;
            change.last_position = info.first_position;
            change.metric_sum = info.metric_sum;
            change.metric_sum_sq = info.metric_sum_sq;
            change.abandon = false;
        }

        // do a batch of current work
        double ts = timestamp();
        size_t stop;
        if (update_style == INCREMENTAL || change.last_position == 0) {
            stop = change.last_position + batch_size;
            if (update_style == LOGARITHMIC)
                batch_size *= 2;
            }
        else
            stop = info.last_position;
        while (change.last_position != stop
               && change.last_position != info.last_position) {
            double log_prob = executor.evaluate(change.last_position);
            change.metric_sum += log_prob;
            change.metric_sum_sq += log_prob * log_prob;
            ++change.last_position;
            ++change.step_count;
        }

        if (usleep_time > 0)
            usleep(usleep_time);
        change.abandon = change.last_position == info.last_position;
        change.step_time = timestamp() - ts;
        change.worker_wait_time = wait_time;
        reqs.push_back(world.isend(0, MSG_W2M_UPDATE, change));

        // check whether messages sent successfully
        while (!reqs.empty() && reqs.front().test())
            reqs.pop_front();

        // check if master told us to abandon our work
        boost::optional<boost::mpi::status> status;
        while ((status = world.iprobe(0, MSG_M2W_ABANDON))) {
            world.recv(0, MSG_M2W_ABANDON, abandon);
            if (abandon.id == change.id
                && abandon.generation == info.generation)
                change.id = 0;
        }
    }
}

#endif
