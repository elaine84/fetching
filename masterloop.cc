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
#ifndef MWMCMC_MASTERLOOP_HH
#define MWMCMC_MASTERLOOP_HH
#include "master.hh"
#include "time.hh"
#include <deque>

template <typename Communicator, typename Request, typename Status>
void JobTree::run_master(Communicator& world, const run_behavior& rb) {
    int nquit = 0;
    size_t iprobe_status = 0, iprobe_no_status = 0;
    bool want_to_quit = false;
    double wait_time = 0.0, work_time = 0.0;

    std::deque<int> waiting_for_work;
    std::deque<Request> reqs;
    std::vector<size_t> work(world.size(), 0);

    int rank = 1;
    double t0 = timestamp();
    while (1) {
        boost::optional<Status> status = world.iprobe(rank, MPI_ANY_TAG);
        if (status) {
            double t1 = timestamp();
            wait_time += t1 - t0;
            ++iprobe_status;

            if (status.get().tag() == MSG_W2M_WANTWORK) {
                if (work[rank])
                    abandon(work[rank], rank);

                (void) world.recv(rank, status.get().tag());
                if (want_to_quit) {
                    JobInfo info_stop;
                    info_stop.id = 0;
                    reqs.push_back(world.isend(rank, MSG_M2W_HAVEWORK, info_stop));
                } else {
                    if (rb.print_sending_work())
                        fprintf(stderr, "%f sending work to: %d\n", timestamp(), rank);
                    waiting_for_work.push_back(rank);
                }

            } else if (status.get().tag() == MSG_W2M_SETPROPOSAL) {
                ProposalInfo proposal;
                (void) world.recv(rank, status.get().tag(), proposal);
                set_proposal(proposal, rank);
                if (node_type* n = node(proposal.id))
                    create_children(n);

            } else if (status.get().tag() == MSG_W2M_UPDATE) {
                MetricChange update;
                (void) world.recv(rank, status.get().tag(), update);
                set_metric(update, rb, wait_time, work_time, rank,
                           iprobe_status, iprobe_no_status);

                const node_type* node = this->node(update.id);
                // reassign work if further work on this node is irrelevant
                if (update.abandon
                    || !node
                    || node->is_dead()
                    || node->complete())
                    abandon(update.id, rank);
                else if (reassign_ratio_ > 0
                         && node->utility() * reassign_ratio_ <= high_utility_pending(update.id).second) {
                    // also reassign work if much better work could be done
                    // (based on reassign_ratio_)
                    AbandonInfo ab(update.id, node ? node->generation() : 0);
                    reqs.push_back(world.isend(rank, MSG_M2W_ABANDON, ab));
                }

            } else if (status.get().tag() == MSG_W2M_REGISTER) {
                size_t data_size;
                (void) world.recv(rank, status.get().tag(), data_size);
                set_data_size(data_size);

            } else if (status.get().tag() == MSG_W2M_QUIT) {
                (void) world.recv(rank, status.get().tag());
                ++nquit;
                if (nquit == world.size() - 1) {
                    std::stringstream sbuf;
                    sbuf << "master wait time: " << wait_time << "\n"
                         << "master work time: " << work_time << "\n";
                    std::cerr << sbuf.str();
                    return;
                }

            } else if (status.get().tag() == MSG_PAUSE) {
                (void) world.recv(rank, status.get().tag());

            } else {
                std::cerr << "Unknown message tag " << status.get().tag() << "\n";
                assert(0);
            }

            std::pair<size_t, double> pending;
            while (!waiting_for_work.empty()
                   && (pending = high_utility_pending()).first != 0) {
                int source = waiting_for_work.front();
                JobInfo info = execute(pending.first, source);
                assert(info.id != 0);
                work[source] = info.id;
                reqs.push_back(world.isend(source, MSG_M2W_HAVEWORK, info));
                waiting_for_work.pop_front();
            }

            while (!reqs.empty() && reqs.front().test())
                reqs.pop_front();

            if (max_depth_ && completed_depth() > max_depth_)
                want_to_quit = true;

            t0 = timestamp();
            work_time += t0 - t1;
        } else
            ++iprobe_no_status;

        ++rank;
        if (rank == world.size())
            rank = 1;
    }
}

#endif
