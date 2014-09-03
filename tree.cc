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
#include "master.hh"
#include "worker.hh"
#include <boost/mpi/environment.hpp>
#include <boost/mpi/communicator.hpp>
#include <signal.h>
#include <unistd.h>
#include "time.hh"
#include "stype.hh"


/* Print usage statement.  */
void usage(void)
{
    char format[] = "`fetching` is a prefetching Metropolis-Hastings sampler\n\n"
        "Usage: fetching [-l CHAINLENGTH]\n\n"
        "Options:\n"
        " -h              Print this help statement and exit.\n"
        " -S NAME         Draw samples from NAME (default is normal).\n"
        " -l SAMPLES      Draw SAMPLES number of samples.\n"
        " -r REASSIGNRATIO Reassign if REASSIGNRATIO* better work exists [1.1].\n"
        " -s              Threshold, don't split.\n"
        "Example:                                          \n"
        " fetching -l 1024 normal\n";

    if (write(1, format, sizeof(format)) < 0)
        raise(SIGTRAP);
}

int main(int argc, char** argv) {
    const char *sampler_name = 0;
    int rank, N, chain_length = 100, batch_size = 1000, usleep_time = 0, dimension = 1;
    double tic, toc, reassign_ratio = 1.1, correlation = 0.99;
    bool threshold = false, print_sending_work = false, quiet = true;
    int update_style = INCREMENTAL, scheduling_policy = PREDICTIVE;

    /* Start timer.  */
    tic = timestamp();

    /* Setup the MPI environment.  */
    boost::mpi::environment env(argc, argv);
    boost::mpi::communicator world;

    /* Determine our MPI rank and the number of processes.  */
    rank = world.rank();
    N = world.size();

    /* Parse command line options.  */
    for (int optc; (optc = getopt(argc, argv, "htfpS:l:D:b:d:r:u:y:s:c:")) != -1; ) {
        switch (optc) {
        case 'h':              /* Print help statement.  */
            if (rank == 0) {
                SamplerType::print_list();
                usage();
            }
            return 0;
        case 'S':              /* Sampler name.  */
            sampler_name = optarg;
            break;
        case 'D':              /* Number of samples to draw.  */
            if (rank == 0)
                std::cerr << "Prefer -l CHAINLENGTH argument.\n";
            /* fallthru */
        case 'l':              /* Chain length.  */
            chain_length = strtol(optarg, 0, 10);
            break;
        case 'b':              /* Batch size.  */
            batch_size = strtol(optarg, 0, 10);
            break;
        case 'd':              /* Dimension of theta.  */
            dimension = strtol(optarg, 0, 10);
            break;
        case 'r':
            reassign_ratio = strtod(optarg, 0);
            break;
        case 'u':
            usleep_time = strtol(optarg, 0, 10);
            break;
        case 'y':
            update_style = strtol(optarg, 0, 10);
            break;
        case 's':
            scheduling_policy = strtol(optarg, 0, 10);
            break;
        case 'c':
            correlation = strtod(optarg, 0);
            break;
        case 't':
            threshold = true;
            break;
        case 'p':
            quiet = false;
            break;
        default:
            return 1;
        }
    }

    // Additional arguments are passed to the sampler.
    std::vector<const char*> sampler_args;
    if (!sampler_name && optind < argc)
        sampler_name = argv[optind++];
    while (optind < argc)
        sampler_args.push_back(argv[optind++]);
    if (!sampler_name)
        sampler_name = "normal";

    /* Lookup sampler in list of registered samplers.  */
    SamplerType* sampler;
    if (!(sampler = SamplerType::find(sampler_name))) {
        if (rank == 0) {
            std::cerr << "No such sampler \"" << sampler_name << "\"\n";
            SamplerType::print_list();
        }
        return 1;
    }
    if (sampler->set_arguments(sampler_args, rank == 0) < 0)
        return 1;

	std::stringstream sbuf;
    {
        char buf[1000];
        gethostname(buf, sizeof(buf));
        sbuf << "# " << rank << " = " << buf << "\n";
    }
	std::cerr << sbuf.str();
    world.barrier();

	std::stringstream ss;
    /* Kick off the sampling run.  */
    if (rank == 0) {
        ss << "# mpirun -np " << N << " -S " << sampler_name << " -l " << chain_length << " -b " << batch_size << " -d " << dimension << " -r " << reassign_ratio << " -c " << correlation << " -u " << usleep_time << (threshold ? " -t" : "") << " -y " << update_style << " -s " << scheduling_policy << "\n";
		std::cerr << ss.str();
        JobTree tree(N, chain_length, sampler, scheduling_policy, correlation, dimension);
        tree.set_reassign_ratio(reassign_ratio);
        tree.set_require_comparison_theta(false);
        tree.set_threshold(threshold);
        sampler->notify(tree, 0, 0, std::string(), false);
        tree.run_master(world, JobTree::run_behavior(print_sending_work, quiet));
    } else {
        sampler->work(world, batch_size, update_style, usleep_time);
        ss << "rank: " << world.rank() << " wait_time: " << sampler->wait_time << "\n";
		std::cerr << ss.str();
    }

    /* Stop timer.  */
    toc = timestamp();

    /* Exit with success code.  */
    if (rank == 0) {
        fflush(stdout);
        fprintf(stderr,
                "OK sample/s = %-7.3f samples = %-6d seconds = %-7.3f N = %d\n",
                chain_length / (toc - tic), chain_length, toc - tic, N);
    }
    return 0;
}
