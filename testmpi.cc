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
#include <boost/mpi/environment.hpp>
#include <boost/mpi/communicator.hpp>
#include <boost/mpi/request.hpp>
#include <deque>
#include "time.hh"

void master(boost::mpi::communicator& world) {
  while (1) {
    boost::mpi::status status = world.probe(MPI_ANY_SOURCE, MPI_ANY_TAG);
    (void) world.recv(status.source(), status.tag());
    if (status.tag() == 1) {
      fprintf(stderr, "%f sending work to rank: %d\n", timestamp(), status.source());
      world.isend(status.source(), 1);
    }
  }
}

void iprobe_master(boost::mpi::communicator& world) {
    int n = world.size();
    std::deque<boost::mpi::request> sent;

    int rr_rank = 1;
    while (1) {
        boost::optional<boost::mpi::status> s = world.iprobe(rr_rank, MPI_ANY_TAG);
        if (s) {
            (void) world.recv(rr_rank, MPI_ANY_TAG);
            if (s.get().tag() == 1) {
                fprintf(stderr, "%f sending work to rank: %d\n",
                        timestamp(), rr_rank);
                sent.push_back(world.isend(rr_rank, 1));
            }
        }

        ++rr_rank;
        if (rr_rank == n)
            rr_rank = 1;

        while (!sent.empty() && sent.front().test())
            sent.pop_front();
    }
}

void imaster(boost::mpi::communicator& world) {
    int n = world.size();
    boost::mpi::request* requests = new boost::mpi::request[n];
    for (int rank = 1; rank < n; ++rank)
        requests[rank] = world.irecv(rank, MPI_ANY_TAG);
    std::deque<boost::mpi::request> sent;

    int rr_rank = 1;
    while (1) {
        boost::optional<boost::mpi::status> s = requests[rr_rank].test();
        if (s) {
            if (s.get().tag() == 1) {
                fprintf(stderr, "%f sending work to rank: %d\n",
                        timestamp(), rr_rank);
                sent.push_back(world.isend(rr_rank, 1));
            }
            requests[rr_rank] = world.irecv(rr_rank, MPI_ANY_TAG);
        }

        ++rr_rank;
        if (rr_rank == n)
            rr_rank = 1;

        while (!sent.empty() && sent.front().test())
            sent.pop_front();
    }
}

void worker(boost::mpi::communicator& world) {
    std::deque<boost::mpi::request> sent;
    while (1) {
        fprintf(stderr, "%f want work on rank: %d\n", timestamp(), world.rank());
        world.send(0, 1);
        (void) world.recv(0, MPI_ANY_TAG);
        for (int i = 0; i < 1000; ++i) {
            sent.push_back(world.isend(0, 2));
        }
        while (!sent.empty() && sent.front().test())
            sent.pop_front();
    }
}

int main(int argc, char** argv) {
    boost::mpi::environment env(argc, argv);
    boost::mpi::communicator world;
    char master_type = 'i';

    for (int optc; (optc = getopt(argc, argv, "spi")) != -1; ) {
        switch (optc) {
        case 's':
        case 'p':
        case 'i':
            master_type = (char) optc;
            break;
        default:
            fprintf(stderr, "usage: testmpi [-s|-p|-i]\n");
            return 1;
        }
    }

    if (world.rank() == 0) {
        if (master_type == 'i')
            imaster(world);
        else if (master_type == 'p')
            iprobe_master(world);
        else
            master(world);
    } else
        worker(world);
    return 0;
}
