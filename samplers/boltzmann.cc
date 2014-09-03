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
#include "protocol.hh"
#include "worker.hh"
#include "stype.hh"
#include "randomsequence.hh"
#include "boltzmannsampler.hh"

namespace {
class BoltzmannSamplerType : public SamplerType {
  public:
    BoltzmannSamplerType();

    int work(boost::mpi::communicator& world, size_t batch_size, int usleep_time);
    std::string unparse(const std::string& theta) const;
};

BoltzmannSamplerType::BoltzmannSamplerType() {
}

int BoltzmannSamplerType::work(boost::mpi::communicator& world, size_t batch_size,
                               int usleep_time) {
    BoltzmannSampler bs;
    worker(bs, world, 1, batch_size, usleep_time);
    return 0;
}

std::string BoltzmannSamplerType::unparse(const std::string& theta) const {
    return load_unparse_theta<double>(theta);
}

__attribute__ ((constructor))
static void constructor()
{
    SamplerType::signup("boltzmann", new BoltzmannSamplerType);
}
}
