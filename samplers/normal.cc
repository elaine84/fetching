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
#include "normalsampler.hh"
#include "worker.hh"
#include "stype.hh"

namespace {
class NormalSamplerType : public SamplerType {
  public:
    enum { NDATA = 1000000 };

    NormalSamplerType() {
    }

    size_t ndata();
    int work(boost::mpi::communicator& world, size_t batch_size,
             size_t update_style, int usleep_time);
    std::string unparse(const std::string& theta) const;
};

size_t NormalSamplerType::ndata() {
    return NDATA;
}

int NormalSamplerType::work(boost::mpi::communicator& world, size_t batch_size,
                            size_t update_style, int usleep_time) {
    typedef NormalSampler<2> executor_type;
    boost::random::mt19937 engine;
    boost::random::normal_distribution<> fart;
    double* data = new double[NDATA * 2];

    for (int i = 0; i != NDATA * 2; ++i)
        data[i] = fart(engine) + (i % 2 ? 100 : 200);

    executor_type executor(data, NDATA);
    wait_time = worker(executor, world, NDATA, batch_size, update_style, usleep_time);

    delete[] data;
    return 0;
}

std::string NormalSamplerType::unparse(const std::string& theta) const {
    return load_unparse_theta<NormalTheta<2> >(theta);
}

__attribute__ ((constructor))
static void constructor() {
    SamplerType::signup("normal", new NormalSamplerType);
}
}
