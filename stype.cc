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
#include "stype.hh"
#include <vector>
#include <string.h>

typedef std::pair<const char*, SamplerType*> named_sampler;
static std::vector<named_sampler>* named_samplers;

SamplerType::SamplerType()
    : wait_time(0.0) {
}

SamplerType::~SamplerType() {
}

int SamplerType::signup(const char* name, SamplerType* sampler) {
    if (find(name))
        return -1;
    named_samplers->push_back(named_sampler(name, sampler));
    return 0;
}

SamplerType* SamplerType::find(const char* name) {
    if (!named_samplers)
        named_samplers = new std::vector<named_sampler>;
    for (size_t i = 0; i != named_samplers->size(); ++i)
        if (strcmp(name, (*named_samplers)[i].first) == 0)
            return (*named_samplers)[i].second;
    return 0;
}

void SamplerType::print_list() {
    if (named_samplers)
        for (size_t i = 0; i != named_samplers->size(); ++i)
            std::cout << (*named_samplers)[i].first << "\n";
}

int SamplerType::set_arguments(std::vector<const char*>&, bool) {
    return 0;
}

int SamplerType::work(boost::mpi::communicator&) {
    return -1;
}

int SamplerType::work(boost::mpi::communicator&, size_t, size_t, int) {
    return -1;
}

std::string SamplerType::unparse(const std::string&) const {
    return "<proposal>";
}

void SamplerType::notify(const JobTree&, size_t depth, const double metric_sum,
                         const std::string& theta, bool accept) {
    if (depth == 0 && !accept) {
        t0_ = timestamp();
        std::cout << "time\t"
                  << "depth\t"
                  << "accept\t"
                  << "metric_sum\t"
                  << "theta\n";
        return;
    }
    std::cout.precision(5);
    std::cout << unparse_timestamp(timestamp() - t0_) << "\t"
              << depth << "\t"
              << accept << "\t" << std::fixed;
    std::cout << metric_sum  << "\t"
              << unparse(theta) << "\n";
}
