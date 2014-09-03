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
#ifndef MWMCMC_STYPE_HH
#define MWMCMC_STYPE_HH
#include <boost/mpi/communicator.hpp>
#include <iostream>
#include <vector>
#include "time.hh"
#include "protocol.hh"

class SamplerType {
  public:
    SamplerType();
    virtual ~SamplerType();

    static int signup(const char* name, SamplerType* s);
    static SamplerType* find(const char* name);
    static void print_list();

    virtual int set_arguments(std::vector<const char*>& argv, bool complain);
    virtual int work(boost::mpi::communicator& world);
    virtual int work(boost::mpi::communicator& world, size_t batch_size,
                     size_t update_style, int usleep_time);
    virtual std::string unparse(const std::string& theta) const;
    virtual void notify(const JobTree& tree, size_t depth,
                        const double metric_sum, const std::string& theta,
                        bool accept);

    double wait_time;

  protected:
    double t0_;
};

class EmptySamplerType : public SamplerType {
  public:
    EmptySamplerType() {}
    void notify(const JobTree&, size_t, const double, const std::string&, bool) {}
};

template <typename T>
class ThetaUnparser : public SamplerType {
  public:
    ThetaUnparser() {}

    std::string unparse(const std::string& theta) const;
    void notify(const JobTree& tree, size_t depth, const double, const std::string& theta, bool accept);
};

template <typename T>
std::string ThetaUnparser<T>::unparse(const std::string& theta) const {
    return load_unparse_theta<T>(theta);
}

template <typename T>
void ThetaUnparser<T>::notify(const JobTree&, size_t, const double, const std::string&, bool) {
}

#endif
