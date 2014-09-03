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
#ifndef MW_MCMC_BOLTZMANNSAMPLER_HH
#define MW_MCMC_BOLTZMANNSAMPLER_HH

class BoltzmannSampler {
  public:
    inline BoltzmannSampler();

    inline std::string load(const JobInfo& job);
    inline void prepare();
    inline double log_prior();
    inline size_t random_position() const;
    inline double evaluate(size_t position);
    inline void evaluate_many(size_t position, size_t count,
                              double& metric_sum, double& metric_sum_sq);

  private:
    double proposal_;
    double temperature_;
    RandomSequence<boost::random::mt19937> engine_;
};

inline BoltzmannSampler::BoltzmannSampler()
    : temperature_(1.0) {
}

inline std::string BoltzmannSampler::load(const JobInfo& job) {
    if (job.need_initial())
        // initial proposal: 0.25
        proposal_ = 0.25;
    else {
        proposal_ = load_theta<double>(job.theta);
        if (job.need_proposal()) {
            engine_.set_position(job.random_position);
            // new proposal: normally distributed around previous proposal
            // with standard deviation 0.1, constrained (by repeated draws)
            // to be non negative
            boost::random::normal_distribution<> distribution(proposal_, 0.1);
            double tentative_proposal = distribution(engine_);
            if (tentative_proposal >= 0)
                proposal_ = tentative_proposal;
        }
    }
    return save_theta(proposal_);
}

inline size_t BoltzmannSampler::random_position() const {
    return engine_.position();
}

inline void BoltzmannSampler::prepare() {
}

inline double BoltzmannSampler::log_prior() {
    return 0.0;
}

inline double BoltzmannSampler::evaluate(size_t) {
    // return log prob (i.e. log(exp(-proposal_ / temperature_)))
    return -proposal_ / temperature_;
}

inline void BoltzmannSampler::evaluate_many(size_t, size_t count,
                                            double& metric_sum,
                                            double& metric_sum_sq) {
    double x = evaluate(0);
    metric_sum += x * count;
    metric_sum_sq += x * x * count;
}

#endif
