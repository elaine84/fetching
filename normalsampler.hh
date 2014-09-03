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
#ifndef MWMCMC_NORMALSAMPLER_HH
#define MWMCMC_NORMALSAMPLER_HH
#include "normaltheta.hh"
#include "strideindex.hh"
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include "randomsequence.hh"

#ifndef M_LOG_2PI
#define M_LOG_2PI 1.837877066409345483560659472810
#endif

template <int D>
class NormalSampler {
  public:
    NormalSampler(const double* data, size_t data_size);

    std::string load(const JobInfo& job);
    void prepare();
    double log_prior();
    size_t random_position() const;
    double evaluate(size_t position);
    inline void evaluate_many(size_t position, size_t count,
                              double& metric_sum, double& metric_sum_sq);

  private:
    const double* data_;
    size_t data_size_;
    StrideIndex permuter_;

    bool identity_covariance_;
    gsl_vector* mean_;
    gsl_matrix* decomposed_variance_inverse_;
    double scale_;
    gsl_vector* work1_;
    gsl_vector* work2_;

    RandomSequence<boost::random::mt19937> engine_;

    NormalTheta<D> proposal_;
};

template <int D>
NormalSampler<D>::NormalSampler(const double* data, size_t data_size)
    : data_(data), data_size_(data_size), mean_(gsl_vector_alloc(D)),
      decomposed_variance_inverse_(gsl_matrix_alloc(D, D)),
      work1_(gsl_vector_alloc(D)), work2_(gsl_vector_alloc(D)) {
}

template <int D>
std::string NormalSampler<D>::load(const JobInfo& job) {
    if (job.need_initial())
        proposal_ = NormalTheta<D>();
    else {
        engine_.set_position(job.random_position);
        proposal_ = load_theta<NormalTheta<D> >(job.theta);
        if (job.need_proposal()) {
            boost::random::normal_distribution<> distribution(0, 0.1);
            for (int i = 0; i != D; ++i)
                proposal_.add_element(i, distribution(engine_));
        }
    }
    return save_theta(proposal_);
}

template <int D>
size_t NormalSampler<D>::random_position() const {
    return engine_.position();
}

template <int D>
void NormalSampler<D>::prepare() {
    identity_covariance_ = proposal_.has_identity_covariance();
    if (identity_covariance_)
        scale_ = -D * M_LOG_2PI / 2;
    else {
        for (int i = 0; i != D; ++i)
            gsl_vector_set(mean_, i, proposal_[i]);
        gsl_matrix* work = gsl_matrix_alloc(D, D);
        for (int i = 0, idx = 0; i != D; ++i)
            for (int j = 0; j <= i; ++j, ++idx) {
                gsl_matrix_set(work, i, j, proposal_[D + idx]);
                gsl_matrix_set(work, j, i, proposal_[D + idx]);
            }
        gsl_permutation* p = gsl_permutation_alloc(D);
        int signum;
        gsl_linalg_LU_decomp(work, p, &signum);
        gsl_linalg_LU_invert(work, p, decomposed_variance_inverse_);
        double determinant = gsl_linalg_LU_det(work, signum);
        gsl_matrix_free(work);
        gsl_permutation_free(p);

        // scale_ = -log(sqrt(pow(2 * M_PI, D) * determinant))
        //        = -0.5 * log(pow(2 * M_PI, D) * determinant)
        //        = -0.5 * (log(pow(2 * M_PI, D)) + log(determinant))
        //        = -0.5 * (D * log(2 * M_PI) + log(determinant))
        scale_ = -(D * M_LOG_2PI + log(determinant)) / 2;
    }

    boost::random::uniform_int_distribution<> dist(0, data_size_ - 1);
    size_t stride = dist(engine_), offset = dist(engine_);
    permuter_ = StrideIndex(data_size_,
                            StrideIndex::valid_stride(data_size_, stride),
                            offset);
}

template <int D>
double NormalSampler<D>::log_prior() {
    return 0.0;
}

template <int D>
double NormalSampler<D>::evaluate(size_t position) {
    const double* data = &data_[D * permuter_(position)];

    double dot_product;
    if (identity_covariance_) {
        dot_product = 0.0;
        for (int i = 0; i != D; ++i) {
            double x = proposal_[i] - data[i];
            dot_product += x * x;
        }
    } else {
        // work1_ := mean - data
        gsl_vector_const_view data_vector = gsl_vector_const_view_array(data, D);
        gsl_vector_memcpy(work1_, mean_);
        gsl_vector_sub(work1_, &data_vector.vector);

        // work2_ := winv * (mean - data)
        gsl_blas_dsymv(CblasUpper, 1.0, decomposed_variance_inverse_,
                       work1_, 0.0, work2_);

        // dot_product := (mean - data) [dot] (winv * (mean - data))
        gsl_blas_ddot(work1_, work2_, &dot_product);
    }

    double log_dmvnorm = -0.5 * dot_product + scale_;
    return log_dmvnorm;
}

template <int D>
inline void NormalSampler<D>::evaluate_many(size_t position, size_t count,
                                            double& metric_sum,
                                            double& metric_sum_sq) {
    for (size_t p = position; p != position + count; ++p) {
        double x = evaluate(p);
        metric_sum += x;
        metric_sum_sq += x * x;
    }
}

#endif
