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
#ifndef MW_MCMC_STRIDEINDEX_HH
#define MW_MCMC_STRIDEINDEX_HH
#include <boost/math/common_factor_rt.hpp>

class StrideIndex {
  public:
    inline StrideIndex();
    inline StrideIndex(size_t n, size_t stride, size_t start = 0);

    inline size_t operator()(size_t i) const;

    static inline size_t valid_stride(size_t n, size_t stride);

  private:
    size_t n_;
    size_t stride_;
    size_t start_;
};

inline StrideIndex::StrideIndex()
    : n_(0), stride_(1), start_(0) {
}

inline StrideIndex::StrideIndex(size_t n, size_t stride, size_t start)
    : n_(n), stride_(stride), start_(start) {
    assert(stride_ != 0);
    assert(n_ == 0 || stride_ < n_);
    assert(n_ ? boost::math::gcd(n_, stride_) == 1 : stride_ & 1);
    assert(n_ == 0 || start_ < n_);
}

inline size_t StrideIndex::operator()(size_t i) const {
    return (start_ + i * stride_) % n_;
}

inline size_t StrideIndex::valid_stride(size_t n, size_t stride) {
    if (stride <= 1)
        return 1;
    size_t g;
    while ((g = boost::math::gcd(n, stride)) > 1)
        stride = (stride + g - 1) % n;
    return stride <= 1 ? 1 : stride;
}

#endif
