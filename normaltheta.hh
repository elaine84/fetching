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
#ifndef MW_MCMC_NORMALTHETA_HH
#define MW_MCMC_NORMALTHETA_HH
#include <boost/serialization/access.hpp>
#include <boost/mpi/datatype.hpp>
#include <iostream>
#include "compiler.hh"

template <int D> class NormalTheta {
  public:
    enum { size = D + D*(D+1)/2 };
    typedef NormalTheta<D> delta_type;

    inline NormalTheta();
    explicit inline NormalTheta(double data[size]);
    explicit inline NormalTheta(const double data[size]);

    inline double mean(int i) const;
    inline double variance(int i, int j) const;
    inline bool has_identity_covariance() const;

    inline double operator[](int i) const;
    NormalTheta<D> operator+(const delta_type& delta) const;

    inline void set_element(int i, double x);
    inline void add_element(int i, double x);

    friend class boost::serialization::access;
    template <typename Archive>
    void serialize(Archive& ar, const unsigned int) {
        for (int i = 0; i != size; ++i)
            ar & v_[i];
    }

  private:
    double v_[size];

    inline void set_mean(int i, double x);
    inline void set_variance(int i, int j, double x);
};

template <int D>
inline NormalTheta<D>::NormalTheta() {
    // mean = 0, variance = identity
    for (int i = 0; i != size; ++i)
        v_[i] = 0;
    for (int i = 0; i != D; ++i)
        set_variance(i, i, 1);
}

template <int N>
inline NormalTheta<N>::NormalTheta(double data[size]) {
    memcpy(v_, data, sizeof(double) * size);
}

template <int N>
inline NormalTheta<N>::NormalTheta(const double data[size]) {
    memcpy(v_, data, sizeof(double) * size);
}

template <int D>
inline double NormalTheta<D>::mean(int i) const {
    assert(i >= 0 && i < D);
    return v_[i];
}

template <int D>
inline double NormalTheta<D>::variance(int i, int j) const {
    assert(i >= 0 && i < D && j >= 0 && j < D);
    if (i > j)
        return v_[D + i*(i+1)/2 + j];
    else
        return v_[D + j*(j+1)/2 + i];
}

template <int D>
inline void NormalTheta<D>::set_mean(int i, double x) {
    assert(i >= 0 && i < D);
    v_[i] = x;
}

template <int D>
inline void NormalTheta<D>::set_variance(int i, int j, double x) {
    assert(i >= 0 && i < D && j >= 0 && j < D);
    if (i > j)
        v_[D + i*(i+1)/2 + j] = x;
    else
        v_[D + j*(j+1)/2 + i] = x;
}

template <int D>
inline void NormalTheta<D>::set_element(int i, double x) {
    assert(i >= 0 && i < size);
    v_[i] = x;
}

template <int D>
inline void NormalTheta<D>::add_element(int i, double x) {
    assert(i >= 0 && i < size);
    v_[i] = v_[i] + x;
}

template <int D>
inline bool NormalTheta<D>::has_identity_covariance() const {
    ALLOW_FLOAT_EQUALITY;
    for (int i = 0, x = 0; i != D; ++i)
        for (int j = 0; j <= i; ++j, ++x)
            if (v_[D + x] != (j == i ? 1.0 : 0.0))
                return false;
    return true;
    DISALLOW_FLOAT_EQUALITY;
}

template <int D>
inline double NormalTheta<D>::operator[](int i) const {
    return v_[i];
}

template <int D>
inline NormalTheta<D> NormalTheta<D>::operator+(const NormalTheta<D>& x) const {
    NormalTheta<D> result;
    for (int i = 0; i != size; ++i)
        result.v_[i] = v_[i] + x.v_[i];
    for (int i = 0; i != D; ++i)
        for (int j = 0; j <= i; ++j)
            if (result.variance(i, j) <= 0) {
                if (i != j)
                    result.set_variance(i, j, 0);
                else
                    result.set_variance(i, j, variance(i, j));
            }
    return result;
}

namespace boost { namespace mpi {
template <int D> struct is_mpi_datatype<NormalTheta<D> > : public mpl::true_ {};
} }

template <int D>
std::ostream& operator<<(std::ostream& s, const NormalTheta<D>& theta) {
    const char* sep = "[";
    for (int i = 0; i != D; ++i) {
        s << sep << theta.mean(i);
        sep = ", ";
    }
    sep = "; ";
    for (int i = 0; i != D; ++i)
        for (int j = 0; j <= i; ++j) {
            s << sep << theta.variance(i, j);
            sep = ", ";
        }
    s << "]";
    return s;
}

#endif
