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
#ifndef MWMCMC_RANDOMSEQUENCE_HH
#define MWMCMC_RANDOMSEQUENCE_HH
#include <boost/random.hpp>

template <typename Generator = boost::random::mt19937,
          typename Allocator = std::allocator<typename Generator::result_type> >
class RandomSequence : private Allocator {
    enum { capacity = (size_t) 4096 };

  public:
    typedef Generator generator_type;
    typedef typename generator_type::result_type result_type;

    inline RandomSequence()
        : position_(0), first_position_(0), last_position_(0),
          results_(Allocator::allocate(capacity)),
          state_(), initial_state_(state_) {
    }
    template <typename Seed>
    explicit inline RandomSequence(Seed seed)
        : position_(0), first_position_(0), last_position_(0),
          results_(Allocator::allocate(capacity)),
          state_(seed), initial_state_(state_) {
    }
    RandomSequence(const RandomSequence<Generator, Allocator>& x)
        : Allocator(x),
          position_(x.position_), first_position_(x.first_position_),
          last_position_(x.last_position_),
          results_(Allocator::allocate(capacity)),
          state_(x.state_), initial_state_(x.initial_state_) {
        for (size_t p = first_position_; p != last_position_; ++p)
            Allocator::construct(&results_[p & (capacity - 1)],
                                 x.results_[p & (capacity - 1)]);
    }
    template <typename Allocator2>
    RandomSequence(const RandomSequence<Generator, Allocator2>& x)
        : position_(x.position_), first_position_(x.first_position_),
          last_position_(x.last_position_),
          results_(Allocator::allocate(capacity)),
          state_(x.state_), initial_state_(x.initial_state_) {
        for (size_t p = first_position_; p != last_position_; ++p)
            Allocator::construct(&results_[p & (capacity - 1)],
                                 x.results_[p & (capacity - 1)]);
    }
    inline ~RandomSequence() {
        for (size_t p = first_position_; p != last_position_; ++p)
            Allocator::destroy(&results_[p & (capacity - 1)]);
        Allocator::deallocate(results_, capacity);
    }

    inline size_t position() const {
        return position_;
    }
    void set_position(size_t position);

    inline result_type operator()() {
        assert(position_ >= first_position_ && position_ <= last_position_);
        assert(last_position_ - first_position_ <= capacity);
        if (position_ == last_position_) {
            if (last_position_ - first_position_ == capacity) {
                for (size_t p = first_position_;
                     p != first_position_ + capacity / 2; ++p)
                    Allocator::destroy(&results_[p & (capacity - 1)]);
                first_position_ += capacity / 2;
            }
            Allocator::construct(&results_[position_ & (capacity - 1)],
                                 state_());
            ++last_position_;
        }
        ++position_;
        return results_[(position_ - 1) & (capacity - 1)];
    }

    static inline result_type min() {
        return generator_type::min();
    }
    static inline result_type max() {
        return generator_type::max();
    }

    template <typename Allocator2>
    RandomSequence<Generator, Allocator>& operator=(const RandomSequence<Generator, Allocator2>& x) {
        if (&x != this) {
            for (size_t p = first_position_; p != last_position_; ++p)
                Allocator::destroy(&results_[p & (capacity - 1)]);
            position_ = x.position_;
            first_position_ = x.first_position_;
            last_position_ = x.last_position_;
            state_ = x.state_;
            initial_state_ = x.initial_state_;
            for (size_t p = first_position_; p != last_position_; ++p)
                Allocator::construct(&results_[p & (capacity - 1)],
                                     x.results_[p & (capacity - 1)]);
        }
        return *this;
    }

  private:
    size_t position_;
    size_t first_position_;
    size_t last_position_;
    result_type* results_;
    generator_type state_;
    generator_type initial_state_;
};

template <typename Generator, typename Allocator>
void RandomSequence<Generator, Allocator>::set_position(size_t position) {
    if (position < first_position_) {
        for (size_t p = first_position_; p != last_position_; ++p)
            Allocator::destroy(&results_[p & (capacity - 1)]);
        first_position_ = last_position_ = position_ = 0;
        state_ = initial_state_;
    }
    if (position <= last_position_)
        position_ = position;
    else {
        position_ = last_position_;
        while (position_ != position)
            (*this)();
    }
}

#endif
