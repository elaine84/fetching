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
#ifndef MW_MCMC_SERVER_HH
#define MW_MCMC_SERVER_HH
#include <stddef.h>
#include <boost/random.hpp>
#include <boost/serialization/access.hpp>
#include <boost/mpi/datatype.hpp>
#include <sstream>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>


enum {
    FAST_SLOW = 0,
    INCREMENTAL = 1,                  // default
    LOGARITHMIC = 2
};

enum {
    MSG_PAUSE = 1,                  // either way
    MSG_W2M_REGISTER,               // worker->master: type size_t (ndata)
    MSG_W2M_WANTWORK,               // worker->master: type void
    MSG_M2W_HAVEWORK,               // master->worker: type JobInfo
    MSG_W2M_SETPROPOSAL,            // worker->master: type ProposalInfo
    MSG_W2M_UPDATE,                 // worker->master: type MetricChange
    MSG_M2W_ABANDON,                // master->worker: type size_t (id)
    MSG_W2M_QUIT                    // worker->master: type void
};

struct MetricChange {
  public:
    size_t id;
    size_t first_position;
    size_t last_position;
    double metric_sum;
    double metric_sum_sq;
    double step_time;
    size_t step_count;
    bool abandon;
    double worker_wait_time;

    MetricChange() {
    }
    MetricChange(size_t id, size_t first_position, size_t last_position,
                 double metric_sum, double metric_sum_sq, double step_time,
                 size_t step_count, bool abandon, double worker_wait_time)
        : id(id), first_position(first_position), last_position(last_position),
          metric_sum(metric_sum), metric_sum_sq(metric_sum_sq),
          step_time(step_time), step_count(step_count), abandon(abandon),
          worker_wait_time(worker_wait_time) {
    }

    friend class boost::serialization::access;
    template <typename Archive>
    void serialize(Archive& ar, const unsigned int) {
        ar & id;
        ar & first_position;
        ar & last_position;
        ar & metric_sum;
        ar & metric_sum_sq;
        ar & step_time;
        ar & step_count;
        ar & abandon;
        ar & worker_wait_time;
    }

    friend std::ostream& operator<<(std::ostream& s, const MetricChange& c) {
        s << c.id << "[" << c.first_position << ", " << c.last_position << "): "
          << "sum=" << c.metric_sum << ", sum_sq=" << c.metric_sum_sq;
        if (c.abandon)
            s << ", abandon";
        return s;
    }
};


struct JobInfo {
    enum { root_id = (size_t) 1 };

    size_t id;
    size_t generation;
    double log_u;
    size_t random_position;
    size_t first_position;
    size_t last_position;
    double utility;
    double metric_sum;
    double metric_sum_sq;
    double scale;
    enum { s_need_initial = 0, s_need_proposal = 1, s_has_proposal = 2 };
    int state;
    std::string theta;
    std::string path;

    inline bool need_initial() const {
        return state == s_need_initial;
    }
    inline bool need_proposal() const {
        return state == s_need_proposal;
    }
    inline bool has_proposal() const {
        return state == s_has_proposal;
    }

    friend class boost::serialization::access;
    template <typename Archive>
    void serialize(Archive& ar, const unsigned int) {
        ar & id;
        ar & generation;
        ar & log_u;
        ar & random_position;
        ar & first_position;
        ar & last_position;
        ar & utility;
        ar & metric_sum;
        ar & metric_sum_sq;
        ar & scale;
        ar & state;
        if (state)
            ar & theta;
        if (state != s_has_proposal)
            ar & path;
    }
};


struct ProposalInfo {
  public:
    size_t id;
    std::string proposal;
    double log_prior;
    size_t random_position;
    size_t next_random_position;
    double proposal_time;

    ProposalInfo() {
    }
    explicit ProposalInfo(size_t node_id, const std::string& proposal,
                          double log_prior,
                          size_t random_position, size_t next_random_position,
                          double proposal_time = 0.0)
        : id(node_id), proposal(proposal), log_prior(log_prior),
          random_position(random_position),
          next_random_position(next_random_position),
          proposal_time(proposal_time) {
    }

    friend class boost::serialization::access;
    template <typename Archive>
    void serialize(Archive& ar, const unsigned int) {
        ar & id;
        ar & proposal;
        ar & log_prior;
        ar & random_position;
        ar & next_random_position;
        ar & proposal_time;
    }
};


struct AbandonInfo {
  public:
    size_t id;
    size_t generation;

    AbandonInfo() {
    }
    explicit AbandonInfo(size_t node_id, size_t generation)
        : id(node_id), generation(generation) {
    }

    friend class boost::serialization::access;
    template <typename Archive>
    void serialize(Archive& ar, const unsigned int) {
        ar & id;
        ar & generation;
    }
};


namespace boost { namespace mpi {
template <> struct is_mpi_datatype<MetricChange> : public mpl::true_ {};
template <> struct is_mpi_datatype<AbandonInfo> : public mpl::true_ {};
} }


template <typename T>
std::string save_theta(const T& theta) {
    std::stringbuf sbuf;
    boost::archive::binary_oarchive archive(sbuf, boost::archive::no_header);
    archive & theta;
    return sbuf.str();
}

template <typename T>
T load_theta(const std::string& str) {
    std::stringbuf sbuf(str);
    boost::archive::binary_iarchive archive(sbuf, boost::archive::no_header);
    T theta;
    archive & theta;
    return theta;
}

template <typename T>
std::string load_unparse_theta(const std::string& str) {
    std::stringstream sbuf;
    sbuf << load_theta<T>(str);
    return sbuf.str();
}

class JobTree;
std::ostream& operator<<(std::ostream&, const JobTree&);

#endif
