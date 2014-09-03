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
#ifndef MW_MCMC_MASTER_HH
#define MW_MCMC_MASTER_HH
#include "compiler.hh"
#include "protocol.hh"
#include <deque>
#include <limits>
#include <algorithm>
#include <boost/math/special_functions/erf.hpp>
#include <float.h>
#include <boost/mpi/communicator.hpp>
class JobTree;
class JobNode;
class SamplerType;

enum {                              // scheduling policies
    NAIVE = 0,                      // pursue all paths
    CONSTANT = 1,                   // acceptance ratio
    PREDICTIVE = 2                  // default: use approximation
};

class JobLevel {
  public:
    double log_u;

    template <typename Engine> inline JobLevel(Engine& generator) {
        boost::random::uniform_01<> distribution;
        log_u = log(distribution(generator));
    }
};


class JobTree {
  public:
    class run_behavior;
    typedef std::string theta_type;
    typedef JobNode node_type;
    typedef JobInfo job_info_type;
    typedef JobLevel job_level_type;
    typedef size_t id_type;
    enum { root_id = (id_type) JobInfo::root_id };

    JobTree(int world_size,
            size_t max_depth = std::numeric_limits<std::size_t>::max(),
            SamplerType* completion_hook = 0, size_t scheduling_policy = 2,
            double correlation = 0.99, size_t dimension = 1);
    ~JobTree();

    void set_data_size(size_t data_size);
    void set_require_comparison_theta(bool x);
    void set_reassign_ratio(double ratio);
    void set_threshold(bool x);

    inline const node_type* root() const;
    inline const node_type* cnode(id_type id) const;
    inline node_type* node(id_type id);
    inline const node_type* node(id_type id) const;
    inline SamplerType* sampler() const;

    inline bool has_comparison_theta(id_type id) const;
    inline const theta_type& comparison_theta(id_type id) const;
    inline bool has_theta(id_type id) const;
    inline const theta_type& theta(id_type id) const;

    inline size_t completed_depth() const;
    inline size_t allocated_nodes() const;
    inline size_t recycled_nodes() const;

    inline double local_accept_rate() const;

    std::pair<size_t, double> high_utility_pending(size_t stopper = 0) const;
    JobInfo execute(size_t id, int executor);

    void set_proposal(const ProposalInfo& proposal, int executor);
    void set_metric(const MetricChange& update, const run_behavior& rb,
                    const double master_wait_time, const double master_work_time,
                    const int rank, const size_t iprobe_status,
                    const size_t iprobe_no_status);
    void abandon(id_type id, int executor);

    class run_behavior {
      public:
        run_behavior()
            : print_sending_work_(false), quiet_(true) {
        }
        run_behavior(bool print_sending_work, bool quiet)
            : print_sending_work_(print_sending_work), quiet_(quiet) {
        }
        bool print_sending_work() const {
            return print_sending_work_;
        }
        run_behavior& print_sending_work(bool x) {
            print_sending_work_ = x;
            return *this;
        }
        bool quiet() const {
            return quiet_;
        }
        run_behavior& quiet(bool x) {
            quiet_ = x;
            return *this;
        }
      private:
        bool print_sending_work_;
        bool quiet_;
    };

    void run_master(boost::mpi::communicator& world, const run_behavior& rb);
    template <typename Communicator, typename Request, typename Status>
    void run_master(Communicator& world, const run_behavior& rb);

    friend std::ostream& operator<<(std::ostream&, const JobTree&);

  private:
    node_type* root_;
    mutable boost::random::mt11213b utility_generator_;

    boost::random::mt19937 u_generator_;
    std::deque<job_level_type> levels_;
    size_t first_depth_;

    std::deque<node_type*> nodes_by_id_;
    id_type first_id_;
    size_t ndata_;

    bool require_comparison_theta_;
    bool threshold_;
    size_t allocated_nodes_;
    size_t recycled_nodes_;
    size_t max_depth_;
    double reassign_ratio_;
    SamplerType* completion_hook_;

    size_t scheduling_policy_;
    double correlation_;
    double accept_count_;
    std::deque<bool> accept_buffer_;
    size_t dimension_;

    double wasted_proposal_time_;
    double wasted_step_time_;
    double useful_proposal_time_;
    double useful_step_time_;
    std::vector<double> wait_time_;

    inline id_type next_id() const;
    void create_child(node_type* parent, bool isaccept);
    inline void create_children(node_type* parent);
    void delete_subtree(node_type* n);
    void delete_loser_child(node_type* n);
    void trim(const run_behavior& rb,
              const double master_wait_time, const double master_work_time,
              const size_t iprobe_status, const size_t iprobe_no_status);
};


class JobNode {
  public:
    typedef std::string theta_type;
    typedef JobNode node_type;

    explicit JobNode(size_t ndata, size_t scheduling_policy, double correlation, size_t dimension);

    inline bool is_root() const;
    inline bool is_accept() const;
    inline size_t id() const;
    inline size_t generation() const;
    inline size_t true_depth() const;
    inline bool has_theta() const;
    inline const theta_type& theta() const;
    inline size_t position() const;
    inline size_t data_size() const;
    inline size_t random_position() const;
    inline size_t next_random_position() const;
    inline double metric_sum() const;
    inline double metric_sum_sq() const;
    inline double metric_mean() const;
    inline double metric_var() const;
    inline double metric_stderr() const;
    inline size_t find_closest(size_t pos) const;
    inline double metric_sum(size_t ind) const;
    inline double metric_sum_sq(size_t ind) const;
    inline double metric_mean(size_t ind) const;
    inline double metric_var(size_t ind) const;
    inline double metric_stderr(size_t ind) const;
    double utility() const;
    inline bool is_dead() const;
    inline bool pending() const;
    inline bool complete() const;
    inline double branch_probability() const;
    inline bool has_final_branch_probability() const;
    inline double log_u() const;
    inline node_type* child(bool isright) const;

  private:
    size_t id_;
    size_t generation_;
    size_t depth_;
    size_t true_depth_;
    JobNode* child_[2];
    JobNode* parent_;
    JobNode* comparison_parent_;
    bool has_theta_;
    theta_type theta_;
    size_t random_position_;
    size_t next_random_position_;
    size_t ndata_;
    size_t position_;
    double scale_;
    double log_prior_;
    double metric_sum_;
    double metric_sum_sq_;
    std::deque<size_t> trace_position_;
    std::deque<double> trace_metric_sum_;
    std::deque<double> trace_metric_sum_sq_;

    double log_u_;
    double branch_probability_;

    int executor_;
    int last_executor_;
    double proposal_time_;
    double step_time_;
    size_t step_count_;
    int flip_count_;
    int threshold_position_;
    size_t scheduling_policy_;
    double correlation_;

    JobNode(size_t id, JobNode* parent, bool isright, double u, size_t ndata,
            double local_accept_rate, size_t scheduling_policy, double correlation);
    void set_branch_probability(bool threshold, double local_accept_rate);
    void set_metric(const MetricChange& update, bool threshold, double local_accept_rate);

    friend class JobTree;
};


inline bool JobNode::is_root() const {
    return !parent_;
}

inline bool JobNode::is_accept() const {
    return parent_ && this == parent_->child_[1];
}

inline size_t JobNode::id() const {
    return id_;
}

inline size_t JobNode::generation() const {
    return generation_;
}

inline size_t JobNode::true_depth() const {
    return true_depth_;
}

inline bool JobNode::has_theta() const {
    return has_theta_;
}

inline const JobNode::theta_type& JobNode::theta() const {
    assert(has_theta_);
    return theta_;
}

inline size_t JobNode::position() const {
    return position_;
}

inline size_t JobNode::data_size() const {
    return ndata_;
}

inline size_t JobNode::random_position() const {
    return random_position_;
}

inline size_t JobNode::next_random_position() const {
    return next_random_position_;
}

inline double JobNode::metric_sum() const {
    return metric_sum_;
}

inline double JobNode::metric_sum_sq() const {
    return metric_sum_sq_;
}

inline double JobNode::metric_mean() const {
    return metric_sum_ / position_;
}

inline double JobNode::metric_var() const {
    // XXX define metric_stderr_ instead of computing every time
    double variance = (metric_sum_sq_ - (metric_sum_ * metric_sum_) / position_) / position_;
    assert(variance > -0.1);
    if (variance < 0.0)
        variance = 0.0;
    return variance;
}

inline double JobNode::metric_stderr() const {
    // XXX define metric_stderr_ instead of computing every time
    double variance = (metric_sum_sq_ - (metric_sum_ * metric_sum_) / position_) / position_;
    assert(variance > -0.1);
    if (variance < 0.0)
        variance = 0.0;
    return sqrt(variance / (position_ - 1));
}

inline size_t JobNode::find_closest(size_t pos) const {
    size_t i;
    int diff, diff_prev;
    diff_prev = abs(trace_position_[0] - pos);
    for (i = 1; i < trace_position_.size(); i++) {  // obviously not optimal
        diff = abs(trace_position_[i] - pos);
        if (diff >= diff_prev)
            return (i - 1);
        diff_prev = diff;
    }
    return trace_position_.size() - 1;
}

inline double JobNode::metric_sum(size_t ind) const {
    return trace_metric_sum_[ind];
}

inline double JobNode::metric_sum_sq(size_t ind) const {
    return trace_metric_sum_sq_[ind];
}

inline double JobNode::metric_mean(size_t ind) const {
    return trace_metric_sum_[ind] / trace_position_[ind];
}

inline double JobNode::metric_var(size_t ind) const {
    // XXX define metric_stderr_ instead of computing every time
    double variance = (trace_metric_sum_sq_[ind] - (trace_metric_sum_[ind] * trace_metric_sum_[ind]) / trace_position_[ind]) / trace_position_[ind];
    assert(variance > -0.1);
    if (variance < 0.0)
        variance = 0.0;
    return variance;
}

inline double JobNode::metric_stderr(size_t ind) const {
    // XXX define metric_stderr_ instead of computing every time
    double variance = (trace_metric_sum_sq_[ind] - (trace_metric_sum_[ind] * trace_metric_sum_[ind]) / trace_position_[ind]) / trace_position_[ind];
    assert(variance > -0.1);
    if (variance < 0.0)
        variance = 0.0;
    return sqrt(variance / (trace_position_[ind] - 1));
}

inline bool JobNode::is_dead() const {
    ALLOW_FLOAT_EQUALITY;
    return utility() == 0.0;
    DISALLOW_FLOAT_EQUALITY;
}

inline bool JobNode::pending() const {
    return executor_ == 0;
}

inline bool JobNode::complete() const {
    return position_ == ndata_;
}

inline double JobNode::branch_probability() const {
    return branch_probability_;
}

inline bool JobNode::has_final_branch_probability() const {
    ALLOW_FLOAT_EQUALITY;
    return branch_probability_ == 0.0 || branch_probability_ == 1.0;
    DISALLOW_FLOAT_EQUALITY;
}

inline double JobNode::log_u() const {
    return log_u_;
}

inline JobNode* JobNode::child(bool isright) const {
    return child_[isright];
}


inline JobTree::id_type JobTree::next_id() const {
    return first_id_ + nodes_by_id_.size();
}

inline const JobTree::node_type* JobTree::root() const {
    return root_;
}

inline const JobNode* JobTree::cnode(id_type id) const {
    if (id - first_id_ < nodes_by_id_.size())
        return nodes_by_id_[id - first_id_];
    else if (id == 1)           // id == 1 always returns current root
        return root_;
    else
        return 0;
}

inline const JobNode* JobTree::node(id_type id) const {
    return cnode(id);
}

inline JobNode* JobTree::node(id_type id) {
    return const_cast<node_type*>(cnode(id));
}

inline SamplerType* JobTree::sampler() const {
    return completion_hook_;
}

inline bool JobTree::has_comparison_theta(id_type id) const {
    const node_type* n = node(id);
    return n && n->comparison_parent_ && n->comparison_parent_->has_theta();
}

inline const JobTree::theta_type& JobTree::comparison_theta(id_type id) const {
    const node_type* n = node(id);
    assert(n && n->comparison_parent_->has_theta());
    return n->comparison_parent_->theta();
}

inline bool JobTree::has_theta(id_type id) const {
    const node_type* n = node(id);
    return n && n->has_theta();
}

inline const JobTree::theta_type& JobTree::theta(id_type id) const {
    const node_type* n = node(id);
    assert(n && n->has_theta());
    return n->theta();
}

inline size_t JobTree::completed_depth() const {
    return first_depth_;
}

inline size_t JobTree::allocated_nodes() const {
    return allocated_nodes_;
}

inline size_t JobTree::recycled_nodes() const {
    return recycled_nodes_;
}

inline double JobTree::local_accept_rate() const {
    return accept_count_ / 100.;
}

inline void JobTree::create_children(node_type* parent) {
    create_child(parent, false);
    create_child(parent, true);
}

#endif
