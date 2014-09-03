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
#include "master.hh"
#include "heap.hh"
#include "stype.hh"

JobNode::JobNode(size_t ndata, size_t scheduling_policy, double correlation,
                 size_t dimension)
    : id_(1), generation_(0),
      depth_(0), true_depth_(0), parent_(0), comparison_parent_(0),
      has_theta_(false), random_position_(0), next_random_position_(0),
      ndata_(ndata), position_(0),
      scale_(log(2.38 * 2.38 / (double) dimension)),
      log_prior_(0.0), metric_sum_(0.0), metric_sum_sq_(0.0),
      log_u_(0.0), branch_probability_(1.0),
      executor_(-1), proposal_time_(0.0), step_time_(0.0),
      step_count_(0), flip_count_(0), threshold_position_(0),
      scheduling_policy_(scheduling_policy), correlation_(correlation) {
    assert(ndata_ != 0);
    child_[0] = child_[1] = 0;
}

JobNode::JobNode(size_t id, JobNode* parent, bool isright,
                 double log_u, size_t ndata, double local_accept_rate,
                 size_t scheduling_policy, double correlation)
    : id_(id), generation_(0),
      depth_(1 + parent->depth_), true_depth_(1 + parent->depth_),
      parent_(parent),
      comparison_parent_(isright ? parent : parent->comparison_parent_),
      has_theta_(false), theta_(comparison_parent_->theta_),
      ndata_(ndata), position_(0),
      scale_(parent->scale_ + 1.0 / sqrt(parent->depth_ + 2.0) * ((double) isright - 0.234)),
      log_prior_(0.0), metric_sum_(0.0), metric_sum_sq_(0.0),
      log_u_(log_u), branch_probability_(local_accept_rate),
      executor_(-1), last_executor_(0), proposal_time_(0.0), step_time_(0.0),
      step_count_(0), flip_count_(0), threshold_position_(0),
      scheduling_policy_(scheduling_policy), correlation_(correlation) {
    assert(ndata_ == parent->ndata_);
    child_[0] = child_[1] = 0;
    parent_->child_[isright] = this;
}

double JobNode::utility() const {
    double u = 1.0;
    for (const JobNode* n = this; n->parent_; n = n->parent_) {
        double p_me = n->parent_->branch_probability_;
        if (!n->is_accept())
            p_me = 1.0 - p_me;
        if (n->parent_->complete()) {
            if (p_me < 0.5)
                return 0.0;
        } else {
            u *= p_me;
            // Make sure that utility == 0.0 only if this node is truly dead.

            // NB: It doesn't suffice to compare utility == 0.0, because of
            // Intel x86 extended precision!!
            if (u <= DBL_MIN)
                u = DBL_MIN;
        }
    }
    return u;
}

void JobNode::set_branch_probability(bool threshold, double local_accept_rate) {
    double old_branch_probability = branch_probability_;
    double delta_prior = 0.0;

    if (comparison_parent_)
        delta_prior = log_prior_ - comparison_parent_->log_prior_;

    if (!comparison_parent_)
        // root has no meaningful branch probability
        branch_probability_ = 1.0;
    else if (complete() && comparison_parent_->complete()) {
        // either accept or reject, given that the parent is on the true path
        if (delta_prior + metric_sum() - comparison_parent_->metric_sum() > log_u_)
            branch_probability_ = 1.0;
        else
            branch_probability_ = 0.0;
    } else if (scheduling_policy_ == NAIVE)
        branch_probability_ = 0.5;
    else if (scheduling_policy_ == CONSTANT)
        branch_probability_ = local_accept_rate;
    else if (scheduling_policy_ == PREDICTIVE) {
        if (threshold) {
            if (position_ == 0 || comparison_parent_->position_ == 0)
                branch_probability_ = local_accept_rate;
            else if (delta_prior + ndata_ * (metric_mean() - comparison_parent_->metric_mean()) > log_u_)
                branch_probability_ = 0.999999;
            else
                branch_probability_ = 0.000001;
        } else if (position_ < 2 || comparison_parent_->position_ < 2) {
            // `position_ == 0` means we have no data for `n`
            // and thus can't compute a meaningful branch probability.
            // `position_ == 1` means we don't yet have a variance.
            // default to recent acceptance rate.
            branch_probability_ = local_accept_rate;
        } else {
            double log_ratio_estimate, error, mm, mv, cp_mm, cp_mv;
            size_t ind;
            if (comparison_parent_->position() == position_) {
                mm = metric_mean();
                cp_mm = comparison_parent_->metric_mean();
                mv = metric_var();
                cp_mv = comparison_parent_->metric_var();
            } else if (comparison_parent_->position() > position_) {
                ind = comparison_parent_->find_closest(position_);
                mm = metric_mean();
                cp_mm = comparison_parent_->metric_mean(ind);
                mv = metric_var();
                cp_mv = comparison_parent_->metric_var(ind);
            } else {
                ind = find_closest(comparison_parent_->position());
                mm = metric_mean(ind);
                cp_mm = comparison_parent_->metric_mean();
                mv = metric_var(ind);
                cp_mv = comparison_parent_->metric_var();
            }
            log_ratio_estimate = delta_prior + ndata_ * (mm - cp_mm);
            // approximation based on corr(metric, cp->metric) ~ 0.99
            error = (ndata_ - position_) * (1. / sqrt(position_) + 1. / sqrt(ndata_ - position_)) * sqrt(mv + cp_mv - 2 * correlation_ * sqrt(mv * cp_mv));
            branch_probability_ = 0.5 + 0.5 * erf((log_ratio_estimate - log_u_) / sqrt(2) / error);
            // ensure that has_final_branch_probability() iff complete()
            if (branch_probability_ < 0.000001)
                branch_probability_ = 0.000001;
            else if (branch_probability_ > 0.999999)
                branch_probability_ = 0.999999;
            //std::cerr << "cpos: " << comparison_parent_->position() << " pos: " << position() << " bp: " << branch_probability_ << " dp: " << delta_prior << " lr: " << log_ratio_estimate << " er: " << error << "\n";
        }
    }

    if ((old_branch_probability < 0.5 && branch_probability_ > 0.5)
        || (old_branch_probability > 0.5 && branch_probability_ < 0.5))
        flip_count_ += 1;
    double delta = 0.1;
    if (((old_branch_probability < 1.0 - delta) && (branch_probability_ > 1.0 - delta))
        || ((old_branch_probability > delta) && (branch_probability_ < delta))) {
        threshold_position_ = position_;
    }
}

void JobNode::set_metric(const MetricChange& u, bool threshold, double local_accept_rate) {
    assert(!complete() || u.last_position <= position_);
    step_time_ += u.step_time;
    step_count_ += u.step_count;
    if (u.last_position > position_) {
        // update metric
        if (u.first_position == 0) {
            metric_sum_ = u.metric_sum;
            metric_sum_sq_ = u.metric_sum_sq;
        } else {
            assert(u.last_position == position_);
            metric_sum_ += u.metric_sum;
            metric_sum_sq_ += u.metric_sum_sq;
        }
        position_ = u.last_position;
        trace_metric_sum_.push_back(u.metric_sum);
        trace_metric_sum_sq_.push_back(u.metric_sum_sq);
        trace_position_.push_back(u.last_position);

        // calculate our own branch probability
        set_branch_probability(threshold, local_accept_rate);

        // calculate branch probabilities on nodes for which we are the
        // comparison parent
        for (JobNode* n = child_[1]; n; n = n->child_[0])
            n->set_branch_probability(threshold, local_accept_rate);
    }
}


void JobTree::create_child(node_type* parent, bool isright) {
    // the true root never has a left child
    if (parent->id() == root_id && !isright)
        return;
    // don't re-create children
    if (parent->child_[isright])
        return;

    // create level information if needed
    if (levels_.size() + first_depth_ <= parent->depth_ + 1) {
        assert(levels_.size() + first_depth_ == parent->depth_ + 1);
        levels_.push_back(job_level_type(u_generator_));
    }

    const job_level_type& level = levels_[parent->depth_ - first_depth_];
    node_type* n = new node_type(next_id(), parent, isright, level.log_u, ndata_,
                                 local_accept_rate(), scheduling_policy_, correlation_);
    nodes_by_id_.push_back(n);
    ++allocated_nodes_;
    n->executor_ = 0;

    // protect against wraparound (though it will never happen):
    // - no real node has id 0
    // - only the true root has id 1
    while (next_id() == 0 || next_id() == 1)
        nodes_by_id_.push_back(0);
}

JobTree::JobTree(int world_size, size_t max_depth, SamplerType* hook,
                 size_t scheduling_policy, double correlation, size_t dimension)
    : root_(0),
      first_depth_(0), first_id_(1), ndata_(0),
      require_comparison_theta_(true), threshold_(true),
      allocated_nodes_(0), recycled_nodes_(0),
      max_depth_(max_depth), reassign_ratio_(0), completion_hook_(hook),
      scheduling_policy_(scheduling_policy), correlation_(correlation),
      accept_count_(50), dimension_(dimension),
      wasted_proposal_time_(0.0), wasted_step_time_(0.0),
      useful_proposal_time_(0.0), useful_step_time_(0.0),
      wait_time_(world_size, 0.0) {
    levels_.push_back(job_level_type(u_generator_));
    if (!completion_hook_)
        completion_hook_ = new EmptySamplerType; // XXX memory leak
    while (accept_buffer_.size() < 100) {
        accept_buffer_.push_back(true);
        accept_buffer_.push_back(false);
    }
}

JobTree::~JobTree() {
    for (std::deque<JobNode*>::iterator it = nodes_by_id_.begin();
         it != nodes_by_id_.end(); ++it)
        if (*it)
            delete *it;
}

void JobTree::set_data_size(size_t ndata) {
    assert(ndata_ == 0 || ndata_ == ndata);
    assert(ndata != 0);

    if (ndata_ == 0) {
        assert(!root_ && first_id_ == 1 && nodes_by_id_.empty());
        ndata_ = ndata;
        root_ = new node_type(ndata_, scheduling_policy_, correlation_, dimension_);
        nodes_by_id_.push_back(root_);
        ++allocated_nodes_;
        root_->executor_ = 0;
    }
}

void JobTree::set_require_comparison_theta(bool x) {
    assert(!root_ && first_id_ == 1 && nodes_by_id_.empty());
    require_comparison_theta_ = x;
}

void JobTree::set_proposal(const ProposalInfo& proposal, int) {
    if (node_type* n = node(proposal.id)) {
        if (n->has_theta_ && n->theta_ != proposal.proposal)
            std::cerr << "node " << n->id_ << ": proposal conflict, " << completion_hook_->unparse(n->theta_) << " vs. " << completion_hook_->unparse(proposal.proposal) << "\n";
        assert(!n->has_theta_ || n->theta_ == proposal.proposal);
        n->has_theta_ = true;
        n->theta_ = proposal.proposal;
        n->log_prior_ = proposal.log_prior;
        n->random_position_ = proposal.random_position;
        n->next_random_position_ = proposal.next_random_position;
        n->proposal_time_ = proposal.proposal_time;
    }
}

void JobTree::set_reassign_ratio(double ratio) {
    reassign_ratio_ = ratio;
}

void JobTree::set_threshold(bool x) {
    threshold_ = x;
}

void JobTree::delete_subtree(node_type* n) {
    for (int i = 0; i != 2; ++i)
        if (n->child_[i])
            delete_subtree(n->child_[i]);
    nodes_by_id_[n->id_ - first_id_] = 0;
    wasted_proposal_time_ += n->proposal_time_;
    wasted_step_time_ += n->step_time_;
    delete n;
    --allocated_nodes_;
    ++recycled_nodes_;
}

void JobTree::delete_loser_child(node_type* n) {
    assert(n && n->complete() && n->comparison_parent_
           && n->comparison_parent_->complete()
	   && n->has_final_branch_probability());
    bool loser = n->branch_probability_ < 0.5;
    if (n->child_[loser])
        delete_subtree(n->child_[loser]);
    n->child_[loser] = 0;
}

void JobTree::trim(const run_behavior& rb,
                   const double master_wait_time, const double master_work_time,
                   const size_t iprobe_status, const size_t iprobe_no_status) {
    assert(root_->depth_ == first_depth_);
    assert(root_->complete());
    assert(!root_->child_[0] && root_->child_[1]);
    node_type* child = root_->child_[1];
    assert(child->complete());
    assert(child->depth_ == child->true_depth_);
    assert(child->has_final_branch_probability());

    bool accept = child->branch_probability_ > 0.5;
    if (!(child->child_[accept] && !child->child_[!accept])) {
        fprintf(stderr, "fuck @%zu: accept %d, brp %.8f, left %p, right %p\n",
                child->depth_, accept, child->branch_probability_,
                child->child_[0], child->child_[1]);
        std::cerr << *this;
    }
    assert(child->child_[accept] && !child->child_[!accept]);

    if (!rb.quiet()) {
        if (child->depth_ == 1) {
            std::cerr << "#depth " << "last_executor " << "allocated_nodes "
                      << "useful_proposal_time " << "useful_step_time "
                      << "wasted_proposal_time " << "wasted_step_time "
                      << "master_wait_time " << "master_work_time "
                      << "worker_wait_time "
                      << "iprobe_status " << "iprobe_no_status "
                      << "flip_count " << "threshold_position "
                      << "local_accept_rate " << "scale\n";
        }
        std::stringstream sbuf;
        sbuf << child->depth_ << " " << child->last_executor_ << " " << allocated_nodes() << " "
             << useful_proposal_time_ << " " << useful_step_time_ << " "
             << wasted_proposal_time_ << " " << wasted_step_time_ << " "
             << master_wait_time << " " << master_work_time << " "
             << std::accumulate(wait_time_.begin(), wait_time_.end(), 0.0) << " "
             << iprobe_status << " " << iprobe_no_status << " "
             << child->flip_count_ << " " << child->threshold_position_ << " "
             << local_accept_rate() << " " << child->scale_ << "\n";
        std::cerr << sbuf.str();
    }

    if (root_->depth_ == 0) {
        completion_hook_->notify(*this, root_->depth_, root_->metric_sum(),
                                 root_->theta(), true);
    }
    if (child->depth_ <= max_depth_) {
        const node_type* th_node = accept ? child : child->comparison_parent_;
        completion_hook_->notify(*this, child->depth_, th_node->metric_sum(),
                                 th_node->theta(), accept);
    }

    // remove root from lookup structure
    nodes_by_id_[root_->id_ - first_id_] = 0;
    // shift to new root
    if (accept) {
        // `child` was accepted
        useful_proposal_time_ += root_->proposal_time_;
        useful_step_time_ += root_->step_time_;
        delete root_;
        root_ = child;
        child->parent_ = child->comparison_parent_ = 0;
    } else {
        // `child` was rejected; reparent, updating root ID
        node_type* gchild = child->child_[0];
        root_->id_ = child->id_;
        root_->child_[1] = gchild;
        ++root_->depth_;
        root_->random_position_ = child->random_position_;
        root_->next_random_position_ = child->next_random_position_;
        nodes_by_id_[root_->id_ - first_id_] = root_;
        useful_proposal_time_ += child->proposal_time_;
        useful_step_time_ += child->step_time_;
        delete child;
        gchild->parent_ = root_;
    }
    --allocated_nodes_;
    ++recycled_nodes_;
    ++first_depth_;
    levels_.pop_front();

    // clean nodes_by_id_
    while (nodes_by_id_.front() == 0) {
        ++first_id_;
        nodes_by_id_.pop_front();
    }

    // update accept_count_
    accept_count_ += (int) accept;
    if (accept_buffer_.size() == 100) {
        accept_count_ -= (int) accept_buffer_.front();
        accept_buffer_.pop_front();
    }
    accept_buffer_.push_back(accept);
}

void JobTree::set_metric(const MetricChange& u, const run_behavior& rb,
                         const double master_work_time, const double master_wait_time,
                         const int rank, const size_t iprobe_status,
                         const size_t iprobe_no_status) {
    wait_time_[rank] = u.worker_wait_time;
    if (node_type* n = node(u.id)) {
        n->set_metric(u, threshold_, local_accept_rate());

        if (n->complete()) {
            n->last_executor_ = n->executor_;
            n->executor_ = -1;

            if (n->comparison_parent_ && n->comparison_parent_->complete())
                delete_loser_child(n);
            for (node_type* nn = n->child_[1]; nn; nn = nn->child_[0])
                if (nn->complete())
                    delete_loser_child(nn);

            // maybe we just completed a level
            while (root_->complete() && root_->child_[1]->complete())
                trim(rb, master_work_time, master_wait_time,
                     iprobe_status, iprobe_no_status);
        }
    } else
        wasted_step_time_ += u.step_time;
}

void JobTree::abandon(size_t id, int executor) {
    if (node_type* n = node(id)) {
        if (!n->complete() && n->executor_ == executor)
            n->executor_ = 0;
    }
}

std::pair<size_t, double> JobTree::high_utility_pending(size_t stopper) const {
    const JobNode* n = root_;
    double u = 1.0;
    boost::random::uniform_01<> dist;
    while (n && !n->pending()) {
        if (n->id() == stopper)
            return std::make_pair(size_t(0), 0.0);
        else if (n->complete())
            n = n->child_[n->branch_probability_ > 0.5];
        else {
            bool side = n->branch_probability_ >= dist(utility_generator_);
            u *= side ? n->branch_probability_ : 1 - n->branch_probability_;
            n = n->child_[side];
        }
    }
    if (n)
        return std::make_pair(n->id(), u);
    else
        return std::make_pair(size_t(0), 0.0);
}

JobInfo JobTree::execute(size_t id, int executor) {
    JobInfo info;
    node_type* node = this->node(id);
    if (node) {
        assert(node->pending());
        node->executor_ = executor;
        ++node->generation_;

        if (!require_comparison_theta_ || node->has_theta())
            create_children(node);

        info.id = node->id_;
        info.generation = node->generation_;
        info.utility = node->utility();
        if (node->has_theta()) {
            info.state = JobInfo::s_has_proposal;
            info.theta = node->theta();
            info.random_position = node->random_position_;
        } else {
            JobNode* cp = node->comparison_parent_;
            while (cp && !cp->has_theta() && cp->comparison_parent_)
                cp = cp->comparison_parent_;
            JobNode* n = node, *p;
            while ((p = n->parent_) && p != cp
                   && (!p->has_theta()
                       || p->comparison_parent_ != cp
                       || !cp->has_theta())) {
                info.path += n->is_accept() ? "a" : "r";
                n = p;
            }
            if (cp && cp->has_theta()) {
                info.state = JobInfo::s_need_proposal;
                info.theta = cp->theta();
            } else {
                info.state = JobInfo::s_need_initial;
                if (node->parent_)
                    info.path += "a"; // "accept" the initial proposal
            }
            if (p)
                info.random_position = p->next_random_position_;
            else
                info.random_position = 0;
        }
        const job_level_type& level = levels_[node->depth_ - first_depth_];
        info.log_u = level.log_u;
        info.first_position = node->position_;
        info.last_position = node->ndata_;
        info.metric_sum = node->metric_sum_;
        info.metric_sum_sq = node->metric_sum_sq_;
        info.scale = node->scale_;
    } else
        info.id = 0;
    return info;
}


static void print_node(std::ostream& s, const JobNode* node,
                       const SamplerType* sampler, int indent) {
    s << std::setw(indent) << "" << "#" << node->id()
      << " D" << node->true_depth();
    if (node->complete())
        s << (node->branch_probability() > 0.5 ? " (accept): " : " (reject): ");
    else
        s << " (" << node->branch_probability() << "): ";
    assert(node->complete() == node->has_final_branch_probability());
    if (node->has_theta())
        s << sampler->unparse(node->theta());
    else
        s << "<pending>";
    s << " (u " << node->utility();
    if (node->position()) {
        std::streamsize precision = s.precision(3);
        s << ", "
          << (100 * double(node->position()) / node->data_size()) << "%";
        s.precision(precision);
    }
    if (node->pending())
        s << ", pending";
    s << ")\n";
    for (int i = 0; i != 2; ++i)
        if (node->child(i))
            print_node(s, node->child(i), sampler, indent + 2);
}

std::ostream& operator<<(std::ostream& s, const JobTree& tree) {
    if (tree.root())
        print_node(s, tree.root(), tree.sampler(), 0);
    else
        s << "EMPTY\n";
    return s;
}


#include "masterloop.cc"

void JobTree::run_master(boost::mpi::communicator& world,
                         const run_behavior& rb) {
    run_master<boost::mpi::communicator,
               boost::mpi::request,
               boost::mpi::status>(world, rb);
}
