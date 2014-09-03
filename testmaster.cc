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
#include "masterloop.cc"
#include "randomsequence.hh"
#include "time.hh"
#include "strideindex.hh"
#include "boltzmannsampler.hh"
#include "normalsampler.hh"
#include "stype.hh"
#include "worker.hh"
#include <getopt.h>
#include <string.h>
#include <iomanip>

static JobTree *tree;
static size_t max_levels = 0;

class MTRequest {
  public:
    bool test() const {
        return true;
    }
};

class MTStatus {
  public:
    MTStatus() {
    }
    MTStatus(int source, int tag)
        : source_(source), tag_(tag) {
    }
    int source() const {
        return source_;
    }
    int tag() const {
        return tag_;
    }
  private:
    int source_;
    int tag_;
};


template <typename S>
class MasterTester {
  public:
    explicit MasterTester(size_t size, size_t duration, size_t data_size,
                          size_t batch_size, const S& sampler);

    int size() const;
    void run(JobTree& tree);

    boost::optional<MTStatus> iprobe(int source, int tag);

    void recv(int source, int tag);
    void recv(int source, int tag, size_t& value);
    void recv(int source, int tag, ProposalInfo& value);
    void recv(int source, int tag, MetricChange& value);

    MTRequest isend(int source, int tag, const JobInfo& job);
    MTRequest isend(int source, int tag, const AbandonInfo& abandon);

  private:
    std::deque<MTStatus> m_types_;
    std::deque<ProposalInfo> m_proposals_;
    std::deque<MetricChange> m_metric_changes_;

    struct Worker {
        size_t id;
        size_t generation;
        size_t position;
    };
    std::vector<Worker> nodes_;
    int update_source_;
    S sampler_;

    size_t worker1_steps_;
    double serial_steps_;
    size_t probes_;
    size_t reboots_;
    size_t duration_;
    size_t data_size_;
    size_t batch_size_;
    RandomSequence<> engine_;
    double t0_;

    void print(std::ostream& s);
    void check_update();
};

template <typename S>
MasterTester<S>::MasterTester(size_t size, size_t duration, size_t data_size,
                              size_t batch_size, const S& sampler)
    : nodes_(size, Worker()), update_source_(0), sampler_(sampler),
      worker1_steps_(0), probes_(0), reboots_(0), duration_(duration),
      data_size_(data_size), batch_size_(batch_size), t0_(timestamp()) {
    assert(size >= 2 && data_size_ > 0 && batch_size_ > 0);
    for (size_t i = 1; i != size; ++i) {
        m_types_.push_back(MTStatus(i, MSG_W2M_REGISTER));
        m_types_.push_back(MTStatus(i, MSG_W2M_WANTWORK));
    }

    // serial_steps_ is the expected # of "worker1_steps_" to calculate a
    // depth
    serial_steps_ = 1 /* to fetch the work */
        + 1 /* to do the proposal */
        + ceil(double(data_size) / batch_size) /* to do the work */;
}

template <typename S>
int MasterTester<S>::size() const {
    return nodes_.size();
}

template <typename S>
void MasterTester<S>::run(JobTree& tree) {
    tree.run_master<MasterTester<S>, MTRequest, MTStatus>(*this, JobTree::run_behavior().quiet(true));
}

template <typename S>
void MasterTester<S>::print(std::ostream& s) {
    s << tree;
    for (size_t i = 0; i != nodes_.size(); ++i)
        s << "worker " << i << ": node " << nodes_[i].id << "\n";
}

static size_t count_min_utility(const JobNode* n) {
    if (!n)
        return 0;
    else
        return (n->utility() <= DBL_MIN)
            + count_min_utility(n->child(0))
            + count_min_utility(n->child(1));
}

template <typename S>
void MasterTester<S>::check_update() {
    int source = update_source_;
    if (source == 0)
        return;
    update_source_ = 0;

    Worker& node = nodes_[source];
    if (node.position == data_size_)
        m_types_.push_back(MTStatus(source, MSG_W2M_WANTWORK));
    else {
        // evaluate proposal
        size_t first = node.position;
        node.position += batch_size_;
        if (node.position > data_size_)
            node.position = data_size_;

        JobNode* jn = tree->node(node.id);
        MetricChange update(node.id, 0, node.position,
                            0, 0, 0, 0, node.position == data_size_, 0);

        if (jn && !jn->is_dead()) {
            update.metric_sum = jn->metric_sum();
            update.metric_sum_sq = jn->metric_sum_sq();

            JobInfo job;
            job.state = JobInfo::s_has_proposal;
            job.theta = jn->theta();
            job.random_position = jn->random_position();

            (void) sampler_.load(job);
            sampler_.prepare();

            sampler_.evaluate_many(first, node.position - first,
                                   update.metric_sum, update.metric_sum_sq);
        }

        m_metric_changes_.push_back(update);
        m_types_.push_back(MTStatus(source, MSG_W2M_UPDATE));
    }
}

template <typename S>
boost::optional<MTStatus> MasterTester<S>::iprobe(int rank, int tag) {
    assert(tag == MPI_ANY_TAG);
    check_update();
    if ((probes_ == duration_ && duration_ != 0)
        || (max_levels != 0 && tree->completed_depth() >= max_levels)) {
        m_types_.clear();
        for (size_t i = 0; i != nodes_.size(); ++i)
            m_types_.push_back(MTStatus(i, MSG_W2M_QUIT));
    }
    if (m_types_.front().source() != rank)
        return boost::optional<MTStatus>();
    ++probes_;
    if (probes_ % (1 << 20) == 0) {
        double delta_t = timestamp() - t0_;
        std::cerr << unparse_timestamp(delta_t) << ": "
                  << tree->completed_depth() << " completed levels, "
                  << probes_ << " messages, "
                  << (worker1_steps_ / double(serial_steps_)) << " serial levels, "
                  << tree->allocated_nodes() << "/"
                  << tree->recycled_nodes() << " nodes, "
                  << reboots_ << " reboots, "
                  << count_min_utility(tree->root()) << " min utility\n";
        std::cerr << "  " << (tree->completed_depth() * serial_steps_ / (double) worker1_steps_) << " hypothetical speedup, "
                  << (tree->completed_depth() / delta_t) << " levels/sec\n";
    }
    assert(!m_types_.empty());
    if (m_types_.front().source() == 1)
        ++worker1_steps_;
    return m_types_.front();
}

template <typename S>
void MasterTester<S>::recv(int, int) {
    check_update();
    if (!m_types_.empty())
        m_types_.pop_front();
}

template <typename S>
void MasterTester<S>::recv(int, int, size_t& value) {
    check_update();
    value = data_size_;
    m_types_.pop_front();
}

template <typename S>
void MasterTester<S>::recv(int source, int, ProposalInfo& value) {
    check_update();
    using std::swap;
    assert(!m_proposals_.empty());
    swap(value, m_proposals_.front());
    m_proposals_.pop_front();
    m_types_.pop_front();
    update_source_ = source;
}

template <typename S>
void MasterTester<S>::recv(int source, int, MetricChange& value) {
    check_update();
    using std::swap;
    swap(value, m_metric_changes_.front());
    m_metric_changes_.pop_front();
    m_types_.pop_front();
    update_source_ = source;
}

template <typename S>
MTRequest MasterTester<S>::isend(int source, int tag, const AbandonInfo& abandon) {
    assert(tag == MSG_M2W_ABANDON);
    assert(source == update_source_);
    Worker& node = nodes_[source];
    assert(node.id == abandon.id && node.generation == abandon.generation);
    if (node.position != data_size_)
        ++reboots_;
    m_types_.push_back(MTStatus(source, MSG_W2M_WANTWORK));
    update_source_ = 0;
    return MTRequest();
}

template <typename S>
MTRequest MasterTester<S>::isend(int source, int, const JobInfo& job_input) {
    JobInfo job(job_input);
    Worker& node = nodes_[source];

    ProposalInfo pi;
    ::Worker::setup_job(sampler_, job, pi);

    node.position = job.first_position;
    node.id = job.id;
    node.generation = job.generation;
    assert(job.last_position == data_size_);

    if (job.state != JobInfo::s_has_proposal) {
        tree->set_proposal(pi, source);
        m_proposals_.push_back(pi);
        m_types_.push_back(MTStatus(source, MSG_W2M_SETPROPOSAL));
    } else
        // restarting a job: send next batch
        update_source_ = source;

    return MTRequest();
}

void test_random_sequence() {
    RandomSequence<> r(18340);
    boost::random::mt19937 rm(18340);

    uint32_t r0 = r();
    uint32_t rm0 = rm();
    assert(r0 == rm0);
    assert(r.position() == 1);

    uint32_t r1 = r();
    uint32_t rm1 = rm();
    assert(r1 == rm1);
    assert(r.position() == 2);

    r.set_position(0);
    uint32_t r0b = r();
    assert(r0 == r0b);

    uint32_t r1b = r();
    assert(r1b == r1);

    for (int i = 0; i != 10384; ++i) {
        uint32_t a = r(), b = rm();
        assert(a == b);
    }

    uint32_t r10386 = r();

    r.set_position(1);
    uint32_t r1c = r();
    if (r1c != r1)
        std::cerr << r1 << " " << r1c << "/" << r0 << "\n";
    assert(r1c == r1);

    r.set_position(10386);
    uint32_t r10386b = r();
    assert(r10386b == r10386);
}

void test_proposals() {
    JobTree tree(0, 0);
    tree.set_require_comparison_theta(false);
    tree.set_data_size(2);

    JobInfo info = tree.execute(1, 1);
    assert(info.need_initial());
    assert(info.theta == "");
    assert(info.path == "");

    info = tree.execute(2, 1);
    assert(info.need_initial());
    assert(info.theta == "");
    assert(info.path == "a");
    tree.abandon(2, 1);

    info = tree.execute(3, 1);
    assert(info.need_initial());
    assert(info.theta == "");
    assert(info.path == "ra");
    tree.abandon(3, 1);

    info = tree.execute(4, 1);
    assert(info.need_initial());
    assert(info.theta == "");
    assert(info.path == "aa");
    tree.abandon(4, 1);

    tree.set_proposal(ProposalInfo(1, "R", 0, 1, 1), 1);

    info = tree.execute(2, 1);
    assert(info.need_proposal());
    assert(info.theta == "R");
    assert(info.path == "");
    tree.abandon(2, 1);

    info = tree.execute(3, 1);
    assert(info.need_proposal());
    assert(info.theta == "R");
    assert(info.path == "r");
    tree.abandon(3, 1);

    info = tree.execute(4, 1);
    assert(info.need_proposal());
    assert(info.theta == "R");
    assert(info.path == "a");
    tree.abandon(4, 1);
}

void test_proposals2() {
    JobTree tree(0, 0);
    tree.set_require_comparison_theta(false);
    tree.set_data_size(2);

    JobInfo info = tree.execute(1, 1);
    assert(info.need_initial());
    assert(info.theta == "");
    assert(info.path == "");

    info = tree.execute(2, 1);
    tree.set_proposal(ProposalInfo(2, "X", 0, 1, 1), 1);
    tree.abandon(2, 1);

    info = tree.execute(2, 1);
    assert(info.has_proposal());
    assert(info.theta == "X");
    assert(info.path == "");
    tree.abandon(2, 1);

    info = tree.execute(3, 1);
    assert(info.need_initial());
    assert(info.theta == "");
    assert(info.path == "ra");
    tree.abandon(3, 1);

    info = tree.execute(4, 1);
    assert(info.need_proposal());
    assert(info.theta == "X");
    assert(info.path == "");
    tree.abandon(4, 1);

    info = tree.execute(7, 1);
    assert(info.need_proposal());
    assert(info.theta == "X");
    assert(info.path == "r");
    tree.abandon(7, 1);

    info = tree.execute(8, 1);
    assert(info.need_proposal());
    assert(info.theta == "X");
    assert(info.path == "a");
    tree.abandon(8, 1);
}

int main(int argc, char** argv) {
    size_t cores = 1000000, messages = 1000000000, data_size = 10000,
        batch_size = 1000;
    double reassign_ratio = 1.1;
    bool threshold = false;
    //test_random_sequence();
    //test_proposals();
    //test_proposals2();

    for (int optc; (optc = getopt(argc, argv, "htd:n:m:r:l:b:")) != -1; )
        switch (optc) {
        case 'd':
            data_size = strtoull(optarg, 0, 0);
            break;
        case 'm':
            messages = strtoull(optarg, 0, 0);
            break;
        case 'n':
            cores = strtoull(optarg, 0, 0);
            break;
        case 'r':
            reassign_ratio = strtod(optarg, 0);
            break;
        case 'l':
            max_levels = strtoull(optarg, 0, 0);
            break;
        case 'b':
            batch_size = strtoull(optarg, 0, 0);
            break;
        case 't':
            threshold = true;
            break;
        case 'h':
        default:
            fprintf(stderr, "Usage: testmaster [-d DATASIZE] [-n CORES] [-m MESSAGES] [-r REASSIGNRATIO]\n                  [-l CHAINLENGTH] [-b BATCHSIZE] [-t] [boltzmann|normal]\n");
            return 1;
        }

    const char* sampler_name = "boltzmann";
    if (optind < argc && strcmp(argv[optind], "normal") == 0)
        sampler_name = "normal";

    std::cerr << "# ./testmaster -n " << cores << " -d " << data_size << " -b " << batch_size << " -r " << reassign_ratio << (threshold ? " -t " : " ") << sampler_name << "\n";

    SamplerType* stype;
    if (strcmp(sampler_name, "boltzmann") == 0)
        stype = new ThetaUnparser<double>;
    else
        stype = new ThetaUnparser<NormalTheta<2> >;

    tree = new JobTree(cores, 0, stype);
    tree->set_reassign_ratio(reassign_ratio);
    tree->set_require_comparison_theta(false);
    tree->set_threshold(threshold);

    if (strcmp(sampler_name, "boltzmann") == 0) {
        MasterTester<BoltzmannSampler> mtester(cores, messages, data_size, batch_size, BoltzmannSampler());
        mtester.run(*tree);
    } else {
        boost::random::mt19937 engine;
        boost::random::normal_distribution<> fart;
        double* data = new double[data_size * 2];

        for (size_t i = 0; i != data_size * 2; ++i)
            data[i] = fart(engine) + (i % 2 ? 100 : 200);

        NormalSampler<2> sampler(data, data_size);
        MasterTester<NormalSampler<2> > mtester(cores, messages, data_size, batch_size, sampler);
        mtester.run(*tree);
    }

    return 0;
}
