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
#include "protocol.hh"
#include "worker.hh"
#include "stype.hh"
#include "randomsequence.hh"
#include <fstream>
#include <sys/wait.h>
#include <boost/iostreams/stream.hpp>
#include <boost/iostreams/device/file_descriptor.hpp>
#include <boost/iostreams/stream_buffer.hpp>
#include <ctype.h>
#include <errno.h>

namespace {
class ShellSampler : public SamplerType {
  public:
    ShellSampler();

    int set_arguments(std::vector<const char*>& argv, bool complain);
    int work(boost::mpi::communicator& world);
    int work(boost::mpi::communicator& world, size_t batch_size,
             size_t update_style, int usleep_time);

    std::string load(const JobInfo& job);
    void prepare();
    double log_prior();
    size_t random_position() const;
    double evaluate(size_t position);

    std::string unparse(const std::string& theta) const;

  private:
    const char* initial_cmd_;
    const char* proposal_cmd_;
    const char* eval_cmd_;
    std::string proposal_;
    //RandomSequence<boost::random::mt19937> engine_;
};

ShellSampler::ShellSampler() {
}

int ShellSampler::set_arguments(std::vector<const char*>& argv, bool complain) {
    if (argv.size() < 3) {
        if (complain)
            std::cerr << "too few arguments, expected 'fetching sh INITIAL PROPOSAL EVALUATOR\n";
        return -1;
    }

    initial_cmd_ = argv[0];
    proposal_cmd_ = argv[1];
    eval_cmd_ = argv[2];
    return 0;
}

std::string run(const char* command, const std::string& input) {
    namespace io = boost::iostreams;

    /* Set up file descriptors for a bidirectional pipe.  */
    int in[2], out[2];
    int r = pipe(in);
    if (r != 0)
        return "";
    r = pipe(out);
    if (r != 0)
        return "";

    /* Create a child process.  */
    pid_t child = fork();
    if (child == 0) {
        /* Shut down file descriptors.  */
        if (out[0] != 0) {
            dup2(out[0], 0);
            close(out[0]);
        }
        if (in[1] != 0) {
            dup2(in[1], 1);
            close(in[1]);
        }
        close(out[1]);
        close(in[0]);

        /* Replace child process image with log likelihood program.  */
        const char* argv[] = { command, 0 };
        execvp(argv[0], (char**) argv);

        /* We should not be here.  */
        perror("mcmc: execvp");
        exit(1);
    } else if (child < 0)
        return "";

    /* Shut down file descriptors.  */
    close(out[0]);
    close(in[1]);

    /* Optionally send a message to child process.  */
    if (!input.empty()) {
        io::file_descriptor_sink outsink(out[1], io::never_close_handle);
        io::stream<io::file_descriptor_sink> outf(outsink);
        outf << input;
    }
    close(out[1]);

    io::file_descriptor_source insource(in[0], io::never_close_handle);
    io::stream<io::file_descriptor_source> inf(insource);
    std::stringstream insstr;
    insstr << inf.rdbuf();
    close(in[0]);

    pid_t w;
    do {
        w = waitpid(child, 0, 0);
    } while (w == -1 && errno == EINTR);
    // XXX check status

    return insstr.str();
}

std::string ShellSampler::load(const JobInfo& job) {
    if (job.need_initial())
        // initial proposal: 0.25
        proposal_ = run(initial_cmd_, "");
    else {
        proposal_ = job.theta;
        if (job.need_proposal()) {
            // XXX random seed suckery
            proposal_ = run(proposal_cmd_, proposal_);
        }
    }

    return proposal_;
    // XXX proposal.random_position = engine_.position();
}

void ShellSampler::prepare() {
}

double ShellSampler::log_prior() {
    return 0.0;
}

size_t ShellSampler::random_position() const {
    return 0;
}

double ShellSampler::evaluate(size_t) {
    // return log prob (i.e. log(exp(-proposal_ / temperature_)))
    std::string str = run(eval_cmd_, proposal_);
    const char* s = str.c_str();
    char* end;
    double prob = strtod(s, &end);
    while (end && end != s && *end && isspace((unsigned char) *end))
        ++end;
    assert(s != end && end && !*end);

    if (prob > 0)               // assume we got real probability
        prob = log(prob);
    return prob;
}

int ShellSampler::work(boost::mpi::communicator& world) {
    worker(*this, world, 1);
    return 0;
}

int ShellSampler::work(boost::mpi::communicator& world, size_t batch_size,
                       size_t update_style, int usleep_time) {
    worker(*this, world, 1, batch_size, update_style, usleep_time);
    return 0;
}

std::string ShellSampler::unparse(const std::string& theta) const {
    size_t len = theta.length();
    while (len != 0 && isspace((unsigned char) theta[len - 1]))
        --len;
    return theta.substr(0, len);
}

__attribute__ ((constructor))
static void constructor()
{
    SamplerType::signup("sh", new ShellSampler);
    SamplerType::signup("shell", new ShellSampler);
}
}
