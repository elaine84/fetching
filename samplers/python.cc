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
#include <boost/python.hpp>
#include <ctype.h>
#include <sstream>

namespace {
using namespace boost::python;
typedef boost::python::object pyobject;
typedef boost::python::dict pydict;

// how dare boost::python not provide this
bool hasattr(const boost::python::object& obj, const char* attr) {
    return PyObject_HasAttrString(obj.ptr(), attr);
}

class Rand {
  public:
    inline Rand(RandomSequence<boost::random::mt19937>& engine)
        : engine_(engine) {
    }
    inline double random() {
        boost::random::uniform_01<> dist;
        return dist(engine_);
    }
    inline double uniform(double a, double b) {
        boost::random::uniform_real_distribution<> dist(a, b);
        return dist(engine_);
    }
    inline double normalvariate(double mu, double sigma) {
        boost::random::normal_distribution<> dist(mu, sigma);
        return dist(engine_);
    }
    inline double lognormvariate(double mu, double sigma) {
        boost::random::lognormal_distribution<> dist(mu, sigma);
        return dist(engine_);
    }
  private:
    RandomSequence<boost::random::mt19937>& engine_;
};

class PythonSampler : public SamplerType {
  public:
    PythonSampler();

    int set_arguments(std::vector<const char*>& argv, bool complain);
    int work(boost::mpi::communicator& world, size_t batch_size,
             size_t update_style, int usleep_time);

    std::string load(const JobInfo& job);
    void prepare();
    double log_prior();

    size_t random_position() const;
    double evaluate(size_t position);

    std::string unparse(const std::string& theta) const;

  private:
    pydict main_namespace_;
    pyobject module_;
    pydict module_namespace_;
    pyobject pickle_module_;

    pyobject sampler_;
    pyobject sampler_eval_;
    pyobject theta_;

    bool has_prepare_;
    bool has_log_prior_;
    bool has_data_size_;
    RandomSequence<boost::random::mt19937> engine_;
    Rand rand_;
};

PythonSampler::PythonSampler()
    : rand_(engine_) {
}

int PythonSampler::set_arguments(std::vector<const char*>& argv, bool complain) {
    if (argv.size() < 1) {
        if (complain)
            std::cerr << "too few arguments, expected 'fetching py MODULE\n";
        return -1;
    }

    if (!Py_IsInitialized()) {
        Py_Initialize();
        // XXX add `pysamplers` directory to path
        PyObject* sysPath = PySys_GetObject((char*) "path");
        PyList_Insert(sysPath, 0, PyString_FromString("pysamplers"));
    }

    try {
        main_namespace_ = extract<pydict>(boost::python::import("__main__").attr("__dict__"));
        module_ = boost::python::import(argv[0]);
        module_namespace_ = extract<dict>(module_.attr("__dict__"));

        module_namespace_["RandClass"] = class_<Rand>("Rand", no_init)
            .def("random", &Rand::random)
            .def("uniform", &Rand::uniform)
            .def("normalvariate", &Rand::normalvariate)
            .def("lognormvariate", &Rand::lognormvariate);
        module_namespace_["rand"] = ptr(&rand_);

        pickle_module_ = boost::python::import("cPickle");

        boost::python::list args;
        for (size_t i = 1; i != argv.size(); ++i)
            args.append(argv[i]);
        sampler_ = module_namespace_["sampler"](args);

        sampler_eval_ = sampler_.attr("evaluate");
        has_prepare_ = hasattr(sampler_, "prepare");
        has_log_prior_ = hasattr(sampler_, "log_prior");
        has_data_size_ = hasattr(sampler_, "data_size");
    } catch (error_already_set const&) {
        PyErr_Print();
        return -1;
    }

    return 0;
}

std::string PythonSampler::load(const JobInfo& job) {
    engine_.set_position(job.random_position);
    std::string proposal;
    try {
        if (job.need_initial())
            theta_ = sampler_.attr("first_proposal")();
        else {
            theta_ = pickle_module_.attr("loads")(job.theta);
            if (job.need_proposal())
                theta_ = sampler_.attr("next_proposal")(theta_, job.scale);
        }
        proposal = extract<std::string>(pickle_module_.attr("dumps")(theta_));
    } catch (const error_already_set&) {
        PyErr_Print();
    }
    return proposal;
}

void PythonSampler::prepare() {
    if (has_prepare_)
        try {
            sampler_.attr("prepare")();
        } catch (const error_already_set&) {
            PyErr_Print();
        }
}

double PythonSampler::log_prior() {
    if (has_log_prior_)
        try {
            return extract<double>(sampler_.attr("log_prior")(theta_));
        } catch (const error_already_set&) {
            PyErr_Print();
        }
    return 0.0;
}

size_t PythonSampler::random_position() const {
    return engine_.position();
}

double PythonSampler::evaluate(size_t position) {
    try {
        if (has_data_size_)
            return extract<double>(sampler_eval_(theta_, position));
        else
            return extract<double>(sampler_eval_(theta_));
    } catch (const error_already_set&) {
        PyErr_Print();
        return 0;
    }
}

int PythonSampler::work(boost::mpi::communicator& world, size_t batch_size,
                        size_t update_style, int usleep_time) {
    size_t data_size = 1;
    if (has_data_size_)
        try {
            data_size = extract<size_t>(sampler_.attr("data_size")());
        } catch (const error_already_set&) {
            PyErr_Print();
            return -1;
        }
    worker(*this, world, data_size, batch_size, update_style, usleep_time);
    return 0;
}

std::string PythonSampler::unparse(const std::string& theta) const {
    try {
        pyobject proposal = pickle_module_.attr("loads")(theta);
        if (hasattr(sampler_, "unparse_proposal"))
            proposal = sampler_.attr("unparse_proposal")(proposal);
        return extract<std::string>(boost::python::str(proposal));
    } catch (const error_already_set&) {
        PyErr_Print();
        return "<error>";
    }
}

__attribute__ ((constructor))
static void constructor()
{
    SamplerType::signup("py", new PythonSampler);
    SamplerType::signup("python", new PythonSampler);
}
}
