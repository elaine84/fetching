"""
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

Bayesian LASSO

With OPV data.

"""

import cPickle
import os
import sys

import numpy as np


class BLasso:

    def __init__(self, train_size=600000, test_size=10000, tau=5.0, quiet=True, seed=0, jump=None):

        self.seed = seed
        self.jump = jump

        if os.getcwd().startswith('/Users/'):
            fr = '../../../data/opv/opv-response.npy'
            ff = '../../../data/opv/opv-features.npy'
        else:
            fr = '../../data/opv/opv-response.npy'
            ff = '../../data/opv/opv-features.npy'

        y = np.load(fr)
        X = np.load(ff)
        try:
            seed = int(rand.random() * 10**12)
            np.random.seed(seed)
        except:
            print 'Did not use speculative rand.random()'
        ind = np.random.permutation(train_size + test_size)
        y = y[ind]
        X = X[ind]

        cols = [j for j in range(X.shape[1]) if X[:train_size, j].std() > 0.001]
        if not quiet:
            print 'Keeping %d / %d columns' % (len(cols), X.shape[1])
        X = X[:, cols]

        X = (X - X.mean(axis=0)) / X.std(axis=0)
        y = y - y.mean()

        self.data = X[:train_size]
        self.test_data = X[train_size:(train_size + test_size)]
        self.response = y[:train_size]
        self.test_response = y[train_size:(train_size + test_size)]

        (self.data_size_, self.ndim) = self.data.shape
        self.data_size_ = self.data_size_ / 18000
        self.tau = tau
        self.permutation = None
        self.theta_dim = self.ndim + 1
        self.jump_cov = np.diag(np.ones(self.theta_dim) * 2.0e-08)

        self.mle = None
        self.lr = None


    def data_size(self):
        return self.data_size_

    def first_proposal(self):
        theta = np.array([  2.00000000e-01,
         1.40579036e+00,   3.60878266e-01,   4.94455574e-01,
        -2.74122043e-02,   1.68881900e-01,   9.20233386e-01,
         3.87164277e-01,  -8.06430469e-02,   2.78239208e-01,
        -1.54895307e+00,   1.43484507e-01,   8.16043012e-04,
         3.14120798e-02,  -1.53496686e-01,  -9.92468016e-03,
        -1.25028385e+00,  -8.51259409e-01,  -3.71299715e-01,
        -1.68059888e+00,  -8.43351644e-03,  -7.21228079e-01,
         3.10229542e-02,   3.90164125e-01,  -3.71647066e-01,
        -1.02620741e+00,  -2.23799427e-03,   4.10134792e-02,
        -3.47502563e-02,  -5.92303867e-03,  -3.44586726e-02,
         9.20233382e-01,   7.38275802e-02,  -9.92467997e-03,
         5.84726376e-01,   1.05192495e+00,   9.02679365e-03,
        -1.51799484e+00,  -2.65142512e-02,   7.68202217e-03,
         4.69544161e+00,  -7.20872443e-01,   7.73013464e-01,
        -3.18937138e+00,  -6.13193416e-02,  -6.19170226e-02,
        -1.06369593e-03,  -2.23766532e+00,   1.84532898e-02,
         2.51643690e-01,   8.52457449e-01,   8.52457449e-01,
        -5.43505858e-01,   5.10997534e-01,  -2.85750869e-01,
        -4.08767576e-02,  -1.03239989e-01])
        #theta = np.zeros(self.theta_dim)
        np.random.seed(self.seed)
        #theta = np.random.normal(0, 0.5, self.theta_dim)
        if (self.jump is not None):
            theta += np.random.normal(0, self.jump, self.theta_dim)
            theta[0] = np.abs(theta[0])
        return theta

    def next_proposal(self, theta, scale):
        try:
            seed = int(rand.random() * 10**12)
            np.random.seed(seed)
        except:
            print 'Did not use speculative rand.random()'
        while (1):
            jump = np.random.normal(0.0, np.sqrt(np.exp(scale)), self.theta_dim)
            proposal = theta + jump
            if (proposal[0] > 0.0):
                return proposal

    def log_prior(self, theta):
        sigma = theta[0]
        beta = theta[1:]

        # p(sigma^2) = 1 / sigma^2
        log_prior_sigma = -2.0 * np.log(sigma)

        # p(beta | sigma^2) = (tau / 2 / sigma)^d * exp(-tau |beta|_1 / sigma)
        log_laplace = (self.ndim * np.log(self.tau / 2 / sigma)
                       - (self.tau * np.abs(beta).sum()) / sigma)
        return log_prior_sigma + log_laplace

    def evaluate(self, theta, position):
        X = self.data[(position * 18000):((position + 1) * 18000)]
        y = self.response[(position * 18000):((position + 1) * 18000)]
        sigma = theta[0]
        beta = theta[1:]
        lp = (np.log(1 / sigma / np.sqrt(2 * np.pi)) - (np.dot(beta, X.T) - y)**2 / 2 / sigma**2).sum()
        return lp
    """
    def evaluate(self, theta, position):
        X = self.data[position]
        y = self.response[position]
        sigma = theta[0]
        beta = theta[1:]
        lp = (np.log(1 / sigma) - (np.dot(beta, X) - y)**2 / 2 / sigma**2).sum()
        return lp
    """

    def log_contributions(self, theta):
        sigma = theta[0]
        beta = theta[1:(self.ndim + 1)]
        return np.log(1 / sigma / np.sqrt(2 * np.pi)) - (np.dot(beta, self.data.T) - self.response)**2 / 2 / sigma**2

    def log_posterior(self, theta):
        return self.log_prior(theta) + self.log_contributions(theta).sum()

    def unparse_proposal(self, theta):
        return theta.tolist().__repr__()

    def check(self, MAKE_PLOT=False):
        import sklearn.linear_model as lm
        """
        r = np.zeros((self.data_size_, self.ndim * self.ndim))
        t = np.zeros((len(self.test_data), self.ndim * self.ndim))

        for i in range(self.ndim):
            for j in range(self.ndim):
                r[:,(i+j)] = self.data[:,i] * self.data[:,j]
                t[:,(i+j)] = self.test_data[:,i] * self.test_data[:,j]
        """     
        lr = lm.Ridge(alpha = 1.0)
        lr.fit(self.data, self.response)
        q1 = lr.predict(self.data)
        q2 = lr.predict(self.test_data)
        self.lr = lr
        if MAKE_PLOT:
            import matplotlib
            import matplotlib.pyplot as plt
            plt.ion()
            plt.figure(4)
            plt.clf()
            plt.plot(q1, self.response, '.')
            plt.plot(q2, self.test_response, '.')

    def sequential(self, nsteps=100, PLOT_TRACES=False):
        if PLOT_TRACES:
            import matplotlib
            import matplotlib.pyplot as plt
            plt.ion()
            plt.figure(2)
            plt.clf()
            plt.plot([1, self.data_size()], [0, 0], 'k-', linewidth=2)

        theta = np.zeros(self.ndim + 1)
        theta[0] = 1.0 #np.random.normal(0.0, 10, 1)
        theta[1:] = np.random.normal(0.0, 10, self.ndim)
        self.check(PLOT_TRACES)
        theta[1:] = self.lr.coef_
        #theta[1:] = 0.0
        #self.beta_scale = 0.0008
        log_theta_density = self.log_prior(theta)
        theta_trace = np.zeros(self.data_size())
        for position in range(self.data_size()):
            log_theta_density += self.evaluate(theta, position)
            theta_trace[position] = log_theta_density / (position + 1) * self.data_size()
        accept = 0
        accept_vec = np.zeros(100)
        mle = theta
        log_mle_density = log_theta_density
        for i in range(1, nsteps + 1):
            seed = int(np.random.random() * 10**12)
            np.random.seed(seed)
            log_u = np.log(np.random.random())
            proposal = self.next_proposal(theta)
            log_proposal_density = self.log_prior(theta)
            trace = np.zeros(self.data_size())
            # perm = np.random.permutation(self.data_size())
            perm = np.arange(self.data_size())

            for position in range(self.data_size()):
                log_proposal_density += self.evaluate(proposal, perm[position])
                trace[position] = log_proposal_density / (position + 1) * self.data_size()
            if PLOT_TRACES and (i > 1):
                plt.figure(2)
                estimate = trace - theta_trace - log_u
                guess = estimate > 0
                changes = guess[1:] != guess[:-1]
                nz = changes.sum() - changes.cumsum()
                plt.figure(2)
                plt.semilogx(estimate, '-', alpha=0.2)
                #plt.axis([0, self.data_size(), -10000, 10000])
                plt.draw()
            if (log_proposal_density - log_theta_density > log_u):
                theta = proposal.copy()
                log_theta_density = log_proposal_density
                theta_trace = trace.copy()
                accept += 1
                accept_vec[i % 100] = 1
                print i, 'accept', log_theta_density, float(accept) / i, accept_vec.sum() / 100.
            else:
                accept_vec[i % 100] = 0
                print i, 'reject', log_theta_density, float(accept) / i, accept_vec.sum() / 100.
            if (log_theta_density > log_mle_density):
                mle = theta.copy()
                log_mle_density = log_theta_density
            print theta
        self.mle = mle
        return mle


def sampler(args):
    train_size = 100000
    test_size = 10000
    seed = 0
    a = dict([arg.split('=') for arg in args])
    if ('train_size') in a.keys():
        train_size = int(a['train_size'])
    if ('test_size') in a.keys():
        test_size = int(a['test_size'])
    if 'seed' in a.keys():
        seed = int(a['seed'])
    if 'jump' in a.keys():
        jump = float(a['jump'])
    return BLasso(train_size=train_size, test_size=test_size, seed=seed, jump=jump)

if __name__ == '__main__':
    nsteps = 1000
    PLOT_TRACES = False
    args = sys.argv[1:]
    a = dict([arg.split('=') for arg in args])
    if 'nsteps' in a.keys():
        nsteps = int(a['nsteps'])
    if 'PLOT_TRACES' in a.keys():
        PLOT_TRACES = bool(a['PLOT_TRACES'])
    s = sampler(args)
    mle = s.sequential(nsteps=nsteps, PLOT_TRACES=PLOT_TRACES)
    try:
        import matplotlib
        import matplotlib.pyplot as plt
        plt.ion()
        plt.figure(3)
        plt.clf()
        plt.plot(np.dot(s.data, mle[1:]), s.response, '.')
        plt.plot(np.dot(s.test_data, mle[1:]), s.test_response, '.')
    except:
        print 'Failed to execute test'
