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

Based on an 8-dimensional Gaussian mixture model from Robert Nishihara.

From http://arxiv.org/abs/1210.7477:

"An eight-component mixture of Gaussians in eight dimensions. Each component is
a spherical Gaussian with unit variance. The components are distributed
uniformly at random within a hypercube of edge length four."

"""

import sys

import numpy as np

class Mixture:
    def __init__(self, train_size=1000, seed=0):
        np.random.seed(0)
        self.seed = seed
        self.train_size_ = train_size
        self.scale = 0.01
        self.length = 4.0
        pos = np.array([[ 0.2456,  0.8211,  0.3065,  0.9171,
                          0.9674,  0.5055,  0.535 ,  0.7781],
                        [ 0.1852,  0.774 ,  0.9248,  0.8285,
                          0.7948,  0.460 ,  0.9904,  0.6430],
                        [ 0.7135,  0.8969,  0.7882,  0.7179,
                          0.8707,  0.1549,  0.364 ,  0.7309],
                        [ 0.3507,  0.8099,  0.0669,  0.2366,
                          0.7635,  0.5878,  0.5188,  0.7846],
                        [ 0.186 ,  0.3913,  0.7746,  0.3846,
                          0.1483,  0.4110,  0.5936,  0.5528],
                        [ 0.2550,  0.7924,  0.5779,  0.5291,
                          0.2643,  0.7684,  0.3859,  0.9556],
                        [ 0.3698,  0.1247,  0.1504,  0.8657,
                          0.9061,  0.2281,  0.9170,  0.9552],
                        [ 0.354 ,  0.3176,  0.2076,  0.0267,
                          0.6507,  0.0931,  0.2434,  0.2387]])
        self.positions = self.length * pos - (self.length / 2.0)
        (self.ncomponents, self.ndimensions) = pos.shape
        self.cov = np.diag(np.ones(self.ndimensions))
        self.data = np.zeros((self.train_size_, self.ndimensions))
        components = np.cast[int](np.random.random(self.train_size_) * self.ncomponents)
        for i in range(self.ncomponents):
            ind = (components == i)
            self.data[ind, :] = np.random.multivariate_normal(self.positions[i],
                                                              self.cov, int(ind.sum()))
        self.train_size_ = train_size / 10000

    def data_size(self):
        return self.train_size_

    def first_proposal(self):
        #return np.zeros((self.ncomponents, self.ndimensions))
        np.random.seed(self.seed)
        pos = np.random.random((self.ncomponents, self.ndimensions))
        return self.length * pos - (self.length / 2.0)

    def next_proposal(self, theta, scale):
        try:
            seed = int(rand.random() * 10**12)
            np.random.seed(seed)
        except:
            print 'Did not use speculative rand.random()'
        std = np.sqrt(np.exp(scale))
        theta = theta[np.random.permutation(self.ncomponents)] # permute labels
        jump = np.random.normal(0.0, std, (self.ncomponents, self.ndimensions))
        return theta + jump

    def evaluate(self, theta, position):
        """
        datum = self.data[position].repeat(self.ncomponents).reshape(self.ncomponents, self.ndimensions).T
        td = theta - datum
        log_density = np.log(np.exp(-(td * td).sum(axis=1) / 2).sum())
        return log_density
        """
        lp = np.zeros(10000)
        for i in range(self.ncomponents):
            td = theta[i] - self.data[(position * 10000):((position + 1) * 10000)]
            lp += np.exp(-(td**2) / 2.).prod(axis=1)
        return np.log(lp).sum()

    def log_contributions(self, theta):
        lp = np.zeros(len(self.data))
        for i in range(self.ncomponents):
            td = theta[i] - self.data
            lp += np.exp(-(td**2) / 2.).prod(axis=1)
        return np.log(lp)

    def log_posterior(self, theta):
        lp = 0
        for i in range(self.ncomponents):
            td = theta[i] - self.data
            lp += np.exp(-(td**2) / 2.).prod(axis=1)
        return np.log(lp).sum()

    def unparse_proposal(self, theta):
        return theta.tolist().__repr__()

    def sequential(self, nsteps=100, PLOT_TRACES=False):
        if PLOT_TRACES:
            import matplotlib
            import matplotlib.pyplot as plt
            plt.ion()
            plt.figure(2)
            plt.clf()
            plt.semilogx([1, self.train_size()], [0, 0], 'k-', linewidth=2)

        theta = self.first_proposal()
        log_theta_density = 0
        for position in range(self.train_size()):
            log_theta_density += self.evaluate(theta, position)
        accept = 0
        accept_vec = np.cast[int](np.random.rand(100) * 2)
        samples = np.zeros((nsteps, self.ncomponents * self.ndimensions))
        for i in range(1, nsteps + 1):
            seed = int(np.random.random() * 10**12)
            np.random.seed(seed)
            log_u = np.log(np.random.random())
            mu = np.zeros(self.ndimensions)
            cov = self.scale * np.diag(np.ones(self.ndimensions))
            proposal = theta.copy()
            proposal += np.random.multivariate_normal(mu, cov, self.ncomponents)
            log_proposal_density = 0
            trace = np.zeros(self.train_size())
            # perm = np.random.permutation(self.train_size())
            for position in range(self.train_size()):
                log_proposal_density += self.evaluate(proposal, position)
                trace[position] = log_proposal_density / (position + 1) * self.train_size()
            if PLOT_TRACES:
                plt.figure(2)
                plt.semilogx(trace - log_theta_density - log_u, '-', alpha=0.2)
                plt.axis([0, self.train_size(), -100000, 100000])
                plt.draw()
            if (log_proposal_density - log_theta_density > log_u):
                theta = proposal.copy()
                log_theta_density = log_proposal_density
                accept += 1
                accept_vec[i % 100] = 1
                print i, 'accept', log_theta_density, float(accept) / i, accept_vec.sum() / 100.
            else:
                accept_vec[i % 100] = 0
                print i, 'reject', log_theta_density, float(accept) / i, accept_vec.sum() / 100.
            samples[i-1] = theta.flatten()

        try:
            import matplotlib
        except:
            print "Couldn't import matplotlib, not plotting."
            return

        import time
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt

        plt.figure(1)
        plt.clf()

        for j in range(1, 5):
            plt.figure(1)
            plt.subplot(2, 2, j)
            x = (j - 1) * 2
            y = (j * 2) - 1
            plt.plot(self.data[:,x][:1000], self.data[:,y][:1000], 'b.', alpha=0.1)
            plt.plot(self.positions[:,x], self.positions[:,y], 'ko')
            for k in range(self.ncomponents):
                plt.plot(samples[:, k * self.ncomponents + x],
                         samples[:, k * self.ncomponents + y])
            plt.xlabel(x)
            plt.ylabel(y)

        figname = ('../results/%s-%s=%d-%s=%d-%d.png' %
                    ('pymixture', 'd', self.train_size(), 'D', nsteps, time.time()))
        plt.savefig(figname)
        try:
            import os
            os.system('open %s' % figname)
        except:
            pass

def sampler(args):
    train_size = 500000
    seed = 0
    a = dict([arg.split('=') for arg in args])
    if 'train_size' in a.keys():
        train_size = int(a['train_size'])
    if 'seed' in a.keys():
        seed = int(a['seed'])
    return Mixture(train_size=train_size, seed=seed)

if __name__ == '__main__':
    nsteps = 10
    PLOT_TRACES = True
    args = sys.argv[1:]
    a = dict([arg.split('=') for arg in args])
    if 'nsteps' in a.keys():
        nsteps = int(a['nsteps'])
    if 'PLOT_TRACES' in a.keys():
        PLOT_TRACES = bool(a['PLOT_TRACES'])
    sampler(args).sequential(nsteps=nsteps, PLOT_TRACES=PLOT_TRACES)
