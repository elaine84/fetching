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

mpirun -np 2 fetching -p 1> out 2> err; python summary.py o:out e:err

depth
last_executor
allocated_nodes
useful_proposal_time
useful_step_time
wasted_proposal_time
wasted_step_time
flip_count
threshold_position

"""

import os
import re
import sys

import numpy
import tabular as tb

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def make_plot(x, xname, yname, fs=14, marker='.'):
    plt.plot(x[xname], x[yname], marker)
    plt.xlabel(xname.replace('_', ' '), fontsize=fs)
    plt.ylabel(yname.replace('_', ' '), fontsize=fs)
    plt.xticks(fontsize=fs)
    plt.yticks(fontsize=fs)
    s = set(x[yname])
    if (len(s) == 2):
        d = max(s) - min(s)
        a = list(plt.axis())
        a[2] = min(s) - 0.1 * d
        a[3] = max(s) + 0.1 * d
        plt.axis(a)
        plt.title('%s = %d: %2.2f' % (yname, max(s), (x[yname]==max(s)).sum() / float(len(x))))

if __name__ == '__main__':

    dfig = '../results/'
    fout = 'out'
    ferr = 'err'
    arg_str = ''
    cmd = ''
    fs = 14

    args = sys.argv[1:]
    a = dict([arg.split(':') for arg in args])
    if ('d') in a.keys():
        dfig = a['d']
    if ('o') in a.keys():
        fout = a['o']
    if ('e') in a.keys():
        ferr = a['e']

    if os.path.exists(ferr):
        comments = []
        data = []
        f = open(ferr, 'rU')
        for line in f:
            if line.startswith('#depth'):
                header = [line]
            elif line.startswith('bash'):
                comments += ['# ' + line]
            elif line.startswith('#'):
                comments += [line]
            else:
                data += [line]
            line = line.strip('#').strip()
            if line.startswith('mpirun'):
                cmd = line
                arg_list = [arg.strip() for arg in re.split('-*', line)[1:]]
                arg_dict = dict([arg.split() if len(arg.split()) > 1 else (arg, '') for arg in arg_list])
                arg_str = '_'.join([arg.replace(' ', '=') for arg in arg_list])
        f.close()
        f = open(ferr, 'w')
        f.write(''.join(comments + header + data))
        f.close()

        assert 'np' in arg_dict.keys()
        np = int(arg_dict['np'])
        assert 'l' in arg_dict.keys()
        l = int(arg_dict['l'])

        (nr, nc) = (4, 4)

        z = tb.tabarray(SVfile=ferr, delimiter=" ", uselines=(0, l))

        plt.figure(1)
        plt.clf()
        plt.suptitle(cmd, fontsize=fs+2)

        znames = [n for n in z.dtype.names if n != 'depth']
        for (i, name) in enumerate(znames):
            plt.figure(1)
            plt.subplot(nr, nc, i + 1)
            make_plot(z, xname='depth', yname=name, fs=fs)

        plt.figure(1)
        plt.subplots_adjust(wspace=0.8, hspace=0.5)
        fig = plt.gcf()
        fig.set_size_inches(14, 10)

        if ('n=' in ferr.split('/')[-1]):
            figname = '%s.png' % (ferr.split('/')[-1])
        else:
            figname = '%s-%s-%s.png' % (ferr.split('/')[-1], arg_str, os.path.getmtime(ferr))
        ffig = os.path.join(dfig, figname)
        plt.savefig(ffig)

        ##############
        plt.figure(10)
        plt.clf()
        xname = 'depth'
        yname = 'threshold_position'
        fs = 18
        plt.plot(z[xname], z[yname], '.')
        plt.xlabel('iteration', fontsize=fs)
        plt.ylabel(yname.replace('_', ' '), fontsize=fs)
        plt.xticks(fontsize=fs)
        plt.yticks(fontsize=fs)
        s = set(z[yname])
        figname = '%s-%s.png' % (ferr.split('/')[-1], yname)
        plt.savefig(os.path.join(dfig, figname))

        ##############
        plt.figure(10)
        plt.clf()
        xname = 'depth'
        yname = 'local_accept_rate'
        fs = 18
        plt.plot(z[xname][100:], z[yname][100:], '-', linewidth=4)
        plt.xlabel('iteration', fontsize=fs)
        plt.ylabel(yname.replace('_', ' '), fontsize=fs)
        a = list(plt.axis())
        plt.plot([a[0], a[1]], [0.234, 0.234], 'k:', linewidth=2)
        plt.xticks(fontsize=fs)
        plt.yticks(fontsize=fs)
        s = set(z[yname])
        figname = '%s-%s.png' % (ferr.split('/')[-1], yname)
        plt.savefig(os.path.join(dfig, figname))

        ##############
        plt.figure(10)
        plt.clf()
        xname = 'depth'
        yname = 'scale'
        fs = 18
        plt.plot(z[xname][100:], z[yname][100:], '-', linewidth=4)
        plt.xlabel('iteration', fontsize=fs)
        plt.ylabel(yname.replace('_', ' '), fontsize=fs)
        plt.xticks(fontsize=fs)
        plt.yticks(fontsize=fs)
        s = set(z[yname])
        figname = '%s-%s.png' % (ferr.split('/')[-1], yname)
        plt.savefig(os.path.join(dfig, figname))

    if os.path.exists(fout):

        (nr, nc) = (2, 2)

        y = tb.tabarray(SVfile=fout, delimiter="\t", uselines=(0, l+1))

        plt.figure(2)
        plt.clf()
        plt.suptitle(cmd, fontsize=fs+2)

        ynames = [n for n in y.dtype.names if n not in ['depth', 'theta']]
        for (i, name) in enumerate(ynames):
            plt.subplot(nr, nc, i + 1)
            make_plot(y, xname='depth', yname=name, fs=fs)

        plt.subplots_adjust(wspace=0.5)

        if ('n=' in fout.split('/')[-1]):
            figname = '%s.png' % (fout.split('/')[-1])
        else:
            figname = '%s-%s-%s.png' % (fout.split('/')[-1], arg_str, os.path.getmtime(fout))
        ffig = os.path.join(dfig, figname)
        plt.savefig(ffig)
        #if (sys.platform == 'darwin'):
        #    os.system('open %s' % ffig)

    if os.path.exists(ferr) and os.path.exists(fout):
        print len(z), len(y)
        assert (z['depth'] == y['depth'][1:]).all()
        x = z.colstack(y[ynames][1:])

        plt.figure(3)
        plt.clf()
        plt.suptitle(cmd, fontsize=fs+2)

        names = ['useful_proposal_time', 'wasted_proposal_time',
                 'master_work_time', 'master_wait_time',
                 'worker_wait_time', 'useful_step_time', 'wasted_step_time']

        a = numpy.zeros(len(x))
        for n in names:
            a += x[n]
            plt.plot(x['depth'], a, '.')

        plt.plot(x['depth'], x['time'] * np, 'k-')

        plt.legend([n.replace('_', ' ') for n in names] + ['total time'], loc='upper left')

        fname = figname.replace('.png', '-stacked.png')
        ffig = os.path.join(dfig, fname)
        plt.savefig(ffig)
        #if (sys.platform == 'darwin'):
        #    os.system('open %s' % ffig)

        ##############
        plt.figure(11)
        plt.clf()
        useful_fraction = x['useful_step_time'] /  x['time'] / np
        ideal = x['useful_step_time'] / (x['useful_step_time'] + x['wasted_step_time'])
        plt.plot(x['time'][20:], ideal[20:], 'g:', linewidth=4)
        plt.plot(x['time'][20:], useful_fraction[20:], 'b-', linewidth=4)
        plt.legend(('ideal', 'actual'), loc='upper right')
        plt.xlabel('time (s)', fontsize=fs)
        plt.ylabel('efficiency', fontsize=fs)
        a = (0, x['time'][-1], 0.0, 1.0)
        plt.axis(a)
        plt.text(a[1] * 0.85, ideal[-100] * 1.1, 'N=%d' % np, fontsize=fs)
        plt.xticks(fontsize=fs)
        plt.yticks(fontsize=fs)
        fname = figname.replace('.png', '-efficiency-time.png')
        ffig = os.path.join(dfig, fname)
        plt.savefig(ffig)

        ##############
        plt.figure(11)
        plt.clf()
        plt.plot(x['depth'][20:], ideal[20:], 'g:', linewidth=4)
        plt.plot(x['depth'][20:], useful_fraction[20:], 'b-', linewidth=4)
        plt.legend(('ideal', 'actual'), loc='upper right')
        plt.xlabel('iteration', fontsize=fs)
        plt.ylabel('efficiency', fontsize=fs)
        plt.xticks(fontsize=fs)
        plt.yticks(fontsize=fs)
        a = (0, x['depth'][-1], 0.0, 1.0)
        plt.axis(a)
        plt.text(a[1] * 0.85, ideal[-100] * 1.1, 'N=%d' % np, fontsize=fs)
        fname = figname.replace('.png', '-efficiency-depth.png')
        ffig = os.path.join(dfig, fname)
        plt.savefig(ffig)

        ##############
        plt.figure(11)
        plt.clf()
        useful_fraction = x['useful_step_time'] /  x['time']
        ideal = x['useful_step_time'] / (x['useful_step_time'] + x['wasted_step_time']) * np
        plt.plot(x['time'][20:], ideal[20:], 'g:', linewidth=4)
        plt.plot(x['time'][20:], useful_fraction[20:], 'b-', linewidth=4)
        plt.legend(('ideal', 'actual'), loc='upper right')
        plt.xlabel('time (s)', fontsize=fs)
        plt.ylabel('speedup', fontsize=fs)
        a = list(plt.axis())
        a[:3] = (0, x['time'][-1], 0)
        plt.axis(a)
        plt.text(a[1] * 0.85, ideal[-100] * 1.1, 'N=%d' % np, fontsize=fs)
        plt.xticks(fontsize=fs)
        plt.yticks(fontsize=fs)
        fname = figname.replace('.png', '-speedup-time.png')
        ffig = os.path.join(dfig, fname)
        plt.savefig(ffig)

        ##############
        plt.figure(11)
        plt.clf()
        plt.plot(x['depth'][20:], ideal[20:], 'g:', linewidth=4)
        plt.plot(x['depth'][20:], useful_fraction[20:], 'b-', linewidth=4)
        plt.legend(('ideal', 'actual'), loc='upper right')
        plt.xlabel('iteration', fontsize=fs)
        plt.ylabel('speedup', fontsize=fs)
        plt.xticks(fontsize=fs)
        plt.yticks(fontsize=fs)
        a = list(plt.axis())
        a[:3] = (0, x['depth'][-1], 0)
        plt.axis(a)
        plt.text(a[1] * 0.85, ideal[-100] * 1.1, 'N=%d' % np, fontsize=fs)
        fname = figname.replace('.png', '-speedup-depth.png')
        ffig = os.path.join(dfig, fname)
        plt.savefig(ffig)
