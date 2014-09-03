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
 
 """
import os
import re
import sys

import tabular as tb
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def overview(fbase):
    ferr = fbase + '.err'
    fout = fbase + '.out'
    os.system('echo "%s" ' \
    '"$(grep \'cloud\' %s | awk \'{ print $4 }\' | sort -u | wc -l)" ' \
    '"$(%s | awk -F\'\t\' \'{ print $5 }\' | md5sum)" >> overview' % (fbase, ferr, fout))

if __name__ == '__main__':

	dfig = '../results/'

	args = sys.argv[1:]

	if (len(args) > 0):
		dir = args[0]
	else:
		dir = '.'

	flist = [os.path.join(dir, f[:-4]) for f in os.listdir(dir) if f.endswith('.out')]
	os.system('echo "" > overview')
	for f in flist:
		overview(f)

	flist = [f + '.err' for f in flist]

	recs = []

	for f in flist:
		x = open(f, 'rU').read().strip().split('\n')

		m = f.split('-')[1:]
		meta = [i.split('=') if len(i.split('='))==2 else (i, '') for i in m[:-2]]
		id = ' '.join([i for i in m[:-2] if not i.startswith('n=')])
		meta = meta + [['rev', m[-2]], ['tstamp', m[-1]]]
		mdict = dict(meta)
		meta = [(i, float(j)) if j.replace('.', '').isdigit() else (i, j) for (i, j) in meta]
		print f

		kvp = []
		keys = []
		vals = []
		record = [('id', id)] + meta
		for line in x:
			if line.startswith('OK'):
				line = line.split('OK')[1].strip()
				keys = [i for i in line.strip().split() if re.match('\D+', i) and i != '=']
				vals = [float(i) for i in line.strip().split() if re.match('\d+', i)]
				record += zip(keys, vals)
				break

		g =f.replace('err', 'out')
		y = tb.tabarray(SVfile=g, uselines=(0, int(mdict['l'])+1), delimiter='\t')
		assert len(record) > 1, f
		record += [('time', y[-1]['time'])]
		recs += [record]

	x = tb.tabarray(kvpairs=recs)
	x.sort(order=['s', 'y', 'b'])
	x = x[::-1]
	id_vec = tb.utils.uniqify(x['id'])
	bvec = list(set([int(i[2:]) for i in ' '.join(id_vec).split() if i.startswith('b=')]))
	bvec.sort()

	code_vec = [['s=0', 'y=0'], ['s=0', 'y=1'], ['s=2', 'y=0'], ['s=2', 'y=1']]

	plt.figure(1)
	plt.clf()
	plt.figure(2)
	plt.clf()
	plt.figure(3)

	cvec = ['orange', 'b', 'm', 'g', 'c', 'k']
	svec = ['solid', 'dashed', 'dashdot', 'dotted']
	legend = []

	for id in id_vec:
		b = [int(i[2:]) for i in id.split(' ') if i.startswith('b=')][0]
		s = svec[bvec.index(b)]
		i = [j for j in range(len(code_vec)) if all([c in id for c in code_vec[j]])][0]
		y = x[x['id'] == id]
		y.sort(order=['n'])
		plt.figure(1)
		baseline = y[y['n']==y['n'].min()]['time'].mean()
		print id, baseline, zip(y['n'], baseline / y['time'])
		plt.plot(y['n'], baseline / y['time'], 'o', linestyle=s, color=cvec[i], markersize=10, linewidth=3)
		plt.figure(2)
		plt.plot(y['n'], y['time'],  'o', linestyle=s, color=cvec[i], markersize=10, linewidth=3)
		legend += [id.replace('r=0', 'r=$\infty$')]

	fs = 16

	nvec = list(set(x['n']))
	nvec.sort()

	plt.figure(1)
	plt.xlabel('Number of cores', fontsize=fs)
	plt.ylabel('Speedup relative to baseline', fontsize=fs)
	plt.xticks(nvec, fontsize=fs)
	plt.yticks(fontsize=fs)
	plt.axis([0, nvec[-1] + 2, 0, 100])
	plt.legend(legend, loc='upper left')
	plt.title('Speedup as a function of the number of cores', fontsize=fs)
	plt.savefig(os.path.join(dfig, 'speedup.png'))

	plt.figure(2)
	plt.xlabel('Number of cores', fontsize=fs)
	plt.ylabel('Time (s)', fontsize=fs)
	plt.xticks(nvec, fontsize=fs)
	plt.yticks(fontsize=fs)
	a = list(plt.axis())
	a[0] = 0
	a[1] = nvec[-1] + 2
	plt.axis(a)
	plt.legend(legend, loc='upper right')
	plt.title('Time as a function of the number of cores', fontsize=fs)
	plt.savefig(os.path.join(dfig, 'time.png'))

