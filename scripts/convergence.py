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
import numpy as np

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def convergence_stats(x):
	# x is a scalar estimand 

	(m, n) = x.shape
	
	xm = x.mean()
	xmj = x.mean(axis=1)

	# between-sequence variance
	B = ((xmj - xm)**2).sum() * n / (m - 1)

	# within-sequence variance
	W = (((x - xmj.reshape((m, 1)))**2).sum(axis=1) / (n - 1)).sum() / m

	# var(x | y) = marginal posterior variance of the estimand
	varhat = W * (n - 1) / n + B / n

	# "for most examples, values below 1.1 are acceptable"
	Rhat = np.sqrt(varhat / W)
	
	# effective number of independent draws
	neff = m * n * varhat / B
	
	print 'B:', B
	print 'W:', W
	print 'varhat:', varhat
	print 'Rhat:', Rhat
	print 'neff:', neff
	return (B, W, varhat, Rhat, neff)

x = np.random.random((4, 250))
(B, W, varhat, Rhat, neff) = convergence_stats(x)
