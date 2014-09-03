Fetching is an implementation of parallel predictive prefetching for Metropolis-Hastings. Prefetching algorithms use speculative execution to parallelize MCMC.

References
==========

[Accelerating MCMC via parallel predictive prefetching][1]
Elaine Angelino, Eddie Kohler, Amos Waterland, Margo Seltzer, Ryan P. Adams
UAI 2014: 30th Conference on Uncertainty in Artificial Intelligence

PhD thesis: [Accelerating Markov chain Monte Carlo via parallel predictive prefetching][2]
Author: Elaine Angelino
Advisors: Margo Seltzer and Ryan P. Adams
	

[1]: http://auai.org/uai2014/proceedings/individuals/286.pdf
[2]: http://www.eecs.harvard.edu/~elaine/draft.pdf


Cloning
=======

	$ git clone git@github.com:elaine84/repulsive.git


If you are using MacPorts
=========================

	$ port install boost +openmpi
	$ port install gsl


If you use things in scripts/
=============================

	$ git clone git@github.com:yamins81/tabular.git


Examples
========

	$ openmpi -np 2 ./fetching

	$ openmpi -np 2 ./fetching -py pyblasso


