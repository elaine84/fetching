Fetching is an implementation of parallel predictive prefetching for Metropolis-Hastings. Prefetching algorithms use speculative execution to parallelize MCMC.

References
----------

Elaine Angelino, Eddie Kohler, Amos Waterland, Margo Seltzer, Ryan P. Adams. [Accelerating MCMC via parallel predictive prefetching][1]. In *30th Conference on Uncertainty in Artificial Intelligence*, UAI ’14, 2014. 

Elaine Angelino. Accelerating Markov chain Monte Carlo via parallel predictive prefetching. PhD thesis, School of Engineering and Applied Sciences, Harvard University, 2014. [Harvard version][2]. [Living version][3].


[1]: http://auai.org/uai2014/proceedings/individuals/286.pdf
[2]: http://www.eecs.harvard.edu/~elaine/thesis-harvard.pdf
[3]: http://www.eecs.harvard.edu/~elaine/thesis-living.pdf


Cloning
-------

	$ git clone git@github.com:elaine84/fetching.git


If you are using MacPorts
-------------------------

	$ port install boost +openmpi
	$ port install gsl


If you use things in scripts/
-----------------------------

	$ git clone git@github.com:yamins81/tabular.git


Examples
--------

	$ mpirun -np 2 ./fetching

	$ mpirun -np 2 ./fetching -py pyblasso


