#include <iostream>
#include <cmath>
#include <vector>
#include <string>
#include "action.hpp"
// #include "../auxiliary_files/vector.hpp"
// #include "../auxiliary_files/read_data.hpp"
// #include "../auxiliary_files/write_scripts.hpp"
#include "../auxiliary_files/hmc.hpp"
#include <fstream>
#include "array3d.hpp"


int main ()
{
	typedef double number_type;
	typedef std::size_t size_type;

	


	// Set up Lattice action
	size_type N = 8; // number of lattice points per dimension
	size_type d = 3; // number of dimensions
	

	number_type L = 5.; // lattice size per dimension, L = N*a, with lattice spacing "a"

	number_type m = 1.;
    number_type lambda = 0.0;

    // Initialize action
	LatticeAction<number_type> RealScalar(N, d, L, m, lambda);




	//initialize HMC object:
	/*
		- from action initialized before
		- leapfrog step size (2th argument)
		- number of leapfrog steps for a given MC trajectory is drawn from uniform distribution
		  between 3. and 4. argument. This can help in case autocorrelation time for some lattice
		  sites blow up (usually not needed and thus these arguments are the same)
	*/
	HMC<number_type> sampler(RealScalar, 3e-1, 5, 5);
	
	// Define grid sites
	Vector<number_type> sites(pow(N,d));

	// Lattice will be initialized by drawing random numbers between
	// range_min and range_max (here: for each site a number between -1 and 1)
	Vector<number_type> range_min(sites.size(), -1.);
	Vector<number_type> range_max(sites.size(), 1.);



	/* Draw random vector from search region as starting point and generate field configurations */
	fill_from_region(sites, range_min, range_max);
	sampler.walk(1e3, 60*40, sites, 100, "output/correlator.txt");
	/*
		arguments:
			- number of configurations generated
			- maximal computation time in minutes
			  (camputation stops when one of the two criteria is met)
			- array of lattice sites
			- Number of times the progress of computation is updated in the console
			- Output file name for the operator expectation value(s) that one wanted to compute
			  and which is defined in "action.hpp"
	*/

	// sampler.autocorrelation(10, 120, 2,"output/autocorrel_times.txt");
	

	

	return 0;
}