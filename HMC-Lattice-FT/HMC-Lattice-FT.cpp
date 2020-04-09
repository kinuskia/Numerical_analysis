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
#include "array4d.hpp"


int main ()
{
	typedef double number_type;
	typedef std::size_t size_type;

	


	// Set up Lattice action
	size_type N = 7;
	size_type d = 4;
	number_type T = 10.;
	number_type m = 1;
	LatticeAction<number_type> RealScalar(N, d, T, m);


	// Define grid sites
	Vector<number_type> sites(pow(N,d));


	// Initialize lattice with random entries
	Vector<number_type> range_min(sites.size(), -1.);
	Vector<number_type> range_max(sites.size(), 1.);




	//initialize HMC object
	HMC<number_type> sampler(RealScalar, range_min, range_max, 5e-2, 20, 40);
	


	// /* ACTUAL RUN */



	/* Draw random vector from search region as starting point and generate Markov chain */
	fill_from_region(sites, range_min, range_max);
	sampler.walk(1e2, 60*30, sites, 100, "output/correlator.txt");

	// sampler.autocorrelation(10, 120, 2,"output/autocorrel_times.txt");
	

	

	return 0;
}