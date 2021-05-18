#include <iostream>
#include <cmath>
#include <vector>
#include <string>
#include "model.hpp"
#include "../auxiliary_files/read_data.hpp"
#include "../auxiliary_files/vector.hpp"
#include "../auxiliary_files/write_scripts.hpp"
#include "../auxiliary_files/hmc_fit.hpp"
#include <fstream>


int main (int argc, char* argv[])
{
	typedef double number_type;
	typedef std::size_t size_type;

	// Initialize empty vectors for experimental input
	std::vector<Vector<number_type>> experiment(3);
	

	// Read in measured data
	read_data("experimental_data/measurements.txt", experiment);
	Vector<number_type> t(experiment[0]);
	Vector<number_type> y(experiment[1]);
	Vector<number_type> dy(experiment[2]);


	// Set up fitting model
	Model<number_type> gaussian(
	t,
	y,
	dy
	);


	// Vector for fitting parameters
	Vector<number_type> popt(3);


	// Estimated search region
	Vector<number_type> range_min(popt.size());
	Vector<number_type> range_max(popt.size());
	// Characteristic length scales for the parameters // default 1
	Vector<number_type> c_lengths(popt.size(), 1);
	// characteristic length scales are here relative to range_max-range-min ...
	range_min[0] = 4;
	range_max[0] = 6;
	c_lengths[0] = 1;
	range_min[1] = 3;
	range_max[1] = 5;
	c_lengths[1] = 1;
	range_min[2] = 0.2;
	range_max[2] = 2;
	c_lengths[2] = 1.0;


	

	
	
	// ... and are now made absolute
	c_lengths *= (range_max-range_min);


	//initialize HMC object
	HMC<number_type> sampler(gaussian, range_min, range_max, c_lengths, 3e-3, 90, 130, 1e1);
	
	

	/* PRELIMINARY RUN TOOLS */
	
	sampler.tune_parameters();


	/* 	
		Routine for computer clusters: Produce Markov chain of specific length and save to output file
		of the following format dataX.txt, where X is given in the console when executing, e.g. ./fitter 5
	*/
	std::string filenumber;
	if (argc > 1)
	{
		filenumber = argv[1];
	}
	//sampler.walk_silently(2e3, "data", filenumber);
	//sampler.walk_silently_disregarding(1e4, 40.225, "data", filenumber);

	

	


	/* ACTUAL RUN */



	// Draw random vector from search region as starting point and generate Markov chain
	fill_from_region(popt, range_min, range_max);
	//sampler.walk(1e4, 60*30, popt, 10);

	

	

	return 0;
}