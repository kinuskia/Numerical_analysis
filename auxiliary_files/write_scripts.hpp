#ifndef WRITE_SCRIPTS_HPP
#define WRITE_SCRIPTS_HPP

#include <fstream>
#include "vector.hpp"

typedef std::size_t size_type;

/* free functions to write n bash scripts */

std::string int_to_string(size_type index, size_type n_scripts)
{
	assert(n_scripts > 0);
	// Step 1: Determine the number of digits to display string number
	int n_digits = log(1.0*n_scripts)/log(10.) + 1;

	// Step 2: Determine the number of leading zeros needed
	int n_digits_index;
	if (index == 0)
	{
		n_digits_index = 1;
	}
	else
	{
		n_digits_index = log(1.0*index)/log(10.) + 1;
	}
	int leading_zeros = n_digits - n_digits_index;

	// Step 3: Create file number
	std::string filenumber = "";
	for (size_type i = 0; i < leading_zeros; ++i)
	{
		filenumber += "0";
	}
	filenumber += std::to_string(index);

	return filenumber;

}


void write_scripts(size_type n_scripts, std::string filename_begin)
{
	for (size_type i = 0; i < n_scripts; ++i)
	{
		// Create correct filename
		std::string filename = filename_begin;
		std::string filenumber = int_to_string(i, n_scripts);
		filename += filenumber;
		std::string format = ".sh";
		filename += format;
		std::ofstream outfile(filename);

		// Write header
		outfile << "#!/usr/bin/tcsh" << "\n";
		outfile << "#PBS -l walltime=15:00:00" << "\n";
		outfile << "#PBS -l mem=128mb" << "\n";

		// Set directory
		outfile << "cd /lpt/xrxjktfa/HMC-Fitter" << "\n";

		// Define name of executable file
		outfile << "./fitter" << "\n";

		// Move output file in a specific directory
		std::string dataname = "data";
		dataname += filenumber;
		outfile << "mv /lpt/xrxjktfa/HMC-Fitter/data.txt /lpt/xrxjktfa/HMC-Fitter/LPT-Cluster/" << dataname << ".txt" << "\n";


	}


}

#endif