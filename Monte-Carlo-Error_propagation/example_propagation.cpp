#include <iostream>


#include "../auxiliary_files/error_propagation.hpp"
#include "model.hpp"

int main ()
{
	typedef double number_type;
	typedef std::size_t size_type;

	Vector<number_type> input(2);
	Vector<number_type> input_err(2);
	input[0] = 1;
	input_err[0] = 0.5;
	input[1] = 5;
	input_err[1] = 1;


	input_err *= 1.0;

	Model<number_type> model;

	number_type y = model.f(input);
	number_type y_err_gauss = y * sqrt(pow(input_err[0]/input[0],2) + pow(input_err[1]/input[1],2) );

	ErrorPropagation<number_type> propagator;
	number_type y_err_mc;
	number_type y_err_err_mc;

	propagator.get_propagated_error(model, input, input_err, 1e4, y_err_mc, y_err_err_mc, "outfile.txt");

	std::cout << "Uncertainty using standard formula: " << y << " +/- " << y_err_gauss << "\n";
	std::cout << "Uncertainty using Monte Carlo     : " << y << " +/- " << y_err_mc  << "\n";
	std::cout << "Relative uncertainty of MC uncertainty: " << y_err_err_mc/y_err_mc << "\n";

	return 0;
}