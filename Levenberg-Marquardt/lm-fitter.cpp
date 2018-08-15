#include <iostream>
#include "../auxiliary_files/read_data.hpp"
#include "../auxiliary_files/vector.hpp"
#include "../auxiliary_files/matrix.hpp"
#include "../auxiliary_files/matrix.hpp"
#include "model.hpp" // load before lm.hpp !!
#include "../auxiliary_files/lm.hpp"
#include <iomanip>



int main ()
{
	typedef double number_type;
	typedef std::size_t size_type;


	// Read in experimental data
	Vector<number_type> x(0);
	Vector<number_type> y(0);
	Vector<number_type> dy(0);

	read_data("measurements.txt", x, y, dy);

	// Set up fitting model
	Model<number_type> gaussian;


	// Set up Levenberg-Marquardt solver
	LM<number_type> minimizer(gaussian, x, y, dy);

	// lower and upper bounds for search region
	Vector<number_type> range_min(gaussian.n_parameters());
	Vector<number_type> range_max(gaussian.n_parameters());
	range_min[0] = 5;
	range_max[0] = 7;
	range_min[1] = 3;
	range_max[1] = 5;
	range_min[2] = 0.1;
	range_max[2] = 2;

	// Set up initial guess
	Vector<number_type> popt(gaussian.n_parameters());


	
	for (size_type i = 0; i < 1e3; ++i)
	{
		fill_from_region(popt, range_min, range_max);
		minimizer.solve(popt);
		std::cout << minimizer.chi2_red(popt) << "\n";
	}

	if (!minimizer.converged())
	{
		std::cout << "No convergence!\n";
	}
	for (size_type i = 0; i < popt.size(); ++i)
	{
		std::cout << i << " : " << popt[i] << "\n";
	}
	std::cout << "chi2/d.o.f = " << minimizer.chi2_red(popt) << "\n";


	fill_from_region(popt, range_min, range_max);
	std::cout << std::setprecision(14) << minimizer.first_derivative_unbiased(2, popt, 2, 1e-1) << "\n";
	std::cout << std::setprecision(14) << minimizer.first_derivative(2, popt, 2, 1e-1) << "\n";
	std::cout << std::setprecision(14) << minimizer.first_derivative_unbiased(2, popt, 2, 1e-3) << "\n";
	std::cout << std::setprecision(14) << minimizer.first_derivative(2, popt, 2, 1e-3) << "\n";



	return 0;
}