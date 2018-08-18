#include <iostream>
#include "../auxiliary_files/read_data.hpp"
#include "../auxiliary_files/vector.hpp"
#include "../auxiliary_files/matrix.hpp"
#include "model.hpp" // load before lm.hpp !!
#include "../auxiliary_files/lm.hpp"
#include <iomanip>



int main ()
{
	typedef double number_type;
	typedef std::size_t size_type;


	// Read in experimental data
	std::vector<Vector<number_type>> data(3); // Read in three columns
	read_data("measurements.txt", data);
	Vector<number_type> x(data[0]);
	Vector<number_type> y(data[1]);
	Vector<number_type> dy(data[2]);


	// Set up fitting model
	Model<number_type> gaussian;


	// Set up Levenberg-Marquardt solver
	std::vector<Vector<number_type>> X(1);
	X[0] = x;
	LM<number_type> minimizer(gaussian, X, y, dy);


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

	popt[0] = 5.56108;
	popt[1] = 5.84062;
	popt[2] = 0.0538888;
	fill_from_region(popt, range_min, range_max);
	minimizer.solve(popt);
	if (!minimizer.converged())
	{
		std::cout << "No convergence!\n";
	}
	for (size_type i = 0; i < popt.size(); ++i)
	{
		std::cout << i << " : " << popt[i] << "\n";
	}
	std::cout << "chi2/d.o.f = " << minimizer.chi2_red(popt) << "\n";




	
	// for (size_type i = 0; i < 1e3; ++i)
	// {
	// 	fill_from_region(popt, range_min, range_max);
	// 	Vector<number_type> popt_copy = popt;
	// 	minimizer.solve(popt);
	// 	if (!minimizer.converged())
	// 	{
	// 		std::cout << minimizer.chi2_red(popt) << "\n";
	// 		for (size_type i = 0; i < popt.size(); ++i)
	// 		{
	// 			std::cout << i << " : " << popt_copy[i] << "\n";
	// 		}
	// 	}
	// }

	// if (!minimizer.converged())
	// {
	// 	std::cout << "No convergence!\n";
	// }
	// for (size_type i = 0; i < popt.size(); ++i)
	// {
	// 	std::cout << i << " : " << popt[i] << "\n";
	// }
	// std::cout << "chi2/d.o.f = " << minimizer.chi2_red(popt) << "\n";


	



	return 0;
}