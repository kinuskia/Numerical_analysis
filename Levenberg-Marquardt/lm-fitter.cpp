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
	std::vector<Vector<number_type>> data(4); // Read in columns
	read_data("measurements.txt", data);
	Vector<number_type> x(data[0]);
	Vector<number_type> dx(data[1]);
	Vector<number_type> y(data[2]);
	Vector<number_type> dy(data[3]);


	// Set up fitting model
	Model<number_type> fitting_model;


	// Set up Levenberg-Marquardt solver
	std::vector<Vector<number_type>> X(1);
	X[0] = x;
	std::vector<Vector<number_type>> dX(1);
	dX[0] = dx;
	LM<number_type> minimizer(fitting_model, X, dX, y, dy);


	// lower and upper bounds for search region
	Vector<number_type> range_min(fitting_model.n_parameters());
	Vector<number_type> range_max(fitting_model.n_parameters());
	range_min[0] = 5;
	range_max[0] = 7;
	range_min[1] = 0.8;
	range_max[1] = 1.4;
	range_min[2] = 4;
	range_max[2] = 7;
	range_min[3] = 0.4;
	range_max[3] = 0.6;
	range_min[4] = 0.9;
	range_max[4] = 1.1;

	number_type success_ratio;
	Vector<number_type> popt(fitting_model.n_parameters());

	// Calculate fitting result
	success_ratio = minimizer.find_best_fit(popt, range_min, range_max, 1e2);
	std::cout << "Success ratio: " << success_ratio << "\n";

	Vector<number_type> uncertainty(popt.size());
	Vector<number_type> uncertainty_err(popt.size());
	minimizer.get_fit_uncertainty("fit_samples.txt", 5e3, uncertainty, uncertainty_err);
	std::cout << "Fitting result: \n";
	for (size_type i = 0; i < popt.size(); ++i)
	{
		std::cout << popt[i] << " +/- " << uncertainty[i]  << "\n";
	}
	uncertainty_err = uncertainty_err / uncertainty;
	std::cout << "Maximal relative uncertainty error: " << *std::max_element(uncertainty_err.begin(), uncertainty_err.end()) << "\n";

	std::cout << "chi2/d.o.f = " << minimizer.chi2_red(popt) << "\n";







	



	return 0;
}