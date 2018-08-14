#include <iostream>
#include "../auxiliary_files/read_data.hpp"
#include "../auxiliary_files/vector.hpp"
#include "../auxiliary_files/matrix.hpp"
#include "../auxiliary_files/matrix.hpp"
#include "model.hpp" // load before lm.hpp !!
#include "../auxiliary_files/lm.hpp"



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

	// Set up initial guess
	Vector<number_type> popt(gaussian.n_parameters());
	popt[0] = 6;
	popt[1] = 4;
	popt[2] = 0.78;
	minimizer.set_J(popt);

	std::cout << minimizer.J(2, 1) << "\n";





	return 0;
}