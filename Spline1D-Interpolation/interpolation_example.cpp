#include <iostream>
#include <cmath>


#include "../auxiliary_files/vector.hpp"
#include "../auxiliary_files/spline1d.hpp"
#include "../auxiliary_files/to_file.hpp"

int main ()
{
	typedef std::size_t size_type;
	typedef double number_type;

	// Data points for interpolation
	size_type N = 14;
	Vector<number_type> x(N);
	Vector<number_type> y(N);
	number_type x_min = 0; 
	number_type x_max = 3.1415926*2;
	for (size_type i = 0; i < N; ++i)
	{
		x[i] = x_min+(x_max-x_min)*i/(N-1);
		y[i] = sin(x[i])+0.4*cos(2.*x[i]);
	}

	Spline1D<number_type> spline;

	spline.initialize(x, y, "natural");
	//spline.print_data();

	// Evaluate spline
	size_type n_points = 200;
	Vector<number_type> x_spline(n_points);
	Vector<number_type> y_spline(n_points);
	for (size_type i = 0; i < n_points; ++i)
	{
		x_spline[i] = x_min+(x_max-x_min)*i/n_points;
		y_spline[i] = spline.evaluate(x_spline[i]);
	}
	// for (size_type i = 0; i < N; ++i)
	// {
	// 	x_spline[i] = x[i];
	// 	y_spline[i] = spline.evaluate(x_spline[i]);
	// }


	std::vector<Vector<number_type>> input(2);
	input[0] = x;
	input[1] = y; 
	to_file("data_points.txt", input);

	std::vector<Vector<number_type>> output(2);
	output[0] = x_spline;
	output[1] = y_spline; 
	to_file("spline.txt", output);



	return 0;
}