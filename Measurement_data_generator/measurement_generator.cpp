#include <iostream>
#include <cmath>
#include <random>
#include <fstream>
#include "../auxiliary_files/vector.hpp"

/* model fitting function */
template<typename number_type>
number_type f(Vector<number_type> x, Vector<number_type> popt)
{
	return popt[0] + sqrt((x[0]-popt[1])*(x[0]-popt[1]) + (x[1]-popt[2])*(x[1]-popt[2]));
}


int main ()
{
	std::size_t n_data = 50;
	Vector<double> popt(3);
	popt[0] = 0;
	popt[1] = 0;
	popt[2] = 0;


	Vector<double> lower_x(2, -150);
	Vector<double> upper_x(2, 150);


	std::random_device rd;
	std::mt19937 gen(rd());
	
	std::vector<Vector<double>> x(n_data);
	std::vector<Vector<double>> dx(n_data);
	Vector<double> y(n_data);
	Vector<double> dy(n_data);

	double rel = 0.02*0;
	double abs = 0.02*0+0.1;

	double rel_x = 0.02*0;
	double abs_x = 0.02*0+0.1;

	// Create x and dy values
	for (std::size_t i = 0; i < n_data; ++i)
	{
		std::uniform_real_distribution<> dis_unif_x(lower_x[0], upper_x[0]);
		std::uniform_real_distribution<> dis_unif_y(lower_x[1], upper_x[1]);
		(x[i]).push_back(dis_unif_x(gen));
		(x[i]).push_back(dis_unif_y(gen));
		dx[i] = rel_x * x[i] + abs_x;
	}
	//std::sort(x.begin(), x.end());
	

	// Create y and dy values
	for (std::size_t i = 0; i<x.size(); ++i)
	{
		double ytheo = f(x[i], popt);
		double error = rel*ytheo + abs;
		std::normal_distribution<> dis_norm(0, error);
		y[i] = ytheo + dis_norm(gen);
		dy[i] = error;
	}

	// save output to text file
	std::ofstream outfile("measurements.txt");
	for (std::size_t i = 0; i < n_data; ++i)
	{
		for (std::size_t j = 0; j < x[0].size(); ++j)
		{
			outfile << (x[i])[j] << " ";
			outfile << (dx[i])[j] << " ";
		}
		outfile << y[i] << " " << dy[i] << "\n";
	}





	return 0;
}