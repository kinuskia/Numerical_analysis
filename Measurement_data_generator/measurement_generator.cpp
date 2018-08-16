#include <iostream>
#include <cmath>
#include <random>
#include <fstream>
#include "../auxiliary_files/vector.hpp"

double f(double x, Vector<double> popt)
{
	return popt[0] * exp(-popt[1]*x) + popt[2];
}


int main ()
{
	std::size_t n_data = 1e2;
	Vector<double> popt(3);
	popt[0] = 6;
	popt[1] = 0.6;
	popt[2] = 1;


	
	Vector<double> pstd(popt.size());
	double rel = 0.02;
	double abs = 0.02;
	pstd = rel*popt + abs;

	double x_min = 0;
	double x_max = 10;

	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_real_distribution<> dis_unif(x_min, x_max);

	Vector<double> x(n_data);
	Vector<double> y(n_data);
	Vector<double> dy(n_data);

	// Create x values
	for (std::size_t i = 0; i < n_data; ++i)
	{
		x[i] = dis_unif(gen);
	}
	std::sort(x.begin(), x.end());
	

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
		outfile << x[i] << "  " << y[i] << " " << dy[i] << "\n";
	}





	return 0;
}