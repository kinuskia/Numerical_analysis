#include <iostream>
#include <cmath>
#include <random>
#include <fstream>
#include "../auxiliary_files/vector.hpp"

/* model fitting function */
template<typename number_type>
number_type f(Vector<number_type> x, Vector<number_type> popt)
{
	return popt[0] * exp(-(x[0]-popt[1])*(x[0]-popt[1])/2/popt[2]/popt[2]-(x[1]-popt[3])*(x[1]-popt[3])/2/popt[4]/popt[4]);
}


int main ()
{
	std::size_t n_data = 170;
	std::size_t x_dim = 2;
	Vector<double> popt(5);
	popt[0] = 1;
	popt[1] = 2;
	popt[2] = 0.5;
	popt[3] = 4;
	popt[4] = 1.7;

	// x-uncertainty?
	bool x_uncertainty = false;


	Vector<double> lower_x(2);
	Vector<double> upper_x(2);
	lower_x[0] = 0.5;
	upper_x[0] = 3.5;
	lower_x[1] = -1.1;
	upper_x[1] = 9.1;


	std::random_device rd;
	std::mt19937 gen(rd());
	
	std::vector<Vector<double>> x(n_data);
	std::vector<Vector<double>> dx(n_data);
	Vector<double> y(n_data);
	Vector<double> dy(n_data);

	double rel = 0.02;
	double abs = 0.02;

	double rel_x = 0.002;
	double abs_x = 0.002;

	// Create (theoretical) x-values
	for (std::size_t i = 0; i < n_data; ++i)
	{
		for (std::size_t j = 0; j < x_dim; ++j)
		{
			std::uniform_real_distribution<> dis_unif_x(lower_x[j], upper_x[j]);
			(x[i]).push_back(dis_unif_x(gen));
		}	
	}
	
	// Sort by x-values in case of 1D fit
	if (x[0].size() == 1)
	{
		std::sort(x.begin(), x.end(),
		[] (Vector<double> a, Vector<double> b)
		{
			return a[0] < b[0];
		}
		);
	}

	// calculate dx-values
	for (std::size_t i = 0; i < n_data; ++i)
	{
		dx[i] = rel_x * x[i] + abs_x;	
	}
	

	// Create y- and dy-values
	for (std::size_t i = 0; i<x.size(); ++i)
	{
		double ytheo = f(x[i], popt);
		double error = rel*ytheo + abs;
		std::normal_distribution<> dis_norm(0, error);
		y[i] = ytheo + dis_norm(gen);
		dy[i] = error;
	}

	// Randomly shift x-values
	std::vector<Vector<double>> x_theo(x);
	for (std::size_t i = 0; i < n_data; ++i)
	{
		for (std::size_t j = 0; j < x_dim; ++j)
		{
			std::normal_distribution<> dis_norm_x((x_theo[i])[j], (dx[i])[j]);
			(x[i])[j] = dis_norm_x(gen);
		}
	}

	// save output to text file
	std::ofstream outfile("measurements.txt");
	for (std::size_t i = 0; i < n_data; ++i)
	{
		for (std::size_t j = 0; j < x[0].size(); ++j)
		{
			outfile << (x[i])[j] << " ";
			if (x_uncertainty)
			{
				outfile << (dx[i])[j] << " ";
			}
		}
		outfile << y[i] << " " << dy[i] << "\n";
	}





	return 0;
}