#ifndef TO_FILE_HPP
#define TO_FILE_HPP
#include <vector>
#include <fstream>
#include <iomanip>

#include "vector.hpp"
#include "matrix.hpp"


// Function which takes an GSL vector of GSL vector columns and writes them to a text file
template<typename number_type>
void to_file(std::string filename, std::vector<std::vector<number_type>> data_columns)
{
	typedef std::size_t size_type;
	std::ofstream outfile(filename);
	for (size_type i = 0; i < data_columns[0].size(); ++i)
	{
		for (size_type j = 0; j < data_columns.size(); ++j)
		{
			number_type value = (data_columns[j])[i] ;
			if (std::isnan(value))
			{
				outfile << 0;
			}
			else 
			{
				outfile << std::setprecision(14) << value;
			}
			if (j < data_columns.size()-1)
			{
				outfile << " ";
			}
		}
		outfile << "\n";
	}	
}

// overload for my own vectors
template<typename number_type>
void to_file(std::string filename, std::vector<Vector<number_type>> data_columns)
{
	typedef std::size_t size_type;
	std::ofstream outfile(filename);
	for (size_type i = 0; i < data_columns[0].size(); ++i)
	{
		for (size_type j = 0; j < data_columns.size(); ++j)
		{
			number_type value = (data_columns[j])[i] ;
			if (isnan(value))
			{
				outfile << 0;
			}
			else 
			{
				outfile << std::setprecision(14) << value;
			}
			if (j < data_columns.size()-1)
			{
				outfile << " ";
			}
		}
		outfile << "\n";
	}	
}

// overload for my own matrices
template<typename number_type>
void to_file(std::string filename, Matrix<number_type> data_columns)
{
	typedef std::size_t size_type;
	std::ofstream outfile(filename);
	for (size_type i = 0; i < data_columns.rowsize(); ++i)
	{
		for (size_type j = 0; j < data_columns.colsize(); ++j)
		{
			number_type value = data_columns(i, j);
			if (isnan(value))
			{
				outfile << 0;
			}
			else 
			{
				outfile << std::setprecision(14) << value;
			}
			if (j < data_columns.colsize()-1)
			{
				outfile << " ";
			}
		}
		outfile << "\n";
	}	
}


#endif