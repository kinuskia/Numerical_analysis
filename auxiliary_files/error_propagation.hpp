#ifndef ERROR_PROPAGATION_HPP
#define ERROR_PROPAGATION_HPP

#include "vector.hpp"
#include <cmath>
#include <algorithm>
#include "storage.hpp"
#include <random>



template <class REAL>
class ErrorPropagation
{
public:
	typedef std::size_t size_type;
	typedef REAL number_type;
	// Constructor with x-uncertainties
	ErrorPropagation(
	)
	{
	}


	// Monte-Carlo analysis of error propagation by generating n experimental data sets (Gaussian statistics)
		template<class Model>
		void get_propagated_error(const Model & model, const Vector<number_type> & x, const Vector<number_type> & dx, size_type n, number_type & uncertainty, number_type & uncertainty_error, std::string filename)
		{
			assert(model.n_parameters() == x.size());
			Storage<number_type> data(x, 1); 

			for (size_type i = 0; i < n; ++i)
			{
				// Generate x-sample
				Vector<number_type> x_current(x.size());
				fill_gaussian(x_current, x, dx);

				// compute y-value
				number_type y_current = model.f(x_current);

				// save to storage file
				data.read_in(x_current);
				data.read_in(y_current);
			}

			// compute resulting uncertainty
			data.std_deviation(x.size(), uncertainty, uncertainty_error);

			// save in output file
			data.write(filename, false);
		}

		// Overload for no text file output
		template<class Model>
		void get_propagated_error(const Model & model, const Vector<number_type> & x, const Vector<number_type> & dx, size_type n, number_type & uncertainty, number_type & uncertainty_error)
		{
			assert(model.n_parameters() == x.size());
			Storage<number_type> data(x, 1); 

			for (size_type i = 0; i < n; ++i)
			{
				// Generate x-sample
				Vector<number_type> x_current(x.size());
				fill_gaussian(x_current, x, dx);

				// compute y-value
				number_type y_current = model.f(x_current);

				// save to storage file
				data.read_in(x_current);
				data.read_in(y_current);
			}

			// compute resulting uncertainty
			data.std_deviation(x.size(), uncertainty, uncertainty_error);

		}

			// Monte-Carlo analysis of error propagation by generating n experimental data sets (uniform distribution)
		template<class Model>
		void get_propagated_error_uniform(const Model & model, const Vector<number_type> & x, const Vector<number_type> & dx, size_type n, number_type & uncertainty, number_type & uncertainty_error, std::string filename)
		{
			assert(model.n_parameters() == x.size());
			Storage<number_type> data(x, 1); 


			for (size_type i = 0; i < n; ++i)
			{
				// Generate x-sample
				Vector<number_type> x_current(x.size());
				
				fill_from_region(x_current, x-dx, x+dx);

				// compute y-value
				number_type y_current = model.f(x_current);

				// save to storage file
				data.read_in(x_current);
				data.read_in(y_current);
			}

			// compute resulting uncertainty
			data.std_deviation(x.size(), uncertainty, uncertainty_error);

			// save in output file
			data.write(filename, false);
		}

private:




};

#endif