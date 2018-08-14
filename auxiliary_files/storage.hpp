#ifndef STORAGE_HPP
#define STORAGE_HPP

#include "vector.hpp"
#include <string>
#include <fstream>
#include <iomanip>
#include <algorithm>

template<typename REAL>
class Storage
{
public:
	typedef std::size_t size_type;
	typedef REAL number_type;
	Storage (Vector<number_type> popt, size_type n_of_extra = 0)
	: n_popt_(popt.size())
	, n_extra_(n_of_extra) 
	, n_variables_(popt.size() + n_of_extra)
	, data_()
	, thermalization_(0)
	, entries_per_variable_(0)
	{
	}

	Storage (size_type n_popt, size_type n_of_extra = 0)
	: n_popt_(n_popt)
	, n_extra_(n_of_extra) 
	, n_variables_(n_popt + n_of_extra)
	, data_()
	, thermalization_(0)
	, entries_per_variable_(0)
	{
	}

	/* Clear storage but keep number of popt, extras etc */
	void clear ()
	{
		data_.clear();
		thermalization_ = 0;
		entries_per_variable_ = 0;
	}

	/* Getter methods */
	size_type n_variables() const
	{
		return n_variables_;
	}

	size_type entries_per_variable() const
	{
		return entries_per_variable_;
	}



	/* receive and save input value */
	void read_in (number_type value)
	{
		data_.push_back(value);
		entries_per_variable_ = (data_.size()/n_variables_ - thermalization_); // discarding non thermalized data
	}
	/* receive and save input vector */
	void read_in (Vector<number_type> popt)
	{
		for (size_type i = 0; i < popt.size(); ++i)
		{
			data_.push_back(popt[i]);
		}
		entries_per_variable_ = data_.size()/n_variables_;
	}

	/* Determine median of parameter i */
	number_type median(size_type index)
	{
		// Fill a vector with parameter values
		Vector<number_type> values(data_.size()/n_variables_);
		for (size_type i = 0; i < values.size(); ++i)
		{
			values[i] = data_[index + i*n_variables_];
		}
		// Place the median value in the middle of the vector
		std::nth_element(values.begin(), values.begin()+values.size()/2, values.end());
		return values[values.size()/2];
	}

	/* 
	Return a vector containing, for each parameter, the percentile of all values for which chi2 smaller than potential_max
	 e.g min_max_percentile(minimals, 20, 0.05) returns the value below which 5 % of the values with chi2red < 20
	 are found and above which 5 % are found.
	*/
	void min_max_percentile(Vector<number_type> & minimals, Vector<number_type> & maximals, number_type potential_max, number_type percentile ) const
	{
		assert(minimals.size() == n_popt_);
		assert(maximals.size() == n_popt_);

		// fill vector of vectors with relevant data
		std::vector<Vector<number_type>> values_cutoff(n_popt_);
		Vector<number_type> current_values(n_popt_);
		for (size_type i = 0; i < data_.size(); ++i)
		{
			if (i%n_variables_>n_popt_)
			{
				continue;
			}
			// Save values to a given chi2 in a vector and add it to cutoff vector if chi2 < potential_max of if no maximal limit has been set
			if (i%n_variables_ < n_popt_)
			{
				current_values[i%n_variables_] = data_[i];
				continue;
			}
			if (i%n_variables_ == n_popt_)
			{
				if ((data_[i] < potential_max))
				{
					for (size_type j = 0; j < values_cutoff.size(); ++j)
					{
						values_cutoff[j].push_back(current_values[j]);
					}
				}
				continue;
			}
		}

		size_type entries_per_variable = values_cutoff[0].size();

		// Find the lower and upper percentile for each parameter
		if (entries_per_variable == 0)
		{

		}
		for (size_type i = 0; i < values_cutoff.size(); ++i)
		{
			std::sort(values_cutoff[i].begin(), values_cutoff[i].end());
			size_type minimal_index = entries_per_variable * percentile; // type conversion wanted !
			size_type maximal_index = entries_per_variable -1 - minimal_index;
			if (entries_per_variable == 0)
			{
				maximal_index = 0;
			}
			// // Change parameter range if they constrain the parameter region
			// number_type minimal_current = minimals[i];
			// number_type maximal_current = maximals[i];
			// number_type minimal_new = values_cutoff[i][minimal_index];
			// number_type maximal_new = values_cutoff[i][maximal_index];
			// if (minimal_new > minimal_current && maximal_new < maximal_current)
			// {
			// 	minimals[i] = values_cutoff[i][minimal_index];
			// 	maximals[i] = values_cutoff[i][maximal_index];
			// }
			minimals[i] = values_cutoff[i][minimal_index];
			maximals[i] = values_cutoff[i][maximal_index];
		}

	}

	/* Determine burn_in length */
	void determine_burn_in_length()
	{
		number_type burn_in_length = 0;
		for (size_type i = 0; i < n_popt_; ++i)
		{
			number_type median = this->median(i);
			size_type counter = 0;
			// Count steps until median is passed
			bool sign = (data_[i] >= median);
			bool sign_old = sign;
			for (size_type j = 0; j < data_.size()/n_variables_; ++j)
			{
				sign = (data_[i + n_variables_*j] >= median);
				if (j > 0)
				{
					if (sign != sign_old)
					{
						burn_in_length = counter + 2;
						if (thermalization_ < burn_in_length)
						{
							thermalization_ = burn_in_length;
						}
						break; // break when median passed
					}
					counter++;
				}
				sign_old = sign;
			}
		}

		// Recalculate entries per variable
		entries_per_variable_ = (data_.size()/n_variables_ - thermalization_); // discarding non thermalized data

	}

	number_type get_burn_in_length() const
	{
		return thermalization_;
	}

	/* calculate average of variable i */
	number_type mean(size_type index)
	{
		assert((index >= 0) && (index < n_variables_)); //check if index within valid bounds
		number_type result = 0;
		for (size_type i = index + n_variables_*thermalization_; i < data_.size(); i+=n_variables_)
		{
			result += data_[i];
		}
		result /= entries_per_variable_;
		return result;

	}

	/* calculate average of variable i as well as its error */
	void mean(size_type index, number_type & result, number_type & result_err)
	{
		// mean
		assert((index >= 0) && (index < n_variables_)); //check if index within valid bounds
		result = 0;
		for (size_type i = index + n_variables_*thermalization_; i < data_.size(); i+=n_variables_)
		{
			result += data_[i];
		}
		result /= entries_per_variable_;

		// and its error
		number_type autocorrel_time;
		number_type autocorrel_time_err;
		autocorr_time(index, autocorrel_time, autocorrel_time_err);
		number_type ESS = entries_per_variable_/(2.*autocorrel_time); //effective sample size
		number_type variance = this->variance(index);
		result_err = sqrt(variance/ESS);
	}

	/* 
	calculate 
	- mean of each variable at once
	- its fitting uncertainty
	- the statistical uncertainty of the mean
	- the statistical uncertainty of the fitting uncertainty
	 */
	void mean(Vector<number_type> & average_vector, Vector<number_type> & err_vector, Vector<number_type> & average_err, Vector<number_type> & err_err, size_type d_of_freedom, number_type beta)
	{
		assert(average_vector.size() == err_vector.size());
		assert(average_vector.size() == err_err.size());
		assert(average_vector.size() == n_popt_);

		// Determine burn-in length
		determine_burn_in_length();


		// Calculation of the mean vector
		average_vector = 0;
		for (size_type i = n_variables_*thermalization_; i < data_.size(); ++i)
		{
			if ((i%n_variables_) >= n_popt_)
			{
				continue;
			}
			average_vector[i%n_variables_] += data_[i];
		}
		average_vector = average_vector / entries_per_variable_;

		// Calculation of the fitting error
		err_vector = 0;
		for (size_type i = n_variables_*thermalization_; i < data_.size(); ++i)
		{
			if ((i%n_variables_) >= n_popt_)
			{
				continue;
			}
			err_vector[i%n_variables_] += (data_[i]-average_vector[i%n_variables_])*(data_[i]-average_vector[i%n_variables_]); // calculate unbiased variance
		}
		err_vector = err_vector / (entries_per_variable_ - 1);
		sqrt(err_vector); // standard deviation
		Vector<number_type> standard_deviation(err_vector.size());
		standard_deviation = err_vector;
		err_vector *= sqrt(2./d_of_freedom*beta); // fitting vector


		// Calculation of the statistical uncertainties
		Vector<number_type> autocorr_times(average_vector.size());
		Vector<number_type> autocorr_times_err(err_err.size());
		autocorr_time(autocorr_times, autocorr_times_err); // calculate integrated autocorrelation times
		std::cout << "Integrated autocorrelation times: \n";
		for (size_type i = 0; i < autocorr_times.size(); ++i)
		{
			std::cout << "Parameter " << i << " : " << autocorr_times[i] << " + - " << autocorr_times_err[i] << "\n";
		}
	

		for (size_type i = 0; i<n_popt_; ++i)
		{
			number_type ESS = entries_per_variable_/(2.*autocorr_times[i]); //effective sample size
			average_err[i] = standard_deviation[i]/sqrt(ESS);
			err_err[i] = err_vector[i]*sqrt(2./(4.*ESS-1));
		}
	}

	/* calculate mean of each variable at once and save it into a vector */
	/* as well as uncertitude vector, taking into account autocorrelation */
	void mean(Vector<number_type> & average_vector, Vector<number_type> & err_vector)
	{
		assert(average_vector.size() == err_vector.size());
		assert(average_vector.size() == n_variables_);

		// Determine burn-in length
		determine_burn_in_length();


		// Calculation of the mean vector
		average_vector = 0;
		for (size_type i = n_variables_*thermalization_; i < data_.size(); ++i)
		{
			average_vector[i%n_variables_] += data_[i];
		}
		average_vector = average_vector / entries_per_variable_;


		// Calculation of its uncertainty
		Vector<number_type> autocorr_times(average_vector.size());
		Vector<number_type> autocorr_times_err(err_vector.size());
		autocorr_time(autocorr_times, autocorr_times_err); // calculate integrated autocorrelation times
		std::cout << "Integrated autocorrelation times: \n";
		for (size_type i = 0; i < autocorr_times.size(); ++i)
		{
			std::cout << "Parameter " << i << " : " << autocorr_times[i] << " + - " << autocorr_times_err[i] << "\n";
		}
	

		for (size_type i = 0; i<n_variables_; ++i)
		{
			number_type ESS = entries_per_variable_/(2.*autocorr_times[i]); //effective sample size
			number_type variance = this->variance(i);
			err_vector[i] = sqrt(variance/ESS);
		}
	}





	/* 
	Autocorrelation function \Gamma_{\alpha, \beta}(n)
	for correlation between parameters with indices alpha and beta with time lag n 
	*/
	number_type gamma(size_type alpha, size_type beta, size_type lag)
	{
		number_type mean_alpha = this->mean(alpha);
		number_type mean_beta = this->mean(beta);
		assert(entries_per_variable_ > lag); // condition for validity of this estimator
		number_type result = 0;
		for (size_type i = thermalization_; i < (entries_per_variable_ - lag); ++i)
		{
			result += (data_[alpha + (i)*n_variables_] - mean_alpha) * (data_[beta + (i+lag)*n_variables_] - mean_beta);
		}
		result /= (entries_per_variable_ - lag);
		return result;
	}

	/* Autocorrelation function for chi2_red as derived quantity */
	number_type gamma_chi2red(Vector<number_type> & derivatives, size_type lag)
	{
		assert(derivatives.size() == n_popt_);
		number_type result = 0;
		for (size_type alpha = 0; alpha < n_popt_; ++alpha)
		{
			for (size_type beta = 0; beta < n_popt_; ++beta)
			{
				result += derivatives[alpha] * derivatives[beta] * gamma(alpha, beta, lag);
			}
		}
		
		return result;
	}


	/* calculate a integrated autocorrelation time vector (one entry for each fitting parameter) */
	void autocorr_time(Vector<number_type> & times, Vector<number_type> & times_err)
	{
		assert(times.size() == times_err.size()); // check for correct dimensions
		assert(times_err.size() > 0);

		// // Print cross-correlations
		// for (size_type i = 0; i < times.size(); ++i)
		// {
		// 	for (size_type j = 0; j < times.size(); ++j)
		// 	{
		// 		std::cout << this->gamma(i, j, 0) << " ";
		// 	}
		// 	std::cout << "\n";
		// }

		for (size_type i = 0; i < times.size(); ++i)
		{
			autocorr_time(i, times[i], times_err[i]);	
		}
	}

	/* autocorrelation time for index i */
	void autocorr_time(size_type i, number_type & time, number_type & time_error)
	{
		//std::cout << "Calculating autocorrelation time for parameter " << i+1  << " / " << times.size() << "\n";
		// Calculate integrated autocorrelation time for parameter i
		number_type C = abs(this->gamma(i, i, 0));
		size_type W = 0;
		number_type gamma_old;
		number_type gamma_current;
		for (size_type t = 1; t < entries_per_variable_; ++t)
		{
			gamma_current = abs(this->gamma(i, i, t));
			//if (i == 0)
				//std::cout << this->gamma(i, i, t-1) << "\n";
			C+= 2.*gamma_current;

			/* 
			procedure to find optimal summation window W: 
			break as soon as module of autocorrelation function increases 
			*/
			if (t>1)
			{
				if (gamma_old < gamma_current)
				{
					W = t;
					break;
				}
			}

			gamma_old = abs(gamma_current);
		}
			number_type v = abs(this->gamma(i, i, 0));
			time = C/2.0/v;
			time_error = time*sqrt(4./entries_per_variable_*(W+1./2-time)) + C*exp(-1.*W/time);

	}

	/* calculate (unbiased) variance of variable i */
	number_type variance(size_type index)
	{
		number_type mean = this->mean(index);
		number_type result = 0;
		for (size_type i = index + n_variables_*thermalization_; i < data_.size(); i+=n_variables_)
		{
			result += (data_[i] - mean)*(data_[i] - mean);
		}
		result /= (entries_per_variable_ -1);
		return result;
	}

	/* Calculate variance as well as its statistical uncertainty */
	void variance(size_type index, number_type & result, number_type & result_err)
	{
		// variance
		result = variance(index);
		//autocorrelation time
		number_type autocorrel_time;
		number_type autocorrel_time_err;
		autocorr_time(index, autocorrel_time, autocorrel_time_err);
		// uncertainty
		number_type ESS = entries_per_variable_/(2.*autocorrel_time);

		result_err = result*sqrt(2./(ESS-1));

	}

	/* Calculate standard deviation of variable i */
	number_type std_deviation(size_type index)
	{
		number_type result = sqrt(variance(index));
		return result;
	}


	/* Determine a vector containing the period of each signal */
	void period(Vector<size_type> & period_vector)
	{
		assert(period_vector.size() == n_popt_);

		for (size_type i = 0; i < n_popt_; ++i)
		{
			/* 
			To calculate the period of the signal, we determine the
			number of steps for two sign changes
			*/
			size_type counter_period = 0;
			size_type sign_changes = 0;
			bool sign = true; // true = '+', false = '-'
			bool sign_old = sign;
			for (size_type j = 0; j < entries_per_variable_; ++j)
			{
				number_type gamma = this->gamma(i, i, j);
				sign = (gamma >= 0);
				if (j > 0)
				{
					if (sign_changes == 2)
					{
						break; // stop calculating once period found
					}
					else if (sign != sign_old)
					{
						sign_changes++;
					}
					counter_period++;
				}
				sign_old = sign;
			}
			period_vector[i] = counter_period;
		}

	}

	/* Write data in an output file */
	void write(std::string filename, bool add_counter = true)
	{
		std::ofstream outfile(filename);
		size_type N = data_.size();
		for (size_type i = 0; i < N; ++i)
		{
			if ((i%n_variables_) == 0 && add_counter) // add counter at beginning of line
			{
				outfile << i/n_variables_ + 1;
				outfile << " ";
			}

			outfile << std::setprecision(14) << data_[i];
			
			if((i%n_variables_) == (n_variables_ - 1)) // new line
			{
				outfile << "\n";
				continue;
			}
			outfile << " ";
		}

	}

	/* Overload: Write data in an output file while discarding data with chi2red larger than a limit */
	void write(std::string filename, number_type maxvalue)
	{
		std::ofstream outfile(filename);
		size_type N = data_.size();
		for (size_type i = 0; i < N; ++i)
		{
			if ((i%n_variables_) == 0) // add counter at beginning of line
			{
				// check if chi2 is below maxvalue
				if (data_[i+n_popt_] > maxvalue)
				{
					i += n_variables_-1;
					continue;
				}

				outfile << i/n_variables_ + 1;
				outfile << " ";
			}

			outfile << data_[i];
			
			if((i%n_variables_) == (n_variables_ - 1)) // new line
			{
				outfile << "\n";
				continue;
			}
			outfile << " ";
		}

	}


private:
	Vector<number_type> data_; // data vector
	size_type n_variables_; // number of variables stored in the vector
	size_type n_popt_; // number of fitting parameters
	size_type n_extra_; // number of additional variables
	size_type thermalization_; // number of steps to reach thermalization (values up to this index are discarded)
	size_type entries_per_variable_; // to be calculated once data has been read in
};


#endif 