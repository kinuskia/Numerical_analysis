#ifndef HMC_HPP
#define HMC_HPP

#include "vector.hpp"
#include "matrix.hpp"
#include <assert.h>
#include <random>
#include <ctime>
#include "storage.hpp"
//#include "../model.hpp"
#include <algorithm>
#include <fstream>
#include <string>
#include <stdexcept>

#include "to_file.hpp"

template<typename REAL>
class HMC
{
public:
	typedef std::size_t size_type;
	typedef REAL number_type;
	/* Standard constructor, fed with leapfrog data */
	HMC(LatticeAction<number_type> action, number_type stepsize, size_type n_steps_min, size_type n_steps_max)
	: stepsize_(stepsize)
	, n_steps_min_(n_steps_min)
	, n_steps_max_(n_steps_max)
	, counter_(0)
	, counter_accepted_(0)
	, action_(action)
	{
	}


	/* getter for action */
	number_type action(const Vector<number_type> & q)
	{
		return action_.action(q);

	}

	// getter for number of lattice sites per space-time dimension
	size_type n_sites_per_dim()
	{
		return action_.n_sites_per_dim();
	}

	/* Getter for operator */
	number_type Operator(const Vector<number_type> & q, size_type t)
	{
		return action_.Operator(q, t);
	}



	// /* check if given position is within bounds */
	// bool is_in_range(Vector<number_type> position) const
	// {
	// 	bool inbounds = true;
	// 	for (size_type i = 0; i < position.size(); ++i)
	// 	{
	// 		inbounds = (position[i] < range_max_[i]) && (position[i] > range_min_[i]);
	// 		if (!inbounds)
	// 		{
	// 			break;
	// 		}
	// 	}
	// 	return inbounds;
	// }

	/* gradient of the action at a position in parameter space */
	void grad_action(Vector<number_type> & output, const Vector<number_type> & position)
	{
		if (output.size() != position.size())
			throw std::runtime_error("Wrong dimensions.\n");
		// second order formula
		for (size_type i = 0; i < output.size(); ++i)
		{
			// estimation of partial derivative with respect to q_i
			number_type h = stepsize_;
			Vector<number_type> position_forward = position;
			position_forward[i] += h;
			Vector<number_type> position_backward = position;
			position_backward[i] -= h;
			output[i] = (action(position_forward)-action(position_backward))/2.0/h;
		}
	}

	

	/* kinetic energy function needed for leapfrog dynamics */
	number_type kinetic (const Vector<number_type> & p)
	{
		number_type energy = 0;
		for (size_type i = 0; i < p.size(); ++i)
		{
			energy += p[i]*p[i]/2.0;
		}
		return energy;
	}

	/* Leapfrog integrator */
	void leapfrog(Vector<number_type> & q, Vector<number_type> & p, size_type n_steps)
	{
		// Introducing a step size vector: 
		Vector<number_type> stepsizes(q.size(), stepsize_);
		// Make a half step for momentum at the beginning
		Vector<number_type> grad_S(q.size());
		grad_action(grad_S, q);
		p -= stepsizes * grad_S / 2;
		//std::cout << "p[0] = " << p[0]<<"\n";

		// Alternate full steps for position and momentum
		for (size_type i = 0; i < n_steps; ++i)
		{
			// Make a full step for the position, update gradient of action
			q += stepsizes * p;
			
			grad_action(grad_S, q);
			// Make a full step for the momentum, except at the end of the trajectory
			if (i != n_steps - 1)
			{
				p -= stepsizes * grad_S;
			}
		}
		// make a half step for momentum at the end
		p -= stepsizes * grad_S / 2;
	}




	/* 
		Leapfrog integrator copying position to external vector after each step
		and offering adjustable number of iterations
	 */

	void leapfrog(Vector<number_type> & q, Vector<number_type> & p, Storage<number_type> & q_copy, size_type nb_iterations)
	{
		// Introducing a step size vector: 
		Vector<number_type> stepsizes(q.size(), stepsize_);

		// Make a half step for momentum at the beginning
		Vector<number_type> grad_S(q.size());
		grad_action(grad_S, q);
		p -= stepsizes * grad_S / 2;

		// Alternate full steps for position and momentum
		for (size_type i = 0; i < nb_iterations; ++i)
		{
			// Make a full step for the position, update gradient of action
			q += stepsizes * p;

			// Read values of q into the storage vector
			q_copy.read_in(q);
			

			grad_action(grad_S, q);
			// Make a full step for the momentum, except at the end of the trajectory
			if (i != nb_iterations - 1)
			{
				p -= stepsizes * grad_S;
			}
		}
		// make a half step for momentum at the end
		p -= stepsizes * grad_S / 2;
	}

	



	/* do one Metropolis step */
	void step_forward(Vector<number_type> & current_q)
	{
		counter_++;
		Vector<number_type> q = current_q;

		// generate momenta from gaussian distribution
		std::normal_distribution<> dis_norm(0, 1); // mean 0, std deviation 1
		Vector<number_type> p(q.size());
		fill_random(p, dis_norm);
		Vector<number_type> current_p = p;

		// Draw random number of leapfrog steps
		std::random_device rd;
		std::mt19937 gen(rd());
		std::uniform_int_distribution<> dis_unif_int(n_steps_min_, n_steps_max_);
		size_type n_steps = dis_unif_int(gen);
		// Compute trajectory using the Leapfrog method
		leapfrog(q, p, n_steps);


		/* Negate momentum at the end of the trajectory to make the proposal symmetric
		(doesn't change the outcome of the algorithm, but is mathematically nicer)
		*/
		p = -p;

		// Evaluate potential (U) and kinetic (K) energies at start and end of trajectory
		number_type current_U = action(current_q);
		number_type current_K = kinetic(current_p);
		number_type proposed_U = action(q);
		number_type proposed_K = kinetic(p);


		/*
		Accept or reject the state at end of trajectory, returning either
		the position at the end of the trajectory or the initial position
		*/
		number_type H_change = -current_U-current_K+proposed_U+proposed_K;
		std::uniform_real_distribution<> dis_unif(0, 1);
		bool accepted = dis_unif(gen) < exp(-H_change);
		
		
		//accepted = accepted && model_.constraints_respected(q); // automatical rejection if contraints not respected
		if (accepted)
		{
			current_q = q; // accept
			counter_accepted_++;
		}
		// otherwise q is refused.

	}


	

	/* do n steps */
	void walk (size_type nb_steps, number_type max_duration, Vector<number_type> & initial, size_type progress_steps, std::string filename)
	{
		// reinitialize internal counters
		counter_ = 0;
		counter_accepted_ = 0;


		// Set timer
		std::time_t start = std::time(nullptr);
		number_type expected_duration;
		bool estimate_given = false;

		// Set up storage vector
		Storage<number_type> data(initial, n_sites_per_dim()+ 2); // records values of "initial" plus additional things
		
		// Start walking
		size_type counter = 0;
		while (counter < nb_steps)
		{
			step_forward(initial);

			// save updated data to storage
			data.read_in(initial); // save current lattice sites
		
			for (size_type i = 0; i < n_sites_per_dim(); ++i)
			{
				data.read_in(Operator(initial, i));
			}
			data.read_in(action(initial)); // save action
			data.read_in(acceptance_rate()); // save acceptance rate

			// After 10s: Estimate duration of the walk
			std::time_t end = std::time(nullptr);
			number_type diff = end - start;
			if ( (diff > 10) && !estimate_given && (counter > 0))
			{
				expected_duration = nb_steps * diff / counter / 60;
				std::cout << "Computation time in min: " << expected_duration << "\n";
				estimate_given = true;
			}

			if (diff > max_duration) // Cancel calculation if time exceeds max_duration
			{
				std::cout << "Computation aborted as time ran out. " << "\n";
				break;
			}

			// Print progress
			for (size_type i = 1; i < progress_steps; ++i)
			{
				if (counter == size_type(i*nb_steps/progress_steps))
				{
					std::cout << "Progress: " << size_type(i*100./progress_steps) << "%" << "\n";
				}
			}

			counter++;
		}

		// write data to output file
		data.write("output/lattice_dynamics.txt");

		// print autocorrelation times
		Vector<Vector<number_type>> autocorrel_times(2); 
		Vector<number_type> times(action_.n_sites());
		Vector<number_type> times_err(action_.n_sites());
		data.autocorr_time(times, times_err);
		autocorrel_times[0] = times;
		autocorrel_times[1] = times_err;
		to_file("output/autocorrel_times.txt", autocorrel_times);


		// print operator result
		Vector<Vector<number_type>> correlator(3);
		Vector<number_type> correlator_time(n_sites_per_dim());
		Vector<number_type> correlator_mean(n_sites_per_dim());
		Vector<number_type> correlator_error(n_sites_per_dim());
		for (size_type i = 0; i < correlator_mean.size(); ++i)
		{
			correlator_time[i] = i;
			data.mean(n_sites_per_dim()+i, correlator_mean[i], correlator_error[i]);
		}
		correlator[0] = correlator_time;
		correlator[1] = correlator_mean;
		correlator[2] = correlator_error;
		to_file(filename, correlator);
		
		

		// report total calculation time
		std::time_t end = std::time(nullptr);
		number_type diff = (end - start)/60;
		std::cout << "Total calculation time in min : " << diff << "\n";

		
	}

	// compute autocorrelation of the lattice sites for different number of leapfrog steps
	void autocorrelation(const Vector<number_type> & range_min, const Vector<number_type> & range_max, size_type Lmin, size_type Lmax, size_type Lstep, std::string filename)
	{
		Storage<number_type> autocorrel_times(action_.n_sites(),1);
		for (size_type L = Lmin; L <= Lmax; L += Lstep) // Loop over different numbers of Leapfrog steps
		{
			std::cout << "L = " << L << "\n";
			n_steps_min_ = L;
			n_steps_max_ = L;


			// Find a starting configuration
			Vector<number_type> sites(action_.n_sites());
			find_start(sites, range_min, range_max);


			Storage<number_type> data(action_.n_sites(), 0);

			size_type N = 1e3;

			for (size_type n = 0; n < N; ++n)
			{
				step_forward(sites);

				// save updated data to storage
				data.read_in(sites); // save current lattice sites
			}


			data.determine_burn_in_length();


			// get autocorrelation times and save to storage
			Vector<number_type> times(sites.size());
			Vector<number_type> times_err(sites.size());
			data.autocorr_time(times, times_err);

			autocorrel_times.read_in(times);
			autocorrel_times.read_in(L);

		}

		autocorrel_times.write(filename);
	}

	/* Find a starting point that respects contraints */
	void find_start(Vector<number_type> & initial, const Vector<number_type> & range_min, const Vector<number_type> & range_max)
	{
		if (initial.size() != action_.n_sites())
		{
			throw std::runtime_error("Wrong dimensions\n");
		}
		
		fill_from_region(initial, range_min, range_max);
	}

	// /* Create a Markov chain without analysis and without output to console */
	// void walk_silently(size_type nb_steps, std::string filename, std::string filenumber = "")
	// {
	// 	// Set up storage vector
	// 	Storage<number_type> data(model_.n_parameters(), 2); // records values of fitting vector and two additional things
		
	// 	// reinitialize internal counters
	// 	counter_ = 0;
	// 	counter_accepted_ = 0;
					
	// 	// choose random starting point that respects constraints
	// 	Vector<number_type> initial(action_.n_parameters());
	// 	find_start(initial);
			
			
	// 	// Start walking
	// 	size_type counter = 0;
	// 	while (counter < nb_steps)
	// 	{
				
	// 		step_forward(initial);

	// 		// save updated data to storage
	// 		data.read_in(initial); // save lattice sites
	// 		data.read_in(action(initial)); // save potential
	// 		data.read_in(acceptance_rate()); // save acceptance rate
						
	// 		counter++;

	// 		// start over if a bad initial condition was chosen
	// 		if (counter == 10 && acceptance_rate() < 0.01) 
	// 		{
	// 			find_start(initial);
	// 			counter = 0;
	// 			counter_ = 0;
	// 			counter_accepted_ = 0;
	// 			data.clear();
	// 		}				
	// 	}

	// 	// Save to data file without counter in the first column
	// 	filename += filenumber;
	// 	filename += ".txt";
	// 	data.write(filename, false);

				
	// }

	// /* Create a Markov chain without analysis and without output to console disregarding data with chi2 over a threshold*/
	// void walk_silently_disregarding(size_type nb_steps, number_type threshold_chi2, std::string filename, std::string filenumber = "")
	// {
	// 	// Set up storage vector
	// 	Storage<number_type> data(model_.n_parameters(), 2); // records values of fitting vector and two additional things
		
	// 	// reinitialize internal counters
	// 	counter_ = 0;
	// 	counter_accepted_ = 0;
					
	// 	// choose random starting point that respects constraints
	// 	Vector<number_type> initial(model_.n_parameters());
	// 	find_start(initial);
			
			
	// 	// Start walking
	// 	size_type counter = 0;
	// 	while (counter < nb_steps)
	// 	{
				
	// 		step_forward(initial);

	// 		// save updated data to storage
	// 		data.read_in(initial); // save fitting parameters
	// 		data.read_in(action(initial)/beta_); // save chi2_red
	// 		data.read_in(acceptance_rate()); // save acceptance rate
						
	// 		counter++;

	// 		// start over if a bad initial condition was chosen
	// 		if (counter == 10 && acceptance_rate() < 0.01) 
	// 		{
	// 			find_start(initial);
	// 			counter = 0;
	// 			counter_ = 0;
	// 			counter_accepted_ = 0;
	// 			data.clear();
	// 			continue;
	// 		}

	// 		// start over if chi2 is not below threshold after 200 steps
			
	// 		if (counter == 100 && action(initial)/beta_ > threshold_chi2)
	// 		{
	// 			find_start(initial);
	// 			counter = 0;
	// 			counter_ = 0;
	// 			counter_accepted_ = 0;
	// 			data.clear();
	// 			continue;
	// 		}
							
	// 	}

	// 	// Save to data file without counter in the first column
	// 	filename += filenumber;
	// 	filename += ".txt";
	// 	data.write(filename, false);

				
	// }

	// /* partially automatic walk to find to global minimum region */
	// void tune_parameters()
	// {
	// 	std::string user_input = "";
	// 	std::string quit = "quit";
	// 	std::string next = "next";
	// 	std::string back = "back";
	// 	std::string next_chain = "continue";
	// 	std::string repeat = "repeat";
	// 	std::string yes = "yes";
	// 	std::string no = "no";
	// 	size_type next_step = 1;
	// 	while(true)
	// 	{
	// 		// Set up storage vector
	// 		Storage<number_type> data(model_.n_parameters(), 2); // records values of fitting vector and two additional things
			
	// 		// Step 1: Set new inverse temperature beta
	// 		if (next_step == 1)
	// 		{
	// 			number_type mean;
	// 			number_type dev;
	// 			get_H(1e3, mean, dev);
	// 			std::cout << "Logarithmic mean of the action: " << mean << " +/- " << dev << "\n";
	// 			std::cout << "Set the value of beta of the system (currently ";
	// 			std::cout << beta_ << "): ";
	// 			std::cin >> beta_;
			
	// 			next_step = 2;
	// 		}



	// 		// Step 2: Set step size
	// 		while(next_step == 2)
	// 		{
	// 			std::cout << "Set a step size (currently: ";
	// 			std::cout << stepsize_ << "): ";
	// 			std::cin >> stepsize_;
	// 			n_steps_min_ = 4;
	// 			n_steps_max_ = 5;
	// 			get_acceptance_rates(range_min_, range_max_, 50, 50, "preliminary_tools/acceptrates.txt");
	// 			std::cout << "Enter \"repeat\" to set a new step size or go to next step (\"next\"), or back (\"back\"): \n";
	// 			std::cin >> user_input;
	// 			if (user_input == next)
	// 			{
	// 				next_step = 3;
	// 			}
	// 			else if (user_input == back)
	// 			{
	// 				next_step = 1;
	// 			}
	// 			else
	// 			{
	// 				next_step = 2;
	// 			}
	// 		}
	
	// 		// Step 3: Set number of leapfrog steps
	// 		while (next_step == 3)
	// 		{
	// 			std::cout << "Set minimal number of leapfrog steps: ";
	// 			std::cin >> n_steps_min_;
	// 			std::cout << "Set maximal number of leapfrog steps: ";
	// 			std::cin >> n_steps_max_;
	// 			std::cout << "Set number of chains: ";
	// 			size_type n_chains;
	// 			std::cin >> n_chains;
	// 			get_acceptance_rates(range_min_, range_max_, n_chains, 50, "preliminary_tools/acceptrates.txt");
	// 			std::cout << "Enter \"repeat\" to try a different range of Leapfrog steps, \"next\" if you want to go the next step, or \"back\" to get to the previous step ";
	// 			std::cin >> user_input;
	// 			if (user_input == repeat)
	// 			{
	// 				next_step = 3;
	// 			}
	// 			else if (user_input == next)
	// 			{
	// 				next_step = 4;
	// 			}
	// 			else if (user_input == back)
	// 			{
	// 				next_step = 2;
	// 			}
				
	// 		}


			
	// 		// Step 4: Generate Markov chains
	// 		while (next_step == 4)
	// 		{
	// 			// Set timer
	// 			std::time_t start = std::time(nullptr);
	// 			number_type expected_duration;
	// 			bool estimate_given = false;

	// 			size_type nb_initials = 1;
	// 			size_type nb_steps = 2e2;

	// 			bool data_saved = false;

	// 			for (size_type i = 1; i <= nb_initials; ++i)
	// 			{
	// 				bool bad_starting_point = false;
	// 				// reinitialize internal counters
	// 				counter_ = 0;
	// 				counter_accepted_ = 0;
	// 				if (i == 1)
	// 				{
	// 					std::cout << "Generating the first Markov chain... \n";
	// 				}
	// 				// choose random starting point that respects constraints
	// 				Vector<number_type> initial(model_.n_parameters());
	// 				find_start(initial);
			
			
	// 				// Start walking
	// 				size_type counter = 0;
	// 				bool region_change = false;
	// 				while (counter < nb_steps)
	// 				{
	// 					if (counter == 10 && acceptance_rate() < 0.01) // break if a bad initial condition was chosen
	// 					{
	// 						--i;
	// 						std::cout << "Bad starting point. Recalculating current chain...\n";
	// 						bad_starting_point = true;
	// 						data.clear();
	// 						break;
	// 					}
				
	// 					if (acceptance_rate() < 0.95 && !region_change)
	// 					{
	// 						region_change = true;
						
	// 					}
	// 					step_forward(initial);

	// 					// save updated data to storage
	// 					data.read_in(initial); // save fitting parameters
	// 					data.read_in(action(initial)/beta_); // save chi2_red
	// 					data.read_in(acceptance_rate()); // save acceptance rate

	// 					// Calculate time
	// 					std::time_t end = std::time(nullptr);
	// 					number_type diff = end - start;

	// 					// Save data after each five minutes
	// 					if ((size_type(diff) % (5*60)) == 0 && !data_saved && diff > 1)
	// 					{
	// 						// write data to output file
	// 						data.write("data.txt");
	// 						data_saved = true;
	// 						std::cout << "Data saved to file.\n";
	// 					}
	// 					// after data has been written, no data can be rewritten for a minute
	// 					if ((size_type(diff) % (5*60)) == 1) 
	// 					{
	// 						data_saved = false;
	// 					}

	// 					counter++;

	// 					if (counter == nb_steps && i == 1) // set length of chain and number of chain
	// 					{
	// 						data.write("data.txt");

	// 						// autocorrelation time analysis
	// 						Vector<number_type> autocorrelation_times(initial.size());
	// 						Vector<number_type> autocorrelation_times_err(initial.size());

	// 						data.determine_burn_in_length();
	// 						data.autocorr_time(autocorrelation_times, autocorrelation_times_err);

	// 						std::cout << "Autocorrelation times: \n";
	// 						for (size_type i = 0; i < autocorrelation_times.size(); ++i)
	// 						{
	// 							std::cout << "Parameter " << i << " : " << autocorrelation_times[i] << " +/- " << autocorrelation_times_err[i]  << "\n";
	// 						}




	// 						std::cout << "End of the first chain reached after " << diff/60 << "min. \n";
	// 						std::cout << "Set chain length to a specific value or continue with the other chains (\"continue\").\n";
	// 						std::cout << "If chain is supposed to be recalculated, enter \"repeat\".\n";
	// 						std::cin >> user_input;
	// 						if (user_input != next_chain && user_input != repeat)
	// 						{
	// 							nb_steps = size_type(std::stod(user_input));
	// 							std::cout << "Change beta to (currently ";
	// 							std::cout << beta_ << "): ";
	// 							std::cin >> beta_;
	// 						}
	// 						if (user_input == repeat)
	// 						{
	// 							i = 0;
	// 							start = std::time(nullptr);
	// 							region_change = false;
	// 							data.clear();
	// 							break;
	// 						}

	// 						if (user_input == next_chain)
	// 						{
	// 							std::cout << "How many chains are supposed to be generated? ";
	// 							std::cin >> nb_initials;
	// 						}
	// 					}
	// 				}
	// 				if (!bad_starting_point)
	// 				{
	// 					std::cout << "Chain " << i << " / " << nb_initials << " generated.\n";
	// 				}

	// 				if (i == nb_initials)
	// 				{
	// 					data.write("data.txt");
	// 					std::cout << "Requested chains have been generated. Increase number of chains, exit (\"next\") or go back one step (\"back\"): ";
	// 					std::cin >> user_input;
	// 					if (user_input != next && user_input != back)
	// 					{
	// 						nb_initials = size_type(std::stod(user_input));
	// 						std::cout << "Change chain length to (currently ";
	// 						std::cout << nb_steps << "): ";
	// 						std::cin >> nb_steps;
	// 						std::cout << "Change minimal number of leapfrog steps to (currently ";
	// 						std::cout << n_steps_min_ << "): ";
	// 						std::cin >> n_steps_min_;
	// 						std::cout << "Change maximal number of leapfrog steps to (currently ";
	// 						std::cout << n_steps_max_ << "): ";
	// 						std::cin >> n_steps_max_;
	// 					}
	// 					else if (user_input == next)
	// 					{
	// 						next_step = 5;
	// 					}
	// 					else
	// 					{
	// 						next_step = 3;
	// 					}
	// 				}
					
	// 			}
			
	// 		}

	// 		if (next_step == 5)
	// 		{
	// 			break;
	// 		}

	// 	}
			
	// }


	number_type acceptance_rate()
	{
		return 1.0 * counter_accepted_ / counter_;
	}




	// 	/* preliminary run that returns the mean and standard deviation of the order of the action in the search region
	// 		=> action ~ 10 to the power of order
	// 	*/
	// void get_H(size_type nb_positions, number_type & pot_mean, number_type & pot_std)
	// {
	// 	Vector<number_type> H_values(0);
		
	// 	for (size_type i = 0; i<nb_positions; ++i)
	// 	{
	// 		// Draw random but permitted position from the search region
	// 		Vector<number_type> popt(range_min_.size());
	// 		find_start(popt);

	// 		// Save H for that position
	// 		H_values.push_back(log(model_.action(popt))/log(10.));	
	// 	}

	// 	pot_mean = H_values.mean();
	// 	pot_std = H_values.std();

	// }
	

	// /* preliminary run to estimate correct step size */
	// void get_acceptance_rates(const Vector<number_type> & range_min, const Vector<number_type> & range_max, size_type nb_positions, size_type nb_steps_forward, std::string filename)
	// {
	// 	std::vector<number_type> acceptance_rates(0);

	// 	for (size_type i = 0; i<nb_positions; ++i)
	// 	{
	// 		// Draw random but permitted position from the search region
	// 		Vector<number_type> popt(range_min.size());
	// 		find_start(popt);
			
	// 		// do some HMC steps and save acceptance rate
	// 		counter_ = 0;
	// 		counter_accepted_ = 0;
		
	// 		for (size_type j = 0; j < nb_steps_forward; ++j)
	// 		{			
	// 			step_forward(popt);
	// 		}
			
			
	// 		acceptance_rates.push_back(acceptance_rate());
	// 	}

	// 	// Save acceptance rate vector in a data file
	// 	std::ofstream outfile(filename);
	// 	for (size_type i = 0; i < acceptance_rates.size(); ++i)
	// 	{
	// 		outfile << acceptance_rates[i] << "\n";
	// 	}

		
	// }

	// /* preliminary run to estimate number of leapfrog steps */
	// void get_optimal_number_of_steps(const Vector<number_type> & range_min, const Vector<number_type> & range_max, size_type nb_positions, size_type nb_iterations, std::string filename)
	// {
	// 	std::vector<number_type> number_steps(0);

	// 	for (size_type i = 0; i<nb_positions; ++i)
	// 	{
	// 		// Draw random position from the search region
	// 		Vector<number_type> popt(range_min.size());
	// 		fill_from_region(popt, range_min, range_max);
			
	// 		// do one leapfrog step and save position values for each iteration
	// 		Storage<number_type> q_values(popt);
	// 		std::random_device rd;
	// 		std::mt19937 gen(rd());
	// 		std::normal_distribution<> dis_norm(0, 1);
	// 		Vector<number_type> p(popt.size());
	// 		fill_random(p, dis_norm);
	// 		leapfrog(popt, p, q_values, nb_iterations);
	// 		//forest_ruth(popt, p, q_values, nb_iterations);

	// 		// calculate period of each position component using autocorrelation
	// 		Vector<size_type> period_vector(popt.size());
	// 		q_values.period(period_vector);
	

	// 		// Half a period corresponds to the optimal number of steps
	// 		for (size_type j = 0; j < period_vector.size(); ++j)
	// 		{
	// 			number_steps.push_back(period_vector[j]/2);
	// 		}
			
	// 	}

	// 	// Save array of optimal number of steps in a data file
	// 	std::ofstream outfile(filename);
	// 	for (size_type i = 0; i < number_steps.size(); ++i)
	// 	{
	// 		outfile << number_steps[i] << "\n";
	// 	}	
	// }


private:
	/* measured data */
	LatticeAction<number_type> action_;


	/* parameters for the leapfrog integrator */
	number_type stepsize_;
	size_type n_steps_min_;
	size_type n_steps_max_;

	/* some statistics */
	size_type counter_;
	size_type counter_accepted_;
	

};

#endif