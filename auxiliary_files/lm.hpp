#ifndef LM_HPP
#define LM_HPP

#include "vector.hpp"
#include "matrix.hpp"
#include <cmath>
#include <algorithm>
#include "storage.hpp"
#include "to_file.hpp"



template <class REAL>
class LM
{
public:
	typedef std::size_t size_type;
	typedef REAL number_type;
	// Constructor with x-uncertainties
	LM(
	Model<number_type> model,
	std::vector<Vector<number_type>> x,
	std::vector<Vector<number_type>> dx,
	Vector<number_type> y,
	Vector<number_type> dy
	)
	: model_(model)
	, x_transposed_(x)
	, x_(x[0].size())
	, dx_transposed_(dx)
	, dx_(dx[0].size())
	, x_error_(true)
	, y_(y)
	, dy_(dy)
	, sigma2_(dy.size())
	, J_(x[0].size(), model.n_parameters(), 0)
	, best_fit_(model.n_parameters())
	, solution_found_(false)
	, n_param_(model.n_parameters())
	, n_data_(x[0].size())
	, d_of_freedom_(x[0].size()-model.n_parameters())
	, maxit_(10*model.n_parameters()+10)
	, lambda0_(1.e-2)
	, scale_up_(11.)
	, scale_down_(9.)
	, limit_gradient_(1.e-12)
	, limit_param_(1.e-12)
	, eps_lm_(1.e-1)
	, converged_(false)
	{
		transpose_x_data();
	}
	// Constructor without x-uncertainties
	LM(
	Model<number_type> model,
	std::vector<Vector<number_type>> x,
	Vector<number_type> y,
	Vector<number_type> dy
	)
	: model_(model)
	, x_transposed_(x)
	, x_(x[0].size())
	, dx_transposed_(x)
	, dx_(x[0].size())
	, x_error_(false)
	, y_(y)
	, dy_(dy)
	, sigma2_(dy.size())
	, J_(x[0].size(), model.n_parameters(), 0)
	, best_fit_(model.n_parameters())
	, solution_found_(false)
	, n_param_(model.n_parameters())
	, n_data_(x[0].size())
	, d_of_freedom_(x[0].size()-model.n_parameters())
	, maxit_(5*model.n_parameters()+10) // Recommended: 10*model.n_parameters()+10
	//, maxit_(model.n_parameters())
	, lambda0_(1.e-2) // Recommended: 1.e-2
	, scale_up_(11.) // Recommended: 11
	, scale_down_(9.) // Recommended: 9
	, limit_gradient_(1.e-12) // Recommended: 1e-12
	, limit_param_(1.e-12) // Recommended: 1e-12
	, eps_lm_(1.e-4) // Recommended: 1e-1
	, converged_(false)
	{
		transpose_x_data();
		set_sigma();
	}

	// Reset LM object
	void clear()
	{
		converged_ = false;
	}

	/* 
	Experimental x-data is read in and saved in a vector x_transposed of the form:
			x1	y1	z1
			x2,	y2,	z2
			x3	y3	z3
			etc.
	However, we need the following form:
			x1	x2	x3	
			y1,	y2,	y3	etc.
			z1	z2	z3

	*/
	void transpose_x_data()
	{
		for (size_type i = 0; i < x_.size(); ++i)
		{
			for (size_type j = 0; j < x_transposed_.size(); ++j)
			{
				(x_[i]).push_back((x_transposed_[j])[i]);
				if (x_error_)
				{
					(dx_[i]).push_back((dx_transposed_[j])[i]);
				}
			}
		}
	}



	void print_x_data() const
	{
		for (size_type i = 0; i < x_.size(); ++i)
		{
			for (size_type j = 0; j < x_[0].size(); ++j)
			{
				std::cout << (x_[i])[j] << " +/- " << (dx_[i])[j] << " ";
			}
			std::cout << "\n";
		}
	}

	// Approximation of partial derivative of f(x,q) with respect to a parameter at position x 
	number_type first_derivative(Vector<number_type> x, const Vector<number_type> & q, size_type parameter_index, number_type stepsize)
	{
		assert(q.size()==n_param_);
		assert(parameter_index >= 0 && parameter_index < q.size());
		Vector<number_type> position_forward = q;
		position_forward[parameter_index] += stepsize;
		Vector<number_type> position_backward = q;
		position_backward[parameter_index] -= stepsize;
		// std::cout << model_.f(x, position_forward) << " " << model_.f(x, position_backward) << "\n";
		return (model_.f(x, position_forward)-model_.f(x, position_backward))/2.0/stepsize;
	}

	// unbiased estimator of first derivative: weighted sum of evaluations with two different step sizes
	number_type first_derivative_unbiased(Vector<number_type> x, const Vector<number_type> & q, size_type parameter_index, number_type stepsize)
	{
		number_type value_h = first_derivative(x, q, parameter_index, stepsize);
		number_type value_h2 = first_derivative(x, q, parameter_index, stepsize/2);
		return (4.*value_h2-value_h)/3.;
	}

	// Approximation of positional derivative of f(x,q) with respect to x with parameters q 
	number_type derivative_x(const Vector<number_type> & x, Vector<number_type>  q, size_type parameter_index, number_type stepsize)
	{
		assert(x.size()==dx_[0].size());
		assert(parameter_index >= 0 && parameter_index < x_[0].size());
		Vector<number_type> position_forward = x;
		position_forward[parameter_index] += stepsize;
		Vector<number_type> position_backward = x;
		position_backward[parameter_index] -= stepsize;
		return (model_.f(position_forward, q)-model_.f(position_backward, q))/2.0/stepsize;
	}

	// unbiased estimator of first derivative: weighted sum of evaluations with two different step sizes
	number_type derivative_x_unbiased(const Vector<number_type> & x, Vector<number_type> q, size_type parameter_index, number_type stepsize)
	{
		number_type value_h = derivative_x(x, q, parameter_index, stepsize);
		number_type value_h2 = derivative_x(x, q, parameter_index, stepsize/2);
		return (4.*value_h2-value_h)/3.;
	}

	// Getter for jacobian
	number_type J(size_type row_index, size_type column_index) const
	{
		assert(row_index < n_data_ && row_index >= 0);
		assert(column_index < n_param_ && column_index >= 0);
		return J_(row_index, column_index);
	}

	// Getter for convergence boolean
	bool converged() const
	{
		return converged_;
	}

	// compute Jacobian for set of parameter values q
	void set_J(const Vector<number_type> & q)
	{
		//std::cout << "Computing Jacobian...\n";
		assert(q.size() == n_param_);
		for (size_type i = 0; i < n_data_; ++i)
		{
			for (size_type j = 0; j < n_param_; ++j)
			{
				number_type stepsize = 1.e-5;
				if ((x_[i])[j]*1.e-5 > stepsize)
				{
					stepsize *= (x_[i])[j];
				}
				J_(i, j) = first_derivative_unbiased(x_[i], q, j, stepsize);


				//if ((i%(n_data_/10) == 0)&&(j==0))
				//	std::cout << "Progress: " << round(100.*i/n_data_) << "%\n";

			}
		}
		//std::cout << "Jacobian has been computed.\n";
	}

	// compute entries of sigma2 vector where sigma^2 = dy^2 + (df/dx)^2*dx^2
	void set_sigma(const Vector<number_type> & q)
	{
		for (size_type i = 0; i < n_data_; ++i)
		{
			sigma2_[i] = dy_[i]*dy_[i];
			if (x_error_)
			{
				for (size_type j = 0; j < dx_[0].size(); ++j)
				{
					sigma2_[i] += ( (dx_[i])[j] * (dx_[i])[j] * derivative_x_unbiased(x_[i], q, j, 1.e-5) * derivative_x_unbiased(x_[i], q, j, 1.e-5) );
				}
			}
		}
	}

	// compute entries of sigma2 vector where sigma^2 = dy^2 (overload)
	void set_sigma()
	{
		for (size_type i = 0; i < n_data_; ++i)
		{
			sigma2_[i] = dy_[i]*dy_[i];
		}
	}

	// fill left-hand side of lm equation
	void fill_left(Matrix<number_type> & A, const Vector<number_type> & q, number_type lambda, Vector<number_type> & diagonals)
	{
		// Note: Jacobian must have been computed before calling this method!
		A = 0;
		for (size_type i = 0; i < A.rowsize(); ++i)
		{
			for (size_type j = 0; j <=i; ++j) //matrix symmetrical
			{
				// Compute Jt*W*J
				for (size_type k = 0; k < n_data_; ++k)
				{
					A(i,j) += J_(k,i)*J_(k,j)/sigma2_[k];
				}

				// Add Levenberg-Marquardt parameter
				if (i == j)
				{
					diagonals[i] = A(i,j);
					A(i,j) *= (lambda + 1.);
				}

				// matrix symmetrical
				if (i != j)
				{
					A(j,i) = A(i,j);
				}
			}
		}
	}

	// fill right-hand side of lm equation
	void fill_right(Vector<number_type> & r, const Vector<number_type> & q)
	{
		r = 0;
		for (size_type i = 0; i < r.size(); ++i)
		{
			for (size_type k = 0; k < n_data_; ++k)
			{
				r[i] += J_(k,i)/sigma2_[k]*(y_[k]-model_.f(x_[k], q));
			}
		}
	}

	// compute current chi2 value
	number_type chi2(const Vector<number_type> & q)
	{
		number_type result = 0;
		for (size_type i = 0; i < n_data_; ++i)
		{
			result += pow((y_[i] - model_.f(x_[i], q)), 2)/sigma2_[i];
		}

		return result;
	}

	// compute reduced chi2
	number_type chi2_red(const Vector<number_type> & q)
	{
		return chi2(q)/d_of_freedom_;
	}

	// measure for how lambda is to be changed
	number_type check_improvement(const Vector<number_type> & q, const Vector<number_type> & h, number_type lambda, const Vector<number_type> & diagonals, const Vector<number_type> & residual)
	{
		number_type result(0.);

		Vector<number_type> y(n_param_);
		y = diagonals*h;
		y *= lambda;
		y += residual;
		result = (chi2_red(q) - chi2_red(q+h))/inner_product(h, y);

		std::cout  << "Inner product: " << inner_product(h, y) << " " << "Metric: " << result << "\n";

		return result;
	}

	bool has_converged (const Vector<number_type> & q, const Vector<number_type> & h, const Vector<number_type> & residual)
	{
		bool result = false;

		// convergence in the gradient?
		Vector<number_type> r_abs(residual.size());
		r_abs = residual;
		abs(r_abs);
		//std::cout << "Gradient: " << *std::max_element(r_abs.begin(), r_abs.end()) << "\n";
		if(*std::max_element(r_abs.begin(), r_abs.end()) < limit_gradient_)
		{
			result = true;
			return result;
		}

		// convergence in parameters?
		Vector<number_type> h_scaled(h.size());
		h_scaled = h/q;
		abs(h_scaled);
		//std::cout << "param: " << *std::max_element(h_scaled.begin(), h_scaled.end()) << "\n";
		if(*std::max_element(h_scaled.begin(), h_scaled.end()) < limit_param_)
		{
			result = true;
			return result;
		}

		return result;

	}

	// solve curve-fitting problem
	void solve(Vector<number_type> & q)
	{
		// reinitialize LM object
		clear();
		// right-hand side / residual
		Vector<number_type> r(n_param_);
		// left-hand side
		Matrix<number_type> A(n_param_, n_param_);
		// solution of linear equation
		Vector<number_type> z(n_param_);

		// scaling factors
		Vector<number_type> s(n_param_);
		// row permutations
		Vector<size_type> perm_r(n_param_);
		// column permutations
		Vector<size_type> perm_c(n_param_);

		// Levenberg-Marquardt parameter
		number_type lambda = lambda0_;

		bool non_invertible = false; // boolean for whether coefficient matrix is non-invertible

		// Compute Jacobian at the beginning
		set_J(q);
		set_sigma(q);

		for (size_type i = 0; i < maxit_; ++i)
		{
			// for (size_type i = 0; i < q.size(); ++i)
			// {
			// 	std::cout << i << " : " << q[i] << "\n";
			// }
			std::cout << i << ": " << "chi2/d.o.f = " << chi2_red(q) << " lambda: " << lambda << "\n";
			// std::cout << non_invertible << "\n";
			// Fill entries of coefficient matrix and save diagonal entries of JtWJ
			Vector<number_type> diagonals(n_param_);
			fill_left(A, q, lambda, diagonals);
			// A.print();

			// Fill entries of residual vector and make a copy
			fill_right(r, q);
			Vector<number_type> r_copy(r);

			// solve system for update
			row_equilibrate(A, s, non_invertible); // equilibrate rows
			if (non_invertible) // break in case of a non-invertible matrix
			{
				std::cout << "Not invertible! \n";
				converged_ = false;
				break;
			}
			lu_fullpivot(A, perm_r, perm_c); // LU decomposition of A
			z = number_type(0.0); // clear solution
			apply_equilibrate(s, r); // equilibration of right-hand side
			permute_forward(perm_r, r); // permutation of right-hand side
			solveL(A, r, r); // forward substitution
			solveU(A, z, r); // backward substitution
			permute_backward(perm_c, z); // backward permutation

			// check whether solution sufficiently decreases chi2
			number_type metric = check_improvement(q, z, lambda, diagonals, r_copy);
			//std::cout << "Metric: " <<  metric << "\n";

			if (metric > eps_lm_)
			{
				q += z;
				lambda = std::max(lambda/scale_down_, 1e-7);
			}
			else
			{
				lambda = std::min(lambda*scale_up_, 1e7);
			}
			

			
			

			// Recomputing jacobian
			set_J(q); 
			set_sigma(q);



			// check convergence
			if (has_converged(q, z, r_copy))
			{
				std::cout << "Converged! \n";
				converged_ = true;
				break;
			}

		}

	}



	// find best-fit result by drawing n random starting points from a given region
	// return ratio of number of starting points that led to the best-fit result
	number_type find_best_fit(Vector<number_type> & popt, const Vector<number_type> & lower_range, const Vector<number_type> & upper_range, size_type n)
	{
		assert(popt.size() == n_param_);
		assert(lower_range.size() == n_param_);
		assert(upper_range.size() == n_param_);

		// temporary best-fit
		Vector<number_type> popt_temp(popt.size());
		number_type chi2_temp;
		size_type counter_best_fit = 0;
		number_type tol = 1.e-8;

		Storage<number_type> popts(popt, 1); // Save popt and chi2/dof for each starting point

		for (size_type i = 0; i < n; ++i)
		{
			// Draw random starting point
			fill_from_region(popt, lower_range, upper_range);

			// minimize from that point
			solve(popt);
			
			

			for (size_type i = 0; i < popt.size(); ++i)
			{
				std::cout << i << " : " << popt[i] << "\n";
			}

	

			// compare to minimum found so far
			number_type chi2_current = chi2_red(popt);
			number_type chi2_diff = chi2_current - chi2_temp;
			if ((chi2_diff < -tol && converged_) || (counter_best_fit == 0 && converged_))
			{
				chi2_temp = chi2_current;
				popt_temp = popt;
				counter_best_fit = 1;
				solution_found_ = true;
			}
			else if (abs(chi2_diff) < tol)
			{
				counter_best_fit++;
			}

			// Save in storage
			popts.read_in(popt);
			popts.read_in(chi2_current);



		}

		// save result
		if (solution_found_)
		{
			best_fit_ = popt_temp;
			popt = popt_temp;
		}

		popts.write("popts.txt", false);

		return 1.*counter_best_fit/n;

	}

	// Monte-Carlo uncertainty analysis by generating n experimental data sets
		void get_fit_uncertainty(std::string filename, size_type n, Vector<number_type> & uncertainty, Vector<number_type> & uncertainty_error)
		{
			// Data storage for the best-fit result vectors
			Storage<number_type> popts(best_fit_, 1);

			// For this method, the solution must have been found beforehand
			assert(solution_found_);

			// save original x- and y-values
			std::vector<Vector<number_type>> x(x_);
			Vector<number_type> y(y_);


			// generate experimental data
			for (size_type i = 0; i < n; ++i)
			{
				// Generate x-values
				if (x_error_)
				{
					std::vector<Vector<number_type>> x_values(x_);
					for (size_type j = 0; j < n_data_; ++j)
					{
						fill_gaussian(x_values[j], x[j], dx_[j]);
					}
					x_ = x_values;
				}
			
				// Generate y-values
				Vector<number_type> y_values(y_);
				fill_gaussian(y_values, y, dy_);
				y_ = y_values;

				// solve minimization problem with solution of original
				// problem as starting point
				Vector<number_type> popt(best_fit_.size());
				popt = best_fit_;
				solve(popt);

				// save result to storage
				if (converged_)
				{
					popts.read_in(popt);
					popts.read_in(chi2_red(popt));
				}
			}

			// Restore original x- and y-values
			x_ = x;
			y_ = y;

			// Calculate uncertainty result
			popts.std_deviation(uncertainty, uncertainty_error);

			// print output file
			popts.write(filename, false);
		}

private:
	// theory model and experimental data
	Model<number_type> model_;
	std::vector<Vector<number_type>> x_transposed_;
	std::vector<Vector<number_type>> x_;
	std::vector<Vector<number_type>> dx_transposed_;
	std::vector<Vector<number_type>> dx_;
	bool x_error_;
	Vector<number_type> y_;
	Vector<number_type> dy_;
	Vector<number_type> sigma2_;

	// Jacobian
	Matrix<number_type> J_;

	//	solution
	Vector<number_type> best_fit_;
	bool solution_found_;

	// useful constants
	size_type n_param_;
	size_type n_data_;
	size_type d_of_freedom_;

	// parameters specific to Levenberg-Marquardt algorithm
	size_type maxit_;
	number_type lambda0_;
	number_type scale_up_;
	number_type scale_down_;

	// convergence criteria
	number_type limit_gradient_;
	number_type limit_param_;
	number_type eps_lm_;
	bool converged_;





};

#endif