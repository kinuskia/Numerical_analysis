#ifndef LM_HPP
#define LM_HPP

#include "vector.hpp"
#include "matrix.hpp"
#include <cmath>
#include <algorithm>



template <class REAL>
class LM
{
public:
	typedef std::size_t size_type;
	typedef REAL number_type;
	// Default constructor
	LM(
	Model<number_type> model,
	Vector<number_type> x,
	Vector<number_type> y,
	Vector<number_type> dy
	)
	: model_(model)
	, x_(x)
	, y_(y)
	, dy_(dy)
	, J_(x.size(), model.n_parameters(), 0)
	, n_param_(model.n_parameters())
	, n_data_(x.size())
	, maxit_(10*model.n_parameters()+10)
	, maxit_lambda_(10)
	, lambda0_(1.e-2)
	, scale_up_(11.)
	, scale_down_(9.)
	, limit_gradient_(1.e-3)
	, limit_param_(1.e-3)
	, limit_chi2_(1.e-1)
	, eps_lm_(1.e-1)
	{}

	// Approximation of partial derivative of f(x,q) with respect to a parameter at position x 
	number_type first_derivative(number_type x, const Vector<number_type> & q, size_type parameter_index, number_type stepsize)
	{
		assert(q.size()==n_param_);
		assert(parameter_index >= 0 && parameter_index < q.size());
		Vector<number_type> position_forward = q;
		position_forward[parameter_index] += stepsize;
		Vector<number_type> position_backward = q;
		position_backward[parameter_index] -= stepsize;
		return (model_.f(x, position_forward)-model_.f(x, position_backward))/2.0/stepsize;
	}

	// Getter for jacobian
	number_type J(size_type row_index, size_type column_index) const
	{
		assert(row_index < n_data_ && row_index >= 0);
		assert(column_index < n_param_ && column_index >= 0);
		return J_(row_index, column_index);
	}

	// compute Jacobian for set of parameter values q
	void set_J(const Vector<number_type> q)
	{
		assert(q.size() == n_param_);
		for (size_type i = 0; i < n_data_; ++i)
		{
			for (size_type j = 0; j < n_param_; ++j)
			{
				J_(i, j) = first_derivative(x_[i], q, j, 1e-3);
			}
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
					A(i,j) += J_(k,i)*J_(k,j)/dy_[k]/dy_[k];
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
				r[i] += J_(k,i)/dy_[k]/dy_[k]*(y_[k]-model_.f(x_[k], q));
			}
		}
	}

	number_type chi2(const Vector<number_type> & q)
	{
		number_type result = 0;
		for (size_type i = 0; i < n_data_; ++i)
		{
			result += pow((y_[i] - model_.f(x_[i], q))/dy_[i], 2);
		}

		return result;
	}

	number_type check_improvement(const Vector<number_type> & q, const Vector<number_type> & h, number_type lambda, const Vector<number_type> & diagonals, const Vector<number_type> & residual)
	{
		number_type result(0.);

		Vector<number_type> y(n_param_);
		y = diagonals*h;
		y *= lambda;
		y += residual;
		result = (chi2(q) - chi2(q+h))/inner_product(h, y);


		return result;
	}

	// solve curve-fitting problem
	void solve(Vector<number_type> & q)
	{
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

		for (size_type i = 0; i < maxit_; ++i)
		{
				std::cout << "Iteration: " << i  << ", chi2: " << chi2(q) << "\n";
				// Compute jacobian
				set_J(q); 
				// Fill entries of coefficient matrix and save diagonal entries of JtWJ
				Vector<number_type> diagonals(n_param_);
				fill_left(A, q, lambda, diagonals);

				// Fill entries of residual vector and make a copy
				fill_right(r, q);
				Vector<number_type> r_copy(r);

				// solve system for update
				row_equilibrate(A, s); // equilibrate rows
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
		}

	}

private:
	// theory model and experimental data
	Model<number_type> model_;
	Vector<number_type> x_;
	Vector<number_type> y_;
	Vector<number_type> dy_;

	// Jacobian
	Matrix<number_type> J_;

	// useful constants
	size_type n_param_;
	size_type n_data_;

	// parameters specific to Levenberg-Marquardt algorithm
	size_type maxit_;
	size_type maxit_lambda_;
	number_type lambda0_;
	number_type scale_up_;
	number_type scale_down_;

	// convergence criteria
	number_type limit_gradient_;
	number_type limit_param_;
	number_type limit_chi2_;
	number_type eps_lm_;





};

#endif