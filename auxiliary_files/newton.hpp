#ifndef NEWTON_HPP
#define NEWTON_HPP

#include "vector.hpp"
#include "matrix.hpp"
#include <algorithm>


template<class REAL>
class Newton
{
public:
	typedef std::size_t size_type;
	typedef REAL number_type;
	// Default constructor
	Newton(
	Model<number_type> model
	)
	: model_(model)
	, dim_(model.size())
	, J_(model.size(), model.size(), 0)
	, maxit_(10*model.size()+10)
	, limit_residual_(1.e-12)
	, converged_(false)
	{}

	// Reset Newton object
	void clear()
	{
		converged_ = false;
	}

	// Central difference jacobian column i at position x 
	void first_derivative(const Vector<number_type> & x, Vector<number_type> & func, size_type i, number_type stepsize)
	{
		assert(x.size()==dim_);
		assert(x.size()==func.size());
		assert(i >= 0 && i < x.size());

		Vector<number_type> position_forward = x;
		position_forward[i] += stepsize;
		Vector<number_type> position_backward = x;
		position_backward[i] -= stepsize;
		Vector<number_type> f_forward(x.size());
		Vector<number_type> f_backward(x.size());
		model_.f(position_forward, f_forward);
		model_.f(position_backward, f_backward);
		func = (f_forward-f_backward)/2./stepsize;
	}

	// Getter for jacobian
	number_type J(size_type row_index, size_type column_index) const
	{
		assert(row_index < dim_ && row_index >= 0);
		assert(column_index < dim_ && column_index >= 0);
		return J_(row_index, column_index);
	}

	// compute Jacobian at position x
	void set_J(const Vector<number_type> & x)
	{
		assert(x.size() == dim_);
		for (size_type i = 0; i < dim_; ++i)
		{
			Vector<number_type> column(dim_);
			first_derivative(x, column, i, 1.e-5);
			for (size_type j = 0; j < dim_; ++j)
			{
				J_(j, i) = column[j];
			}
		}
	}

	bool has_converged (const Vector<number_type> & q, const Vector<number_type> & h, const Vector<number_type> & residual)
	{
		bool result = false;

		// convergence in the residual?
		Vector<number_type> r_abs(residual.size());
		r_abs = residual;
		abs(r_abs);
		if(*std::max_element(r_abs.begin(), r_abs.end()) < limit_residual_)
		{
			result = true;
			return result;
		}

		return result;

	}

	// Getter for convergence boolean
	bool converged() const
	{
		return converged_;
	}

	// solve Newton problem
	void solve(Vector<number_type> & q)
	{
		clear();
		// right-hand side / residual
		Vector<number_type> r(dim_);
		// left-hand side
		Matrix<number_type> A(dim_, dim_);
		// solution of linear equation
		Vector<number_type> z(dim_);

		// scaling factors
		Vector<number_type> s(dim_);
		// row permutations
		Vector<size_type> perm_r(dim_);
		// column permutations
		Vector<size_type> perm_c(dim_);

		bool non_invertible = false; // boolean for whether coefficient matrix is non-invertible

		for (size_type i = 0; i < maxit_; ++i)
		{
			// Fill entries of coefficient matrix
			set_J(q);
			A = J_;

			// Fill entries of residual vector and make a copy
			model_.f(q, r);
			r = -r;
			Vector<number_type> r_copy(r);

			// solve system for update
			row_equilibrate(A, s, non_invertible); // equilibrate rows
			if (non_invertible) // break in case of a non-invertible matrix
			{
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

			// update solution vector
			q += z;
				


			// check convergence
			if (has_converged(q, z, r_copy))
			{
				converged_ = true;
				break;
			}

		}

	}



private:
	Model<number_type> model_;
	size_type dim_; // number of parameters of the problem

	// Jacobian
	Matrix<number_type> J_;

	// algorithm-specific parameters
	size_type maxit_;
	number_type limit_residual_;
	bool converged_;
};


#endif