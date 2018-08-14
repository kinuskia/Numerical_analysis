#ifndef LM_HPP
#define LM_HPP

#include "vector.hpp"



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
	Vector<number_type> dy,
	number_type lambda0
	)
	: model_(model)
	, x_(x)
	, y_(y)
	, dy_(dy)
	, J_(x.size(), model.n_parameters(), 0)
	, n_param_(model.n_parameters())
	, n_data_(x.size())
	, lambda_(lambda0)
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
	number_type lambda_;



};

#endif