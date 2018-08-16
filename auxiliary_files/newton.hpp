#ifndef NEWTON_HPP
#define NEWTON_HPP

#include "vector.hpp"
#include "matrix.hpp"


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
	, J_(model.size(), model.size(), 0)
	, maxit_(10*model.size()+10)
	, limit_residual_(1.e-12)
	, limit_param_(1.e-12)
	{}

	// Central difference jacobian column i at position x 
	void first_derivative(const Vector<number_type> & x, Vector<number_type> & func, size_type i, number_type stepsize)
	{
		assert(x.size()==model_.size());
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



private:
	Model<number_type> model_;

	// Jacobian
	Matrix<number_type> J_;

	// algorithm-specific parameters
	size_type maxit_;
	number_type limit_residual_;
	number_type limit_param_;
};


#endif