#ifndef MODEL_HPP
#define MODEL_HPP
#include <cmath>
#include <assert.h>
#include "../auxiliary_files/vector.hpp"

template<typename REAL>
class Model
{
public:
	typedef std::size_t size_type;
	typedef REAL number_type;
	Model(
	Vector<number_type> t,
	Vector<number_type> y,
	Vector<number_type> dy
	) // default constructor
	: t_(t)
	, y_(y)
	, dy_(dy)
	{}

	/* number of fitting parameters */
	size_type n_parameters() const
	{
		return 3;
	}
public: 
	/* degrees of freedom */
	size_type d_of_freedom() const
	{
		return t_.size() - n_parameters();
	}
private:

	/* model fitting function */
	number_type f(number_type x, Vector<number_type> popt)
	{
		return popt[0] * exp(-(x-popt[1])*(x-popt[1])/2/popt[2]/popt[2]);
	}
	

public:
	/* reduced chi2 sum as potential function to be minimized */
	number_type potential(const Vector<number_type> & q)
	{
		/* checks if data arrays have equal lengths */
		size_type N = t_.size();
		assert(y_.size() == N);
		assert(dy_.size() == N);
		
		number_type chi2 = 0;
		for (size_type k = 0; k<t_.size(); ++k)
		{
			chi2 += pow(y_[k]-f(t_[k], q), 2)/pow(dy_[k], 2);
		}

		return chi2/d_of_freedom();

	}

	/* Constraints to the parameters */
	bool constraints_respected(const Vector<number_type> &q)
	{
		bool respected = true;

		if (q[2] < 0)
		{
			respected = false;
			return respected;
		}

		return respected;
	}

private:
	// experimental data
	Vector<number_type> t_;
	Vector<number_type> y_;
	Vector<number_type> dy_;




};


#endif