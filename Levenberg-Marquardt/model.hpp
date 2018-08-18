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
	) // default constructor
	{}

	
public:

	/* number of fitting parameters */
	size_type n_parameters() const
	{
		return 3;
	}

	/* dimension of the x vector */
	size_type x_dim() const
	{
		return 1;
	}

	/* model fitting function */
	number_type f(Vector<number_type> x, Vector<number_type> popt) const
	{
		assert(x.size() == this->x_dim());
		return popt[0] * exp(-(x[0]-popt[1])*(x[0]-popt[1])/2/popt[2]/popt[2]);
	}
	



private:
	




};


#endif