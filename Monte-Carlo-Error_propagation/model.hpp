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

	/* number of  parameters */
	size_type n_parameters() const
	{
		return 2;
	}

	/* model fitting function */
	number_type f(Vector<number_type> x) const
	{
		assert(x.size() == this->n_parameters());
		return x[0]/x[1];
	}
	



private:
	




};


#endif