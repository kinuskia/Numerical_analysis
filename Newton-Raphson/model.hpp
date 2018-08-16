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

	/* dimension of the problem */
	size_type size() const
	{
		return 2;
	}

	

	/* vector function	f(x) = 0 */
	void f(const Vector<number_type> & x, Vector<number_type> & result)
	{
		assert(result.size() == x.size());
		assert(x.size() == this->size());
		result[0] = 2.*x[0]*x[0] - x[1] - 7.;
		result[1] = x[0]*x[1] + 3.*x[1]*x[1]*x[1] - 5.;
	}
	



private:
	




};


#endif