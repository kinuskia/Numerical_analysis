#ifndef IVP_HPP
#define IVP_HPP

#include <stdexcept>
#include "../auxiliary_files/vector.hpp"


template<class N>
class IVP
{
public:
	typedef N number_type;
	typedef std::size_t size_type;
	IVP()
	: dim_(2) // IVP is two-dimensional (position and velocity as functions of time)
	{};


	// return derivative of the problem, i.e. right-hand side of ODE: y' = f(state)
	void F(Vector<number_type> & deriv, const Vector<number_type> & state) const
	{
		/*
		"state" imput:
		state[0]: parameter with respect to which derivatives appear, i.e. here time
		state[1]: position
		state[2]: velocity
		*/
		/*
		"output":
		deriv[0]: derivative of first coordinate (position)
		deriv[1]: derivative of second coordinate (velocity)
		*/
		if (deriv.size() != state.size()-1 || deriv.size() != dim_)
		{
			throw std::runtime_error("State and derivative vector do not match the size of the ODE.");
		}
		deriv[0] = state[2]; 
		deriv[1] = -state[1];
	}


	// getter for member variables
	size_type size() const
	{
		return dim_;
	}


private:
	size_type dim_;

};







#endif