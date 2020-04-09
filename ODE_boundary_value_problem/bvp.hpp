#ifndef BVP_HPP
#define BVP_HPP

#include "../auxiliary_files/vector.hpp"
#include "../auxiliary_files/matrix.hpp"
#include <stdexcept>



template<class N>
class BVP
{
public:
	typedef N number_type;
	typedef std::size_t size_type;
	BVP(number_type cw, size_type dim)
	: cw_(cw)
	, dim_(dim)
	{};


	// return derivative of the problem
	void F(Vector<number_type> & deriv, const Vector<number_type> & state) const
	{
		if (deriv.size() != state.size()-1 || deriv.size() != dim_)
			throw std::runtime_error("State and derivative vector do not match the size of the ODE.");
		deriv[0] = state[2]; 
		deriv[1] = cw_/2*state[2]*sqrt(1+state[4]*state[4]);
		deriv[2] = state[4]; 
		deriv[3] = -state[2]*state[2];
	}


	void dF(Matrix<number_type> & J, const Vector<number_type> & state) const
	{
		bool correct_input = (J.colsize() == J.rowsize() && J.colsize() == dim_ && J.colsize() == state.size() - 1);
		if (!correct_input)
			throw std::runtime_error("State vector and Jacobian don't match.");
		J[0][0] = 0;
		J[0][1] = 1;
		J[0][2] = 0; 
		J[0][3] = 0;
		J[1][0] = 0;
		J[1][1] = cw_/2*sqrt(1+state[4]*state[4]);
		J[1][2] = 0;
		J[1][3] = cw_/2*state[2]*state[4]/sqrt(1+state[4]*state[4]);
		J[2][0] = 0;
		J[2][1] = 0;
		J[2][2] = 0; 
		J[2][3] = 1;
		J[3][0] = 0;
		J[3][1] = -2*state[2];
		J[3][2] = 0;
		J[3][3] = 0;
	}

	// getter for member variables
	size_type size() const
	{
		return dim_;
	}
	// return boundary condition
	//void bound(Vector<number_type> & r, con)





private:
	number_type cw_;
	size_type dim_;

};








#endif