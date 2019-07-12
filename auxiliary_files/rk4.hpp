#ifndef RK4_HPP
#define RK4_HPP

#include "vector.hpp"
#include <stdexcept>



template<class N>
class RK4
{
public:
	typedef std::size_t size_type;
	typedef N number_type;
	// constructors
	RK4 ()
	: stepsize_(number_type(.001))
	{};

	// setter and getter functions for member variable 
	void set_stepsize(number_type stepsize)
	{
		stepsize_ = number_type(stepsize);
	}

	number_type stepsize() const
	{
		return stepsize_;
	}

	template<class M>
	void step_forward(const M & IVP, Vector<typename M::number_type> & state)
	{
		if (state.size() != IVP.size() + 1)
		{
			throw std::runtime_error("Incorrect input of initial conditions.");
		}


		// Runge Kutta step calculation

		//k1's
		Vector<typename M::number_type> state_1 = state;
		Vector<typename M::number_type> k1(IVP.size());
		IVP.F(k1, state_1);

		//k2's
		Vector<typename M::number_type> state_2;
		state_2.push_back(state[0] + stepsize_ / 2);
		for (size_type i = 0; i < IVP.size(); ++i)
		{
			state_2.push_back(state[i+1] + stepsize_ / 2 * k1[i]);
		}
		Vector<typename M::number_type> k2(IVP.size());
		IVP.F(k2, state_2);

		//k3's
		Vector<typename M::number_type> state_3;
		state_3.push_back(state[0] + stepsize_ / 2);
		for (size_type i = 0; i < IVP.size(); ++i)
		{
			state_3.push_back(state[i+1] + stepsize_ / 2 * k2[i]);
		}
		Vector<typename M::number_type> k3(IVP.size());
		IVP.F(k3, state_3);

		//k4's
		Vector<typename M::number_type> state_4;
		state_4.push_back(state[0] + stepsize_);
		for (size_type i = 0; i < IVP.size(); ++i)
		{
			state_4.push_back(state[i+1] + stepsize_ * k3[i]);
		}
		Vector<typename M::number_type> k4(IVP.size());
		IVP.F(k4, state_4);

		// update state
		state[0] += stepsize_;
		for (size_type i = 0; i < IVP.size(); ++i)
		{
			state[i+1] += stepsize_ / 6 * (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]);
		}

	}

private:
	number_type stepsize_;

};




#endif