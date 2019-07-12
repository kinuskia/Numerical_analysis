#ifndef DOPRI5_HPP
#define DOPRI5_HPP

#include <cmath>
#include <stdexcept>
#include "vector.hpp"

template<class N>
class DOPRI5
{
public:
	typedef std::size_t size_type;
	typedef N number_type;
	// constructors
	DOPRI5 (number_type atol, number_type rtol)
	: stepsize_(number_type(.001))
	, alpha_(0.9)
	, alpha_max_(1.5)
	, alpha_min_(0.5)
	, atol_(atol)
	, rtol_(rtol)
	, fsal_(false)
	{};
	DOPRI5 ()
	: stepsize_(number_type(.001))
	, alpha_(0.9)
	, alpha_max_(1.5)
	, alpha_min_(0.5)
	, atol_(1.e-16)
	, rtol_(1.e-16)
	, fsal_(false)
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


		// Dormand/Prince step calculation

		//k1's FSAL
		Vector<typename M::number_type> state_1 = state;
		Vector<typename M::number_type> k1(IVP.size());
		if (fsal_ == false)
		{
			IVP.F(k1, state_1);
			fsal_ = true;
		}
		else
		{
			k1 = k_fsal_;
		}

		//k2's
		Vector<typename M::number_type> state_2;
		state_2.push_back(state[0] + stepsize_ / 5);
		for (size_type i = 0; i < IVP.size(); ++i)
		{
			state_2.push_back(state[i+1] + stepsize_ / 5 * k1[i]);
		}
		Vector<typename M::number_type> k2(IVP.size());
		IVP.F(k2, state_2);

		//k3's
		Vector<typename M::number_type> state_3;
		state_3.push_back(state[0] + 3*stepsize_ / 10);
		for (size_type i = 0; i < IVP.size(); ++i)
		{
			state_3.push_back(state[i+1] + 3*stepsize_ / 40 * k1[i]+9*stepsize_/40*k2[i]);
		}
		Vector<typename M::number_type> k3(IVP.size());
		IVP.F(k3, state_3);

		//k4's
		Vector<typename M::number_type> state_4;
		state_4.push_back(state[0] + 4*stepsize_/5);
		for (size_type i = 0; i < IVP.size(); ++i)
		{
			state_4.push_back(state[i+1] + 44*stepsize_/45 * k1[i] - 56*stepsize_/15*k2[i] + 32*stepsize_/9*k3[i]);
		}
		Vector<typename M::number_type> k4(IVP.size());
		IVP.F(k4, state_4);

		//k5's
		Vector<typename M::number_type> state_5;
		state_5.push_back(state[0] + 8*stepsize_/9);
		for (size_type i = 0; i < IVP.size(); ++i)
		{
			state_5.push_back(state[i+1] + 19372*stepsize_/6561 * k1[i] - 25360*stepsize_/2187*k2[i] + 64448*stepsize_/6561*k3[i] - 212*stepsize_/729*k4[i]);
		}
		Vector<typename M::number_type> k5(IVP.size());
		IVP.F(k5, state_5);

		//k6's
		Vector<typename M::number_type> state_6;
		state_6.push_back(state[0] + stepsize_);
		for (size_type i = 0; i < IVP.size(); ++i)
		{
			state_6.push_back(state[i+1] + 9017*stepsize_/3168 * k1[i] - 355*stepsize_/33*k2[i] + 46732*stepsize_/5247*k3[i] + 49*stepsize_/176*k4[i] - 5103*stepsize_/18656*k5[i]);
		}
		Vector<typename M::number_type> k6(IVP.size());
		IVP.F(k6, state_6);

		//k7's
		Vector<typename M::number_type> state_7;
		state_7.push_back(state[0] + stepsize_);
		for (size_type i = 0; i < IVP.size(); ++i)
		{
			state_7.push_back(state[i+1] + 35*stepsize_/384 * k1[i] + 500*stepsize_/1113*k3[i] + 125*stepsize_/192*k4[i] - 2187*stepsize_/6784*k5[i] + 11*stepsize_/84*k6[i]);
		}
		Vector<typename M::number_type> k7(IVP.size());
		IVP.F(k7, state_7);

		// save result according to FSAL
		k_fsal_ = k7;


		// fourth-order estimation:
		Vector<typename M::number_type> state_four(state);
		state_four[0] += stepsize_;
		for (size_type i = 0; i < IVP.size(); ++i)
		{
			state_four[i+1] += stepsize_*(5179*k1[i]/57600 + 7571*k3[i]/16695 + 393*k4[i]/640 -92097*k5[i]/339200 + 187*k6[i]/2100 + k7[i]/40);
		}

		// update state (using FSAL property)
		state[0] += stepsize_;
		for (size_type i = 0; i < IVP.size(); ++i)
		{
			state[i+1] = state_7[i+1];
		}

		// step size calculation
		Vector<typename M::number_type> sk(IVP.size());
		for (int i = 0; i < sk.size(); i++)
		{
			sk[i] = atol_ + fmax(abs(state_four[i+1]), abs(state[i+1]))*rtol_;
		}

		number_type err(0.);
		for (int i = 0; i < IVP.size(); i++)
		{
			err += (state_four[i] - state[i])*(state_four[i] - state[i])/sk[i]/sk[i];
		}
		err = err/15/IVP.size();
		err = sqrt(err);
		err = fmax(err, 1e-50);

		stepsize_ = stepsize_ * fmin(alpha_max_, fmax(alpha_min_, alpha_*pow(1./err, 1./5)));


	}

private:
	number_type stepsize_;
	number_type alpha_; // prevent frequent step repetitions
	number_type alpha_max_; // upper limit
	number_type alpha_min_; // lower limit
	number_type atol_; // absolute tolerance for error
	number_type rtol_; // relative tolerance for error
	bool fsal_; // has fsal already been taken into account?
	Vector<number_type> k_fsal_; // save for fsal


};




#endif