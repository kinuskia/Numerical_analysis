#include <iostream>
#include <cmath>

#include "ivp.hpp"
#include "../auxiliary_files/rk4.hpp"
#include "../auxiliary_files/dopri5.hpp"
#include "../auxiliary_files/vector.hpp"


using namespace std;




int main ()
{
	typedef std::size_t size_type;
	typedef double number_type;
	
	// Initialisation of the problem: harmonic oscillator
	IVP<number_type> harmonic;

	// Set ODE solver (Runge-Kutta-4 or Dormand-Prince-4-5)
	//RK4<number_type> solver;
	DOPRI5<number_type> solver; // adaptive step size method

	// Set step size
	number_type h = 0.001;
	solver.set_stepsize(h);

	// initial conditions
	Vector<number_type> x(3);
	x[0] = 0; // time
	x[1] = 0; // position
	x[2] = 0.1; // velocity

	/
	for (size_type i = 0; i < 100; ++i)
	{
		std::cout << x[0] << " " << x[1] << " " << x[2] << "\n";
		solver.step_forward(harmonic, x);
	}


	


	return 0;
}