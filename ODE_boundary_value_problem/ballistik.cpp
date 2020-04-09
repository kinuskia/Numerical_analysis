#include <iostream>
#include "bvp.hpp"
#include "../auxiliary_files/rk4.hpp"
#include "../auxiliary_files/dopri5.hpp"
#include "../auxiliary_files/multiple_shooting.hpp"
#include <vector>



int main ()
{
	typedef double Number;
	BVP<Number> Kugel(0.45, 4);
	RK4<Number> solver_RK4;
	DOPRI5<Number> solver_DOPRI5;


	// Initial guess for parameters of grid
	Vector<Number> s0(5);
	s0[0] = 0;
	s0[1] = 0;
	s0[2] = .355;
	s0[3] = 0;
	s0[4] = 1;
	Vector<Number> s1(5);
	s1[0] = 2.25;
	s1[1] = 1.1;
	s1[2] = .7;
	s1[3] = 1.7;
	s1[4] = 0.4;
	Vector<Number> s2(5);
	s2[0] = 3.5;
	s2[1] = 2.1;
	s2[2] = .9;
	s2[3] = 1.8;
	s2[4] = -.4;
	Vector<Number> s3(5);
	s3[0] = 5;
	s3[1] = 3.8;
	s3[2] = 1.5;
	s3[3] = 0;
	s3[4] = -2.37;

	std::vector<Vector<Number>> grid(4);
	grid[0] = s0;
	grid[1] = s1;
	grid[2] = s2;
	grid[3] = s3;



	mshooter<Number> mshooter(grid,1e-28,0);

	Vector<Number> opt(5);
	
	mshooter.optimal_initials(Kugel, opt);
	//opt=s0;

	// calculate evolution with optimal initials:
	solver_RK4.set_stepsize(0.01);
	while(opt[0] < s3[0])
	{
		std::cout << opt[0] << "\t" << opt[1] << "\t" << opt[2] << "\t" << opt[3] << "\t" << opt[4] << "\n";
		//std::cout << opt[1] << "\t" << opt[0] << "\t" << opt[3] << "\n";
		solver_RK4.step_forward(Kugel, opt);
	}

	
	


	return 0;
}