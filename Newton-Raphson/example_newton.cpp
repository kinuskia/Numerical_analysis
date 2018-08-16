#include <iostream>

#include "model.hpp"
#include "../auxiliary_files/newton.hpp"

int main ()
{
	typedef std::size_t size_type;
	typedef double number_type;

	Vector<number_type> q(2);
	q[0] = 2.1;
	q[1] = 0.8;

	// initialize root problem
	Model<number_type> problem;

	// initialize Newton solver
	Newton<number_type> solver(problem);



	return 0;
}