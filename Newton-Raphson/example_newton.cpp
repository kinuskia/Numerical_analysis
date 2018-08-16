#include <iostream>

#include "model.hpp"
#include "../auxiliary_files/newton.hpp"


int main ()
{
	typedef std::size_t size_type;
	typedef double number_type;

	Vector<number_type> q(2);
	Vector<number_type> lower_bounds(q.size());
	Vector<number_type> upper_bounds(q.size());
	lower_bounds[0] = 0;
	upper_bounds[0] = 4;
	lower_bounds[1] = -1;
	upper_bounds[1] = 3;

	// initialize root problem
	Model<number_type> problem;

	// initialize Newton solver
	Newton<number_type> solver(problem);

	for (size_type j = 0; j < 10 ; ++j)
	{
		fill_from_region(q, lower_bounds, upper_bounds);
		solver.solve(q);
		std::cout << "Convergence? (1/0) " << solver.converged() << "\n";
		for (size_type i = 0; i < q.size(); ++i)
		{
			std::cout << q[i] << "\n";
		}
	}




	return 0;
}