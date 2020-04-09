#ifndef MULTIPLE_SHOOTING_HPP
#define MULTIPLE_SHOOTING_HPP

#include "rk4.hpp"
#include "dopri5.hpp"
#include <vector>
#include "matrix.hpp"
#include <stdexcept>


// Solve A' = B(A) per Euler method, where A and B are matrices

template<class N>
class Matrix_Euler
{
public:
	typedef std::size_t size_type;
	typedef N number_type;
	// constructors
	Matrix_Euler ()
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
	void step_forward(const M & IVP, Matrix<typename M::number_type> & state, const Vector<typename M::number_type> & parameters)
	{
		//if (state.size() != IVP.size() + 1)
		//{
		//	HDNUM_ERROR("Incorrect input of initial conditions.");
		//}


		// Euler step

		Matrix<N> diff_state(IVP.size(), IVP.size());
		IVP.F(diff_state, state, parameters);

		// update state
		
		state.update(stepsize_, diff_state);

	}

private:
	number_type stepsize_;

};

// Helper class to define an ODE in terms of a matrix

template<class N>
class Matrix_ODE
{
public:
	typedef N number_type;
	typedef std::size_t size_type;
	Matrix_ODE(BVP<N> problem, size_type dim)
	: problem_(problem)
	, dim_(dim)
	{};

	// return derivative of the problem
	void F(Matrix<number_type> & deriv, const Matrix<number_type> & state, const Vector<number_type> parameters) const
	{
		//if (deriv.size() != state.size()-1 || deriv.size() != dim_)
		//	HDNUM_ERROR("State and derivative vector do not match the size of the ODE.");
		Matrix<number_type> Fx(problem_.size(), problem_.size());
		problem_.dF(Fx, parameters);
		deriv = Fx*state;
	}


	size_type size() const
	{
		return dim_;
	}

private:
	BVP<N> problem_;
	size_type dim_;

};

// multiple shooting method

template<class N>
class mshooter
{
public:
	typedef N number_type;
	typedef std::size_t size_type;
	mshooter(std::vector<Vector<N>> grid, number_type tolerance, size_type verbosity=0)
	: grid_(grid)
	, residual_(1)
	, tolerance_(tolerance)
	, verbosity_(verbosity)
	{};

	// return boundary condition and Jacobians with respect to x and y
	void r(Vector<N> & bound) const
	{
		bound[0] = grid_[0][1]; //t(0)=0
		bound[1] = grid_[0][3]; //y(0)=0
		bound[2] = grid_[0][4]-1; //alpha=45°
		bound[3] = grid_[3][3]; //y(1)=0
	}

	void rx(Matrix<N> & boundx) const
	{
		boundx[0][0] = 1;
		boundx[0][1] = 0;
		boundx[0][2] = 0;
		boundx[0][3] = 0;
		boundx[1][0] = 0;
		boundx[1][1] = 0;
		boundx[1][2] = 1;
		boundx[1][3] = 0;
		boundx[2][0] = 0;
		boundx[2][1] = 0;
		boundx[2][2] = 0;
		boundx[2][3] = 1;
		boundx[3][0] = 0;
		boundx[3][1] = 0;
		boundx[3][2] = 0;
		boundx[3][3] = 0;
	}

	void ry(Matrix<N> & boundy) const
	{
		boundy[0][0] = 0;
		boundy[0][1] = 0;
		boundy[0][2] = 0;
		boundy[0][3] = 0;
		boundy[1][0] = 0;
		boundy[1][1] = 0;
		boundy[1][2] = 0;
		boundy[1][3] = 0;
		boundy[2][0] = 0;
		boundy[2][1] = 0;
		boundy[2][2] = 0;
		boundy[2][3] = 0;
		boundy[3][0] = 0;
		boundy[3][1] = 0;
		boundy[3][2] = 1;
		boundy[3][3] = 0;
	}

	// Initial value problem solver
	void initial_solver(BVP<N> problem, RK4<N> IVP_solver,const Vector<N> initial, Vector<N> & output, Matrix_Euler<N> Matrix_solver, Matrix_ODE<N> problem_Matrix, Matrix<N> & output_Matrix, number_type abort)
	{
		output = initial;
		output_Matrix = Matrix<N>(problem.size(), problem.size());
		identity(output_Matrix);
		while(output[0] < abort)
		{
			number_type step_size = IVP_solver.stepsize();
			Matrix_solver.set_stepsize(step_size);
			Matrix_solver.step_forward(problem_Matrix, output_Matrix, output);
			IVP_solver.step_forward(problem, output);
		}
	}

	// do one iteration step to improve grid parameters
	void iterate(BVP<N> problem, std::vector<Vector<N>> & grid_evolved, std::vector<Matrix<N>> & Matrix_evolved, Matrix<N> & A, Matrix<N> & B)
	{
		if (grid_evolved.size() != grid_.size()-1)
		{
			throw std::runtime_error("Incorrect grid dimensions.");
		}
		RK4<N> IVP_solver;
		Matrix_Euler<N> Matrix_solver;
		for (size_type i = 0; i < grid_evolved.size(); i++)
		{
			Matrix_ODE<N> Matrix_problem(problem, problem.size());

			initial_solver(problem, IVP_solver, grid_[i], grid_evolved[i], Matrix_solver, Matrix_problem, Matrix_evolved[i], grid_[i+1][0]); 
		}
		// Calculation of A and B
		rx(A);
		ry(B);

		// // Calculation of change of grid parameters
	
		// d_s0: Solve linear system of equation L d_s0 = W
		Matrix<N> L(problem.size(), problem.size());
		Vector<N> W(problem.size());
		identity(L);
		for (int i = 0; i<Matrix_evolved.size(); i++)
		{
			L =  Matrix_evolved[i]*L;
		}
		L = A + B * L;

		std::vector<Vector<N>> Fs(grid_.size());
		for (int i = 0; i < grid_.size()-1; i++)
		{
			Fs[i] = grid_evolved[i] - grid_[i+1];
			Fs[i] = Fs[i].sub(1, problem.size());
		}
		size_type R = grid_.size()-1;
		Fs[R] = Fs[0];
		r(Fs[R]);

		W = Fs[R];
		
		W += B * Fs[R-1];


		for (int i = 0; i < R-1; i++)
		{
			Matrix<N> H(problem.size(), problem.size());
			identity(H);
			for (int j = R - 1 - i; j < R; j++)
			{
				H = Matrix_evolved[j]*H;
			}
			W += B * H * Fs[R-2-i];
		}
		W *= -1;


		bool non_invertible;
		std::vector<Vector<N>> grid_change(grid_.size());
		grid_change[0] = Vector<N>(problem.size());
		Vector<N> z(problem.size());              // solution of linear system
      	Vector<N> s(problem.size());              	  // scaling factors
      	Vector<size_type> p(problem.size());            // row permutations
      	Vector<size_type> q(problem.size());            // column permutations
      	row_equilibrate(L,s, non_invertible);            // equilibrate rows
	  	lu_fullpivot(L,p,q);                          // LR decomposition of L
	  	z = N(0.0);                                   // clear solution
	 	apply_equilibrate(s,W);                       // equilibration of right hand side
	  	permute_forward(p,W);                         // permutation of right hand side
	  	solveL(L,W,W);                                // forward substitution
	  	solveU(L,z,W);                                // backward substitution
	  	permute_backward(q,z);                        // backward permutation
	  	grid_change[0] = z;
	  	// calculate other ds_i
	  	for (int i = 0; i < R; i++)
	  	{
	  		grid_change[i+1] = Matrix_evolved[i]*grid_change[i] + Fs[i]; 
	  	}

	  	// calculate residual
	  	residual_ = 0;
	  	for (int i = 0; i < R+1; i++)
	  	{
	  		residual_ += grid_change[i].two_norm_2();
	  	}


	  	// calculate new grid parameters
	  	for (int i = 0; i < grid_.size(); i++)
	  	{
	  		for (int j = 0; j < problem.size(); j++)
	  		{
	  			grid_[i][j+1] += grid_change[i][j];
	  		}
	  	}


	}

	// calculate optimal initial conditions

	void optimal_initials(BVP<N> problem, Vector<N> & initials)
	{
		std::vector<Vector<N>> grid_evolved(grid_.size()-1);
		std::vector<Matrix<N>> Matrix_evolved(grid_.size()-1);
		Matrix<N> A(problem.size(), problem.size());
		Matrix<N> B(problem.size(), problem.size());
		int counter = 0;
		while (counter < 50)
		{
			if (verbosity_ > 0)
			{
				std::cout << residual_ << "\n";
			}
			iterate(problem, grid_evolved, Matrix_evolved, A, B);
			counter++;
			if (residual_ < tolerance_)
			{
				break;
			}
		}
		if (counter == 50)
		{
			throw std::runtime_error("No convergence.");
		}
		grid0(initials);

	}

	// getter function for grid
	std::vector<Vector<N>> grid() const
	{
		return grid_;
	}

	// getter for grid[0]
	void grid0(Vector<N> & v) const
	{
		v = grid_[0];
	}



private:
	std::vector<Vector<N>> grid_;
	number_type residual_;
	number_type tolerance_;
	size_type verbosity_;

};











#endif