#ifndef SPLINE1D_HPP
#define SPLINE1D_HPP

#include "vector.hpp"
#include "matrix.hpp"
#include <stdexcept>
#include <string>


template<typename REAL>
class Spline1D
{
public:
	typedef std::size_t size_type;
	typedef REAL number_type;

	// Constructor
	Spline1D ()
	: initialized_(false)
	{};

	void initialize(Vector<number_type> x, Vector<number_type> y, std::string boundary_condition)
	{
		if (x.size() != y.size())
		{
			throw std::runtime_error("Input vectors must be of same length.");
		}

		x_ = x;
		y_ = y;
		N_ = x_.size()-1; // Number of splines

		// // add a grid point that is only taken into account for periodic boundary conditions
		// x_.push_back(x_[N_-1]+.2*(x_[N_-1]-x_[N_-2]));
		// y_.push_back(y_[0]);

		Vector<number_type> h(N_); // step sizes in x
		for (size_type i = 0; i < h.size(); ++i)
		{
			h[i] = x_[i+1] - x_[i];
		}
		h_ = h;

		// For periodic boundary conditions one temporarily discards the last data point as it is identical to the first one
		if (boundary_condition == "periodic")
		{
			N_ -= 1;
		}
		
		// Constants determined by boundary conditions
		number_type mu_i;
		number_type lam_i;
		number_type mu_f;
		number_type lam_f;
		number_type b_i;
		number_type b_f;

		if (boundary_condition == "natural")
		{
			mu_i = 1.;
			lam_i = 0.;
			mu_f = 1.;
			lam_f = 0.;
			b_i = 0.;
			b_f = 0.;

		}
		else if (boundary_condition == "periodic")
		{
			mu_i = (h_[0]+h_[N_])/3.;
			lam_i = h_[0]/6.;
			mu_f = (h_[N_-1]+h_[N_])/3.;
			lam_f = h_[N_-1]/6.;
			b_i = (y_[1]-y_[0])/h_[0] - (y_[0]-y_[N_])/h_[N_];
			b_f = (y_[0]-y_[N_])/h_[N_] - (y_[N_]-y_[N_-1])/h_[N_-1];
		}
		else
		{
			throw std::runtime_error("Undefined boundary conditions.");
		}

		// Determine spline by solving linear system of equations

		Vector<number_type> M(N_+1); // Moments, i.e. second derivatives at the grid points ->unknown

		Matrix<number_type> LHS(M.size(), M.size(), 0);

		for (size_type row = 0; row < LHS.rowsize(); ++row)
		{
			//boundary conditions
			if (row == 0)
			{
				LHS(row, 0) = mu_i;
				LHS(row, 1) = lam_i;
				if (boundary_condition == "periodic")
				{
					LHS(row, N_) = h_[N_]/6;
				}
			}
			else if (row == N_)
			{
				LHS(row, N_-1) = lam_f;
				LHS(row, N_) = mu_f;
				if (boundary_condition == "periodic")
				{
					LHS(row, 0) = h_[N_]/6;
				}
			}
			else // Main part
			{
				LHS(row, row-1) = h_[row-1]/6.;
				LHS(row, row) = (h_[row-1]+h_[row])/3.;
				LHS(row, row+1) = h_[row]/6.;
			}

		}

		Vector<number_type> RHS(M.size(), 0);
		for (size_type row = 0; row < RHS.size(); ++row)
		{
			//boundary conditions
			if (row == 0)
			{
				RHS[row] = b_i;
			}
			else if (row == N_)
			{
				RHS[row] = b_f;
			}
			else // Main part
			{
				RHS[row] = (y_[row+1]-y[row])/h_[row] - (y_[row]-y[row-1])/h_[row-1]; 
			}
		}

		// scaling factors
		Vector<number_type> s(M.size());
		// row permutations
		Vector<size_type> perm_r(M.size());
		// column permutations
		Vector<size_type> perm_c(M.size());
		bool non_invertible = false; // boolean for whether coefficient matrix is non-invertible


		row_equilibrate(LHS, s, non_invertible); // equilibrate rows
		
		if (non_invertible)
		{
			throw std::runtime_error("Singular coefficient matrix.");
		}

		lu_fullpivot(LHS, perm_r, perm_c); // LU decomposition of A
		M = number_type(0.0); // clear solution
		apply_equilibrate(s, RHS); // equilibration of right-hand side
		permute_forward(perm_r, RHS); // permutation of right-hand side
		solveL(LHS, RHS, RHS); // forward substitution
		solveU(LHS, M, RHS); // backward substitution
		permute_backward(perm_c, M); // backward permutation

		M_ = M;

		// Add the last spline again in case of periodic boundary conditions
		if (boundary_condition == "periodic")
		{
			N_ += 1;
			M_.push_back(M_[0]);
		}

		Vector<number_type> c(h_.size());
		Vector<number_type> d(h_.size());
		for (size_type i = 0; i < c.size(); ++i)
		{
			d[i] = y_[i] - M_[i]*h_[i]*h_[i]/6.;
			c[i] = (y_[i+1]-y_[i])/h_[i]-h_[i]/6.*(M_[i+1]-M_[i]);
		}

		c_ = c;
		d_ = d;


		initialized_ = true;
	}

	number_type evaluate(number_type x) const
	{
		if (!initialized_)
		{
			throw std::runtime_error("Initialize Spline1D object first.");
		}
		if ((x <x_[0]) || (x>x_[x_.size()-1]))
		{
			throw std::runtime_error("Outside interpolation domain.");
		}

		// Step 1:  Determine correct spline
		size_type index = N_-1; 
		for (size_type i = 1; i < x_.size(); ++i)
		{
			if (x < x_[i])
			{
				index = i-1;
				break;
			}
		}

		return 1./6*((x_[index+1]-x)*(x_[index+1]-x)*(x_[index+1]-x)/h_[index]*M_[index]+(x-x_[index])*(x-x_[index])*(x-x_[index])/h_[index]*M_[index+1]) + c_[index]*(x-x_[index]) +d_[index];
	}

	void print_data()
	{
		if (!initialized_)
		{
			throw std::runtime_error("Initialize Spline1D object first.");
		}
		for (size_type i = 0; i < N_; ++i)
		{
			std::cout << x_[i] << " " << y_[i] << "\n";
		}
	}


private:
	// interpolation data
	Vector<number_type> x_;
	Vector<number_type> y_;
	size_type N_;

	// objects for computing the spline
	bool initialized_; 
	Vector<number_type> h_;
	Vector<number_type> M_;
	Vector<number_type> c_;
	Vector<number_type> d_;

};


#endif