#ifndef MODEL_HPP
#define MODEL_HPP
#include <cmath>
#include "../auxiliary_files/vector.hpp"
#include <stdexcept>

#include "array3d.hpp"


// Define here the Lattice action of the theory
template<typename REAL>
class LatticeAction
{
public:
	typedef std::size_t size_type;
	typedef REAL number_type;
	LatticeAction(
	size_type N,
	size_type d,
	number_type L,
	number_type m,
	number_type lambda
	) // default constructor
	: N_(N)
	, d_(d)
	, L_(L)
	, m_(m)
	, lambda_(lambda)
	{}

	/* number of lattice grid points */
	size_type n_sites() const
	{
		return size_type(pow(N_, d_));
	}

	/* number of space-time dimensions */
	size_type n_dim() const
	{
		return d_;
	}

	/* number of lattice grid points per dimension */
	size_type n_sites_per_dim() const
	{
		return N_;
	}

	

public:
	// /* define here lattice action (harmonic oscillator) */
	// number_type action(const Vector<number_type> & q)
	// {
	// 	/* checks if data arrays have equal lengths */
	// 	if (q.size() != n_parameters())
	// 	{
	// 		throw std::runtime_error("Array fed to action has incorrect size.\n");
	// 	}

	// 	number_type a = L_/N_;
	// 	number_type result = 0;
	// 	for (size_type i = 0; i < N_; ++i)
	// 	{
	// 		result += (m_/8./a*(q[index(i+2)]-q[i])*(q[index(i+2)]-q[i]) + a*V(q[index(i+1)]));
	// 	}
		

	// 	return result;
	// }

	/* define here lattice action (real scalar) */
	number_type action(const Vector<number_type> & q)
	{

		/* checks if data arrays have equal lengths */
		if (q.size() != n_sites())
		{
			throw std::runtime_error("Array fed to action has incorrect size.\n");
		}

		Array3d<number_type> grid(q);

		number_type a = L_/N_;
		number_type result = 0;
		// sum over lattice sites
		for (int t = 0; t < n_sites_per_dim(); ++t)
		{
			for (int x = 0; x < n_sites_per_dim(); ++x)
			{
				for (int y = 0; y < n_sites_per_dim(); ++y)
				{
					
					number_type phi = grid(index(t), index(x), index(y));
					// compute kinetic term of the Lagrangian
					number_type kinetic = 0;
					for (size_type mu = 0; mu < n_dim(); ++mu)
					{
						number_type diffp1;
					
						if (mu == 0)
						{
							diffp1 = grid(index(t+1),index(x),index(y))-phi;
						}
						else if (mu == 1)
						{
							diffp1 = grid(index(t),index(x+1),index(y))-phi;
						}
						else if (mu == 2)
						{
							diffp1 = grid(index(t),index(x),index(y+1))-phi;
						}

						kinetic += 1./2*(diffp1/a)*(diffp1/a);	
					}
					// compute mass term
					number_type massive = (m_*m_/2)*phi*phi;

					// compute interaction term
					number_type interact = interaction(phi);

					// compute Lagrangian
					result += a*a*a*a*(kinetic + massive + interact);
					
				}
			}
		}
		

		return result;
	}

	// help function: interaction term
	number_type interaction (number_type phi)
	{
		return lambda_/4./3./2.*phi*phi*phi*phi;
	}

	// method to keep track of correct index with periodic boundary conditions
	size_type index(int i)
	{
		int result = i;
		if (result >= int(N_))
		{
			result -= int(N_);
		}
		else if (result < 0)
		{
			result += int(N_);
		}
		return size_type(result);
	}

	// Define (lattice) operator, of which the ensemble average is to be taken
	number_type Operator (const Vector<number_type> & q, size_type t)
	{
		if (t >= n_sites_per_dim())
		{
			throw std::runtime_error("time index overflow\n");
		}
		Array3d<number_type> grid(q);

		// Compute field position average for all time slices

		// Vector<number_type> Gammas(n_sites_per_dim());
		// for (size_type i = 0; i < Gammas.size(); ++i)
		// {
		// 	Gammas[i] = grid.average_at(i);
		// }
		
		// number_type result = 0;
		// for (int i = 0; i < n_sites_per_dim(); ++i)
		// {
		// 	result += Gammas[i]*Gammas[index(i+t)];
		// }

		// return result/n_sites_per_dim();

		return grid.average_at(t);
	}

private:
	number_type m_; // mass of the particle
	number_type lambda_; // coupling constant

	size_type N_; // number of lattice sites per dimension
	size_type d_; // number of dimensions
	number_type L_; // lattice size per direction
};


#endif