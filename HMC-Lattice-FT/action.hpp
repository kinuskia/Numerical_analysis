#ifndef MODEL_HPP
#define MODEL_HPP
#include <cmath>
#include "../auxiliary_files/vector.hpp"
#include <stdexcept>

#include "array4d.hpp"


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
	number_type T,
	number_type m
	) // default constructor
	: N_(N)
	, d_(d)
	, T_(T)
	, m_(m)
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

	// 	number_type a = T_/N_;
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

		Array4d<number_type> grid(q);

		number_type a = T_/N_;
		number_type result = 0;
		// sum over lattice sites
		for (int t = 0; t < n_sites_per_dim(); ++t)
		{
			for (int x = 0; x < n_sites_per_dim(); ++x)
			{
				for (int y = 0; y < n_sites_per_dim(); ++y)
				{
					for (int z = 0; z < n_sites_per_dim(); ++z)
					{
						number_type phi = grid(index(t), index(x), index(y), index(z));
						// get kinetic term of the Lagrangian
						number_type kinetic = 0;
						for (size_type mu = 0; mu < n_dim(); ++mu)
						{
							number_type diffp1;
							number_type diffp2;
							number_type diffm1;
							number_type diffm2;
							if (mu == 0)
							{
								diffp2 = grid(index(t+2),index(x),index(y),index(z))-phi;
								diffp1 = grid(index(t+1),index(x),index(y),index(z))-phi;
								diffm1 = grid(index(t-1),index(x),index(y),index(z))-phi;
								diffm2 = grid(index(t-2),index(x),index(y),index(z))-phi;
							}
							else if (mu == 1)
							{
								diffp2 = grid(index(t),index(x+2),index(y),index(z))-phi;
								diffp1 = grid(index(t),index(x+1),index(y),index(z))-phi;
								diffm1 = grid(index(t),index(x-1),index(y),index(z))-phi;
								diffm2 = grid(index(t),index(x-2),index(y),index(z))-phi;
							}
							else if (mu == 2)
							{
								diffp2 = grid(index(t),index(x),index(y+2),index(z))-phi;
								diffp1 = grid(index(t),index(x),index(y+1),index(z))-phi;
								diffm1 = grid(index(t),index(x),index(y-1),index(z))-phi;
								diffm2 = grid(index(t),index(x),index(y-2),index(z))-phi;
							}
							else if (mu == 3)
							{
								diffp2 = grid(index(t),index(x),index(y),index(z+2))-phi;
								diffp1 = grid(index(t),index(x),index(y),index(z+1))-phi;
								diffm1 = grid(index(t),index(x),index(y),index(z-1))-phi;
								diffm2 = grid(index(t),index(x),index(y),index(z-2))-phi;
							}
							kinetic += 1./2*diffp1*diffp1/a/a;	
						}
						// get mass term
						number_type massive = (m_*m_/2)*phi*phi;

						// get interaction term
						number_type interact = interaction(phi);

						// compute Lagrangian
						result += a*a*a*a*(kinetic + massive + interact);
					}
				}
			}
		}
		

		return result;
	}

	// help function: interaction term
	number_type interaction (number_type phi)
	{
		return 0;
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
		Array4d<number_type> grid(q);

		Vector<number_type> Gammas(n_sites_per_dim());
		for (size_type i = 0; i < Gammas.size(); ++i)
		{
			Gammas[i] = grid.average_at(i);
		}
		
		number_type result = 0;
		for (int i = 0; i < n_sites_per_dim(); ++i)
		{
			result += Gammas[i]*Gammas[index(i+t)];
		}

		return result/n_sites_per_dim();
	}

private:
	number_type m_; // (bare) mass of the particle

	size_type N_; // number of lattice sites per dimension
	size_type d_; // number of lattice sites
	number_type T_; // lattice time size
};


#endif