#ifndef VECTOR_HPP
#define VECTOR_HPP

#include <vector>
#include <assert.h>
#include <random>
#include <fstream>
#include <cmath>

template<typename REAL>
class Vector : public std::vector<REAL> // inherit from the STL vector
{
public:
	typedef std::size_t size_type; // type used for array indices
	Vector() : std::vector<REAL>() // default constructor, inherited
	{
	}

	/* constructor to create an array of given size and filled by a value */
	Vector( const size_type size, const REAL defaultvalue = 0) 
	: std::vector<REAL>(size, defaultvalue)
	{

	}

	/* fill vector with one value */
	Vector & operator= (const REAL value)
	{
		const size_type s = this->size();
		Vector & self = *this;
		for (size_type i = 0; i<s; ++i)
		{
			self[i] = value;
		}
		return *this;
	}

	/*!
    Returns a new vector that is a subset of the components
    of the given vector.

    i first index of the new vector, m size of the new vector
    */
    Vector sub (size_type i, size_type m)
    {
    	Vector v(m);
    	Vector &self = *this;
    	size_type k=0;
    	for (size_type j=i; j<i+m; j++)
    	{
    		v[k]=self[j];
    		k++;
    	}
    	return v;
    }

	/* Multiplication by a scalar */
	Vector & operator*= (const REAL value)
	{
		Vector &self = *this;
		for (size_type i = 0; i < this->size(); ++i)
		{
			self[i] *= value;
		}
		return *this;
	}

	/* Component-wise multiplication with a vector */
	Vector & operator*= (const Vector & y)
	{
		assert(this->size() == y.size());
		Vector &self = *this;
		for (size_type i = 0; i < this->size(); ++i)
		{
			self[i] *= y[i];
		}
		return *this;
	}

	/* Division by a scalar */
	Vector & operator/= (const REAL value)
	{
		Vector &self = *this;
		for (size_type i = 0; i < this->size(); ++i)
		{
			self[i] /= value;
		}
		return *this;
	}

	/* Add another vector */
	Vector & operator+= (const Vector & y)
	{
		assert(this->size() == y.size());
		Vector &self = *this;
		for (size_type i = 0; i < this->size(); ++i)
		{
			self[i] += y[i];
		}
		return *this;
	}
	/* Subtract another vector */
	Vector & operator-= (const Vector & y)
	{
		assert(this->size() == y.size());
		Vector &self = *this;
		for (size_type i = 0; i < this->size(); ++i)
		{
			self[i] -= y[i];
		}
		return *this;
	}

	/* Adding two vectors y+x */
	Vector operator+ (const Vector & x) const
	{
		assert(x.size() == this->size()); // checks if dimensions of the two vectors match
		Vector sum(*this);
		sum += x;
		return sum;
	}

	/* Adding a constant y+alpha (component-wise) */
	Vector operator+ (const REAL & x) const
	{
		Vector sum(*this);
		for (size_type i = 0; i < sum.size(); ++i)
		{
			sum[i] += x;
		}
		return sum;
	}

	/* Subtraction of two vectors y-x */
	Vector operator- (const Vector & x) const
	{
		assert(x.size() == this->size()); // checks if dimensions of the two vectors match
		Vector sum(*this);
		sum -= x;
		return sum;
	}

	/* Negation of a vector */
	Vector & operator-()
	{
		Vector & neg(*this);
		neg *= -1.0;
		return *this;
	}


	/* Square of the Euclidean norm */
	REAL two_norm_2 () const
	{
		REAL sum(0);
		const Vector & self = *this;
		for (size_type i = 0; i < (size_type) this->size(); ++i)
		{
			sum += self[i] * self[i];
		}
		return sum;
	}

	/* Mean of the entries */
	REAL mean() const
	{
		REAL mean = 0;
		for (size_type i = 0; i< this->size(); ++i)
		{
			mean += (*this)[i];
		}
		mean /= this->size();
		return mean;
	}

	/* (Unbiased) standard deviation of the entries */
	REAL std() const
	{
		REAL result = 0;
		REAL mean = this->mean();
		for (size_type i = 0; i< this->size(); ++i)
		{
			result += ((*this)[i] - mean)*((*this)[i] - mean);
		}
		result /= (this->size()-1.);
		return sqrt(result);
	}

	// Save vector entries in text file
	void save_as (std::string filename) const
	{
		std::ofstream outfile(filename);
		for (size_type i = 0; i < this->size(); ++i)
		{
			outfile << (*this)[i] << "\n";
		}
	}

};

// some free functions

/* Multiplication from the left hand side by scalar */
	template<typename REAL>
	Vector<REAL> operator* (const REAL & alpha, const Vector<REAL> x)
	{
		Vector<REAL> result(x);
		result *= alpha;
		return result;
	}
/* Multiplication from the right hand side by scalar */
	template<typename REAL>
	Vector<REAL> operator* (const Vector<REAL> & x, const REAL & alpha)
	{
		Vector<REAL> result(x);
		result *= alpha;
		return result;
	}

/* Component-wise multiplication of two vectors, resulting in a vector (!), no inner product ! */
	template<typename REAL>
	Vector<REAL> operator* (const Vector<REAL> & x, const Vector<REAL> & y)
	{
		Vector<REAL> result(x);
		for (std::size_t i = 0; i < result.size(); ++i)
		{
			result[i] *= y[i];
		}
		return result;
	}

/* Inner product of two vectors */
	template<typename REAL>
	REAL inner_product(const Vector<REAL> & x, const Vector<REAL> & y)
	{
		assert(x.size() == y.size());
		REAL result(0);
		for (std::size_t i = 0; i < x.size(); ++i)
		{
			result += x[i]*y[i];
		}
		return result;
	}

	/* Minkowski product of two vectors (signature +---) */
	template<typename REAL>
	REAL minkowski_product(const Vector<REAL> & x, const Vector<REAL> & y)
	{
		assert(x.size() == y.size());
		REAL result(x[0]*y[0]);
		for (std::size_t i = 1; i < x.size(); ++i)
		{
			result -= x[i]*y[i];
		}
		return result;
	}

/* Division from the right hand side by scalar */
	template<typename REAL1, typename REAL2>
	Vector<REAL1> operator/ (const Vector<REAL1> & x, const REAL2 & alpha)
	{
		Vector<REAL1> result(x);
		for (std::size_t i = 0; i < result.size(); ++i)
		{
			result[i] = x[i]/alpha;
		}
		return result;
	}

/* Component-wise division of two vectors */
	template<typename REAL1, typename REAL2>
	Vector<REAL1> operator/ (const Vector<REAL1> & x, const Vector<REAL2> & y)
	{
		Vector<REAL1> result(x);
		for (std::size_t i = 0; i < result.size(); ++i)
		{
			result[i] = x[i]/y[i];
		}
		return result;
	}

/* Component-wise square root */
	template<typename REAL>
	void sqrt(Vector<REAL> & x)
	{
		for (std::size_t i = 0; i < x.size(); ++i)
		{
			x[i] = sqrt(x[i]);
		}
	}

/* Component-wise absolute value */
	template<typename REAL>
	void abs(Vector<REAL> & x)
	{
		for (std::size_t i = 0; i < x.size(); ++i)
		{
			x[i] = abs(x[i]);
		}
	}

/* Fill vector with random entries following a probability distribution */
	template<typename REAL, typename Distribution>
	void fill_random(Vector<REAL> & vec, Distribution distribution)
	{
		std::random_device rd;
		std::mt19937 gen(rd());
		for (std::size_t i = 0; i < vec.size(); ++i)
		{
			vec[i] = distribution(gen);
		}

	}

	/* Fill vector with gaussian random entries of mean vector und uncertainty vector */
	template<typename REAL>
	void fill_gaussian(Vector<REAL> & vec, const Vector<REAL> & mu, const Vector<REAL> & sigma)
	{
		std::random_device rd;
		std::mt19937 gen(rd());
		for (std::size_t i = 0; i < vec.size(); ++i)
		{
			std::normal_distribution<> normal(mu[i], sigma[i]);
			vec[i] = normal(gen);
		}

	}

/* Fill vector with random entries from a specific region */
	template<typename REAL>
	void fill_from_region(Vector<REAL> &vec, const Vector<REAL> & region_min, const Vector<REAL> & region_max)
	{
		assert(vec.size() == region_min.size() && vec.size() == region_max.size());
		std::random_device rd;
		std::mt19937 gen(rd());
		for (std::size_t i = 0; i < vec.size(); ++i)
		{
			std::uniform_real_distribution<> dis(region_min[i], region_max[i]);
			//std::normal_distribution<> dis(0.5*(region_min[i]+region_max[i]), (region_max[i]-region_min[i])/3.);
			vec[i] = dis(gen);
		}

	}


#endif