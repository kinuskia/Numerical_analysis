#ifndef ARRAY4D_HPP
#define ARRAY4D_HPP

#include <vector>
#include <assert.h>
#include <random>
#include <fstream>
#include <cmath>
#include <stdexcept>
#include "../auxiliary_files/vector.hpp"

template<typename REAL>
class Array3d : public std::vector<REAL> // inherit from the STL vector
{
public:
	typedef std::size_t size_type; // type used for array indices
	typedef REAL number_type;

	/* constructor to create an array of given size and filled by a value */
	Array3d( const size_type size, const REAL defaultvalue = 0) 
	: data_(size*size*size, defaultvalue)
	, n_per_dim_(size)
	{

	}

	/* constructor to create an array from a Vector */
	Array3d(Vector<number_type> q) 
	: data_(q)
	, n_per_dim_(round(pow(q.size(),1./3)))
	{

	}

	// /* fill vector with one value */
	// Array4d & operator= (const REAL value)
	// {
	// 	const size_type s = data_.size();
	// 	Array4d & self = *this;
	// 	for (size_type i = 0; i<s; ++i)
	// 	{
	// 		self[i] = value;
	// 	}
	// 	return *this;
	// }

	// Access operators
	// A(i, j)
	number_type & operator()(const size_type t, const size_type x, const size_type y)
	{
		assert( t < n_per_dim_ && x < n_per_dim_ && y < n_per_dim_ );
		return data_[get_index(t, x, y)];
	}
	const number_type & operator()(const size_type t, const size_type x, const size_type y) const
	{
		assert( t < n_per_dim_ && x < n_per_dim_ && y < n_per_dim_);
		return data_[get_index(t, x, y)];
	}



	/* Getter */
	size_type get_n_per_dim() const
	{
		return n_per_dim_;
	}

	// fill with increasing indices
	void fill()
	{
		for (size_type i = 0; i < data_.size(); ++i)
		{
			data_[i] = number_type(i);
		}
	}

	/* Convert vector index to space-time triple */
	void get_position(const size_type & index, size_type & t, size_type & x, size_type & y)
	{
		size_type i = index;
		t = i % n_per_dim_;
		i = i / n_per_dim_;
		x = i % n_per_dim_;
		i = i / n_per_dim_;
		y = i % n_per_dim_;
		i = i / n_per_dim_;
	}

	/* Convert space-time triple to vector index */
	size_type get_index(const size_type & t, const size_type & x, const size_type & y)
	{
		return t + x*n_per_dim_ + y*n_per_dim_*n_per_dim_;
	}


	/* Get spatial average value of lattice at time index t */
	number_type average_at(size_type t)
	{
		if (t >= n_per_dim_)
			throw std::runtime_error("time index overflow\n");
		number_type Gamma_t = 0;
		
		for (size_type x = 0; x < n_per_dim_; ++x)
		{
			for (size_type y = 0; y < n_per_dim_; ++y)
			{
				Gamma_t += (*this)(t, x, y);	
			}
		}

		return Gamma_t/n_per_dim_/n_per_dim_;
	}



	/*!
    // Returns a new vector that is a subset of the components
    // of the given vector.

    // i first index of the new vector, m size of the new vector
    // */
    // Vector sub (size_type i, size_type m)
    // {
    // 	Vector v(m);
    // 	Vector &self = *this;
    // 	size_type k=0;
    // 	for (size_type j=i; j<i+m; j++)
    // 	{
    // 		v[k]=self[j];
    // 		k++;
    // 	}
    // 	return v;
    // }

	// /* Multiplication by a scalar */
	// Array4d & operator*= (const REAL value)
	// {
	// 	Array4d &self = *this;
	// 	for (size_type i = 0; i < this->size(); ++i)
	// 	{
	// 		self[i] *= value;
	// 	}
	// 	return *this;
	// }

	// /* Component-wise multiplication with a vector */
	// Vector & operator*= (const Vector & y)
	// {
	// 	assert(this->size() == y.size());
	// 	Vector &self = *this;
	// 	for (size_type i = 0; i < this->size(); ++i)
	// 	{
	// 		self[i] *= y[i];
	// 	}
	// 	return *this;
	// }

	// /* Division by a scalar */
	// Vector & operator/= (const REAL value)
	// {
	// 	Vector &self = *this;
	// 	for (size_type i = 0; i < this->size(); ++i)
	// 	{
	// 		self[i] /= value;
	// 	}
	// 	return *this;
	// }

	// /* Add another vector */
	// Array4d & operator+= (const Array4d & y)
	// {
	// 	assert(this->size() == y.size());
	// 	Array4d &self = *this;
	// 	for (size_type i = 0; i < this->size(); ++i)
	// 	{
	// 		self[i] += y[i];
	// 	}
	// 	return *this;
	// }
	// /* Subtract another vector */
	// Array4d & operator-= (const Array4d & y)
	// {
	// 	assert(this->size() == y.size());
	// 	Array4d &self = *this;
	// 	for (size_type i = 0; i < this->size(); ++i)
	// 	{
	// 		self[i] -= y[i];
	// 	}
	// 	return *this;
	// }

	// /* Adding two vectors y+x */
	// Array4d operator+ (const Array4d & x) const
	// {
	// 	assert(x.size() == this->size()); // checks if dimensions of the two vectors match
	// 	Array4d sum(*this);
	// 	sum += x;
	// 	return sum;
	// }

	// /* Adding a constant y+alpha (component-wise) */
	// Array4d operator+ (const REAL & x) const
	// {
	// 	Array4d sum(*this);
	// 	for (size_type i = 0; i < sum.size(); ++i)
	// 	{
	// 		sum[i] += x;
	// 	}
	// 	return sum;
	// }

	// /* Subtraction of two vectors y-x */
	// Array4d operator- (const Array4d & x) const
	// {
	// 	assert(x.size() == this->size()); // checks if dimensions of the two vectors match
	// 	Array4d sum(*this);
	// 	sum -= x;
	// 	return sum;
	// }

	// /* Negation of a vector */
	// Array4d & operator-()
	// {
	// 	Array4d & neg(*this);
	// 	neg *= -1.0;
	// 	return *this;
	// }


	// /* Square of the Euclidean norm */
	// REAL two_norm_2 () const
	// {
	// 	REAL sum(0);
	// 	const Vector & self = *this;
	// 	for (size_type i = 0; i < (size_type) this->size(); ++i)
	// 	{
	// 		sum += self[i] * self[i];
	// 	}
	// 	return sum;
	// }

	// /* Mean of the entries */
	// REAL mean() const
	// {
	// 	REAL mean = 0;
	// 	for (size_type i = 0; i< this->size(); ++i)
	// 	{
	// 		mean += (*this)[i];
	// 	}
	// 	mean /= this->size();
	// 	return mean;
	// }

	// /* (Unbiased) standard deviation of the entries */
	// REAL std() const
	// {
	// 	REAL result = 0;
	// 	REAL mean = this->mean();
	// 	for (size_type i = 0; i< this->size(); ++i)
	// 	{
	// 		result += ((*this)[i] - mean)*((*this)[i] - mean);
	// 	}
	// 	result /= (this->size()-1.);
	// 	return sqrt(result);
	// }

	// // Save vector entries in text file
	// void save_as (std::string filename) const
	// {
	// 	std::ofstream outfile(filename);
	// 	for (size_type i = 0; i < this->size(); ++i)
	// 	{
	// 		outfile << (*this)[i] << "\n";
	// 	}
	// }

private:
	std::vector<number_type> data_; // Array data is stored in an STL vector
	size_type n_per_dim_;


};

// some free functions

// /* Multiplication from the left hand side by scalar */
// 	template<typename REAL>
// 	Array4d<REAL> operator* (const REAL & alpha, const Array4d<REAL> x)
// 	{
// 		Array4d<REAL> result(x);
// 		result *= alpha;
// 		return result;
// 	}
//  Multiplication from the right hand side by scalar 
// 	template<typename REAL>
// 	Array4d<REAL> operator* (const Array4d<REAL> & x, const REAL & alpha)
// 	{
// 		Array4d<REAL> result(x);
// 		result *= alpha;
// 		return result;
// 	}

// /* Component-wise multiplication of two vectors, resulting in a vector (!), no inner product ! */
// 	template<typename REAL>
// 	Vector<REAL> operator* (const Vector<REAL> & x, const Vector<REAL> & y)
// 	{
// 		Vector<REAL> result(x);
// 		for (std::size_t i = 0; i < result.size(); ++i)
// 		{
// 			result[i] *= y[i];
// 		}
// 		return result;
// 	}

// /* Inner product of two vectors */
// 	template<typename REAL>
// 	REAL inner_product(const Vector<REAL> & x, const Vector<REAL> & y)
// 	{
// 		assert(x.size() == y.size());
// 		REAL result(0);
// 		for (std::size_t i = 0; i < x.size(); ++i)
// 		{
// 			result += x[i]*y[i];
// 		}
// 		return result;
// 	}

// 	 Minkowski product of two vectors (signature +---) 
// 	template<typename REAL>
// 	REAL minkowski_product(const Vector<REAL> & x, const Vector<REAL> & y)
// 	{
// 		assert(x.size() == y.size());
// 		REAL result(x[0]*y[0]);
// 		for (std::size_t i = 1; i < x.size(); ++i)
// 		{
// 			result -= x[i]*y[i];
// 		}
// 		return result;
// 	}

// /* Division from the right hand side by scalar */
// 	template<typename REAL1, typename REAL2>
// 	Array4d<REAL1> operator/ (const Array4d<REAL1> & x, const REAL2 & alpha)
// 	{
// 		Array4d<REAL1> result(x);
// 		for (std::size_t i = 0; i < result.size(); ++i)
// 		{
// 			result[i] = x[i]/alpha;
// 		}
// 		return result;
// 	}

// /* Component-wise division of two vectors */
// 	template<typename REAL1, typename REAL2>
// 	Vector<REAL1> operator/ (const Vector<REAL1> & x, const Vector<REAL2> & y)
// 	{
// 		Vector<REAL1> result(x);
// 		for (std::size_t i = 0; i < result.size(); ++i)
// 		{
// 			result[i] = x[i]/y[i];
// 		}
// 		return result;
// 	}

// /* Component-wise square root */
// 	template<typename REAL>
// 	void sqrt(Vector<REAL> & x)
// 	{
// 		for (std::size_t i = 0; i < x.size(); ++i)
// 		{
// 			x[i] = sqrt(x[i]);
// 		}
// 	}

// /* Component-wise absolute value */
// 	template<typename REAL>
// 	void abs(Vector<REAL> & x)
// 	{
// 		for (std::size_t i = 0; i < x.size(); ++i)
// 		{
// 			x[i] = abs(x[i]);
// 		}
// 	}

// /* Fill vector with random entries following a probability distribution */
// 	template<typename REAL, typename Distribution>
// 	void fill_random(Vector<REAL> & vec, Distribution distribution)
// 	{
// 		std::random_device rd;
// 		std::mt19937 gen(rd());
// 		for (std::size_t i = 0; i < vec.size(); ++i)
// 		{
// 			vec[i] = distribution(gen);
// 		}

// 	}

// 	/* Fill vector with gaussian random entries of mean vector und uncertainty vector */
// 	template<typename REAL>
// 	void fill_gaussian(Vector<REAL> & vec, const Vector<REAL> & mu, const Vector<REAL> & sigma)
// 	{
// 		std::random_device rd;
// 		std::mt19937 gen(rd());
// 		for (std::size_t i = 0; i < vec.size(); ++i)
// 		{
// 			std::normal_distribution<> normal(mu[i], sigma[i]);
// 			vec[i] = normal(gen);
// 		}

// 	}

// /* Fill vector with random entries from a specific region */
// 	template<typename REAL>
// 	void fill_from_region(Vector<REAL> &vec, const Vector<REAL> & region_min, const Vector<REAL> & region_max)
// 	{
// 		assert(vec.size() == region_min.size() && vec.size() == region_max.size());
// 		std::random_device rd;
// 		std::mt19937 gen(rd());
// 		for (std::size_t i = 0; i < vec.size(); ++i)
// 		{
// 			std::uniform_real_distribution<> dis(region_min[i], region_max[i]);
// 			//std::normal_distribution<> dis(0.5*(region_min[i]+region_max[i]), (region_max[i]-region_min[i])/3.);
// 			vec[i] = dis(gen);
// 		}

// 	}


#endif