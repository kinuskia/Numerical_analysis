#ifndef MATRIX_HPP
#define MATRIX_HPP

#include "vector.hpp"
#include <assert.h>
#include <string>
#include <fstream>

template<typename REAL>
class Matrix
{
public:
	typedef std::size_t size_type;
	typedef REAL number_type;
	typedef typename std::vector<REAL> VType;
	typedef typename VType::const_iterator ConstVectorIterator;
	typedef typename VType::iterator VectorIterator;

private:
	VType m_data_; // Matrix data is stored in an STL vector
	size_type m_rows_; // Number of matrix rows
	size_type m_cols_; // Number of matrix columns


	/* Get matrix element for write access */
	inline number_type & at(const size_type row, const size_type col)
	{
		return m_data_[row * m_cols_ + col ];
	}
	/* get matrix element for read-only access */
	inline const number_type & at(const size_type row, const size_type col) const
	{
		return m_data_[row * m_cols_ + col];
	}

public:
	//  constructor for empty matrix
	Matrix()
	: m_data_(0, 0)
	, m_rows_(0)
	, m_cols_(0)
	{}

	// Constructor
	Matrix(const size_type rows, const size_type cols, const number_type def_val = 0)
	: m_data_(rows*cols, def_val)
	, m_rows_(rows)
	, m_cols_(cols)
	{}

	// useful getter variables
	size_type rowsize() const
	{
		return m_rows_;
	}

	size_type colsize() const
	{
		return m_cols_;
	}

	// Overloaded element access operators
	// A(i, j)
	inline number_type & operator()(const size_type row, const size_type col)
	{
		assert( row < m_rows_ && col<m_cols_);
		return at(row, col);
	}
	inline const number_type & operator()(const size_type row, const size_type col) const
	{
		assert( row < m_rows_ && col<m_cols_);
		return at(row, col);
	}

	// A[i][j]
	VectorIterator operator[] (const size_type row)
	{
		assert (row < m_rows_);
		return m_data_.begin() + row * m_cols_;
	} 
	const ConstVectorIterator operator[] (const size_type row) const
	{
		assert (row < m_rows_);
		return m_data_.begin() + row * m_cols_;
	} 

	// Assignment operator
	Matrix & operator= (const Matrix & A)
	{
		m_data_ = A.m_data_;
		m_rows_ = A.m_rows_;
		m_cols_ = A.m_cols_;
		return *this;
	}

	// Assignment operator for scalar
	Matrix & operator= (number_type s)
	{
		for (size_type i = 0; i < rowsize(); ++i)
		{
			for (size_type j = 0; j < colsize(); ++j)
			{
				(*this)(i, j) *= s;
			}
		}
		return *this;
	}

	// Addition assignment
	Matrix & operator+= (const Matrix & B)
	{
		for (size_type i = 0; i < rowsize(); ++i)
		{
			for (size_type j = 0; j < colsize(); ++j)
			{
				(*this)(i, j) += B(i, j);
			}
		}
		return *this;
	}

	// Subtraction assignment
	Matrix & operator-= (const Matrix & B)
	{
		for (size_type i = 0; i < rowsize(); ++i)
		{
			for (size_type j = 0; j < colsize(); ++j)
			{
				(*this)(i, j) -= B(i, j);
			}
		}
		return *this;
	}

	// Multiplication assignment
	Matrix & operator*= (const number_type s)
	{
		for (size_type i = 0; i < rowsize(); ++i)
		{
			for (size_type j = 0; j < colsize(); ++j)
			{
				(*this)(i, j) *= s;
			}
		}
		return *this;
	}

	// Division assignment
	Matrix & operator/= (const number_type s)
	{
		for (size_type i = 0; i < rowsize(); ++i)
		{
			for (size_type j = 0; j < colsize(); ++j)
			{
				(*this)(i, j) /= s;
			}
		}
		return *this;
	}

	// in-place Matrix-vector product: y = A*x
	template<class V>
	void mv (Vector<V> & y, const Vector<V> & x) const
	{
		assert(this->rowsize() == y.size()); // check dimensions of Matrix and vectors
		assert(this-> colsize() == x.size());

		for (size_type i = 0; i < rowsize(); ++i)
		{
			y[i] = 0;
			for (size_type j = 0; j < colsize(); ++j)
			{
				y[i] += (*this)(i, j)*x[j];
			}
		}
	}

	// in-place Matrix-matrix product C.mm(A, B)
	void mm (const Matrix<number_type> & A, const Matrix<number_type> & B)
	{
		// check dimensions
		assert(this->rowsize() == A.rowsize());
		assert(this->colsize() == B.colsize());
		assert(A.colsize() == B.rowsize());
		for (size_type i = 0; i < rowsize(); ++i)
		{
			for (size_type j = 0; j < colsize(); ++j)
			{
				(*this)(i, j) = 0;
				for (size_type k = 0; k < A.colsize(); ++k)
				{
					(*this)(i, j) += A(i, k) * B(k, j);
				}
			}
		}
	}

	// out-of place Matrix-Vector Multiplication y=A*x
	Vector<number_type> operator* (const Vector<number_type> &x)
	{
		assert(x.size() == rowsize());

		Vector<number_type> y (rowsize());
		for (size_type i = 0; i < rowsize(); ++i)
		{
			for (size_type j = 0; j < colsize(); ++j)
			{
				y[i] += at(i, j) * x[j];
			}
		}
		return y;
	}

	// out-of place Matrix-Matrix Multiplication C = (*this) * x
	Matrix operator* (const Matrix & x) const
	{
		assert(colsize() == x.rowsize());

		const size_type out_rows = rowsize();
		const size_type out_cols = x.colsize();
		Matrix y(out_rows, out_cols, 0.0);
		for (size_type i = 0; i < out_rows; ++i)
		{
			for (size_type j = 0; j < out_cols; ++j)
			{
				for (size_type k = 0; k < colsize(); ++k)
				{
					y(i, j) += at(i, k) * x(k, j);
				}
			}
		}

		return y;
	}

	// out-of-place Matrix-matrix addition
	Matrix operator+ (const Matrix & x) const
	{
		assert(colsize() == x.rowsize());

		const size_type out_rows = rowsize();
		const size_type out_cols = x.colsize();
		Matrix y(out_rows, out_cols, 0.0);
		y = *this;
		y += x;

		return y;
	}

	// out-of-place Matrix-matrix subtraction
	Matrix operator- (const Matrix & x) const
	{
		assert(colsize() == x.rowsize());

		const size_type out_rows = rowsize();
		const size_type out_cols = x.colsize();
		Matrix y(out_rows, out_cols, 0.0);
		y = *this;
		y -= x;

		return y;
	}

	// Make identity matrix
	template<class T>
	inline void identity (Matrix<T> & A)
	{
		for (typename Matrix<T>::size_type i = 0; i<A.rowsize(); ++i)
		{
			for (typename Matrix<T>::size_type j = 0; j<A.colsize(); ++j)
			{
				if (i == j)
				{
					A[i][j] = T(1);
				}
				else
				{
					A[i][j] = T(0);
				}
			}
		}
	}


};

// free functions for LU decomposition
typedef std::size_t size_type;


// LU decomposition of A with full pivoting
template<class T>
void lu_fullpivot (Matrix<T> & A, Vector<size_type> & p, Vector<size_type> & q)
{
	typedef T number_type;
	// check dimensions
	assert(A.rowsize() == A.colsize() && A.rowsize() != 0); // square and nonempty matrix
	assert(A.rowsize() == p.size());	// permutation vector compatible with matrix

	// initialize permutation
	for (size_type k = 0; k < A.rowsize(); ++k)
	{
		p[k] = q[k] = k;
	}

	// transformation to upper triangular
	for (size_type k = 0; k < A.rowsize()-1; ++k)
	{
		// find pivot element
		for (size_type r = k; r < A.rowsize(); ++r)
		{
			for (size_type s = k; s < A.colsize(); ++s)
			{
				if (abs(A[r][s]) > abs(A[k][k]))
				{
					p[k] = r;	// store permutation in step k
					q[k] = s;
				}
			}
		}

		if (p[k] > k) // exchange complete row if r !=k
		{
			for (size_type j = 0; j < A.colsize(); ++j)
			{
				number_type temp(A[k][j]);
				A[k][j] = A[p[k]][j];
				A[p[k]][j] = temp;
			}
		}
		if (q[k] > k) // exchange complete column if s != k
		{
			for (size_type i = 0; i < A.rowsize(); ++i)
			{
				number_type temp(A[i][k]);
				A[i][k] = A[i][q[k]];
				A[i][q[k]] = temp;
			}
		}

		assert(A[k][k] != 0); // Verify that A is not singular

		// modification
		for (size_type i = k+1; i < A.rowsize(); ++i)
		{
			number_type qik(A[i][k]/A[k][k]);
			A[i][k] = qik;
			for (size_type j = k+1; j < A.colsize(); ++j)
			{
				A[i][j] -= qik * A[k][j];
			}
		}
	}
}

// apply permutations to a right hand side vector
template <class T>
void permute_forward (const Vector<size_type> & p, Vector<T> & b)
{
	typedef T number_type;
	assert(b.size() == p.size()); // Check dimensions
	for (size_type k = 0; k < b.size()-1; ++k)
	{
		if (p[k] != k)
		{
			number_type temp(b[k]);
			b[k] = b[p[k]];
			b[p[k]] = temp;
		}
	}
}

// apply permutations to a solution vector
template <class T>
void permute_backward (const Vector<size_type> & q, Vector<T> & z)
{
	typedef T number_type;
	assert(z.size() == q.size());

	for (int k = z.size()-2; k >= 0; --k)
	{
		if (q[k] != size_type(k))
		{
			number_type temp(z[k]);
			z[k] = z[q[k]];
			z[q[k]] = temp;
		}
	}
}

// perform a row equilibration of a matrix; return scaling for later use.
// if matrix is non-invertible, it is saved in-place
template <class T>
void row_equilibrate (Matrix<T> & A, Vector<T> & s, bool & non_invertible)
{
	typedef T number_type;
	assert(A.rowsize() != 0 && A.colsize() != 0); // nonempty matrix
	assert(A.rowsize() == s.size()); // scaling vector compatible with matrix

	// equilibrate row sums
	for (size_type k = 0; k < A.rowsize(); ++k)
	{
		s[k] = number_type(0.0);
		for (size_type j = 0; j < A.colsize(); ++j)
		{
			s[k] += abs(A[k][j]);
		}
		if (s[k] == 0) // matrix non-invertible if row sum zero
		{
			non_invertible = true;
		}
		//assert(s[k] != 0); // nonzero row sum
		else
		{
			for (size_type j = 0; j < A.colsize(); ++j)
			{
				A[k][j] /= s[k];
			}
		}
	}
}

// Apply row equilibration to right hand side vector
template<class T>
void apply_equilibrate (Vector<T> & s, Vector<T> & b)
{
	typedef T number_type;

	assert(s.size() == b.size());

	// equilibrate row sums
	for (size_type k = 0; k < b.size(); ++k)
	{
		b[k] /= s[k];
	}
}

// Assume L = lower triangle of A with L_ii = 1, solve Lx = b
template<class T>
void solveL (const Matrix<T> & A, Vector<T> & x, const Vector<T> & b)
{
	typedef T number_type;
	assert(A.rowsize() == A.colsize() && A.rowsize() != 0); // square and nonempty matrix
	assert(A.rowsize() == b.size()); // right hand side compatible with matrix

	for (size_type i = 0; i < A.rowsize(); ++i)
	{
		number_type rhs(b[i]);
		for (size_type j = 0; j < i; ++j)
		{
			rhs -= A[i][j] * x[j];
		}
		x[i] = rhs;
	}
}

// Assume U = upper triangle of A and solve U x = b
template<class T>
void solveU (const Matrix<T> & A, Vector<T> & x, const Vector<T> & b)
{
	typedef T number_type;
	assert(A.rowsize() == A.colsize() && A.rowsize() != 0); // square and nonempty matrix
	assert(A.rowsize() == b.size()); // right hand side compatible with matrix

	for (int i = A.rowsize() - 1; i >=0; --i)
	{
		number_type rhs(b[i]);
		for (size_type j = i+1; j < A.colsize(); ++j)
		{
			rhs -= A[i][j] * x[j];
		}
		x[i] = rhs/A[i][i];
	}
}

/* Use LU decomposition to calculate to inverse of a Matrix */

template<class T>
void MatrixInverse (Matrix<T> & A, Matrix<T> & inverse)
{
	typedef T number_type;
	assert(A.colsize() == inverse.colsize());
	assert(A.rowsize() == inverse.rowsize());

	// Row equilibration of A
	Vector<number_type> s(A.colsize());
	row_equilibrate(A, s);

	// LU decomposition with full pivot
	Vector<size_type> p(A.colsize());
	Vector<size_type> q(A.colsize());
	lu_fullpivot(A, p, q);

	// Compute each column of the inverse matrix separately (each column defines a linear system)
	for (size_type j = 0; j < inverse.colsize(); ++j)
	{
		Vector<number_type> x(A.rowsize()); // solution vector
		Vector<number_type> b(A.rowsize(), 0); // jth column of identity matrix
		b[j] = 1.0;

		x = number_type(0.0);
		apply_equilibrate(s, b);
		permute_forward(p, b);
		solveL(A, b, b);
		solveU(A, x, b);
		permute_backward(q, x);

		for (size_type i = 0; i < inverse.rowsize(); ++i)
		{
			inverse[i][j] = x[i];
		}
	}

	

}

// Cholesky decomposition: H is overwritten with entries of G (and G^T)

template<class T>
void cholesky (Matrix<T> & H)
{
	typedef T number_type;
	for (int i = 0; i < H.rowsize(); ++i) 
	{
		for (int j = 0; j <= i; ++j) // using int instead of std::size_t because index becomes negative
		{
			number_type sum = H(i, j);
			for (int k = 0; k <= j-1; ++k)
			{
				sum -= H(i, k) * H(j, k);
			}
			if ( i > j)
			{
				H(i, j) = sum/H(j, j);
				H(j, i) = sum/H(j, j);
			}
			else 
			{
				assert(sum > 0);
				H(i, i) = sqrt(sum);
			}
		}
	}
}

/* Use Cholesky decomposition to calculate to inverse of a Matrix */

template<class T>
void MatrixInverse_cholesky (Matrix<T> & A, Matrix<T> & inverse)
{
	typedef T number_type;
	assert(A.colsize() == inverse.colsize());
	assert(A.rowsize() == inverse.rowsize());
	inverse = number_type(0.0);
	cholesky(A);

	// Compute each column of the inverse matrix separately (each column defines a linear system)
	for (size_type j = 0; j < inverse.colsize(); ++j)
	{
		Vector<number_type> x(A.rowsize()); // solution vector
		Vector<number_type> b(A.rowsize(), 0); // jth column of identity matrix
		b[j] = 1.0;

		x = number_type(0.0);
		// solve G y = b
		Vector<number_type> y = x;
		for (size_type i = 0; i < A.rowsize(); ++i)
		{
			number_type rhs(b[i]);
			for (size_type k = 0; k < i; ++k)
			{
				rhs -= A[i][k] * y[k];
			}
			y[i] = rhs / A[i][i];
		}
		// solve G^T x = y
		for (int i = A.rowsize() - 1; i >=0; --i)
		{
			number_type rhs(y[i]);
			for (size_type k = i+1; k < A.colsize(); ++k)
			{
				rhs -= A[i][k] * x[k];
			}
			x[i] = rhs/A[i][i];
		}

		// save result to output
		for (size_type i = 0; i < inverse.rowsize(); ++i)
		{
			inverse[i][j] = x[i];
		}


	}
	

	


}



#endif

