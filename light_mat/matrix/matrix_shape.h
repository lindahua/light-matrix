/**
 * @file matrix_shape.h
 *
 * The class that represents the shape of a matrix
 *
 * @author Dahua Lin
 */

#ifdef _MSC_VER
#pragma once
#endif

#ifndef LIGHTMAT_MATRIX_SHAPE_H_
#define LIGHTMAT_MATRIX_SHAPE_H_

#include <light_mat/common/basic_defs.h>

#if LMAT_DIAGNOSIS_LEVEL >= 1
#define LMAT_CHECK_DIM_VALIDITY(ct_dim, d) dim_checker<ct_dim>::check(d);
#else
#define LMAT_CHECK_DIM_VALIDITY(ct_dim, d)
#endif

namespace lmat
{
	/********************************************
	 *
	 *  dim checker
	 *
	 ********************************************/

	template<int N>
	struct dim_checker
	{
		LMAT_ENSURE_INLINE
		static void check(index_t d)
		{
			check_arg(d == N, "Input dimension is incompatible with the compile-time spec.");
		}
	};

	template<>
	struct dim_checker<0>
	{
		LMAT_ENSURE_INLINE
		static void check(index_t d) { }
	};



	/********************************************
	 *
	 *  single dimension
	 *
	 ********************************************/

	template<int N>
	class dimension
	{
	public:
		LMAT_ENSURE_INLINE
		dimension() { }

		LMAT_ENSURE_INLINE
		dimension(const index_t n)
		{
			LMAT_CHECK_DIM_VALIDITY(N, n)
		}

		LMAT_ENSURE_INLINE
		index_t value() const { return N; }
	};


	template<>
	class dimension<0>
	{
	public:
		LMAT_ENSURE_INLINE
		dimension() : m_dim(0) { }

		LMAT_ENSURE_INLINE
		dimension(const index_t n) : m_dim(n)
		{ }

		LMAT_ENSURE_INLINE
		index_t value() const { return m_dim; }

	private:
		index_t m_dim;
	};


	/********************************************
	 *
	 *  matrix_shape
	 *
	 ********************************************/

	template<int M, int N>
	class matrix_shape
	{
	public:
		static const int ct_nrows = M;
		static const int ct_ncols = N;

		LMAT_ENSURE_INLINE matrix_shape() { }

		LMAT_ENSURE_INLINE matrix_shape(index_t m, index_t n)
		{
			LMAT_CHECK_DIM_VALIDITY(M, m)
			LMAT_CHECK_DIM_VALIDITY(N, n)
		}

		LMAT_ENSURE_INLINE index_t nrows() const
		{
			return M;
		}

		LMAT_ENSURE_INLINE index_t ncolumns() const
		{
			return N;
		}

		LMAT_ENSURE_INLINE index_t nelems() const
		{
			return M * N;
		}
	};

	template<int M>
	class matrix_shape<M, 0>
	{
	public:
		static const int ct_nrows = M;
		static const int ct_ncols = 0;

		LMAT_ENSURE_INLINE matrix_shape() : m_ncols(0) { }

		LMAT_ENSURE_INLINE matrix_shape(index_t m, index_t n)
		: m_ncols(n)
		{
			LMAT_CHECK_DIM_VALIDITY(M, m)
		}

		LMAT_ENSURE_INLINE index_t nrows() const
		{
			return M;
		}

		LMAT_ENSURE_INLINE index_t ncolumns() const
		{
			return m_ncols;
		}

		LMAT_ENSURE_INLINE index_t nelems() const
		{
			return M * m_ncols;
		}

	private:
		index_t m_ncols;
	};


	template<int N>
	class matrix_shape<0, N>
	{
	public:
		static const int ct_nrows = 0;
		static const int ct_ncols = N;

		LMAT_ENSURE_INLINE matrix_shape() : m_nrows(0) { }

		LMAT_ENSURE_INLINE matrix_shape(index_t m, index_t n)
		: m_nrows(m)
		{
			LMAT_CHECK_DIM_VALIDITY(N, n)
		}

		LMAT_ENSURE_INLINE index_t nrows() const
		{
			return m_nrows;
		}

		LMAT_ENSURE_INLINE index_t ncolumns() const
		{
			return N;
		}

		LMAT_ENSURE_INLINE index_t nelems() const
		{
			return m_nrows * N;
		}

	private:
		index_t m_nrows;
	};


	template<>
	class matrix_shape<0, 0>
	{
	public:
		static const int ct_nrows = 0;
		static const int ct_ncols = 0;

		LMAT_ENSURE_INLINE matrix_shape() : m_nrows(0), m_ncols(0) { }

		LMAT_ENSURE_INLINE matrix_shape(index_t m, index_t n) : m_nrows(m), m_ncols(n) { }

		LMAT_ENSURE_INLINE index_t nrows() const
		{
			return m_nrows;
		}

		LMAT_ENSURE_INLINE index_t ncolumns() const
		{
			return m_ncols;
		}

		LMAT_ENSURE_INLINE index_t nelems() const
		{
			return m_nrows * m_ncols;
		}

	private:
		index_t m_nrows;
		index_t m_ncols;
	};

}

#endif /* MATRIX_SHAPE_H_ */
