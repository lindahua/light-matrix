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

namespace lmat
{

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
			check_arg(n == N,
					"The input dimension is invalid.");
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
		LMAT_ENSURE_INLINE matrix_shape() { }

		LMAT_ENSURE_INLINE matrix_shape(const index_t m, const index_t n)
		: m_dim0(m), m_dim1(n) { }

		LMAT_ENSURE_INLINE index_t nrows() const
		{
			return m_dim0.value();
		}

		LMAT_ENSURE_INLINE index_t ncolumns() const
		{
			return m_dim1.value();
		}

		LMAT_ENSURE_INLINE index_t nelems() const
		{
			return m_dim0.value() * m_dim1.value();
		}

	private:
		dimension<M> m_dim0;
		dimension<N> m_dim1;
	};


	/********************************************
	 *
	 *  indexing
	 *
	 ********************************************/

	struct column_major_layout { };

	template<int M, int N>
	LMAT_ENSURE_INLINE
	inline index_t sub2offset(column_major_layout,
			const matrix_shape<M, N>& shape, index_t i, index_t j)
	{
		return i + shape.nrows() * j;
	}

	template<int M>
	LMAT_ENSURE_INLINE
	inline index_t sub2offset(column_major_layout,
			const matrix_shape<M, 1>& shape, index_t i, index_t j)
	{
		return i;
	}

	template<int N>
	LMAT_ENSURE_INLINE
	inline index_t sub2offset(column_major_layout,
			const matrix_shape<1, N>& shape, index_t i, index_t j)
	{
		return shape.nrows() * j;
	}

	LMAT_ENSURE_INLINE
	inline index_t sub2offset(column_major_layout,
			const matrix_shape<1, 1>& shape, index_t i, index_t j)
	{
		return 0;
	}


	struct column_major_layout_ex
	{
		const index_t lead_dim;

		LMAT_ENSURE_INLINE
		column_major_layout_ex(index_t ld) : lead_dim(ld) { }
	};

	template<int M, int N>
	LMAT_ENSURE_INLINE
	inline index_t sub2offset(const column_major_layout_ex& layout,
			const matrix_shape<M, N>& shape, index_t i, index_t j)
	{
		return i + layout.lead_dim * j;
	}

	template<int M>
	LMAT_ENSURE_INLINE
	inline index_t sub2offset(const column_major_layout_ex& layout,
			const matrix_shape<M, 1>& shape, index_t i, index_t j)
	{
		return i;
	}

	template<int N>
	LMAT_ENSURE_INLINE
	inline index_t sub2offset(const column_major_layout_ex& layout,
			const matrix_shape<1, N>& shape, index_t i, index_t j)
	{
		return layout.lead_dim * j;
	}

	LMAT_ENSURE_INLINE
	inline index_t sub2offset(const column_major_layout_ex& layout,
			const matrix_shape<1, 1>& shape, index_t i, index_t j)
	{
		return 0;
	}

}

#endif /* MATRIX_SHAPE_H_ */
