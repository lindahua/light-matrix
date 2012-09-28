/*
 * @file dense_matrix_internal.h
 *
 * Internal implementation for dense_matrix
 *
 * @author Dahua Lin
 */

#ifndef _MSC_VER
#pragma once
#endif

#ifndef LIGHTMAT_DENSE_MATRIX_INTERNAL_H_
#define LIGHTMAT_DENSE_MATRIX_INTERNAL_H_

#include <light_mat/matrix/matrix_base.h>
#include <light_mat/common/block.h>

namespace lmat { namespace detail {

	template<typename T, int CTRows, int CTCols> class dense_matrix_internal;

	template<typename T, int CTRows, int CTCols>
	class dense_matrix_internal
	{
	public:
		LMAT_ENSURE_INLINE
		dense_matrix_internal() : m_blk() { }

		LMAT_ENSURE_INLINE
		dense_matrix_internal(index_t m, index_t n)
		: m_blk()
		{
			check_arg(m == CTRows && n == CTCols,
					"Attempted to construct a static dense_matrix with incorrect dimensions.");
		}

		LMAT_ENSURE_INLINE void swap(dense_matrix_internal& s)
		{
			m_blk.swap(s.m_blk);
		}

		LMAT_ENSURE_INLINE index_t nelems() const { return CTRows * CTCols; }

		LMAT_ENSURE_INLINE index_t nrows() const { return CTRows; }

		LMAT_ENSURE_INLINE index_t ncolumns() const { return CTCols; }

		LMAT_ENSURE_INLINE index_t offset(const index_t i, const index_t j) const
		{
			return CTRows * j + i;
		}

		LMAT_ENSURE_INLINE const T *ptr_data() const { return m_blk.ptr_data(); }

		LMAT_ENSURE_INLINE T *ptr_data() { return m_blk.ptr_data(); }

		LMAT_ENSURE_INLINE const T *ptr_col(const index_t j) const
		{
			return ptr_data() + j * CTRows;
		}

		LMAT_ENSURE_INLINE T *ptr_col(const index_t j)
		{
			return ptr_data() + j * CTRows;
		}

		LMAT_ENSURE_INLINE void resize(index_t m, index_t n)
		{
			check_arg(m == CTRows && n == CTCols,
					"Attempted to change a fixed dimension in resize.");
		}

	private:
		LMAT_ALIGN(LMAT_DEFAULT_ALIGNMENT)
		sblock<T, CTRows * CTCols> m_blk;
	};


	template<typename T, int CTRows>
	class dense_matrix_internal<T, CTRows, 0>
	{
	public:
		LMAT_ENSURE_INLINE
		dense_matrix_internal() : m_blk(), m_ncols(0) { }

		LMAT_ENSURE_INLINE
		dense_matrix_internal(index_t m, index_t n)
		: m_blk(n * check_forward(m, m == CTRows,
				"Invalid construction of dense_matrix (m != CTRows)") )
		, m_ncols(n) { }

		inline void swap(dense_matrix_internal& s)
		{
			m_blk.swap(s.m_blk);
			std::swap(m_ncols, s.m_ncols);
		}

		LMAT_ENSURE_INLINE index_t nelems() const { return CTRows * m_ncols; }

		LMAT_ENSURE_INLINE index_t nrows() const { return CTRows; }

		LMAT_ENSURE_INLINE index_t ncolumns() const { return m_ncols; }

		LMAT_ENSURE_INLINE index_t offset(const index_t i, const index_t j) const
		{
			return CTRows * j + i;
		}

		LMAT_ENSURE_INLINE const T *ptr_data() const { return m_blk.ptr_data(); }

		LMAT_ENSURE_INLINE T *ptr_data() { return m_blk.ptr_data(); }

		LMAT_ENSURE_INLINE const T *ptr_col(const index_t j) const
		{
			return ptr_data() + j * CTRows;
		}

		LMAT_ENSURE_INLINE T *ptr_col(const index_t j)
		{
			return ptr_data() + j * CTRows;
		}

		LMAT_ENSURE_INLINE void resize(index_t m, index_t n)
		{
			check_arg(m == CTRows,
					"Attempted to change a fixed dimension in resize.");

			if (n != m_ncols)
			{
				m_blk.resize(m * n);
				m_ncols = n;
			}
		}

	private:
		dblock<T, aligned_allocator<T> > m_blk;
		index_t m_ncols;
	};

	template<typename T, int CTCols>
	class dense_matrix_internal<T, 0, CTCols>
	{
	public:
		LMAT_ENSURE_INLINE
		dense_matrix_internal() : m_blk(), m_nrows(0) { }

		LMAT_ENSURE_INLINE
		dense_matrix_internal(index_t m, index_t n)
		: m_blk(m * check_forward(n, n == CTCols,
				"Invalid construction of dense_matrix (n != CTCols)") )
		, m_nrows(m) { }

		LMAT_ENSURE_INLINE void swap(dense_matrix_internal& s)
		{
			m_blk.swap(s.m_blk);
			std::swap(m_nrows, s.m_nrows);
		}

		LMAT_ENSURE_INLINE index_t nelems() const { return m_nrows * CTCols; }

		LMAT_ENSURE_INLINE index_t nrows() const { return m_nrows; }

		LMAT_ENSURE_INLINE index_t ncolumns() const { return CTCols; }

		LMAT_ENSURE_INLINE index_t offset(const index_t i, const index_t j) const
		{
			return m_nrows * j + i;
		}

		LMAT_ENSURE_INLINE const T *ptr_data() const { return m_blk.ptr_data(); }

		LMAT_ENSURE_INLINE T *ptr_data() { return m_blk.ptr_data(); }

		LMAT_ENSURE_INLINE const T *ptr_col(const index_t j) const
		{
			return ptr_data() + j * m_nrows;
		}

		LMAT_ENSURE_INLINE T *ptr_col(const index_t j)
		{
			return ptr_data() + j * m_nrows;
		}

		LMAT_ENSURE_INLINE void resize(index_t m, index_t n)
		{
			check_arg(n == CTCols,
					"Attempted to change a fixed dimension in resize.");

			if (m != m_nrows)
			{
				m_blk.resize(m * n);
				m_nrows = m;
			}
		}

	private:
		dblock<T, aligned_allocator<T> > m_blk;
		index_t m_nrows;
	};


	template<typename T>
	class dense_matrix_internal<T, 0, 0>
	{
	public:
		LMAT_ENSURE_INLINE
		dense_matrix_internal() : m_blk(), m_nrows(0), m_ncols(0) { }

		LMAT_ENSURE_INLINE
		dense_matrix_internal(index_t m, index_t n)
		: m_blk(m * n), m_nrows(m), m_ncols(n) { }

		LMAT_ENSURE_INLINE void swap(dense_matrix_internal& s)
		{
			m_blk.swap(s.m_blk);
			std::swap(m_nrows, s.m_nrows);
			std::swap(m_ncols, s.m_ncols);
		}

		LMAT_ENSURE_INLINE index_t nelems() const { return m_nrows * m_ncols; }

		LMAT_ENSURE_INLINE index_t nrows() const { return m_nrows; }

		LMAT_ENSURE_INLINE index_t ncolumns() const { return m_ncols; }

		LMAT_ENSURE_INLINE index_t offset(const index_t i, const index_t j) const
		{
			return m_nrows * j + i;
		}

		LMAT_ENSURE_INLINE const T *ptr_data() const { return m_blk.ptr_data(); }

		LMAT_ENSURE_INLINE T *ptr_data() { return m_blk.ptr_data(); }

		LMAT_ENSURE_INLINE const T *ptr_col(const index_t j) const
		{
			return ptr_data() + j * m_nrows;
		}

		LMAT_ENSURE_INLINE T *ptr_col(const index_t j)
		{
			return ptr_data() + j * m_nrows;
		}

		LMAT_ENSURE_INLINE void resize(index_t m, index_t n)
		{
			if (nelems() != m * n)
			{
				m_blk.resize(m * n);
			}

			m_nrows = m;
			m_ncols = n;
		}

	private:
		dblock<T, aligned_allocator<T> > m_blk;
		index_t m_nrows;
		index_t m_ncols;
	};


} }

#endif 
