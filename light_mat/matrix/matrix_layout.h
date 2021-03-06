/**
 * @file matrix_layout.h
 *
 * @brief Layout and indexing
 *
 * @author Dahua Lin
 */

#ifdef _MSC_VER
#pragma once
#endif

#ifndef LIGHTMAT_MATRIX_LAYOUT_H_
#define LIGHTMAT_MATRIX_LAYOUT_H_

#include "internal/matrix_layout_internal.h"

namespace lmat
{

	// forward declaration

	template<index_t M, index_t N> class cont_layout_cm;
	template<index_t M, index_t N> class block_layout_cm;
	template<index_t M, index_t N> class grid_layout;

	/********************************************
	 *
	 *  Layout traits
	 *
	 ********************************************/

	template<class Layout> struct layout_traits;

	template<index_t M, index_t N>
	struct layout_traits<cont_layout_cm<M, N> >
	{
		static const index_t ct_num_rows = M;
		static const index_t ct_num_cols = N;

		static const bool ct_is_contiguous = true;
		static const bool ct_is_percol_contiguous = true;

		typedef matrix_shape<M, N> shape_type;
	};

	template<index_t M, index_t N>
	struct layout_traits<block_layout_cm<M, N> >
	{
		static const index_t ct_num_rows = M;
		static const index_t ct_num_cols = N;

		static const bool ct_is_contiguous = (N == 1);
		static const bool ct_is_percol_contiguous = true;

		typedef matrix_shape<M, N> shape_type;
	};

	template<index_t M, index_t N>
	struct layout_traits<grid_layout<M, N> >
	{
		static const index_t ct_num_rows = M;
		static const index_t ct_num_cols = N;

		static const bool ct_is_contiguous = (M == 1 && N == 1);
		static const bool ct_is_percol_contiguous = M == 1;

		typedef matrix_shape<M, N> shape_type;
	};


	/********************************************
	 *
	 *  general interface
	 *
	 ********************************************/

	template<class Layout> struct layout_traits;

	template<class Derived>
	class IMatrixLayout
	{
	public:
		LMAT_CRTP_REF

		LMAT_ENSURE_INLINE
		index_t nrows() const
		{
			return derived().nrows();
		}

		LMAT_ENSURE_INLINE
		index_t ncolumns() const
		{
			return derived().ncolumns();
		}

		LMAT_ENSURE_INLINE
		index_t nelems() const
		{
			return derived().nelems();
		}

		LMAT_ENSURE_INLINE
		typename layout_traits<Derived>::shape_type
		shape() const
		{
			return derived().shape();
		}

		LMAT_ENSURE_INLINE
		index_t row_stride() const
		{
			return derived().row_stride();
		}

		LMAT_ENSURE_INLINE
		index_t col_stride() const
		{
			return derived().col_stride();
		}

		LMAT_ENSURE_INLINE
		bool is_contiguous() const
		{
			return derived().is_contiguous();
		}

		LMAT_ENSURE_INLINE
		bool is_percol_contiguous() const
		{
			return derived().is_percol_contiguous();
		}

		LMAT_ENSURE_INLINE
		index_t offset(index_t i, index_t j) const
		{
			return derived().offset(i, j);
		}

		LMAT_ENSURE_INLINE
		index_t col_offset(index_t j) const
		{
			return derived().col_offset(j);
		}

		LMAT_ENSURE_INLINE
		index_t row_offset(index_t i) const
		{
			return derived().row_offset(i);
		}

	};





	/********************************************
	 *
	 *  specific layout classes
	 *
	 ********************************************/

	template<index_t M, index_t N>
	class cont_layout_cm : public IMatrixLayout<cont_layout_cm<M, N> >
	{
	public:
		LMAT_ENSURE_INLINE
		cont_layout_cm() : m_shape() { }

		LMAT_ENSURE_INLINE
		cont_layout_cm(index_t m, index_t n) : m_shape(m, n) { };

	public:
		LMAT_ENSURE_INLINE
		index_t nrows() const
		{
			return m_shape.nrows();
		}

		LMAT_ENSURE_INLINE
		index_t ncolumns() const
		{
			return m_shape.ncolumns();
		}

		LMAT_ENSURE_INLINE
		index_t nelems() const
		{
			return m_shape.nelems();
		}

		LMAT_ENSURE_INLINE
		matrix_shape<M, N> shape() const
		{
			return m_shape;
		}

		LMAT_ENSURE_INLINE
		index_t row_stride() const
		{
			return 1;
		}

		LMAT_ENSURE_INLINE
		index_t col_stride() const
		{
			return m_shape.nrows();
		}

		LMAT_ENSURE_INLINE
		bool is_contiguous() const
		{
			return true;
		}

		LMAT_ENSURE_INLINE
		bool is_percol_contiguous() const
		{
			return true;
		}

		LMAT_ENSURE_INLINE
		index_t offset(index_t i, index_t j) const
		{
			return internal::cont_cm_sub2offset(m_shape, i, j);
		}

		LMAT_ENSURE_INLINE
		index_t col_offset(index_t j) const
		{
			return m_shape.nrows() * j;
		}

		LMAT_ENSURE_INLINE
		index_t row_offset(index_t i) const
		{
			return i;
		}

		LMAT_ENSURE_INLINE
		index_t lin_offset(index_t i) const
		{
			return i;
		}

	private:
		matrix_shape<M, N> m_shape;
	};


	template<index_t M, index_t N>
	class block_layout_cm : public IMatrixLayout<block_layout_cm<M, N> >
	{
	public:
		LMAT_ENSURE_INLINE
		block_layout_cm(index_t m, index_t n, index_t ldim) : m_shape(m, n), m_leaddim(ldim) { };

	public:
		LMAT_ENSURE_INLINE
		index_t nrows() const
		{
			return m_shape.nrows();
		}

		LMAT_ENSURE_INLINE
		index_t ncolumns() const
		{
			return m_shape.ncolumns();
		}

		LMAT_ENSURE_INLINE
		index_t nelems() const
		{
			return m_shape.nelems();
		}

		LMAT_ENSURE_INLINE
		matrix_shape<M, N> shape() const
		{
			return m_shape;
		}

		LMAT_ENSURE_INLINE
		index_t row_stride() const
		{
			return 1;
		}

		LMAT_ENSURE_INLINE
		index_t col_stride() const
		{
			return m_leaddim;
		}

		LMAT_ENSURE_INLINE
		bool is_contiguous() const
		{
			return m_leaddim == nrows() || ncolumns() == 1;
		}

		LMAT_ENSURE_INLINE
		bool is_percol_contiguous() const
		{
			return true;
		}

		LMAT_ENSURE_INLINE
		index_t offset(index_t i, index_t j) const
		{
			return internal::block_cm_sub2offset(m_shape, i, j, m_leaddim);
		}

		LMAT_ENSURE_INLINE
		index_t col_offset(index_t j) const
		{
			return m_leaddim * j;
		}

		LMAT_ENSURE_INLINE
		index_t row_offset(index_t i) const
		{
			return i;
		}

		LMAT_ENSURE_INLINE
		index_t lin_offset(index_t i) const
		{
			return internal::block_cm_linoffset(m_shape, i, m_leaddim);
		}

	private:
		matrix_shape<M, N> m_shape;
		index_t m_leaddim;
	};

	template<index_t M, index_t N>
	class grid_layout : public IMatrixLayout<grid_layout<M, N> >
	{
	public:
		LMAT_ENSURE_INLINE
		grid_layout(index_t m, index_t n, index_t rs, index_t cs)
		: m_shape(m, n), m_rowstride(rs), m_colstride(cs) { };

	public:
		LMAT_ENSURE_INLINE
		index_t nrows() const
		{
			return m_shape.nrows();
		}

		LMAT_ENSURE_INLINE
		index_t ncolumns() const
		{
			return m_shape.ncolumns();
		}

		LMAT_ENSURE_INLINE
		index_t nelems() const
		{
			return m_shape.nelems();
		}

		LMAT_ENSURE_INLINE
		matrix_shape<M, N> shape() const
		{
			return m_shape;
		}

		LMAT_ENSURE_INLINE
		index_t row_stride() const
		{
			return m_rowstride;
		}

		LMAT_ENSURE_INLINE
		index_t col_stride() const
		{
			return m_colstride;
		}

		LMAT_ENSURE_INLINE
		bool is_contiguous() const
		{
			return (m_rowstride == 1 || nrows() == 1) &&
					(m_colstride == nrows() || ncolumns() == 1);
		}

		LMAT_ENSURE_INLINE
		bool is_percol_contiguous() const
		{
			return m_rowstride() == 1;
		}

		LMAT_ENSURE_INLINE
		index_t offset(index_t i, index_t j) const
		{
			return internal::grid_sub2offset(m_shape, i, j, m_rowstride, m_colstride);
		}

		LMAT_ENSURE_INLINE
		index_t col_offset(index_t j) const
		{
			return m_colstride * j;
		}

		LMAT_ENSURE_INLINE
		index_t row_offset(index_t i) const
		{
			return i * m_rowstride;
		}

		LMAT_ENSURE_INLINE
		index_t lin_offset(index_t i) const
		{
			return internal::grid_linoffset(m_shape, i, m_rowstride, m_colstride);
		}

	private:
		matrix_shape<M, N> m_shape;
		index_t m_rowstride;
		index_t m_colstride;
	};


}

#endif /* MATRIX_LAYOUT_H_ */




