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

#include <light_mat/matrix/matrix_shape.h>

namespace lmat
{

	template<int M, int N> class continuous_layout_cm;
	template<int M, int N> class block_layout_cm;
	template<int M, int N> class grid_layout_cm;

	/********************************************
	 *
	 *  Layout traits
	 *
	 ********************************************/

	template<class Layout> struct layout_traits;

	template<int M, int N>
	struct layout_traits<continuous_layout_cm<M, N> >
	{
		static const int ct_num_rows = M;
		static const int ct_num_cols = N;

		static const bool ct_is_continuous = true;
		static const bool ct_is_percol_continuous = true;

		typedef matrix_shape<M, N> shape_type;
	};

	template<int M, int N>
	struct layout_traits<block_layout_cm<M, N> >
	{
		static const int ct_num_rows = M;
		static const int ct_num_cols = N;

		static const bool ct_is_continuous = (N == 1);
		static const bool ct_is_percol_continuous = true;

		typedef matrix_shape<M, N> shape_type;
	};

	template<int M, int N>
	struct layout_traits<grid_layout_cm<M, N> >
	{
		static const int ct_num_rows = M;
		static const int ct_num_cols = N;

		static const bool ct_is_continuous = (M == 1 && N == 1);
		static const bool ct_is_percol_continuous = N == 1;

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
		bool is_continuous() const
		{
			return derived().is_continuous();
		}

		LMAT_ENSURE_INLINE
		bool is_percol_continuous() const
		{
			return derived().is_percol_continuous();
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

	};


	inline index_t raise_no_linear_offset()
	{
		throw invalid_operation("Linear offset is only supported for compile-time vectors");
	}


	/********************************************
	 *
	 *  specific layout classes
	 *
	 ********************************************/

	namespace detail
	{
		template<int M, int N>
		LMAT_ENSURE_INLINE
		inline index_t cont_cm_sub2offset(const matrix_shape<M, N>& shape, index_t i, index_t j)
		{
			return i + shape.nrows() * j;
		}

		template<int M>
		LMAT_ENSURE_INLINE
		inline index_t cont_cm_sub2offset(const matrix_shape<M, 1>& shape, index_t i, index_t j)
		{
			return i;
		}

		template<int N>
		LMAT_ENSURE_INLINE
		inline index_t cont_cm_sub2offset(const matrix_shape<1, N>& shape, index_t i, index_t j)
		{
			return shape.nrows() * j;
		}

		LMAT_ENSURE_INLINE
		inline index_t cont_cm_sub2offset(const matrix_shape<1, 1>& shape, index_t i, index_t j)
		{
			return 0;
		}
	}

	template<int M, int N>
	class continuous_layout_cm : public IMatrixLayout<continuous_layout_cm<M, N> >
	{
	public:
		LMAT_ENSURE_INLINE
		continuous_layout_cm() : m_shape() { }

		LMAT_ENSURE_INLINE
		continuous_layout_cm(index_t m, index_t n) : m_shape(m, n) { };

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
		bool is_continuous() const
		{
			return true;
		}

		LMAT_ENSURE_INLINE
		bool is_percol_continuous() const
		{
			return true;
		}

		LMAT_ENSURE_INLINE
		index_t offset(index_t i, index_t j) const
		{
			return detail::cont_cm_sub2offset(m_shape, i, j);
		}

		LMAT_ENSURE_INLINE
		index_t col_offset(index_t j) const
		{
			return m_shape.nrows() * j;
		}

		LMAT_ENSURE_INLINE
		index_t lin_offset(index_t i) const
		{
			return i;
		}

	private:
		matrix_shape<M, N> m_shape;
	};


	namespace detail
	{
		template<int M, int N>
		LMAT_ENSURE_INLINE
		inline index_t block_cm_sub2offset(const matrix_shape<M, N>& shape, index_t i, index_t j, index_t ldim)
		{
			return i + ldim * j;
		}

		template<int M>
		LMAT_ENSURE_INLINE
		inline index_t block_cm_sub2offset(const matrix_shape<M, 1>& shape, index_t i, index_t j, index_t ldim)
		{
			return i;
		}

		template<int N>
		LMAT_ENSURE_INLINE
		inline index_t block_cm_sub2offset(const matrix_shape<1, N>& shape, index_t i, index_t j, index_t ldim)
		{
			return ldim * j;
		}

		LMAT_ENSURE_INLINE
		inline index_t block_cm_sub2offset(const matrix_shape<1, 1>& shape, index_t i, index_t j, index_t ldim)
		{
			return 0;
		}

		template<int M, int N>
		LMAT_ENSURE_INLINE
		inline index_t block_cm_linoffset(const matrix_shape<M, N>& shape, index_t i, index_t ldim)
		{
			return raise_no_linear_offset();
		}

		template<int M>
		LMAT_ENSURE_INLINE
		inline index_t block_cm_linoffset(const matrix_shape<M, 1>& shape, index_t i, index_t ldim)
		{
			return i;
		}

		template<int N>
		LMAT_ENSURE_INLINE
		inline index_t block_cm_linoffset(const matrix_shape<1, N>& shape, index_t i, index_t ldim)
		{
			return ldim * i;
		}

		LMAT_ENSURE_INLINE
		inline index_t block_cm_linoffset(const matrix_shape<1, 1>& shape, index_t i, index_t ldim)
		{
			return 0;
		}
	}


	template<int M, int N>
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
		bool is_continuous() const
		{
			return m_leaddim == nrows();
		}

		LMAT_ENSURE_INLINE
		bool is_percol_continuous() const
		{
			return true;
		}

		LMAT_ENSURE_INLINE
		index_t offset(index_t i, index_t j) const
		{
			return detail::block_cm_sub2offset(m_shape, i, j, m_leaddim);
		}

		LMAT_ENSURE_INLINE
		index_t col_offset(index_t j) const
		{
			return m_leaddim * j;
		}

		LMAT_ENSURE_INLINE
		index_t lin_offset(index_t i) const
		{
			return detail::block_cm_linoffset(m_shape, i, m_leaddim);
		}

	private:
		matrix_shape<M, N> m_shape;
		index_t m_leaddim;
	};



	namespace detail
	{
		template<int M, int N>
		LMAT_ENSURE_INLINE
		inline index_t grid_cm_sub2offset(const matrix_shape<M, N>& shape, index_t i, index_t j, index_t rs, index_t cs)
		{
			return rs * i + cs * j;
		}

		template<int M>
		LMAT_ENSURE_INLINE
		inline index_t grid_cm_sub2offset(const matrix_shape<M, 1>& shape, index_t i, index_t j, index_t rs, index_t cs)
		{
			return rs * i;
		}

		template<int N>
		LMAT_ENSURE_INLINE
		inline index_t grid_cm_sub2offset(const matrix_shape<1, N>& shape, index_t i, index_t j, index_t rs, index_t cs)
		{
			return cs * j;
		}

		LMAT_ENSURE_INLINE
		inline index_t grid_cm_sub2offset(const matrix_shape<1, 1>& shape, index_t i, index_t j, index_t rs, index_t cs)
		{
			return 0;
		}

		template<int M, int N>
		LMAT_ENSURE_INLINE
		inline index_t grid_cm_linoffset(const matrix_shape<M, N>& shape, index_t i, index_t rs, index_t cs)
		{
			return raise_no_linear_offset();
		}

		template<int M>
		LMAT_ENSURE_INLINE
		inline index_t grid_cm_linoffset(const matrix_shape<M, 1>& shape, index_t i, index_t rs, index_t cs)
		{
			return i * rs;
		}

		template<int N>
		LMAT_ENSURE_INLINE
		inline index_t grid_cm_linoffset(const matrix_shape<1, N>& shape, index_t i, index_t rs, index_t cs)
		{
			return i * cs;
		}

		LMAT_ENSURE_INLINE
		inline index_t grid_cm_linoffset(const matrix_shape<1, 1>& shape, index_t i, index_t rs, index_t cs)
		{
			return 0;
		}
	}


	template<int M, int N>
	class grid_layout_cm : public IMatrixLayout<grid_layout_cm<M, N> >
	{
	public:
		LMAT_ENSURE_INLINE
		grid_layout_cm(index_t m, index_t n, index_t rs, index_t cs)
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
		bool is_continuous() const
		{
			return m_rowstride == 1 && m_colstride == nrows();
		}

		LMAT_ENSURE_INLINE
		bool is_percol_continuous() const
		{
			return true;
		}

		LMAT_ENSURE_INLINE
		index_t offset(index_t i, index_t j) const
		{
			return detail::grid_cm_sub2offset(m_shape, i, j, m_rowstride, m_colstride);
		}

		LMAT_ENSURE_INLINE
		index_t col_offset(index_t j) const
		{
			return m_colstride * j;
		}

		LMAT_ENSURE_INLINE
		index_t lin_offset(index_t i) const
		{
			return detail::grid_cm_linoffset(m_shape, i, m_rowstride, m_colstride);
		}

	private:
		matrix_shape<M, N> m_shape;
		index_t m_rowstride;
		index_t m_colstride;
	};


}

#endif /* MATRIX_LAYOUT_H_ */




