/**
 * @file ref_matrix_ex.h
 *
 * ref_matrix_ex classes
 *
 * @author Dahua Lin
 */

#ifdef _MSC_VER
#pragma once
#endif

#ifndef LIGHTMAT_REF_MATRIX_EX_H_
#define LIGHTMAT_REF_MATRIX_EX_H_

#include "bits/ref_matrix_ex_internal.h"

namespace lmat
{
	/********************************************
	 *
	 *  cref_matrix_ex
	 *
	 ********************************************/

	template<typename T, int CTRows, int CTCols, typename Align>
	struct matrix_traits<cref_matrix_ex<T, CTRows, CTCols, Align> >
	{
		static const int num_dimensions = 2;
		static const int compile_time_num_rows = CTRows;
		static const int compile_time_num_cols = CTCols;

		static const bool is_readonly = true;

		typedef T value_type;
	};

	template<typename T, int CTRows, int CTCols, typename Align>
	struct has_continuous_layout<cref_matrix_ex<T, CTRows, CTCols, Align> >
	{
		static const bool value = (CTCols == 1);
	};

	template<typename T, int CTRows, int CTCols, typename Align>
	struct is_base_aligned<cref_matrix_ex<T, CTRows, CTCols, Align> >
	{
		static const bool value = is_base_aligned_from_tag<Align>::value;
	};

	template<typename T, int CTRows, int CTCols, typename Align>
	struct is_percol_aligned<cref_matrix_ex<T, CTRows, CTCols, Align> >
	{
		static const bool value = is_percol_aligned_from_tag<Align>::value;
	};

	template<typename T, int CTRows, int CTCols, typename Align>
	struct is_linear_accessible<cref_matrix_ex<T, CTRows, CTCols, Align> >
	{
		static const bool value = (CTRows == 1 || CTCols == 1);
	};


	template<typename T, int CTRows, int CTCols, typename Align>
	class cref_matrix_ex : public IDenseMatrix<cref_matrix_ex<T, CTRows, CTCols, Align>, T>
	{
	public:
		LMAT_MAT_TRAITS_CDEFS(T)

	public:

		LMAT_ENSURE_INLINE
		cref_matrix_ex(const T* pdata, index_type m, index_type n, index_type ldim)
		: m_internal(pdata, m, n, ldim)
		{
		}

	private:
		cref_matrix_ex& operator = (const cref_matrix_ex& );  // no assignment

	public:
		LMAT_ENSURE_INLINE index_type nelems() const
		{
			return m_internal.nelems();
		}

		LMAT_ENSURE_INLINE size_type size() const
		{
			return static_cast<size_type>(nelems());
		}

		LMAT_ENSURE_INLINE index_type nrows() const
		{
			return m_internal.nrows();
		}

		LMAT_ENSURE_INLINE index_type ncolumns() const
		{
			return m_internal.ncolumns();
		}

		LMAT_ENSURE_INLINE const_pointer ptr_data() const
		{
			return m_internal.ptr_data();
		}

		LMAT_ENSURE_INLINE index_type lead_dim() const
		{
			return m_internal.lead_dim();
		}

		LMAT_ENSURE_INLINE index_type offset(index_type i, index_type j) const
		{
			return matrix_indexer<CTRows, CTCols>::offset(lead_dim(), i, j);
		}

		LMAT_ENSURE_INLINE const_reference elem(index_type i, index_type j) const
		{
			return m_internal.ptr_data()[offset(i, j)];
		}

		LMAT_ENSURE_INLINE const_reference operator[] (index_type i) const
		{
			return m_internal.ptr_data()[detail::ref_ex_linear_offset(*this, i)];
		}

	private:
		detail::ref_matrix_ex_internal<const T, CTRows, CTCols> m_internal;

	}; // end class cref_matrix_ex



	/********************************************
	 *
	 *  ref_matrix_ex
	 *
	 ********************************************/


	template<typename T, int CTRows, int CTCols, typename Align>
	struct matrix_traits<ref_matrix_ex<T, CTRows, CTCols, Align> >
	{
		static const int num_dimensions = 2;
		static const int compile_time_num_rows = CTRows;
		static const int compile_time_num_cols = CTCols;

		static const bool is_readonly = false;

		typedef T value_type;
	};

	template<typename T, int CTRows, int CTCols, typename Align>
	struct has_continuous_layout<ref_matrix_ex<T, CTRows, CTCols, Align> >
	{
		static const bool value = (CTCols == 1);
	};

	template<typename T, int CTRows, int CTCols, typename Align>
	struct is_base_aligned<ref_matrix_ex<T, CTRows, CTCols, Align> >
	{
		static const bool value = is_base_aligned_from_tag<Align>::value;
	};

	template<typename T, int CTRows, int CTCols, typename Align>
	struct is_percol_aligned<ref_matrix_ex<T, CTRows, CTCols, Align> >
	{
		static const bool value = is_percol_aligned_from_tag<Align>::value;
	};

	template<typename T, int CTRows, int CTCols, typename Align>
	struct is_linear_accessible<ref_matrix_ex<T, CTRows, CTCols, Align> >
	{
		static const bool value = (CTRows == 1 || CTCols == 1);
	};


	template<typename T, int CTRows, int CTCols, typename Align>
	class ref_matrix_ex : public IDenseMatrix<ref_matrix_ex<T, CTRows, CTCols, Align>, T>
	{
	public:
		LMAT_MAT_TRAITS_DEFS(T)

	public:
		LMAT_ENSURE_INLINE
		ref_matrix_ex(T* pdata, index_type m, index_type n, index_type ldim)
		: m_internal(pdata, m, n, ldim)
		{
		}

	public:
		LMAT_ENSURE_INLINE ref_matrix_ex& operator = (const ref_matrix_ex& r)
		{
			if (this != &r)
			{
				assign_to(r, *this);
			}
			return *this;
		}

		template<class Other>
		LMAT_ENSURE_INLINE ref_matrix_ex& operator = (const IMatrixView<Other, T>& r)
		{
			assign_to(r.derived(), *this);
			return *this;
		}

		template<class Expr>
		LMAT_ENSURE_INLINE ref_matrix_ex& operator = (const IMatrixXpr<Expr, T>& r)
		{
			assign_to(r.derived(), *this);
			return *this;
		}

	public:
		LMAT_ENSURE_INLINE index_type nelems() const
		{
			return m_internal.nelems();
		}

		LMAT_ENSURE_INLINE size_type size() const
		{
			return static_cast<size_type>(nelems());
		}

		LMAT_ENSURE_INLINE index_type nrows() const
		{
			return m_internal.nrows();
		}

		LMAT_ENSURE_INLINE index_type ncolumns() const
		{
			return m_internal.ncolumns();
		}

		LMAT_ENSURE_INLINE const_pointer ptr_data() const
		{
			return m_internal.ptr_data();
		}

		LMAT_ENSURE_INLINE pointer ptr_data()
		{
			return m_internal.ptr_data();
		}

		LMAT_ENSURE_INLINE index_type lead_dim() const
		{
			return m_internal.lead_dim();
		}

		LMAT_ENSURE_INLINE index_type offset(index_type i, index_type j) const
		{
			return matrix_indexer<CTRows, CTCols>::offset(lead_dim(), i, j);
		}

		LMAT_ENSURE_INLINE const_reference elem(index_type i, index_type j) const
		{
			return m_internal.ptr_data()[offset(i, j)];
		}

		LMAT_ENSURE_INLINE reference elem(index_type i, index_type j)
		{
			return m_internal.ptr_data()[offset(i, j)];
		}

		LMAT_ENSURE_INLINE const_reference operator[] (index_type i) const
		{
			return m_internal.ptr_data()[detail::ref_ex_linear_offset(*this, i)];
		}

		LMAT_ENSURE_INLINE reference operator[] (index_type i)
		{
			return m_internal.ptr_data()[detail::ref_ex_linear_offset(*this, i)];
		}

	private:
		detail::ref_matrix_ex_internal<T, CTRows, CTCols> m_internal;

	}; // end ref_matrix_ex

}

#endif /* REF_MATRIX_EX_H_ */
