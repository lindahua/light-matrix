/**
 * @file ref_block.h
 *
 * Classes : cref_block and ref_block
 *
 * @author Dahua Lin
 */

#ifdef _MSC_VER
#pragma once
#endif

#ifndef LIGHTMAT_REF_MATRIX_EX_H_
#define LIGHTMAT_REF_MATRIX_EX_H_

#include <light_mat/matrix/regular_mat_base.h>

namespace lmat
{

	template<typename T, int CTRows, int CTCols>
	struct matrix_traits<cref_block<T, CTRows, CTCols> >
	{
		static const int num_dimensions = 2;
		static const int ct_num_rows = CTRows;
		static const int ct_num_cols = CTCols;

		static const bool is_readonly = true;

		typedef T value_type;
		typedef block_layout_cm<CTRows, CTCols> layout_type;
		typedef cpu_domain domain;
	};


	template<typename T, int CTRows, int CTCols>
	class cref_block : public regular_mat_base<cref_block<T, CTRows, CTCols>, T>
	{
		static_assert( meta::is_supported_matrix_value_type<T>::value,
				"T must be a supported matrix value type");

	public:
		LMAT_DEFINE_REGMAT_CTYPES(T)
		typedef block_layout_cm<CTRows, CTCols> layout_type;

	public:

		LMAT_ENSURE_INLINE
		cref_block(const T* pdata, index_t m, index_t n, index_t ldim)
		: m_data(pdata), m_layout(m, n, ldim)
		{
		}

	private:
		cref_block& operator = (const cref_block& );  // no assignment

	public:
		LMAT_ENSURE_INLINE const layout_type& layout() const
		{
			return m_layout;
		}

		LMAT_ENSURE_INLINE const_pointer ptr_data() const
		{
			return m_data;
		}

		LMAT_DEFINE_NO_RESIZE( cref_block )

	private:
		const T *m_data;
		layout_type m_layout;

	}; // end class cref_block


	template<typename T, int CTRows, int CTCols>
	struct matrix_traits<ref_block<T, CTRows, CTCols> >
	{
		static const int num_dimensions = 2;
		static const int ct_num_rows = CTRows;
		static const int ct_num_cols = CTCols;

		static const bool is_readonly = false;

		typedef T value_type;
		typedef block_layout_cm<CTRows, CTCols> layout_type;
		typedef cpu_domain domain;
	};


	template<typename T, int CTRows, int CTCols>
	class ref_block : public regular_mat_base<ref_block<T, CTRows, CTCols>, T>
	{
		static_assert( meta::is_supported_matrix_value_type<T>::value,
				"T must be a supported matrix value type");

	public:
		LMAT_DEFINE_REGMAT_TYPES(T)
		typedef block_layout_cm<CTRows, CTCols> layout_type;

	public:
		LMAT_ENSURE_INLINE
		ref_block(T* pdata, index_t m, index_t n, index_t ldim)
		: m_data(pdata), m_layout(m, n, ldim)
		{
		}

	public:
		LMAT_ENSURE_INLINE ref_block& operator = (const ref_block& r)
		{
			if (this != &r)
			{
				copy(r, *this);
			}
			return *this;
		}

		template<class Expr>
		LMAT_ENSURE_INLINE ref_block& operator = (const IMatrixXpr<Expr, T>& r)
		{
			assign(r);
			return *this;
		}

	public:
		LMAT_ENSURE_INLINE const layout_type& layout() const
		{
			return m_layout;
		}

		LMAT_ENSURE_INLINE const_pointer ptr_data() const
		{
			return m_data;
		}

		LMAT_ENSURE_INLINE pointer ptr_data()
		{
			return m_data;
		}

		LMAT_DEFINE_NO_RESIZE( ref_block )

	private:
		template<class Expr>
		LMAT_ENSURE_INLINE
		void assign(const IMatrixXpr<Expr, T>& r)
		{
			evaluate(r.derived(), *this);
		}

	private:
		T *m_data;
		layout_type m_layout;

	}; // end ref_block

}

#endif
