/**
 * @file ref_matrix.h
 *
 * ref_matrix classes
 *
 * @author Dahua Lin
 */

#ifdef _MSC_VER
#pragma once
#endif

#ifndef LIGHTMAT_REF_MATRIX_H_
#define LIGHTMAT_REF_MATRIX_H_

#include <light_mat/matrix/regular_mat_base.h>

namespace lmat
{

	/********************************************
	 *
	 *  matrix traits
	 *
	 ********************************************/


	template<typename T, int CM, int CN>
	struct matrix_traits<cref_matrix<T, CM, CN> >
	: public regular_matrix_traits_base<const T, CM, CN, cpu_domain>
	{
		typedef cont_layout_cm<CM, CN> layout_type;
	};

	template<typename T, int CM, int CN>
	struct matrix_traits<ref_matrix<T, CM, CN> >
	: public regular_matrix_traits_base<T, CM, CN, cpu_domain>
	{
		typedef cont_layout_cm<CM, CN> layout_type;
	};


	/********************************************
	 *
	 *  matrix classes
	 *
	 ********************************************/

	template<typename T, int CM, int CN>
	class cref_matrix : public regular_mat_base<cref_matrix<T, CM, CN> >
	{
	public:
		LMAT_DEFINE_REGMAT_TYPES(const T)
		typedef cont_layout_cm<CM, CN> layout_type;

	public:
		LMAT_ENSURE_INLINE
		cref_matrix(const T* pdata, index_t m, index_t n)
		: m_data(pdata), m_layout(m, n)
		{
		}

	private:
		cref_matrix& operator = (const cref_matrix& );  // no assignment

	public:
		LMAT_ENSURE_INLINE const layout_type& layout() const
		{
			return m_layout;
		}

		LMAT_ENSURE_INLINE const_pointer ptr_data() const
		{
			return m_data;
		}

		LMAT_DEFINE_NO_RESIZE( cref_matrix )

	private:
		const T *m_data;
		layout_type m_layout;

	}; // end class cref_matrix


	template<typename T, int CM, int CN>
	class ref_matrix : public regular_mat_base<ref_matrix<T, CM, CN> >
	{
	public:
		LMAT_DEFINE_REGMAT_TYPES(T)
		typedef cont_layout_cm<CM, CN> layout_type;

	public:
		LMAT_ENSURE_INLINE
		ref_matrix(T* pdata, index_t m, index_t n)
		: m_data(pdata), m_layout(m, n)
		{
		}

	public:
		LMAT_ENSURE_INLINE ref_matrix& operator = (const ref_matrix& r)
		{
			if (this != &r)
			{
				copy(r, *this);
			}
			return *this;
		}

		template<class Expr>
		LMAT_ENSURE_INLINE ref_matrix& operator = (const IMatrixXpr<Expr, T>& r)
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

		LMAT_DEFINE_NO_RESIZE( ref_matrix )

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

	}; // end ref_matrix


	/********************************************
	 *
	 *  vectors derived from (c)ref_matrix
	 *
	 ********************************************/

	template<typename T, int CM>
	class cref_col: public cref_matrix<T, CM, 1>
	{
		typedef cref_matrix<T, CM, 1> base_mat_t;

	public:
		LMAT_ENSURE_INLINE
		cref_col(const T* pdata, index_t m)
		: base_mat_t(pdata, m, 1) { }

		LMAT_ENSURE_INLINE
		cref_col(const base_mat_t& s)
		: base_mat_t(s) { }

	};

	template<typename T, int CM>
	class ref_col: public ref_matrix<T, CM, 1>
	{
		typedef ref_matrix<T, CM, 1> base_mat_t;

	public:
		LMAT_ENSURE_INLINE
		ref_col(T* pdata, index_t m)
		: base_mat_t(pdata, m, 1) { }

		LMAT_ENSURE_INLINE
		ref_col(const base_mat_t& s)
		: base_mat_t(s) { }

		template<class Expr>
		LMAT_ENSURE_INLINE ref_col& operator = (const IMatrixXpr<Expr, T>& r)
		{
			base_mat_t::operator = (r);
			return *this;
		}

	};


	template<typename T, int CN>
	class cref_row: public cref_matrix<T, 1, CN>
	{
		typedef cref_matrix<T, 1, CN> base_mat_t;

	public:
		LMAT_ENSURE_INLINE
		cref_row(const T* pdata, index_t n)
		: base_mat_t(pdata, 1, n) { }

		LMAT_ENSURE_INLINE
		cref_row(const base_mat_t& s)
		: base_mat_t(s) { }
	};

	template<typename T, int CN>
	class ref_row: public ref_matrix<T, 1, CN>
	{
		typedef ref_matrix<T, 1, CN> base_mat_t;

	public:
		LMAT_ENSURE_INLINE
		ref_row(T* pdata, index_t n)
		: base_mat_t(pdata, 1, n) { }

		LMAT_ENSURE_INLINE
		ref_row(const base_mat_t& s)
		: base_mat_t(s) { }

		LMAT_ENSURE_INLINE ref_row& operator = (const base_mat_t& r)
		{
			base_mat_t::operator = (r);
			return *this;
		}

		template<class Expr>
		LMAT_ENSURE_INLINE ref_row& operator = (const IMatrixXpr<Expr, T>& r)
		{
			base_mat_t::operator = (r);
			return *this;
		}
	};

}

#endif /* REF_MATRIX_H_ */
