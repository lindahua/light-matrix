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

#ifndef LIGHTMAT_REF_BLOCK_H_
#define LIGHTMAT_REF_BLOCK_H_

#include <light_mat/matrix/regular_mat_base.h>

namespace lmat
{
	/********************************************
	 *
	 *  matrix traits
	 *
	 ********************************************/

	template<typename T, index_t CM, index_t CN>
	struct matrix_traits<cref_block<T, CM, CN> >
	: public regular_matrix_traits_base<const T, CM, CN, cpu_domain>
	{
		typedef block_layout_cm<CM, CN> layout_type;
	};


	template<typename T, index_t CM, index_t CN>
	struct matrix_traits<ref_block<T, CM, CN> >
	: public regular_matrix_traits_base<T, CM, CN, cpu_domain>
	{
		typedef block_layout_cm<CM, CN> layout_type;
	};


	/********************************************
	 *
	 *  matrix classes
	 *
	 ********************************************/

	template<typename T, index_t CM, index_t CN>
	class cref_block : public regular_mat_base<cref_block<T, CM, CN> >
	{
	public:
		LMAT_DEFINE_REGMAT_TYPES(const T)
		typedef block_layout_cm<CM, CN> layout_type;

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




	template<typename T, index_t CM, index_t CN>
	class ref_block : public regular_mat_base<ref_block<T, CM, CN> >
	{
	public:
		LMAT_DEFINE_REGMAT_TYPES(T)
		typedef block_layout_cm<CM, CN> layout_type;

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
