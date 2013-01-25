/**
 * @file ref_grid.h
 *
 * @brief cref_grid and ref_grid classes
 *
 * @author Dahua Lin
 */

#ifdef _MSC_VER
#pragma once
#endif

#ifndef LIGHTMAT_REF_GRID_H_
#define LIGHTMAT_REF_GRID_H_

#include <light_mat/matrix/regular_mat_base.h>

namespace lmat
{

	/********************************************
	 *
	 *  matrix traits
	 *
	 ********************************************/

	template<typename T, index_t CM, index_t CN>
	struct matrix_traits<cref_grid<T, CM, CN> >
	: public regular_matrix_traits_base<const T, CM, CN, cpu_domain>
	{
		typedef grid_layout<CM, CN> layout_type;
	};

	template<typename T, index_t CM, index_t CN>
	struct matrix_traits<ref_grid<T, CM, CN> >
	: public regular_matrix_traits_base<T, CM, CN, cpu_domain>
	{
		typedef grid_layout<CM, CN> layout_type;
	};


	/********************************************
	 *
	 *  matrix classes
	 *
	 ********************************************/

	template<typename T, index_t CM, index_t CN>
	class cref_grid : public regular_mat_base<cref_grid<T, CM, CN> >
	{
	public:
		LMAT_DEFINE_REGMAT_TYPES(const T)
		typedef grid_layout<CM, CN> layout_type;

	public:

		LMAT_ENSURE_INLINE
		cref_grid(const T* pdata, index_t m, index_t n, index_t rs, index_t cs)
		: m_data(pdata), m_layout(m, n, rs, cs)
		{
		}

	private:
		cref_grid& operator = (const cref_grid& );  // no assignment

	public:
		LMAT_ENSURE_INLINE const layout_type& layout() const
		{
			return m_layout;
		}

		LMAT_ENSURE_INLINE const_pointer ptr_data() const
		{
			return m_data;
		}

		LMAT_DEFINE_NO_RESIZE( cref_grid )

	private:
		const T *m_data;
		layout_type m_layout;

	}; // end class cref_grid


	template<typename T, index_t CM, index_t CN>
	class ref_grid : public regular_mat_base<ref_grid<T, CM, CN> >
	{
	public:
		LMAT_DEFINE_REGMAT_TYPES(T)
		typedef grid_layout<CM, CN> layout_type;

	public:
		LMAT_ENSURE_INLINE
		ref_grid(T* pdata, index_t m, index_t n, index_t rs, index_t cs)
		: m_data(pdata), m_layout(m, n, rs, cs)
		{
		}

	public:
		LMAT_ENSURE_INLINE ref_grid& operator = (const ref_grid& r)
		{
			if (this != &r)
			{
				copy(r, *this);
			}
			return *this;
		}

		template<class Expr>
		LMAT_ENSURE_INLINE ref_grid& operator = (const IMatrixXpr<Expr, T>& r)
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

		LMAT_DEFINE_NO_RESIZE( ref_grid )

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

	}; // end ref_grid

}

#endif /* REF_GRID_H_ */
