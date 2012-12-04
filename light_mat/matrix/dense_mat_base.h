/**
 * @file dense_mat_base.h
 *
 * @brief The implementation base of all dense matrix classes
 *
 * @author Dahua Lin
 */


#ifdef _MSC_VER
#pragma once
#endif

#ifndef LIGHTMAT_DENSE_MAT_BASE_H_
#define LIGHTMAT_DENSE_MAT_BASE_H_

#include <light_mat/matrix/matrix_base.h>
#include <light_mat/matrix/matrix_layout.h>

namespace lmat
{
	template<class Derived, typename T>
	class dense_mat_base : public IDenseMatrix<Derived, T>
	{
	public:
		LMAT_MAT_TRAITS_DEFS_FOR_BASE(Derived, T)
		typedef typename meta::shape<Derived>::type shape_type;
		typedef typename matrix_traits<Derived>::layout_type layout_type;

#ifdef LMAT_USE_STATIC_ASSERT
		static_assert(is_same<typename layout_traits<layout_type>::shape_type, shape_type>::value,
				"Inconsistent shape type");
#endif

		LMAT_CRTP_REF

	public:

		LMAT_ENSURE_INLINE const layout_type& layout() const
		{
			return derived().layout();
		}

		LMAT_ENSURE_INLINE const_pointer ptr_data() const
		{
			return derived().ptr_data();
		}

		LMAT_ENSURE_INLINE pointer ptr_data()
		{
			return derived().ptr_data();
		}

	public:

		LMAT_ENSURE_INLINE index_type nelems() const
		{
			return layout().nelems();
		}

		LMAT_ENSURE_INLINE index_type nrows() const
		{
			return layout().nrows();
		}

		LMAT_ENSURE_INLINE index_type ncolumns() const
		{
			return layout().ncolumns();
		}

		LMAT_ENSURE_INLINE index_type row_stride() const
		{
			return layout().row_stride();
		}

		LMAT_ENSURE_INLINE index_type col_stride() const
		{
			return layout().col_stride();
		}

		LMAT_ENSURE_INLINE const_pointer ptr_col(const index_type j) const
		{
			return ptr_data() + layout().col_offset(j);
		}

		LMAT_ENSURE_INLINE pointer ptr_col(const index_type j)
		{
			return ptr_data() + layout().col_offset(j);
		}

		LMAT_ENSURE_INLINE const_pointer ptr_row(const index_type i) const
		{
			return ptr_data() + layout().row_offset(i);
		}

		LMAT_ENSURE_INLINE pointer ptr_row(const index_type i)
		{
			return ptr_data() + layout().row_offset(i);
		}

		LMAT_ENSURE_INLINE const_reference elem(const index_type i, const index_type j) const
		{
			return ptr_data()[layout().offset(i, j)];
		}

		LMAT_ENSURE_INLINE reference elem(const index_type i, const index_type j)
		{
			return ptr_data()[layout().offset(i, j)];
		}

		LMAT_ENSURE_INLINE const_reference operator[] (const index_type i) const
		{
			return ptr_data()[layout().lin_offset(i)];
		}

		LMAT_ENSURE_INLINE reference operator[] (const index_type i)
		{
			return ptr_data()[layout().lin_offset(i)];
		}

	};


}

#endif /* DENSE_MAT_BASE_H_ */
