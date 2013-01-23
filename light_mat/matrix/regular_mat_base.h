/**
 * @file regular_mat_base.h
 *
 * @brief The implementation base of all regular matrix classes
 *
 * @author Dahua Lin
 */


#ifdef _MSC_VER
#pragma once
#endif

#ifndef LIGHTMAT_REGULAR_MAT_BASE_H_
#define LIGHTMAT_REGULAR_MAT_BASE_H_

#include <light_mat/matrix/matrix_base.h>
#include <light_mat/matrix/matrix_layout.h>


#define LMAT_DEFINE_REGMAT_TYPES(QT) \
	typedef typename std::remove_const<QT>::type value_type; \
	typedef const value_type& const_reference; \
	typedef const value_type* const_pointer; \
	typedef QT& reference; \
	typedef QT* pointer;

#define LMAT_DEFINE_NO_RESIZE( classname ) \
		LMAT_ENSURE_INLINE void require_size(index_t m, index_t n) { \
			check_arg(this->nrows() == m && this->ncolumns() == n, \
					"Cannot change the size of an instance of class " #classname); }

namespace lmat
{
	template<class Derived>
	class regular_mat_base
	: public IRegularMatrix<Derived, typename matrix_traits<Derived>::value_type>
	{

	public:
		typedef typename matrix_traits<Derived>::value_type value_type;

		static_assert( meta::is_supported_matrix_value_type<value_type>::value,
				"T must be a supported matrix value type");

		typedef typename meta::pointer_of<Derived>::type pointer;
		typedef typename meta::reference_of<Derived>::type reference;
		typedef typename meta::const_pointer_of<Derived>::type const_pointer;
		typedef typename meta::const_reference_of<Derived>::type const_reference;

		typedef typename meta::shape<Derived>::type shape_type;
		typedef typename matrix_traits<Derived>::layout_type layout_type;

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

		LMAT_ENSURE_INLINE void require_size(index_t m, index_t n)
		{
			derived().require_size(m, n);
		}

	public:

		LMAT_ENSURE_INLINE index_t nelems() const
		{
			return layout().nelems();
		}

		LMAT_ENSURE_INLINE index_t nrows() const
		{
			return layout().nrows();
		}

		LMAT_ENSURE_INLINE index_t ncolumns() const
		{
			return layout().ncolumns();
		}

		LMAT_ENSURE_INLINE shape_type shape() const
		{
			return layout().shape();
		}

		LMAT_ENSURE_INLINE bool is_contiguous() const
		{
			return layout().is_contiguous();
		}

		LMAT_ENSURE_INLINE bool is_percol_contiguous() const
		{
			return layout().is_percol_contiguous();
		}

		LMAT_ENSURE_INLINE index_t row_stride() const
		{
			return layout().row_stride();
		}

		LMAT_ENSURE_INLINE index_t col_stride() const
		{
			return layout().col_stride();
		}

		LMAT_ENSURE_INLINE const_pointer ptr_col(const index_t j) const
		{
			return ptr_data() + layout().col_offset(j);
		}

		LMAT_ENSURE_INLINE pointer ptr_col(const index_t j)
		{
			return ptr_data() + layout().col_offset(j);
		}

		LMAT_ENSURE_INLINE const_pointer ptr_row(const index_t i) const
		{
			return ptr_data() + layout().row_offset(i);
		}

		LMAT_ENSURE_INLINE pointer ptr_row(const index_t i)
		{
			return ptr_data() + layout().row_offset(i);
		}

		LMAT_ENSURE_INLINE const_reference elem(const index_t i, const index_t j) const
		{
			return ptr_data()[layout().offset(i, j)];
		}

		LMAT_ENSURE_INLINE reference elem(const index_t i, const index_t j)
		{
			return ptr_data()[layout().offset(i, j)];
		}

		LMAT_ENSURE_INLINE const_reference operator[] (const index_t i) const
		{
			return ptr_data()[layout().lin_offset(i)];
		}

		LMAT_ENSURE_INLINE reference operator[] (const index_t i)
		{
			return ptr_data()[layout().lin_offset(i)];
		}

	};

}





#endif
