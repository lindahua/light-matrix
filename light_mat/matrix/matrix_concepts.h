/**
 * @file matrix_concepts.h
 *
 * Basic matrix concepts
 *
 * @author Dahua Lin
 */

#ifdef _MSC_VER
#pragma once
#endif

#ifndef LIGHTMAT_MATRIX_CONCEPTS_H_
#define LIGHTMAT_MATRIX_CONCEPTS_H_

#include <light_mat/matrix/matrix_meta.h>


#ifdef LMAT_ENABLE_INDEX_CHECKING
#define LMAT_CHECK_IDX(i, n) check_arg(i >= 0 && i < (n), "Index out of range.");
#define LMAT_CHECK_SUBS(a, i, j) \
		check_arg(i >= 0 && i < (a).nrows() && j >= 0 && j < (a).ncolumns(), "Subscripts out of range.");
#else
#define LMAT_CHECK_IDX(i, n)
#define LMAT_CHECK_SUBS(a, i, j)
#endif


namespace lmat
{

	template<class Mat>
	struct matrix_access_types
	{
		typedef typename matrix_traits<Mat>::value_type value_type;
		typedef typename meta::if_<meta::is_readonly<Mat>,
				const value_type*, value_type*>::type pointer;
		typedef typename meta::if_<meta::is_readonly<Mat>,
				const value_type&, value_type&>::type reference;
	};


	/********************************************
	 *
	 *  IMatrixXpr
	 *
	 ********************************************/

	template<class Derived, typename T>
	class IMatrixXpr
	{
	public:
		typedef T value_type;
		typedef matrix_shape< meta::nrows<Derived>::value, meta::ncols<Derived>::value > shape_type;

		LMAT_CRTP_REF

		LMAT_ENSURE_INLINE index_t nelems() const
		{
			return derived().nelems();
		}

		LMAT_ENSURE_INLINE index_t nrows() const
		{
			return derived().nrows();
		}

		LMAT_ENSURE_INLINE index_t ncolumns() const
		{
			return derived().ncolumns();
		}

		LMAT_ENSURE_INLINE shape_type shape() const
		{
			return derived().shape();
		}

	}; // end class IMatrixXpr


	/********************************************
	 *
	 *  IMatrixXpr
	 *
	 ********************************************/

	template<class Derived, typename T>
	class IEWiseMatrix : public IMatrixXpr<Derived, T>
	{
	public:
		typedef T value_type;
		typedef matrix_shape< meta::nrows<Derived>::value, meta::ncols<Derived>::value > shape_type;

		LMAT_CRTP_REF

		LMAT_ENSURE_INLINE index_t nelems() const
		{
			return derived().nelems();
		}

		LMAT_ENSURE_INLINE index_t nrows() const
		{
			return derived().nrows();
		}

		LMAT_ENSURE_INLINE index_t ncolumns() const
		{
			return derived().ncolumns();
		}

		LMAT_ENSURE_INLINE shape_type shape() const
		{
			return derived().shape();
		}

	}; // end class IEWiseMatrix


	/********************************************
	 *
	 *  IRegularMatrix
	 *
	 ********************************************/

	template<class Derived, typename T>
	class IRegularMatrix : public IEWiseMatrix<Derived, T>
	{
	public:
		typedef typename matrix_access_types<Derived>::value_type value_type;
		typedef const value_type* const_pointer;
		typedef const value_type& const_reference;
		typedef typename matrix_access_types<Derived>::pointer pointer;
		typedef typename matrix_access_types<Derived>::reference reference;

		typedef matrix_shape< meta::nrows<Derived>::value, meta::ncols<Derived>::value > shape_type;

	public:
		LMAT_CRTP_REF

		LMAT_ENSURE_INLINE index_t nelems() const
		{
			return derived().nelems();
		}

		LMAT_ENSURE_INLINE index_t nrows() const
		{
			return derived().nrows();
		}

		LMAT_ENSURE_INLINE index_t ncolumns() const
		{
			return derived().ncolumns();
		}

		LMAT_ENSURE_INLINE shape_type shape() const
		{
			return derived().shape();
		}

		LMAT_ENSURE_INLINE index_t row_stride() const
		{
			return derived().row_stride();
		}

		LMAT_ENSURE_INLINE index_t col_stride() const
		{
			return derived().col_stride();
		}

		LMAT_ENSURE_INLINE const_pointer ptr_data() const
		{
			return derived().ptr_data();
		}

		LMAT_ENSURE_INLINE pointer ptr_data()
		{
			return derived().ptr_data();
		}

		LMAT_ENSURE_INLINE const_pointer ptr_col(const index_t j) const
		{
			LMAT_CHECK_IDX(j, ncolumns())
			return derived().ptr_col(j);
		}

		LMAT_ENSURE_INLINE pointer ptr_col(const index_t j)
		{
			LMAT_CHECK_IDX(j, ncolumns())
			return derived().ptr_col(j);
		}

		LMAT_ENSURE_INLINE const_pointer ptr_row(const index_t i) const
		{
			LMAT_CHECK_IDX(i, nrows())
			return derived().ptr_row(i);
		}

		LMAT_ENSURE_INLINE pointer ptr_row(const index_t i)
		{
			LMAT_CHECK_IDX(i, nrows())
			return derived().ptr_row(i);
		}

		LMAT_ENSURE_INLINE const_reference elem(const index_t i, const index_t j) const
		{
			return derived().elem(i, j);
		}

		LMAT_ENSURE_INLINE reference elem(const index_t i, const index_t j)
		{
			return derived().elem(i, j);
		}

		LMAT_ENSURE_INLINE const_reference operator() (const index_t i, const index_t j) const
		{
			LMAT_CHECK_SUBS(*this, i, j)
			return elem(i, j);
		}

		LMAT_ENSURE_INLINE reference operator() (const index_t i, const index_t j)
		{
			LMAT_CHECK_SUBS(*this, i, j)
			return elem(i, j);
		}

		LMAT_ENSURE_INLINE void require_size(index_t m, index_t n)
		{
			derived().require_size(m, n);
		}

	public:

		// sub-vector views (only apply to compile-time vectors)

		template<class Range>
		LMAT_ENSURE_INLINE
		typename vecview_map<Derived, Range>::const_type
		operator() (const IRange<Range>& rgn) const
		{
			return vecview_map<Derived, Range>::get(derived(), rgn.derived());
		}

		template<class Range>
		LMAT_ENSURE_INLINE
		typename vecview_map<Derived, Range>::type
		operator() (const IRange<Range>& rgn)
		{
			return vecview_map<Derived, Range>::get(derived(), rgn.derived());
		}

		// column views

		LMAT_ENSURE_INLINE
		typename colview_map<Derived, whole>::const_type
		column(const index_t j) const
		{
			LMAT_CHECK_IDX(j, this->ncolumns())
			return colview_map<Derived, whole>::get(derived(), j, whole());
		}

		LMAT_ENSURE_INLINE
		typename colview_map<Derived, whole>::type
		column(const index_t j)
		{
			LMAT_CHECK_IDX(j, this->ncolumns())
			return colview_map<Derived, whole>::get(derived(), j, whole());
		}

		template<class Range>
		LMAT_ENSURE_INLINE
		typename colview_map<Derived, Range>::const_type
		operator()(const IRange<Range>& rgn, const index_t j) const
		{
			LMAT_CHECK_IDX(j, this->ncolumns())
			return colview_map<Derived, Range>::get(derived(), j, rgn.derived());
		}

		template<class Range>
		LMAT_ENSURE_INLINE
		typename colview_map<Derived, Range>::type
		operator()(const IRange<Range>& rgn, const index_t j)
		{
			LMAT_CHECK_IDX(j, this->ncolumns())
			return colview_map<Derived, Range>::get(derived(), j, rgn.derived());
		}

		// row views

		LMAT_ENSURE_INLINE
		typename rowview_map<Derived, whole>::const_type
		row(const index_t i) const
		{
			LMAT_CHECK_IDX(i, this->nrows())
			return rowview_map<Derived, whole>::get(derived(), i, whole());
		}

		LMAT_ENSURE_INLINE
		typename rowview_map<Derived, whole>::type
		row(const index_t i)
		{
			LMAT_CHECK_IDX(i, this->nrows())
			return rowview_map<Derived, whole>::get(derived(), i, whole());
		}

		template<class Range>
		LMAT_ENSURE_INLINE
		typename rowview_map<Derived, Range>::const_type
		operator()(const index_t i, const IRange<Range>& rgn) const
		{
			LMAT_CHECK_IDX(i, this->nrows())
			return rowview_map<Derived, Range>::get(derived(), i, rgn.derived());
		}

		template<class Range>
		LMAT_ENSURE_INLINE
		typename rowview_map<Derived, Range>::type
		operator()(const index_t i, const IRange<Range>& rgn)
		{
			LMAT_CHECK_IDX(i, this->nrows())
			return rowview_map<Derived, Range>::get(derived(), i, rgn.derived());
		}


		// mat-views

		template<class Range0, class Range1>
		LMAT_ENSURE_INLINE
		typename matview_map<Derived, Range0, Range1>::const_type
		operator()(const IRange<Range0>& row_rgn, const IRange<Range1>& col_rgn) const
		{
			return matview_map<Derived, Range0, Range1>::get(derived(),
					row_rgn.derived(), col_rgn.derived());
		}

		template<class Range0, class Range1>
		LMAT_ENSURE_INLINE
		typename matview_map<Derived, Range0, Range1>::type
		operator()(const IRange<Range0>& row_rgn, const IRange<Range1>& col_rgn)
		{
			return matview_map<Derived, Range0, Range1>::get(derived(),
					row_rgn.derived(), col_rgn.derived());
		}

		// diagonal view

		LMAT_ENSURE_INLINE
		typename diagview_map<Derived>::const_type diag() const
		{
			return diagview_map<Derived>::get(derived());
		}

		LMAT_ENSURE_INLINE
		typename diagview_map<Derived>::type diag()
		{
			return diagview_map<Derived>::get(derived());
		}

	}; // end class IRegularMatrix


}

#endif /* MATRIX_CONCEPTS_H_ */


