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

#define LMAT_MAT_TRAITS_DEFS_FOR_BASE(D, T) \
	typedef T value_type; \
	typedef const value_type* const_pointer; \
	typedef const value_type& const_reference; \
	typedef typename mat_access<D>::pointer pointer; \
	typedef typename mat_access<D>::reference reference; \
	typedef size_t size_type; \
	typedef index_t difference_type; \
	typedef index_t index_type;

#define LMAT_MAT_TRAITS_CDEFS(T) \
	typedef T value_type; \
	typedef const value_type* const_pointer; \
	typedef const value_type& const_reference; \
	typedef size_t size_type; \
	typedef index_t difference_type; \
	typedef index_t index_type;

#define LMAT_MAT_TRAITS_DEFS(T) \
	typedef T value_type; \
	typedef const value_type* const_pointer; \
	typedef const value_type& const_reference; \
	typedef value_type* pointer; \
	typedef value_type& reference; \
	typedef size_t size_type; \
	typedef index_t difference_type; \
	typedef index_t index_type;


namespace lmat
{

	/********************************************
	 *
	 *  Concepts
	 *
	 *  Each concept is associated with a
	 *  class as a static polymorphism base
	 *
	 ********************************************/

	template<class Derived, typename T>
	class IMatrixXpr
	{
	public:
		LMAT_MAT_TRAITS_DEFS_FOR_BASE(Derived, T)
		LMAT_CRTP_REF

		LMAT_ENSURE_INLINE index_type nelems() const
		{
			return derived().nelems();
		}

		LMAT_ENSURE_INLINE size_type size() const
		{
			return derived().size();
		}

		LMAT_ENSURE_INLINE index_type nrows() const
		{
			return derived().nrows();
		}

		LMAT_ENSURE_INLINE index_type ncolumns() const
		{
			return derived().ncolumns();
		}

		LMAT_ENSURE_INLINE
		typename unary_expr_map<transpose_t, ref_arg_t, Derived>::type
		trans() const
		{
			return unary_expr_map<transpose_t, ref_arg_t, Derived>::get(
					transpose_t(), ref_arg(derived()));
		}

	}; // end class IMatrixBase


	template<class Mat, typename T>
	LMAT_ENSURE_INLINE
	inline bool is_subscripts_in_range(const IMatrixXpr<Mat, T>& X, index_t i, index_t j)
	{
		return i >= 0 && i < X.nrows() && j >= 0 && j < X.ncolumns();
	}

	template<class Mat, typename T>
	LMAT_ENSURE_INLINE
	inline void check_subscripts_in_range(const IMatrixXpr<Mat, T>& X, index_t i, index_t j)
	{
#ifndef BCSLIB_NO_DEBUG
		check_range(is_subscripts_in_range(X, i, j),
				"Attempted to access element with subscripts out of valid range.");
#endif
	}


	/**
	 * The interfaces for matrix views
	 */
	template<class Derived, typename T>
	class IMatrixView : public IMatrixXpr<Derived, T>
	{
	public:
		LMAT_MAT_TRAITS_DEFS_FOR_BASE(Derived, T)
		LMAT_CRTP_REF

		LMAT_ENSURE_INLINE index_type nelems() const
		{
			return derived().nelems();
		}

		LMAT_ENSURE_INLINE size_type size() const
		{
			return derived().size();
		}

		LMAT_ENSURE_INLINE index_type nrows() const
		{
			return derived().nrows();
		}

		LMAT_ENSURE_INLINE index_type ncolumns() const
		{
			return derived().ncolumns();
		}

		LMAT_ENSURE_INLINE value_type elem(const index_type i, const index_type j) const
		{
			return derived().elem(i, j);
		}

		LMAT_ENSURE_INLINE value_type operator() (const index_type i, const index_type j) const
		{
			check_subscripts_in_range(*this, i, j);
			return elem(i, j);
		}

	}; // end class IDenseMatrixView


	/**
	 * The interfaces for matrix blocks
	 */
	template<class Derived, typename T>
	class IDenseMatrix : public IMatrixView<Derived, T>
	{
	public:
		LMAT_MAT_TRAITS_DEFS_FOR_BASE(Derived, T)
		LMAT_CRTP_REF

	public:
		LMAT_ENSURE_INLINE index_type nelems() const
		{
			return derived().nelems();
		}

		LMAT_ENSURE_INLINE size_type size() const
		{
			return derived().size();
		}

		LMAT_ENSURE_INLINE index_type nrows() const
		{
			return derived().nrows();
		}

		LMAT_ENSURE_INLINE index_type ncolumns() const
		{
			return derived().ncolumns();
		}

		LMAT_ENSURE_INLINE index_t lead_dim() const
		{
			return derived().lead_dim();
		}

		LMAT_ENSURE_INLINE const_pointer ptr_data() const
		{
			return derived().ptr_data();
		}

		LMAT_ENSURE_INLINE pointer ptr_data()
		{
			return derived().ptr_data();
		}

		LMAT_ENSURE_INLINE const_pointer ptr_col(const index_type j) const
		{
			return derived().ptr_col(j);
		}

		LMAT_ENSURE_INLINE pointer ptr_col(const index_type j)
		{
			return derived().ptr_col(j);
		}

		LMAT_ENSURE_INLINE const_reference elem(const index_type i, const index_type j) const
		{
			return derived().elem(i, j);
		}

		LMAT_ENSURE_INLINE reference elem(const index_type i, const index_type j)
		{
			return derived().elem(i, j);
		}

		LMAT_ENSURE_INLINE const_reference operator() (const index_type i, const index_type j) const
		{
			check_subscripts_in_range(derived(), i, j);
			return elem(i, j);
		}

		LMAT_ENSURE_INLINE reference operator() (const index_type i, const index_type j)
		{
			check_subscripts_in_range(derived(), i, j);
			return elem(i, j);
		}


		LMAT_ENSURE_INLINE void require_size(const index_type m, const index_type n)
		{
			derived().require_size(m, n);
		}

	public:

		// column views

		LMAT_ENSURE_INLINE
		typename colview_map<Derived, whole>::const_type
		column(const index_type j) const
		{
			return colview_map<Derived, whole>::get(derived(), j, whole());
		}

		LMAT_ENSURE_INLINE
		typename colview_map<Derived, whole>::type
		column(const index_type j)
		{
			return colview_map<Derived, whole>::get(derived(), j, whole());
		}

		template<class Range>
		LMAT_ENSURE_INLINE
		typename colview_map<Derived, Range>::const_type
		operator()(const IRange<Range>& rgn, const index_t j) const
		{
			return colview_map<Derived, Range>::get(derived(), j, rgn.derived());
		}

		template<class Range>
		LMAT_ENSURE_INLINE
		typename colview_map<Derived, Range>::type
		operator()(const IRange<Range>& rgn, const index_t j)
		{
			return colview_map<Derived, Range>::get(derived(), j, rgn.derived());
		}

		// row views

		LMAT_ENSURE_INLINE
		typename rowview_map<Derived, whole>::const_type
		row(const index_type i) const
		{
			return rowview_map<Derived, whole>::get(derived(), i, whole());
		}

		LMAT_ENSURE_INLINE
		typename rowview_map<Derived, whole>::type
		row(const index_type i)
		{
			return rowview_map<Derived, whole>::get(derived(), i, whole());
		}

		template<class Range>
		LMAT_ENSURE_INLINE
		typename rowview_map<Derived, Range>::const_type
		operator()(const index_t i, const IRange<Range>& rgn) const
		{
			return rowview_map<Derived, Range>::get(derived(), i, rgn.derived());
		}

		template<class Range>
		LMAT_ENSURE_INLINE
		typename rowview_map<Derived, Range>::type
		operator()(const index_t i, const IRange<Range>& rgn)
		{
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


	}; // end class IDenseMatrixBlock


	/********************************************
	 *
	 *  matrix assignment
	 *
	 ********************************************/

	template<class SExpr, class DMat>
	struct matrix_eval_verifier
	{
		static const bool value =
				is_dense_mat<DMat>::value &&
				!is_readonly_mat<DMat>::value &&
				is_mat_xpr<SExpr>::value &&
				has_same_domain<SExpr, DMat>::value &&
				has_same_value_type<SExpr, DMat>::value &&
				has_compatible_ct_size<DMat, SExpr>::value;
	};

	template<typename T, class SExpr, class DMat>
	LMAT_ENSURE_INLINE
	inline typename enable_if<matrix_eval_verifier<SExpr, DMat>, void>::type
	default_evaluate(const IMatrixXpr<SExpr, T>& sexpr, IDenseMatrix<DMat, T>& dmat)
	{
		typedef typename default_matrix_eval_policy<SExpr, DMat>::type policy_t;
		evaluate(sexpr.derived(), dmat.derived(), policy_t());
	}

	template<typename T, class LMat, class RExpr>
	LMAT_ENSURE_INLINE
	inline typename enable_if<matrix_eval_verifier<RExpr, LMat>, void>::type
	default_assign(IDenseMatrix<LMat, T>& lhs, const IMatrixXpr<RExpr, T>& rhs)
	{
		lhs.require_size(rhs.nrows(), rhs.ncolumns());
		default_evaluate(rhs, lhs);
	}

}

#endif /* MATRIX_CONCEPTS_H_ */


