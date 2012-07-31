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
		typename transpose_expr_map<Derived>::type
		trans() const
		{
			return transpose_expr_map<Derived>::get(derived());
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
		typename unary_expr_map<colview_spec<whole>, Derived>::const_type
		column(const index_type j) const
		{
			return unary_expr_map<colview_spec<whole>, Derived>::get(
					colview_spec<whole>(j), derived());
		}

		LMAT_ENSURE_INLINE
		typename unary_expr_map<colview_spec<whole>, Derived>::type
		column(const index_type j)
		{
			return unary_expr_map<colview_spec<whole>, Derived>::get(
					colview_spec<whole>(j), derived());
		}

		template<class Range>
		LMAT_ENSURE_INLINE
		typename unary_expr_map<colview_spec<Range>, Derived>::const_type
		operator()(const IRange<Range>& rgn, const index_t j) const
		{
			return unary_expr_map<colview_spec<Range>, Derived>::get(
					colview_spec<Range>(j, rgn), derived());
		}

		template<class Range>
		LMAT_ENSURE_INLINE
		typename unary_expr_map<colview_spec<Range>, Derived>::type
		operator()(const IRange<Range>& rgn, const index_t j)
		{
			return unary_expr_map<colview_spec<Range>, Derived>::get(
					colview_spec<Range>(j, rgn), derived());
		}

		// row views

		LMAT_ENSURE_INLINE
		typename unary_expr_map<rowview_spec<whole>, Derived>::const_type
		row(const index_type i) const
		{
			return unary_expr_map<rowview_spec<whole>, Derived>::get(
					rowview_spec<whole>(i), derived());
		}

		LMAT_ENSURE_INLINE
		typename unary_expr_map<rowview_spec<whole>, Derived>::type
		row(const index_type i)
		{
			return unary_expr_map<rowview_spec<whole>, Derived>::get(
					rowview_spec<whole>(i), derived());
		}

		template<class Range>
		LMAT_ENSURE_INLINE
		typename unary_expr_map<rowview_spec<Range>, Derived>::const_type
		operator()(const index_t i, const IRange<Range>& rgn) const
		{
			return unary_expr_map<rowview_spec<Range>, Derived>::get(
					rowview_spec<Range>(i, rgn), derived());
		}

		template<class Range>
		LMAT_ENSURE_INLINE
		typename unary_expr_map<rowview_spec<Range>, Derived>::type
		operator()(const index_t i, const IRange<Range>& rgn)
		{
			return unary_expr_map<rowview_spec<Range>, Derived>::get(
					rowview_spec<Range>(i, rgn), derived());
		}

		// sub-views

		template<class Range0, class Range1>
		LMAT_ENSURE_INLINE
		typename unary_expr_map<matview_spec<Range0, Range1>, Derived>::const_type
		operator()(const IRange<Range0>& row_rgn, const IRange<Range1>& col_rgn) const
		{
			return unary_expr_map<matview_spec<Range0, Range1>, Derived>::get(
					matview_spec<Range0, Range1>(row_rgn, col_rgn), derived());
		}

		template<class Range0, class Range1>
		LMAT_ENSURE_INLINE
		typename unary_expr_map<matview_spec<Range0, Range1>, Derived>::const_type
		operator()(const IRange<Range0>& row_rgn, const IRange<Range1>& col_rgn)
		{
			return unary_expr_map<matview_spec<Range0, Range1>, Derived>::get(
					matview_spec<Range0, Range1>(row_rgn, col_rgn), derived());
		}

	}; // end class IDenseMatrixBlock


	/********************************************
	 *
	 *  matrix assignment
	 *
	 ********************************************/

	template<class SExpr, class DMat>
	class mat_eval_verifier
	{
		static const bool value =
				is_dense_mat<DMat>::value &&
				!is_readonly_mat<DMat>::value &&
				is_mat_xpr<SExpr>::value &&
				is_same<
					typename matrix_traits<DMat>::value_type,
					typename matrix_traits<SExpr>::value_type
				>::value &&
				has_compatible_ct_size<DMat, SExpr>::value;
	};

	template<typename T, class SExpr, class DMat>
	LMAT_ENSURE_INLINE
	inline typename enable_if<mat_eval_verifier<SExpr, DMat>, void>::type
	default_evaluate(const IMatrixXpr<SExpr, T>& sexpr, IDenseMatrix<DMat, T>& dmat)
	{
		typedef typename mateval_ctx<default_evaldom, SExpr, DMat>::type ctx_t;
		evaluate(sexpr.derived(), dmat.derived(), ctx_t());
	}

	template<typename T, class LMat, class RExpr>
	LMAT_ENSURE_INLINE
	inline typename enable_if<mat_eval_verifier<RExpr, LMat>, void>::type
	default_assign(IDenseMatrix<LMat, T>& lhs, const IMatrixXpr<RExpr, T>& rhs)
	{
		lhs.require_size(rhs.nrows(), rhs.ncolumns());
		default_evaluate(rhs, lhs);
	}

	template<typename S, typename T, class SExpr, class DMat>
	LMAT_ENSURE_INLINE
	inline typename enable_if<is_implicitly_convertible<S, T>, void>::type
	default_assign_with_implicit_cast(const IMatrixXpr<SExpr, S>& src, IDenseMatrix<DMat, T>& dst)
	{
		default_assign(dst, make_expr(type_converter<S, T>(), src.derived()));
	}

}

#endif /* MATRIX_CONCEPTS_H_ */


