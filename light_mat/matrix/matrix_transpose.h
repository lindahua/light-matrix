/**
 * @file matrix_transpose.h
 *
 * Matrix transposition
 *
 * @author Dahua Lin
 */

#ifdef _MSC_VER
#pragma once
#endif

#ifndef LIGHTMAT_MATRIX_TRANSPOSE_H_
#define LIGHTMAT_MATRIX_TRANSPOSE_H_

#include "bits/matrix_transpose_internal.h"

namespace lmat
{

	// traits

	template<class Expr>
	struct matrix_traits<transpose_expr<Expr> >
	{
		static const int num_dimensions = 2;
		static const int compile_time_num_rows = ct_cols<Expr>::value;
		static const int compile_time_num_cols = ct_rows<Expr>::value;

		static const bool is_readonly = true;
		typedef typename matrix_traits<Expr>::value_type value_type;
	};

	template<class Expr>
	struct ct_has_continuous_layout<transpose_expr<Expr> >
	{
		typedef typename detail::matrix_transpose_base_map<Expr>::type base_t;
		static const bool value = base_t::has_continuous_layout;
	};

	template<class Expr>
	struct is_base_aligned<transpose_expr<Expr> >
	{
		typedef typename detail::matrix_transpose_base_map<Expr>::type base_t;
		static const bool value = base_t::is_base_aligned;
	};

	template<class Expr>
	struct is_percol_aligned<transpose_expr<Expr> >
	{
		typedef typename detail::matrix_transpose_base_map<Expr>::type base_t;
		static const bool value = base_t::is_percol_aligned;
	};

	template<class Expr>
	struct is_linear_accessible<transpose_expr<Expr> >
	{
		typedef typename detail::matrix_transpose_base_map<Expr>::type base_t;
		static const bool value = base_t::is_linear_accessible;
	};

	// class

	template<class Expr>
	class transpose_expr
	: public detail::matrix_transpose_base_map<Expr>::type
	{
	};

	// evaluation

	template<class Expr, class DMat>
	LMAT_ENSURE_INLINE
	void evaluate(const transpose_expr<Expr>& s,
			IDenseMatrix<DMat, typename matrix_traits<Expr>::value_type>& dst)
	{
		s.eval_to(dst.derived());
	}

}

#endif /* MATRIX_TRANSPOSE_H_ */
