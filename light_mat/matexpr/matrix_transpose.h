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

#include <light_mat/matexpr/matrix_ewise_eval.h>
#include "bits/direct_transpose_impl.h"
#include "bits/transpose_expr_internal.h"

namespace lmat
{

	/********************************************
	 *
	 *  Direct transpose
	 *
	 ********************************************/

	template<typename T, class SMat, class DMat>
	LMAT_ENSURE_INLINE
	inline void direct_transpose(const IDenseMatrix<SMat, T>& src, IDenseMatrix<DMat, T>& dst)
	{
		const index_t m = src.nrows();
		const index_t n = src.ncolumns();

		detail::transpose(m, n,
				src.ptr_data(), src.row_stride(), src.col_stride(),
				dst.ptr_data(), dst.row_stride(), dst.col_stride());
	}


	/********************************************
	 *
	 *  Transpose expression
	 *
	 ********************************************/

	template<typename Arg_HP, class Arg> class transpose_expr;

	template<typename Arg_HP, class Arg>
	struct matrix_traits<transpose_expr<Arg_HP, Arg> >
	{
		static const int num_dimensions = 2;
		static const int ct_num_rows = ct_cols<Arg>::value;
		static const int ct_num_cols = ct_rows<Arg>::value;

		static const bool is_readonly = true;

		typedef typename matrix_traits<Arg>::value_type value_type;
		typedef typename matrix_traits<Arg>::domain domain;
	};


	template<typename Arg_HP, class Arg>
	class transpose_expr :
			public unary_expr_base<Arg_HP, Arg>,
			public IMatrixXpr<transpose_expr<Arg_HP, Arg>, typename matrix_traits<Arg>::value_type>
	{
#ifdef LMAT_USE_STATIC_ASSERT
		static_assert(is_mat_xpr<Arg>::value, "Arg must be a matrix expression class.");
#endif

	public:
		typedef unary_expr_base<Arg_HP, Arg> base_t;

		LMAT_ENSURE_INLINE
		transpose_expr(const arg_forwarder<Arg_HP, Arg>& arg_fwd)
		: base_t(arg_fwd) { }

		LMAT_ENSURE_INLINE index_t nelems() const
		{
			return this->arg().nelems();
		}

		LMAT_ENSURE_INLINE index_t nrows() const
		{
			return this->arg().ncolumns();
		}

		LMAT_ENSURE_INLINE index_t ncolumns() const
		{
			return this->arg().nrows();
		}

	};


	/********************************************
	 *
	 *  Expression mapping
	 *
	 ********************************************/

	// specifier

	struct transpose_t { };

	// verifier

	template<class Arg>
	struct unary_expr_verifier<transpose_t, Arg>
	{
		static const bool value = is_mat_xpr<Arg>::value;
	};

	// maps

	template<typename Arg_HP, class Arg>
	struct unary_expr_map<transpose_t, Arg_HP, Arg>
	{
		typedef transpose_expr<Arg_HP, Arg> type;

		LMAT_ENSURE_INLINE
		static type get(transpose_t, const arg_forwarder<Arg_HP, Arg>& arg_fwd)
		{
			return type(arg_fwd);
		}
	};

	template<class Arg, typename T>
	LMAT_ENSURE_INLINE
	inline typename unary_expr_map<transpose_t, ref_arg_t, Arg>::type
	trans(const IMatrixXpr<Arg, T>& arg)
	{
		return unary_expr_map<transpose_t, ref_arg_t, Arg>::get(
				transpose_t(), ref_arg(arg.derived()));
	}


	/********************************************
	 *
	 *  Evaluation
	 *
	 ********************************************/

	template<typename T, typename Arg_HP, class Arg, class DMat>
	LMAT_ENSURE_INLINE
	inline typename detail::transpose_scheme_map<Arg>::type
	get_default_eval_scheme(
			const transpose_expr<Arg_HP, Arg>& sexpr,
			IDenseMatrix<DMat, T>& dmat)
	{
		typedef typename detail::transpose_scheme_map<Arg>::type scheme_t;
		return scheme_t();
	}

}

#endif /* MATRIX_TRANSPOSE_H_ */
