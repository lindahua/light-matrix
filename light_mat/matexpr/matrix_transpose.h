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

		internal::transpose(m, n,
				src.ptr_data(), src.row_stride(), src.col_stride(),
				dst.ptr_data(), dst.row_stride(), dst.col_stride());
	}


	/********************************************
	 *
	 *  Transpose expression
	 *
	 ********************************************/

	template<class QArg> class transpose_expr;

	template<class QArg>
	struct matrix_traits<transpose_expr<QArg> >
	{
		typedef typename QArg::argument arg_t;

		static const int num_dimensions = 2;
		static const int ct_num_rows = meta::ncols<arg_t>::value;
		static const int ct_num_cols = meta::nrows<arg_t>::value;

		static const bool is_readonly = true;

		typedef typename matrix_traits<arg_t>::value_type value_type;
		typedef typename matrix_traits<arg_t>::domain domain;
	};


	template<class QArg>
	class transpose_expr :
		public IMatrixXpr<transpose_expr<QArg>,
			typename matrix_traits<typename QArg::argument>::value_type>
	{
	public:

		LMAT_ENSURE_INLINE
		transpose_expr(const typename QArg::forwarder& arg_fwd)
		: arg_holder(arg_fwd) { }

		LMAT_ENSURE_INLINE
		const typename QArg::argument& arg() const
		{
			return arg_holder.get();
		}

		LMAT_ENSURE_INLINE index_t nelems() const
		{
			return arg().nelems();
		}

		LMAT_ENSURE_INLINE index_t nrows() const
		{
			return arg().ncolumns();
		}

		LMAT_ENSURE_INLINE index_t ncolumns() const
		{
			return arg().nrows();
		}

	private:
		typename QArg::holder arg_holder;

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
	struct expr_verifier<transpose_t, LMAT_TYPELIST_1(Arg) >
	{
		static const bool value = meta::is_mat_xpr<Arg>::value;
	};

	// maps

	template<class QArg>
	struct expr_map<transpose_t, LMAT_TYPELIST_1(QArg) >
	{
		typedef transpose_expr<QArg> type;

		LMAT_ENSURE_INLINE
		static type get(transpose_t,
				const tied_forwarder< LMAT_TYPELIST_1(QArg) >& tfwd)
		{
			return type(tfwd.arg1_fwd);
		}
	};

	template<class Arg, typename T>
	LMAT_ENSURE_INLINE
	inline typename expr_map<transpose_t, LMAT_TYPELIST_1( CRefArg<Arg> ) >::type
	trans(const IMatrixXpr<Arg, T>& arg)
	{
		typedef tied_forwarder< LMAT_TYPELIST_1( CRefArg<Arg> ) > tfwd_t;

		return make_expr(
				transpose_t(), tfwd_t(arg.derived()) );
	}


	/********************************************
	 *
	 *  Evaluation
	 *
	 ********************************************/

	template<class QArg, class DMat>
	LMAT_ENSURE_INLINE
	void evaluate(
			const transpose_expr<QArg>& sexpr,
			IDenseMatrix<DMat, typename matrix_traits< typename QArg::argument >::value_type >& dmat)
	{
		typedef typename internal::transpose_scheme_map<typename QArg::argument>::type scheme_t;
		scheme_t()._eval( sexpr.derived(), dmat.derived() );
	}

}

#endif /* MATRIX_TRANSPOSE_H_ */
