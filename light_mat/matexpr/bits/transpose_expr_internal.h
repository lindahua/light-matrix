/**
 * @file transpose_expr_internal.h
 *
 * Internal implementation of transpose expression
 * 
 * @author Dahua Lin 
 */

#ifdef _MSC_VER
#pragma once
#endif

#ifndef LIGHTMAT_TRANSPOSE_EXPR_INTERNAL_H_
#define LIGHTMAT_TRANSPOSE_EXPR_INTERNAL_H_

#include <light_mat/matrix/matrix_classes.h>

namespace lmat { namespace internal {

	struct dense_col2row_transpose_scheme
	{
		template<class QArg, class DMat>
		LMAT_ENSURE_INLINE
		void _eval(const transpose_expr<QArg>& expr, DMat& dmat)
		{
			const typename QArg::argument& arg = expr.arg();
			const index_t n = dmat.ncolumns();

			internal::vec_transpose(n, arg.ptr_data(), arg.row_stride(),
					dmat.ptr_data(), dmat.col_stride());
		}
	};

	struct dense_row2col_transpose_scheme
	{
		template<class QArg, class DMat>
		LMAT_ENSURE_INLINE
		void _eval(const transpose_expr<QArg>& expr, DMat& dmat)
		{
			const typename QArg::argument& arg = expr.arg();
			const index_t n = dmat.nrows();

			internal::vec_transpose(n, arg.ptr_data(), arg.col_stride(),
					dmat.ptr_data(), dmat.row_stride());
		}
	};

	struct dense_mat_transpose_scheme
	{
		template<class QArg, class DMat>
		LMAT_ENSURE_INLINE
		void _eval(const transpose_expr<QArg>& expr, DMat& dmat)
		{
			direct_transpose(expr.arg(), dmat);
		}
	};

	struct xpr_col2row_transpose_scheme
	{
		template<class QArg, class DMat>
		LMAT_ENSURE_INLINE
		void _eval(const transpose_expr<QArg>& expr, DMat& dmat)
		{
			typedef typename QArg::argument arg_t;
			const int L = meta::common_dim<
					LMAT_INTLIST_2( meta::nrows<arg_t>::value, meta::ncols<DMat>::value ) >::value;
			typedef typename matrix_traits<DMat>::value_type T;

			const index_t n = dmat.ncolumns();
			const index_t step = dmat.col_stride();

			if (step == 1)
			{
				ref_col<T, L> dcol(dmat.ptr_data(), n);
				evaluate(expr.arg(), dcol);
			}
			else
			{
				step_col<T, L> dcol(dmat.ptr_data(), n, step);
				evaluate(expr.arg(), dcol);
			}
		}
	};

	struct xpr_row2col_transpose_scheme
	{
		template<class QArg, class DMat>
		LMAT_ENSURE_INLINE
		void _eval(const transpose_expr<QArg>& expr, DMat& dmat)
		{
			typedef typename QArg::argument arg_t;
			const int L = meta::common_dim<
					LMAT_INTLIST_2( meta::ncols<arg_t>::value, meta::nrows<DMat>::value ) >::value;
			typedef typename matrix_traits<DMat>::value_type T;

			const index_t n = dmat.nrows();
			const index_t step = dmat.row_stride();

			if (step == 1)
			{
				ref_row<T, L> drow(dmat.ptr_data(), n);
				evaluate(expr.arg(), drow);
			}
			else
			{
				step_row<T, L> drow(dmat.ptr_data(), n, step);
				evaluate(expr.arg(), drow);
			}
		}
	};

	struct xpr_transpose_scheme
	{
		template<class QArg, class DMat>
		LMAT_ENSURE_INLINE
		void _eval(const transpose_expr<QArg>& expr, DMat& dmat)
		{
			direct_transpose(eval(expr.arg()), dmat);
		}
	};


	template<class Arg>
	struct transpose_scheme_map
	{
		typedef typename
				meta::if_<meta::is_dense_mat<Arg>,
					// is dense matrix
					typename
					meta::if_<meta::is_col<Arg>,
						dense_col2row_transpose_scheme,
						typename
						meta::if_<meta::is_row<Arg>,
							dense_row2col_transpose_scheme,
							dense_mat_transpose_scheme
						>::type
					>::type,
					// is general expression
					typename
					meta::if_<meta::is_col<Arg>,
						xpr_col2row_transpose_scheme,
						typename
						meta::if_<meta::is_row<Arg>,
							xpr_row2col_transpose_scheme,
							xpr_transpose_scheme
						>::type
					>::type
				>::type type;
	};

} }

#endif 
