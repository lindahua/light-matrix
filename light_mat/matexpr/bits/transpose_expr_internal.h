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

namespace lmat { namespace detail {

	struct dense_col2row_transpose_scheme
	{
		template<typename Arg_HP, class Arg, class DMat>
		LMAT_ENSURE_INLINE
		void evaluate(const transpose_expr<Arg_HP, Arg>& expr, DMat& dmat)
		{
			const Arg& arg = expr.arg();
			const index_t n = dmat.ncolumns();

			detail::vec_transpose(n, arg.ptr_data(), arg.row_stride(),
					dmat.ptr_data(), dmat.col_stride());
		}
	};

	struct dense_row2col_transpose_scheme
	{
		template<typename Arg_HP, class Arg, class DMat>
		LMAT_ENSURE_INLINE
		void evaluate(const transpose_expr<Arg_HP, Arg>& expr, DMat& dmat)
		{
			const Arg& arg = expr.arg();
			const index_t n = dmat.nrows();

			detail::vec_transpose(n, arg.ptr_data(), arg.col_stride(),
					dmat.ptr_data(), dmat.row_stride());
		}
	};

	struct dense_mat_transpose_scheme
	{
		template<typename Arg_HP, class Arg, class DMat>
		LMAT_ENSURE_INLINE
		void evaluate(const transpose_expr<Arg_HP, Arg>& expr, DMat& dmat)
		{
			direct_transpose(expr.arg(), dmat);
		}
	};

	struct xpr_col2row_transpose_scheme
	{
		template<typename Arg_HP, class Arg, class DMat>
		LMAT_ENSURE_INLINE
		void evaluate(const transpose_expr<Arg_HP, Arg>& expr, DMat& dmat)
		{
			const int L = binary_ctdim<ct_rows<Arg>::value, ct_cols<DMat>::value>::value;
			typedef typename matrix_traits<DMat>::value_type T;

			const index_t n = dmat.ncolumns();
			const index_t step = dmat.col_stride();

			if (step == 1)
			{
				ref_col<T, L> dcol(dmat.ptr_data(), n);
				default_evaluate(expr.arg(), dcol);
			}
			else
			{
				step_col<T, L> dcol(dmat.ptr_data(), n, step);
				default_evaluate(expr.arg(), dcol);
			}
		}
	};

	struct xpr_row2col_transpose_scheme
	{
		template<typename Arg_HP, class Arg, class DMat>
		LMAT_ENSURE_INLINE
		void evaluate(const transpose_expr<Arg_HP, Arg>& expr, DMat& dmat)
		{
			const int L = binary_ctdim<ct_cols<Arg>::value, ct_rows<DMat>::value>::value;
			typedef typename matrix_traits<DMat>::value_type T;

			const index_t n = dmat.nrows();
			const index_t step = dmat.row_stride();

			if (step == 1)
			{
				ref_row<T, L> drow(dmat.ptr_data(), n);
				default_evaluate(expr.arg(), drow);
			}
			else
			{
				step_row<T, L> drow(dmat.ptr_data(), n, step);
				default_evaluate(expr.arg(), drow);
			}
		}
	};

	struct xpr_transpose_scheme
	{
		template<typename Arg_HP, class Arg, class DMat>
		LMAT_ENSURE_INLINE
		void evaluate(const transpose_expr<Arg_HP, Arg>& expr, DMat& dmat)
		{
			direct_transpose(eval(expr.arg()), dmat);
		}
	};


	template<class Arg>
	struct transpose_scheme_map
	{
		typedef typename
				if_<is_dense_mat<Arg>,
					// is dense matrix
					typename
					if_<ct_is_col<Arg>,
						dense_col2row_transpose_scheme,
						typename
						if_<ct_is_row<Arg>,
							dense_row2col_transpose_scheme,
							dense_mat_transpose_scheme
						>::type
					>::type,
					// is general expression
					typename
					if_<ct_is_col<Arg>,
						xpr_col2row_transpose_scheme,
						typename
						if_<ct_is_row<Arg>,
							xpr_row2col_transpose_scheme,
							xpr_transpose_scheme
						>::type
					>::type
				>::type type;
	};

} }

#endif 
