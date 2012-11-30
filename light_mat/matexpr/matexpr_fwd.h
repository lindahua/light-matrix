/**
 * @file matexpr_fwd.h
 *
 * @brief Forward declarations for matrix expression
 *
 * @author Dahua Lin
 */

#ifndef LIGHTMAT_MATEXPR_FWD_H_
#define LIGHTMAT_MATEXPR_FWD_H_

#include <light_mat/matrix/matrix_base.h>
#include <light_mat/common/expr_base.h>

namespace lmat
{
	// forward declaration of matrix expressions

	template<typename Op, typename Arg_HP, class Arg> class unary_ewise_expr;

	template<typename Op,
		typename Arg1_HP, class Arg1,
		typename Arg2_HP, class Arg2> class binary_ewise_expr;

	template<typename Arg_HP, class Arg> class transpose_expr;

	template<typename Arg_HP, class Arg, int N> class repeat_col_expr;
	template<typename Arg_HP, class Arg, int M> class repeat_row_expr;

	struct rowwise { };
	struct colwise { };

	template<typename Op, typename Arg_HP, class Arg> class colwise_reduce_expr;
	template<typename Op, typename Arg_HP, class Arg> class rowwise_reduce_expr;

	// scalar expression

	template<typename T>
	struct scalar_expr
	{
		const T value;

		LMAT_ENSURE_INLINE
		scalar_expr(const T& v)
		: value(v) { }
	};

	template<typename T>
	LMAT_ENSURE_INLINE
	inline scalar_expr<T> make_scalar(const T& v)
	{
		return v;
	}

	template<typename T>
	struct matrix_traits<scalar_expr<T> >
	{
		typedef T value_type;
	};

}

#endif /* MATEXPR_FWD_H_ */



