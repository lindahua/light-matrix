/**
 * @file scalar_expr.h
 *
 * @brief Scalar expression
 *
 * @author Dahua Lin
 */

#ifndef LIGHTMAT_SCALAR_EXPR_H_
#define LIGHTMAT_SCALAR_EXPR_H_

#include <light_mat/matrix/matrix_fwd.h>

namespace lmat
{
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
		static const int ct_num_rows = 0;
		static const int ct_num_cols = 0;

		static const bool is_readonly = true;
		typedef cpu_domain domain;
	};
}

#endif /* SCALAR_EXPR_H_ */
