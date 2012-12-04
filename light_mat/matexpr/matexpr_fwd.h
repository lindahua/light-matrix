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

	template<typename Tag, class QArgList> class ewise_expr;

	template<class QArg> class transpose_expr;

	template<class QArg, int N> class repeat_col_expr;
	template<class QArg, int M> class repeat_row_expr;

	struct rowwise { };
	struct colwise { };

	template<typename Tag, typename AlongDim, class QArgList> class partial_reduce_expr;
}

#endif /* MATEXPR_FWD_H_ */



