/**
 * @file mat_pred.h
 *
 * @brief Predicate functions on matrices
 *
 * @author Dahua Lin
 */

#ifndef LIGHTMAT_MAT_PRED_H_
#define LIGHTMAT_MAT_PRED_H_

#include <light_mat/mateval/mat_cast.h>

#define LMAT_DEFINE_MAT_PRED_FUN_G1( FunName, FTag ) \
	template<typename T, class X> \
	LMAT_ENSURE_INLINE \
	inline map_expr<FTag, X> \
	FunName (const IMatrixXpr<X, T>& x) \
	{ return make_map_expr(FTag(), x); }

#define LMAT_DEFINE_MAT_PRED_FUN_G2( FunName, FTag ) \
	template<typename T, class X, class Y> \
	LMAT_ENSURE_INLINE \
	inline map_expr<FTag, X, Y> \
	FunName (const IMatrixXpr<X, T>& x, const IMatrixXpr<Y, T>& y) \
	{ return make_map_expr(FTag(), x, y); } \
	template<typename T, class X> \
	LMAT_ENSURE_INLINE \
	inline map_expr<FTag, X, T> \
	FunName (const IMatrixXpr<X, T>& x, const T& y) \
	{ return make_map_expr_fix2(FTag(), x, y); } \
	template<typename T, class Y> \
	LMAT_ENSURE_INLINE \
	inline map_expr<FTag, T, Y> \
	FunName (const T& x, const IMatrixXpr<Y, T>& y) \
	{ return make_map_expr_fix1(FTag(), x, y); }


namespace lmat
{
	// comparison

	LMAT_DEFINE_MAT_PRED_FUN_G2( operator ==, eq_ )
	LMAT_DEFINE_MAT_PRED_FUN_G2( operator !=, ne_ )
	LMAT_DEFINE_MAT_PRED_FUN_G2( operator >=, ge_ )
	LMAT_DEFINE_MAT_PRED_FUN_G2( operator >,  gt_ )
	LMAT_DEFINE_MAT_PRED_FUN_G2( operator <=, le_ )
	LMAT_DEFINE_MAT_PRED_FUN_G2( operator <,  lt_ )

#ifdef LMAT_HAS_CXX11_MATH

	// classification

	LMAT_DEFINE_MAT_PRED_FUN_G1( signbit,  signbit_ )
	LMAT_DEFINE_MAT_PRED_FUN_G1( isfinite, isfinite_ )
	LMAT_DEFINE_MAT_PRED_FUN_G1( isinf,    isinf_ )
	LMAT_DEFINE_MAT_PRED_FUN_G1( isnan,    isnan_ )

#endif

}

#endif /* MAT_PRED_H_ */
