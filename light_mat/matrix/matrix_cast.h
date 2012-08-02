/**
 * @file matrix_cast.h
 *
 * Value-type casting for matrices
 *
 * @author Dahua Lin
 */

#ifdef _MSC_VER
#pragma once
#endif

#ifndef LIGHTMAT_MATRIX_CAST_H_
#define LIGHTMAT_MATRIX_CAST_H_

#include <light_mat/matrix/matrix_ewise_expr.h>

namespace lmat
{


	template<class SMat, typename S, typename T>
	LMAT_ENSURE_INLINE
	inline typename unary_expr_map<
		ewise_t<type_converter<S, T> >,
		ref_arg_t, SMat >::type
	cast(const IMatrixXpr<SMat, S>& sexpr, type<T> )
	{
		return ewise(type_converter<S, T>(), sexpr.derived());
	}


	template<class SMat, typename S>
	LMAT_ENSURE_INLINE
	inline typename unary_expr_map<
		ewise_t<type_converter<S, bool> >,
		ref_arg_t, SMat >::type
	to_bool(const IMatrixXpr<SMat, S>& sexpr)
	{
		return ewise(type_converter<S, bool>(), sexpr.derived());
	}


	template<class SMat, typename S>
	LMAT_ENSURE_INLINE
	inline typename unary_expr_map<
		ewise_t<type_converter<S, float> >,
		ref_arg_t, SMat >::type
	to_f32(const IMatrixXpr<SMat, S>& sexpr)
	{
		return ewise(type_converter<S, float>(), sexpr.derived());
	}

	template<class SMat, typename S>
	LMAT_ENSURE_INLINE
	inline typename unary_expr_map<
		ewise_t<type_converter<S, double> >,
		ref_arg_t, SMat >::type
	to_f64(const IMatrixXpr<SMat, S>& sexpr)
	{
		return ewise(type_converter<S, double>(), sexpr.derived());
	}

	template<class SMat, typename S>
	LMAT_ENSURE_INLINE
	inline typename unary_expr_map<
		ewise_t<type_converter<S, int32_t> >,
		ref_arg_t, SMat >::type
	to_i32(const IMatrixXpr<SMat, S>& sexpr)
	{
		return ewise(type_converter<S, int32_t>(), sexpr.derived());
	}

	template<class SMat, typename S>
	LMAT_ENSURE_INLINE
	inline typename unary_expr_map<
		ewise_t<type_converter<S, uint32_t> >,
		ref_arg_t, SMat >::type
	to_u32(const IMatrixXpr<SMat, S>& sexpr)
	{
		return ewise(type_converter<S, uint32_t>(), sexpr.derived());
	}
}

#endif /* MATRIX_CAST_H_ */
