/**
 * @file matrix_expr_base.h
 *
 * The basis of all matrix expressions
 *
 * @author Dahua Lin
 */

#ifdef _MSC_VER
#pragma once
#endif

#ifndef LIGHTMAT_MATRIX_EXPR_BASE_H_
#define LIGHTMAT_MATRIX_EXPR_BASE_H_

#include <light_mat/matrix/matrix_properties.h>
#include <light_mat/math/functor_base.h>
#include <light_mat/core/expr_base.h>

// Useful macros


#define LMAT_DEFINE_UNARY_MATFUNCTION( matfun_name, functor_name ) \
	template<typename T, class Arg> \
	LMAT_ENSURE_INLINE \
	inline typename unary_expr_map<ewise_t<functor_name<T> >, ref_arg_t, Arg>::type \
	matfun_name(const IMatrixXpr<Arg, T>& A) { \
		return ewise(functor_name<T>(), A.derived()); }


#define LMAT_DEFINE_BINARY_MATFUNCTION( matfun_name, functor_name ) \
	template<typename T, class LArg, class RArg> \
	LMAT_ENSURE_INLINE \
	inline typename binary_expr_map<ewise_t<functor_name<T> >, ref_arg_t, LArg, ref_arg_t, RArg>::type \
	matfun_name(const IMatrixXpr<LArg, T>& A, const IMatrixXpr<RArg, T>& B) { \
		return ewise(functor_name<T>(), A.derived(), B.derived()); }


#define LMAT_DEFINE_BINARY_MATFUNCTION_EX( matfun_name, functor_name ) \
	LMAT_DEFINE_BINARY_MATFUNCTION( matfun_name, functor_name ) \
	template<typename T, class RArg> \
	LMAT_ENSURE_INLINE \
	inline typename binary_fix1_ewise_expr_map< functor_name<T>, ref_arg_t, RArg>::type \
	matfun_name(const T& a, const IMatrixXpr<RArg, T>& B) { \
		return ewise(functor_name<T>(), a, B.derived()); } \
	template<typename T, class LArg> \
	LMAT_ENSURE_INLINE \
	inline typename binary_fix2_ewise_expr_map< functor_name<T>, ref_arg_t, LArg>::type \
	matfun_name(const IMatrixXpr<LArg, T>& A, const T& b) { \
		return ewise(functor_name<T>(), A.derived(), b); }

#endif






