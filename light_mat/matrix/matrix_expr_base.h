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


// Useful macros

#define LMAT_DECLARE_UNARY_TYPE_MAP( map_prefix, functor_name ) \
	template<class Arg> \
	struct map_prefix##_expr_map { \
		typedef typename unary_ewise_expr_map< \
				functor_name<typename matrix_traits<Arg>::value_type>, Arg>::type type; };

#define LMAT_DECLARE_BINARY_TYPE_MAP( map_prefix, functor_name ) \
	template<class LArg, class RArg> \
	struct map_prefix##_expr_map { \
		typedef typename binary_ewise_expr_map< \
				functor_name<typename binary_value_type<LArg, RArg>::type>, LArg, RArg>::type type; };

#define LMAT_DECLARE_BINARY_TYPE_MAP_EX( map_prefix, functor_name ) \
	template<class LArg, class RArg> \
	struct map_prefix##_expr_map { \
		typedef typename binary_ewise_expr_map< \
				functor_name<typename binary_value_type<LArg, RArg>::type>, LArg, RArg>::type type; }; \
	template<class RArg> \
	struct map_prefix##_fix1_expr_map { \
		typedef typename binary_fix1_ewise_expr_map< \
				functor_name<typename matrix_traits<RArg>::value_type>, RArg>::type type; }; \
	template<class LArg> \
	struct map_prefix##_fix2_expr_map { \
		typedef typename binary_fix2_ewise_expr_map< \
				functor_name<typename matrix_traits<LArg>::value_type>, LArg>::type type; };


#define LMAT_DEFINE_UNARY_MATFUNCTION( map_prefix, matfun_name, functor_name ) \
	template<typename T, class Arg> \
	LMAT_ENSURE_INLINE \
	inline typename map_prefix##_expr_map<Arg>::type \
	matfun_name(const IMatrixXpr<Arg, T>& A) { \
		return ewise(functor_name<T>(), A.derived()); }


#define LMAT_DEFINE_BINARY_MATFUNCTION( map_prefix, matfun_name, functor_name ) \
	template<typename T, class LArg, class RArg> \
	LMAT_ENSURE_INLINE \
	inline typename map_prefix##_expr_map<LArg, RArg>::type \
	matfun_name(const IMatrixXpr<LArg, T>& A, const IMatrixXpr<RArg, T>& B) { \
		return ewise(functor_name<T>(), A.derived(), B.derived()); }


#define LMAT_DEFINE_BINARY_MATFUNCTION_EX( map_prefix, matfun_name, functor_name ) \
	LMAT_DEFINE_BINARY_MATFUNCTION( map_prefix, matfun_name, functor_name ) \
	template<typename T, class RArg> \
	LMAT_ENSURE_INLINE \
	inline typename map_prefix##_fix1_expr_map<RArg>::type \
	matfun_name(const T& a, const IMatrixXpr<RArg, T>& B) { \
		return ewise(functor_name<T>(), a, B.derived()); } \
	template<typename T, class LArg> \
	LMAT_ENSURE_INLINE \
	inline typename map_prefix##_fix2_expr_map<LArg>::type \
	matfun_name(const IMatrixXpr<LArg, T>& A, const T& b) { \
		return ewise(functor_name<T>(), A.derived(), b); }

#endif






