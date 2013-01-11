/**
 * @file math_functors.h
 *
 * Math functors
 * 
 * @author Dahua Lin 
 */

#ifdef _MSC_VER
#pragma once
#endif

#ifndef LIGHTMAT_MATH_FUNCTORS_H_
#define LIGHTMAT_MATH_FUNCTORS_H_

#include <light_mat/math/math_base.h>
#include <light_mat/math/special.h>

#include <light_mat/math/fun_tags.h>
#include <light_mat/math/functor_base.h>
#include <light_mat/math/simd.h>


#define _LMAT_DECL_SIMDIZABLE_S(FTag, FunT) \
		template<typename Kind> \
		struct is_simdizable<FunT<float>, Kind> : public meta::has_simd_support<FTag, float, Kind> { }; \
		template<typename Kind> \
		struct is_simdizable<FunT<double>, Kind> : public meta::has_simd_support<FTag, double, Kind> { };

#define _LMAT_DEFINE_GENERIC_MATH_FUNCTOR( Name, NA, Expr ) \
	template<typename T> \
	struct Name##_fun { \
		typedef T result_type; \
		LMAT_ENSURE_INLINE \
		T operator()(LMAT_REPEAT_ARGS_P##NA(const T& x)) const \
		{ return Expr; } }; \
	_LMAT_DECL_SIMDIZABLE_S( ftags::Name##_, Name##_fun ) \
	LMAT_DEF_TRIVIAL_SIMDIZE_MAP( Name##_fun )

#define _LMAT_DEFINE_MATH_FUNCTOR( Name, NA ) \
	_LMAT_DEFINE_GENERIC_MATH_FUNCTOR( Name, NA, math::Name(LMAT_REPEAT_ARGS_P##NA(x)) )

#define _LMAT_DEFINE_GENERIC_NUMPRED_FUNCTOR( Name, NA, Expr ) \
	template<typename T> \
	struct Name##_fun { \
		typedef typename pred_result<T>::type result_type; \
		LMAT_ENSURE_INLINE \
		result_type operator()(LMAT_REPEAT_ARGS_P##NA(const T& x)) const \
		{ return Expr; } }; \
	_LMAT_DECL_SIMDIZABLE_S( ftags::Name##_, Name##_fun ) \
	LMAT_DEF_TRIVIAL_SIMDIZE_MAP( Name##_fun )

#define _LMAT_DEFINE_NUMPRED_FUNCTOR( Name, NA ) \
	_LMAT_DEFINE_GENERIC_NUMPRED_FUNCTOR( Name, NA, math::Name(LMAT_REPEAT_ARGS_P##NA(x)) )

#define _LMAT_DEFINE_COMPARISON_FUNCTOR( Name, Expr ) \
	_LMAT_DEFINE_GENERIC_NUMPRED_FUNCTOR( Name, 2, Expr )

#define LMAT_DEFINE_LOGICAL_FUNCTOR( Name, NA, MExpr, BExpr ) \
	template<typename T> \
	struct Name##_fun { \
		typedef mask_t<T> result_type; \
		LMAT_ENSURE_INLINE \
		result_type operator()(LMAT_REPEAT_ARGS_P##NA(const mask_t<T>& x)) const { return MExpr; } }; \
	template<typename T, typename Kind> \
	struct Name##_fun<math::simd_pack<T, Kind> > { \
		typedef math::simd_bpack<T, Kind> result_type; \
		LMAT_ENSURE_INLINE \
		math::simd_bpack<T, Kind> operator()(LMAT_REPEAT_ARGS_P##NA(const LMAT_SIMD_BPACK_M(T, Kind)& x)) const { return MExpr; } }; \
	template<> \
	struct Name##_fun<bool> { \
		typedef bool result_type; \
		LMAT_ENSURE_INLINE \
		bool operator()(LMAT_REPEAT_ARGS_P##NA(bool x)) const { return BExpr; } }; \
	_LMAT_DECL_SIMDIZABLE_S( ftags::Name##_, Name##_fun ) \
	LMAT_DEF_TRIVIAL_SIMDIZE_MAP( Name##_fun )


namespace lmat {

	// Arithmetic functors

	_LMAT_DEFINE_GENERIC_MATH_FUNCTOR( add, 2, x1 + x2 )
	_LMAT_DEFINE_GENERIC_MATH_FUNCTOR( sub, 2, x1 - x2 )
	_LMAT_DEFINE_GENERIC_MATH_FUNCTOR( mul, 2, x1 * x2 )
	_LMAT_DEFINE_GENERIC_MATH_FUNCTOR( div, 2, x1 / x2 )
	_LMAT_DEFINE_GENERIC_MATH_FUNCTOR( neg, 1, -x1 )

	_LMAT_DEFINE_MATH_FUNCTOR( abs, 1 )
	_LMAT_DEFINE_MATH_FUNCTOR( sqr, 1 )
	_LMAT_DEFINE_MATH_FUNCTOR( cube, 1 )
	_LMAT_DEFINE_MATH_FUNCTOR( rcp, 1 )
	_LMAT_DEFINE_MATH_FUNCTOR( fma, 3 )

	_LMAT_DEFINE_MATH_FUNCTOR( max, 2 )
	_LMAT_DEFINE_MATH_FUNCTOR( min, 2 )
	_LMAT_DEFINE_MATH_FUNCTOR( clamp, 3 )

	// Comparison functors

	_LMAT_DEFINE_COMPARISON_FUNCTOR( eq, x1 == x2 )
	_LMAT_DEFINE_COMPARISON_FUNCTOR( ne, x1 != x2 )
	_LMAT_DEFINE_COMPARISON_FUNCTOR( ge, x1 >= x2 )
	_LMAT_DEFINE_COMPARISON_FUNCTOR( gt, x1 >  x2 )
	_LMAT_DEFINE_COMPARISON_FUNCTOR( le, x1 <= x2 )
	_LMAT_DEFINE_COMPARISON_FUNCTOR( lt, x1 <  x2 )

	// logical functors

	LMAT_DEFINE_LOGICAL_FUNCTOR( logical_not, 1, ~x1, !x1 )
	LMAT_DEFINE_LOGICAL_FUNCTOR( logical_eq,  2, x1 == x2, x1 == x2 )
	LMAT_DEFINE_LOGICAL_FUNCTOR( logical_ne,  2, x1 != x2, x1 != x2 )
	LMAT_DEFINE_LOGICAL_FUNCTOR( logical_or,  2, x1 | x2,  x1 || x2 )
	LMAT_DEFINE_LOGICAL_FUNCTOR( logical_and, 2, x1 & x2,  x1 && x2 )

	// power

	_LMAT_DEFINE_MATH_FUNCTOR( sqrt, 1 )
	_LMAT_DEFINE_MATH_FUNCTOR( rsqrt, 1 )
	_LMAT_DEFINE_MATH_FUNCTOR( pow, 2 )

	// floor & ceil

	_LMAT_DEFINE_MATH_FUNCTOR( floor, 1 )
	_LMAT_DEFINE_MATH_FUNCTOR( ceil, 1 )

	// exp & log

	_LMAT_DEFINE_MATH_FUNCTOR( exp, 1 )
	_LMAT_DEFINE_MATH_FUNCTOR( log, 1 )
	_LMAT_DEFINE_MATH_FUNCTOR( log10, 1 )
	_LMAT_DEFINE_MATH_FUNCTOR( xlogy, 2 )
	_LMAT_DEFINE_MATH_FUNCTOR( xlogx, 1 )

	// trigonometry

	_LMAT_DEFINE_MATH_FUNCTOR( sin, 1 )
	_LMAT_DEFINE_MATH_FUNCTOR( cos, 1 )
	_LMAT_DEFINE_MATH_FUNCTOR( tan, 1 )

	_LMAT_DEFINE_MATH_FUNCTOR( asin, 1 )
	_LMAT_DEFINE_MATH_FUNCTOR( acos, 1 )
	_LMAT_DEFINE_MATH_FUNCTOR( atan, 1 )
	_LMAT_DEFINE_MATH_FUNCTOR( atan2, 2 )

	// hyperbolic

	_LMAT_DEFINE_MATH_FUNCTOR( sinh, 1 )
	_LMAT_DEFINE_MATH_FUNCTOR( cosh, 1 )
	_LMAT_DEFINE_MATH_FUNCTOR( tanh, 1 )

#ifdef LMAT_HAS_CXX11_MATH

	// cbrt and hypot

	_LMAT_DEFINE_MATH_FUNCTOR( cbrt, 1 )
	_LMAT_DEFINE_MATH_FUNCTOR( hypot, 2 )

	// rounding

	_LMAT_DEFINE_MATH_FUNCTOR( round, 1 )
	_LMAT_DEFINE_MATH_FUNCTOR( trunc, 1 )

	// exp & log

	_LMAT_DEFINE_MATH_FUNCTOR( exp2, 1 )
	_LMAT_DEFINE_MATH_FUNCTOR( log2, 1 )
	_LMAT_DEFINE_MATH_FUNCTOR( expm1, 1 )
	_LMAT_DEFINE_MATH_FUNCTOR( log1p, 1 )

	// inverse hyperbolic

	_LMAT_DEFINE_MATH_FUNCTOR( asinh, 1 )
	_LMAT_DEFINE_MATH_FUNCTOR( acosh, 1 )
	_LMAT_DEFINE_MATH_FUNCTOR( atanh, 1 )

	// numeric predicates

	_LMAT_DEFINE_NUMPRED_FUNCTOR( signbit, 1 )
	_LMAT_DEFINE_NUMPRED_FUNCTOR( isfinite, 1 )
	_LMAT_DEFINE_NUMPRED_FUNCTOR( isinf, 1 )
	_LMAT_DEFINE_NUMPRED_FUNCTOR( isnan, 1 )

	// special functions

	_LMAT_DEFINE_MATH_FUNCTOR( erf, 1 )
	_LMAT_DEFINE_MATH_FUNCTOR( erfc, 1 )
	_LMAT_DEFINE_MATH_FUNCTOR( lgamma, 1 )
	_LMAT_DEFINE_MATH_FUNCTOR( tgamma, 1 )

#endif

	_LMAT_DEFINE_MATH_FUNCTOR( norminv, 1 )

	// conditional selection

	template<typename T>
	struct cond_fun
	{
		typedef T result_type;

		LMAT_ENSURE_INLINE
		T operator() (bool b, const T& x, const T& y) const { return math::cond(b, x, y); }

		LMAT_ENSURE_INLINE
		T operator() (mask_t<T> m, const T& x, const T& y) const { return math::cond(m.bvalue, x, y); }
	};

	template<typename T, typename Kind>
	struct cond_fun<math::simd_pack<T, Kind> >
	{
		typedef math::simd_pack<T, Kind> result_type;

		LMAT_ENSURE_INLINE
		math::simd_pack<T, Kind> operator() (const math::simd_bpack<T, Kind>& b,
				const math::simd_pack<T, Kind>& x, const math::simd_pack<T, Kind>& y) const
		{ return cond(b, x, y); }
	};

	_LMAT_DECL_SIMDIZABLE_S( ftags::cond_, cond_fun )
	LMAT_DEF_TRIVIAL_SIMDIZE_MAP( cond_fun )
}


#endif 


