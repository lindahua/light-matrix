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
#include <light_mat/math/fun_tags.h>
#include <light_mat/math/simd.h>


#define LMAT_DEFINE_GENERIC_MATH_FUNCTOR( Name, NA, Expr ) \
	template<typename T> \
	struct Name##_fun { \
		typedef T result_type; \
		LMAT_ENSURE_INLINE \
		T operator()(LMAT_REPEAT_ARGS_P##NA(const T& x)) const { return Expr; } \
		template<typename Kind> \
		LMAT_ENSURE_INLINE \
		simd_pack<T, Kind> operator()(LMAT_REPEAT_ARGS_P##NA(const LMAT_SIMD_PACK_(T, Kind)& x)) const \
		{ return Expr; } \
	}; \
	template<typename T, typename Kind> \
	struct fun_simd_pack<Name##_fun<T>, Kind> { \
		typedef simd_pack<T, Kind> type; \
	};

#define LMAT_DEFINE_GENERIC_NUMPRED_FUNCTOR( Name, NA, Expr ) \
	template<typename T> \
	struct Name##_fun { \
		typedef mask_t<T> result_type; \
		LMAT_ENSURE_INLINE \
		mask_t<T> operator()(LMAT_REPEAT_ARGS_P##NA(const T& x)) const { return Expr; } \
		template<typename Kind> \
		LMAT_ENSURE_INLINE \
		simd_bpack<T, Kind> operator()(LMAT_REPEAT_ARGS_P##NA(const LMAT_SIMD_PACK_(T, Kind)& x)) const \
		{ return Expr; } \
	}; \
	template<typename T, typename Kind> \
	struct fun_simd_pack<Name##_fun<T>, Kind> { \
		typedef simd_bpack<T, Kind> type; \
	};

#define LMAT_DEFINE_COMPARISON_FUNCTOR( Name, Expr ) \
	template<typename T> \
	struct Name##_fun { \
		typedef mask_t<T> result_type; \
		LMAT_ENSURE_INLINE \
		mask_t<T> operator()(const T& x, const T& y) const { return Expr; } \
		template<typename Kind> \
		LMAT_ENSURE_INLINE \
		simd_bpack<T, Kind> operator()(const simd_pack<T, Kind>& x, const simd_pack<T, Kind>& y) const \
		{ return Expr; } \
	}; \
	template<typename T, typename Kind> \
	struct fun_simd_pack<Name##_fun<T>, Kind> { \
		typedef simd_bpack<T, Kind> type; \
	};

#define LMAT_DEFINE_LOGICAL_FUNCTOR_1( Name, MExpr, BExpr ) \
	template<typename T> \
	struct Name##_fun { \
		typedef mask_t<T> result_type; \
		LMAT_ENSURE_INLINE \
		mask_t<T> operator()(const mask_t<T>& x) const { return MExpr; } \
		template<typename Kind> \
		LMAT_ENSURE_INLINE \
		simd_bpack<T, Kind> operator()(const simd_bpack<T, Kind>& x) const \
		{ return MExpr; } \
	}; \
	template<> \
	struct Name##_fun<bool> { \
		typedef bool result_type; \
		LMAT_ENSURE_INLINE \
		bool operator()(const bool& x) const { return BExpr; } \
	}; \
	template<typename Kind> \
	struct fun_simd_pack<Name##_fun<float>, Kind> { typedef simd_bpack<float, Kind> type; }; \
	template<typename Kind> \
	struct fun_simd_pack<Name##_fun<double>, Kind> { typedef simd_bpack<double, Kind> type; };


#define LMAT_DEFINE_LOGICAL_FUNCTOR_2( Name, MExpr, BExpr ) \
	template<typename T> \
	struct Name##_fun { \
		typedef mask_t<T> result_type; \
		LMAT_ENSURE_INLINE \
		mask_t<T> operator()(const mask_t<T>& x, const mask_t<T>& y) const { return MExpr; } \
		template<typename Kind> \
		LMAT_ENSURE_INLINE \
		simd_bpack<T, Kind> operator()(const simd_bpack<T, Kind>& x, const simd_bpack<T, Kind>& y) const \
		{ return MExpr; } \
	}; \
	template<> \
	struct Name##_fun<bool> { \
		typedef bool result_type; \
		LMAT_ENSURE_INLINE \
		bool operator()(const bool& x, const bool& y) const { return BExpr; } \
	}; \
	template<typename Kind> \
	struct fun_simd_pack<Name##_fun<float>, Kind> { typedef simd_bpack<float, Kind> type; }; \
	template<typename Kind> \
	struct fun_simd_pack<Name##_fun<double>, Kind> { typedef simd_bpack<double, Kind> type; };


#define LMAT_DEFINE_MATH_FUNCTOR( Name, NA ) \
	LMAT_DEFINE_GENERIC_MATH_FUNCTOR( Name, NA, Name(LMAT_REPEAT_ARGS_P##NA(x)) )

#define LMAT_DEFINE_NUMPRED_FUNCTOR( Name, NA ) \
	LMAT_DEFINE_GENERIC_NUMPRED_FUNCTOR( Name, NA, Name(LMAT_REPEAT_ARGS_P##NA(x)) )


namespace lmat { namespace math {

	template<class Fun, typename Kind> struct fun_simd_pack;

	// Arithmetic functors

	LMAT_DEFINE_GENERIC_MATH_FUNCTOR( add, 2, x1 + x2 )
	LMAT_DEFINE_GENERIC_MATH_FUNCTOR( sub, 2, x1 - x2 )
	LMAT_DEFINE_GENERIC_MATH_FUNCTOR( mul, 2, x1 * x2 )
	LMAT_DEFINE_GENERIC_MATH_FUNCTOR( div, 2, x1 / x2 )
	LMAT_DEFINE_GENERIC_MATH_FUNCTOR( neg, 1, -x1 )

	LMAT_DEFINE_MATH_FUNCTOR( abs, 1 )
	LMAT_DEFINE_MATH_FUNCTOR( sqr, 1 )
	LMAT_DEFINE_MATH_FUNCTOR( cube, 1 )
	LMAT_DEFINE_MATH_FUNCTOR( rcp, 1 )
	LMAT_DEFINE_MATH_FUNCTOR( fma, 3 )

	LMAT_DEFINE_MATH_FUNCTOR( max, 2 )
	LMAT_DEFINE_MATH_FUNCTOR( min, 2 )
	LMAT_DEFINE_MATH_FUNCTOR( clamp, 3 )

	// Comparison functors

	LMAT_DEFINE_COMPARISON_FUNCTOR( eq, x == y )
	LMAT_DEFINE_COMPARISON_FUNCTOR( ne, x != y )
	LMAT_DEFINE_COMPARISON_FUNCTOR( ge, x >= y )
	LMAT_DEFINE_COMPARISON_FUNCTOR( gt, x > y )
	LMAT_DEFINE_COMPARISON_FUNCTOR( le, x <= y )
	LMAT_DEFINE_COMPARISON_FUNCTOR( lt, x < y )

	// logical functors

	LMAT_DEFINE_LOGICAL_FUNCTOR_1( logical_not, ~x, !x )
	LMAT_DEFINE_LOGICAL_FUNCTOR_2( logical_eq,  x == y, x == y )
	LMAT_DEFINE_LOGICAL_FUNCTOR_2( logical_ne,  x != y, x != y )
	LMAT_DEFINE_LOGICAL_FUNCTOR_2( logical_or,  x | y,  x || y )
	LMAT_DEFINE_LOGICAL_FUNCTOR_2( logical_and, x & y,  x && y )

	// power

	LMAT_DEFINE_MATH_FUNCTOR( sqrt, 1 )
	LMAT_DEFINE_MATH_FUNCTOR( rsqrt, 1 )
	LMAT_DEFINE_MATH_FUNCTOR( pow, 2 )

	// floor & ceil

	LMAT_DEFINE_MATH_FUNCTOR( floor, 1 )
	LMAT_DEFINE_MATH_FUNCTOR( ceil, 1 )

	// exp & log

	LMAT_DEFINE_MATH_FUNCTOR( exp, 1 )
	LMAT_DEFINE_MATH_FUNCTOR( log, 1 )
	LMAT_DEFINE_MATH_FUNCTOR( log10, 1 )
	LMAT_DEFINE_MATH_FUNCTOR( xlogy, 2 )
	LMAT_DEFINE_MATH_FUNCTOR( xlogx, 1 )

	// trigonometry

	LMAT_DEFINE_MATH_FUNCTOR( sin, 1 )
	LMAT_DEFINE_MATH_FUNCTOR( cos, 1 )
	LMAT_DEFINE_MATH_FUNCTOR( tan, 1 )

	LMAT_DEFINE_MATH_FUNCTOR( asin, 1 )
	LMAT_DEFINE_MATH_FUNCTOR( acos, 1 )
	LMAT_DEFINE_MATH_FUNCTOR( atan, 1 )
	LMAT_DEFINE_MATH_FUNCTOR( atan2, 2 )

	// hyperbolic

	LMAT_DEFINE_MATH_FUNCTOR( sinh, 1 )
	LMAT_DEFINE_MATH_FUNCTOR( cosh, 1 )
	LMAT_DEFINE_MATH_FUNCTOR( tanh, 1 )

#ifdef LMAT_HAS_CXX11_MATH

	// cbrt and hypot

	LMAT_DEFINE_MATH_FUNCTOR( cbrt, 1 )
	LMAT_DEFINE_MATH_FUNCTOR( hypot, 2 )

	// rounding

	LMAT_DEFINE_MATH_FUNCTOR( round, 1 )
	LMAT_DEFINE_MATH_FUNCTOR( trunc, 1 )

	// exp & log

	LMAT_DEFINE_MATH_FUNCTOR( exp2, 1 )
	LMAT_DEFINE_MATH_FUNCTOR( log2, 1 )
	LMAT_DEFINE_MATH_FUNCTOR( expm1, 1 )
	LMAT_DEFINE_MATH_FUNCTOR( log1p, 1 )

	// inverse hyperbolic

	LMAT_DEFINE_MATH_FUNCTOR( asinh, 1 )
	LMAT_DEFINE_MATH_FUNCTOR( acosh, 1 )
	LMAT_DEFINE_MATH_FUNCTOR( atanh, 1 )

	// error and gamma

	LMAT_DEFINE_MATH_FUNCTOR( erf, 1 )
	LMAT_DEFINE_MATH_FUNCTOR( erfc, 1 )
	LMAT_DEFINE_MATH_FUNCTOR( lgamma, 1 )
	LMAT_DEFINE_MATH_FUNCTOR( tgamma, 1 )

	// numeric predicates

	LMAT_DEFINE_NUMPRED_FUNCTOR( signbit, 1 )
	LMAT_DEFINE_NUMPRED_FUNCTOR( isfinite, 1 )
	LMAT_DEFINE_NUMPRED_FUNCTOR( isinf, 1 )
	LMAT_DEFINE_NUMPRED_FUNCTOR( isnan, 1 )

#endif


	// conditional selection

	template<typename T>
	struct cond_fun
	{
		typedef T result_type;

		LMAT_ENSURE_INLINE
		T operator() (bool b, const T& x, const T& y) const { return cond(b, x, y); }

		LMAT_ENSURE_INLINE
		T operator() (mask_t<T> m, const T& x, const T& y) const { return cond(m.bvalue, x, y); }
	};

	template<> struct cond_fun<float>
	{
		typedef float result_type;

		LMAT_ENSURE_INLINE
		float operator() (bool b, const float& x, const float& y) const
		{ return cond(b, x, y); }

		LMAT_ENSURE_INLINE
		float operator() (mask_t<float> b, const float& x, const float& y) const
		{ return cond(b, x, y); }

		template<typename Kind>
		LMAT_ENSURE_INLINE
		simd_pack<float, Kind> operator() (const simd_bpack<float, Kind>& b,
				const simd_pack<float, Kind>& x, const simd_pack<float, Kind>& y) const
		{ return cond(b, x, y); }
	};

	template<> struct cond_fun<double>
	{
		typedef double result_type;

		LMAT_ENSURE_INLINE
		double operator() (bool b, const double& x, const double& y) const
		{ return cond(b, x, y); }

		LMAT_ENSURE_INLINE
		double operator() (mask_t<double> b, const double& x, const double& y) const
		{ return cond(b, x, y); }

		template<typename Kind>
		LMAT_ENSURE_INLINE
		simd_pack<double, Kind> operator() (const simd_bpack<double, Kind>& b,
				const simd_pack<double, Kind>& x, const simd_pack<double, Kind>& y) const
		{ return cond(b, x, y); }
	};

	template<typename Kind>
	struct fun_simd_pack<cond_fun<float>, Kind> { typedef simd_pack<float, Kind> type; };

	template<typename Kind>
	struct fun_simd_pack<cond_fun<double>, Kind> { typedef simd_pack<double, Kind> type; };

} }


#endif 


