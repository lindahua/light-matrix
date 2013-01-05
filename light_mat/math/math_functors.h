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


#define LMAT_DEFINE_GENERIC_MATH_FUNCTOR_1( Name, Expr ) \
	template<typename T> \
	struct Name##_fun { \
		typedef T value_type; \
		LMAT_ENSURE_INLINE \
		T operator()(const T& x) const { return Expr; } \
		template<typename Kind> \
		LMAT_ENSURE_INLINE \
		simd_pack<T, Kind> operator()(const simd_pack<T, Kind>& x) const \
		{ return Expr; } \
	}; \
	template<typename T, typename Kind> \
	struct fun_simd_pack<Name##_fun<T>, Kind> { \
		typedef simd_pack<T, Kind> type; \
	};

#define LMAT_DEFINE_GENERIC_MATH_FUNCTOR_2( Name, Expr ) \
	template<typename T> \
	struct Name##_fun { \
		typedef T value_type; \
		LMAT_ENSURE_INLINE \
		T operator()(const T& x, const T& y) const { return Expr; } \
		template<typename Kind> \
		LMAT_ENSURE_INLINE \
		simd_pack<T, Kind> operator()(const simd_pack<T, Kind>& x, const simd_pack<T, Kind>& y) const \
		{ return Expr; } \
	}; \
	template<typename T, typename Kind> \
	struct fun_simd_pack<Name##_fun<T>, Kind> { \
		typedef simd_pack<T, Kind> type; \
	};

#define LMAT_DEFINE_GENERIC_MATH_FUNCTOR_3( Name, Expr ) \
	template<typename T> \
	struct Name##_fun { \
		typedef T value_type; \
		LMAT_ENSURE_INLINE \
		T operator()(const T& x, const T& y, const T& z) const { return Expr; } \
		template<typename Kind> \
		LMAT_ENSURE_INLINE \
		simd_pack<T, Kind> operator()(const simd_pack<T, Kind>& x, const simd_pack<T, Kind>& y, const simd_pack<T, Kind>& z) const \
		{ return Expr; } \
	}; \
	template<typename T, typename Kind> \
	struct fun_simd_pack<Name##_fun<T>, Kind> { \
		typedef simd_pack<T, Kind> type; \
	};

#define LMAT_DEFINE_NUMPRED_FUNCTOR( Name, Expr ) \
	template<typename T> \
	struct Name##_fun { \
		typedef mask_t<T> value_type; \
		LMAT_ENSURE_INLINE \
		mask_t<T> operator()(const T& x) const { return Expr; } \
		template<typename Kind> \
		LMAT_ENSURE_INLINE \
		simd_bpack<T, Kind> operator()(const simd_pack<T, Kind>& x) const \
		{ return Expr; } \
	}; \
	template<typename T, typename Kind> \
	struct fun_simd_pack<Name##_fun<T>, Kind> { \
		typedef simd_bpack<T, Kind> type; \
	};

#define LMAT_DEFINE_COMPARISON_FUNCTOR( Name, Expr ) \
	template<typename T> \
	struct Name##_fun { \
		typedef mask_t<T> value_type; \
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
		typedef mask_t<T> value_type; \
		LMAT_ENSURE_INLINE \
		mask_t<T> operator()(const mask_t<T>& x) const { return MExpr; } \
		template<typename Kind> \
		LMAT_ENSURE_INLINE \
		simd_bpack<T, Kind> operator()(const simd_bpack<T, Kind>& x) const \
		{ return MExpr; } \
	}; \
	template<> \
	struct Name##_fun<bool> { \
		typedef bool value_type; \
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
		typedef mask_t<T> value_type; \
		LMAT_ENSURE_INLINE \
		mask_t<T> operator()(const mask_t<T>& x, const mask_t<T>& y) const { return MExpr; } \
		template<typename Kind> \
		LMAT_ENSURE_INLINE \
		simd_bpack<T, Kind> operator()(const simd_bpack<T, Kind>& x, const simd_bpack<T, Kind>& y) const \
		{ return MExpr; } \
	}; \
	template<> \
	struct Name##_fun<bool> { \
		typedef bool value_type; \
		LMAT_ENSURE_INLINE \
		bool operator()(const bool& x, const bool& y) const { return BExpr; } \
	}; \
	template<typename Kind> \
	struct fun_simd_pack<Name##_fun<float>, Kind> { typedef simd_bpack<float, Kind> type; }; \
	template<typename Kind> \
	struct fun_simd_pack<Name##_fun<double>, Kind> { typedef simd_bpack<double, Kind> type; };


#define LMAT_DEFINE_MATH_FUNCTOR_1( Name ) \
	LMAT_DEFINE_GENERIC_MATH_FUNCTOR_1( Name, Name(x) )

#define LMAT_DEFINE_MATH_FUNCTOR_2( Name ) \
	LMAT_DEFINE_GENERIC_MATH_FUNCTOR_2( Name, Name(x, y) )

#define LMAT_DEFINE_MATH_FUNCTOR_3( Name ) \
	LMAT_DEFINE_GENERIC_MATH_FUNCTOR_3( Name, Name(x, y, z) )


#define LMAT_DEFINE_NUMPRED_FUNCTOR_1( Name ) \
	LMAT_DEFINE_NUMPRED_FUNCTOR( Name, Name(x) )


namespace lmat { namespace math {

	template<class Fun, typename Kind> struct fun_simd_pack;

	// Arithmetic functors

	LMAT_DEFINE_GENERIC_MATH_FUNCTOR_2( add, x + y )
	LMAT_DEFINE_GENERIC_MATH_FUNCTOR_2( sub, x - y )
	LMAT_DEFINE_GENERIC_MATH_FUNCTOR_2( mul, x * y )
	LMAT_DEFINE_GENERIC_MATH_FUNCTOR_2( div, x / y )
	LMAT_DEFINE_GENERIC_MATH_FUNCTOR_1( neg, -x )

	LMAT_DEFINE_MATH_FUNCTOR_1( abs )
	LMAT_DEFINE_MATH_FUNCTOR_1( sqr )
	LMAT_DEFINE_MATH_FUNCTOR_1( cube )
	LMAT_DEFINE_MATH_FUNCTOR_1( rcp )
	LMAT_DEFINE_MATH_FUNCTOR_3( fma )

	LMAT_DEFINE_MATH_FUNCTOR_2( diff_abs )
	LMAT_DEFINE_MATH_FUNCTOR_2( diff_sqr )

	LMAT_DEFINE_MATH_FUNCTOR_2( max )
	LMAT_DEFINE_MATH_FUNCTOR_2( min )
	LMAT_DEFINE_MATH_FUNCTOR_3( clamp )

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

	LMAT_DEFINE_MATH_FUNCTOR_1( sqrt )
	LMAT_DEFINE_MATH_FUNCTOR_1( rsqrt )
	LMAT_DEFINE_MATH_FUNCTOR_2( pow )

	// floor & ceil

	LMAT_DEFINE_MATH_FUNCTOR_1( floor )
	LMAT_DEFINE_MATH_FUNCTOR_1( ceil )

	// exp & log

	LMAT_DEFINE_MATH_FUNCTOR_1( exp )
	LMAT_DEFINE_MATH_FUNCTOR_1( log )
	LMAT_DEFINE_MATH_FUNCTOR_1( log10 )
	LMAT_DEFINE_MATH_FUNCTOR_2( xlogy )
	LMAT_DEFINE_MATH_FUNCTOR_1( xlogx )

	// trigonometry

	LMAT_DEFINE_MATH_FUNCTOR_1( sin )
	LMAT_DEFINE_MATH_FUNCTOR_1( cos )
	LMAT_DEFINE_MATH_FUNCTOR_1( tan )

	LMAT_DEFINE_MATH_FUNCTOR_1( asin )
	LMAT_DEFINE_MATH_FUNCTOR_1( acos )
	LMAT_DEFINE_MATH_FUNCTOR_1( atan )
	LMAT_DEFINE_MATH_FUNCTOR_2( atan2 )

	// hyperbolic

	LMAT_DEFINE_MATH_FUNCTOR_1( sinh )
	LMAT_DEFINE_MATH_FUNCTOR_1( cosh )
	LMAT_DEFINE_MATH_FUNCTOR_1( tanh )

#ifdef LMAT_HAS_CXX11_MATH

	// cbrt and hypot

	LMAT_DEFINE_MATH_FUNCTOR_1( cbrt )
	LMAT_DEFINE_MATH_FUNCTOR_2( hypot )

	// rounding

	LMAT_DEFINE_MATH_FUNCTOR_1( round )
	LMAT_DEFINE_MATH_FUNCTOR_1( trunc )

	// exp & log

	LMAT_DEFINE_MATH_FUNCTOR_1( exp2 )
	LMAT_DEFINE_MATH_FUNCTOR_1( log2 )
	LMAT_DEFINE_MATH_FUNCTOR_1( expm1 )
	LMAT_DEFINE_MATH_FUNCTOR_1( log1p )

	// inverse hyperbolic

	LMAT_DEFINE_MATH_FUNCTOR_1( asinh )
	LMAT_DEFINE_MATH_FUNCTOR_1( acosh )
	LMAT_DEFINE_MATH_FUNCTOR_1( atanh )

	// error and gamma

	LMAT_DEFINE_MATH_FUNCTOR_1( erf )
	LMAT_DEFINE_MATH_FUNCTOR_1( erfc )
	LMAT_DEFINE_MATH_FUNCTOR_1( lgamma )
	LMAT_DEFINE_MATH_FUNCTOR_1( tgamma )

	// numeric predicates

	LMAT_DEFINE_NUMPRED_FUNCTOR_1( signbit )
	LMAT_DEFINE_NUMPRED_FUNCTOR_1( isfinite )
	LMAT_DEFINE_NUMPRED_FUNCTOR_1( isinf )
	LMAT_DEFINE_NUMPRED_FUNCTOR_1( isnan )

#endif


	// conditional selection

	template<typename T>
	struct cond_fun
	{
		typedef T value_type;

		LMAT_ENSURE_INLINE
		T operator() (bool b, const T& x, const T& y) const { return cond(b, x, y); }

		LMAT_ENSURE_INLINE
		T operator() (mask_t<T> m, const T& x, const T& y) const { return cond(m.bvalue, x, y); }
	};

	template<> struct cond_fun<float>
	{
		typedef float value_type;

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
		typedef double value_type;

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


