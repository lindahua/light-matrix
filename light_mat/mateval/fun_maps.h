/**
 * @file fun_maps.h
 *
 * Functor maps
 * 
 * @author Dahua Lin 
 */

#ifdef _MSC_VER
#pragma once
#endif

#ifndef LIGHTMAT_FUN_MAPS_H_
#define LIGHTMAT_FUN_MAPS_H_

#include <light_mat/mateval/mateval_fwd.h>
#include <light_mat/math/math_functors.h>


#define _LMAT_DEFINE_GENERIC_FUNMAP( Name, NA ) \
	template<typename T> \
	struct fun_map<ftags::Name##_, LMAT_REPEAT_ARGS_##NA(T)> { \
		typedef Name##_fun<T> type; };

#define _LMAT_DEFINE_REAL_FUNMAP( Name, NA ) \
	template<> struct fun_map<ftags::Name##_, LMAT_REPEAT_ARGS_##NA(float)> { \
		typedef Name##_fun<float> type; }; \
	template<> struct fun_map<ftags::Name##_, LMAT_REPEAT_ARGS_##NA(double)> { \
		typedef Name##_fun<double> type; };

#define _LMAT_DEFINE_LOGICAL_FUNMAP_1( Name ) \
	template<> struct fun_map<ftags::Name##_, bool> { \
		typedef Name##_fun<bool> type; }; \
	template<> struct fun_map<ftags::Name##_, mask_t<float> > { \
		typedef Name##_fun<float> type; }; \
	template<> struct fun_map<ftags::Name##_, mask_t<double> > { \
		typedef Name##_fun<double> type; };

#define _LMAT_DEFINE_LOGICAL_FUNMAP_2( Name ) \
	template<> struct fun_map<ftags::Name##_, bool, bool> { \
		typedef Name##_fun<bool> type; }; \
	template<> struct fun_map<ftags::Name##_, bool, mask_t<float> > { \
		typedef Name##_fun<bool> type; }; \
	template<> struct fun_map<ftags::Name##_, mask_t<float>, bool > { \
		typedef Name##_fun<bool> type; }; \
	template<> struct fun_map<ftags::Name##_, bool, mask_t<double> > { \
		typedef Name##_fun<bool> type; }; \
	template<> struct fun_map<ftags::Name##_, mask_t<double>, bool > { \
		typedef Name##_fun<bool> type; }; \
	template<> struct fun_map<ftags::Name##_, mask_t<float>, mask_t<float> > { \
		typedef Name##_fun<float> type; }; \
	template<> struct fun_map<ftags::Name##_, mask_t<double>, mask_t<double> > { \
		typedef Name##_fun<double> type; };


namespace lmat
{
	template<typename FTag, typename... T>
	struct fun_map;

	template<typename FTag, typename... T>
	struct fun_result
	{
		typedef typename std::result_of<typename fun_map<FTag, T...>::type(T...)>::type type;
	};

	// arithmetics

	_LMAT_DEFINE_GENERIC_FUNMAP( add, 2 )
	_LMAT_DEFINE_GENERIC_FUNMAP( sub, 2 )
	_LMAT_DEFINE_GENERIC_FUNMAP( mul, 2 )
	_LMAT_DEFINE_GENERIC_FUNMAP( div, 2 )
	_LMAT_DEFINE_GENERIC_FUNMAP( neg, 1 )

	_LMAT_DEFINE_GENERIC_FUNMAP( abs, 1 )
	_LMAT_DEFINE_GENERIC_FUNMAP( sqr, 1 )
	_LMAT_DEFINE_GENERIC_FUNMAP( cube, 1 )

	_LMAT_DEFINE_GENERIC_FUNMAP( max, 2 )
	_LMAT_DEFINE_GENERIC_FUNMAP( min, 2 )

	// comparison

	_LMAT_DEFINE_GENERIC_FUNMAP( eq, 2 )
	_LMAT_DEFINE_GENERIC_FUNMAP( ne, 2 )
	_LMAT_DEFINE_GENERIC_FUNMAP( ge, 2 )
	_LMAT_DEFINE_GENERIC_FUNMAP( gt, 2 )
	_LMAT_DEFINE_GENERIC_FUNMAP( le, 2 )
	_LMAT_DEFINE_GENERIC_FUNMAP( lt, 2 )

	// logical

	_LMAT_DEFINE_LOGICAL_FUNMAP_1( logical_not )
	_LMAT_DEFINE_LOGICAL_FUNMAP_2( logical_and )
	_LMAT_DEFINE_LOGICAL_FUNMAP_2( logical_or )
	_LMAT_DEFINE_LOGICAL_FUNMAP_2( logical_eq )
	_LMAT_DEFINE_LOGICAL_FUNMAP_2( logical_ne )

	// real math

	_LMAT_DEFINE_REAL_FUNMAP( clamp, 3 )
	_LMAT_DEFINE_REAL_FUNMAP( fma, 3 )

	_LMAT_DEFINE_REAL_FUNMAP( rcp, 1 )
	_LMAT_DEFINE_REAL_FUNMAP( sqrt, 1 )
	_LMAT_DEFINE_REAL_FUNMAP( rsqrt, 1 )
	_LMAT_DEFINE_REAL_FUNMAP( pow, 2 )

	_LMAT_DEFINE_REAL_FUNMAP( floor, 1 )
	_LMAT_DEFINE_REAL_FUNMAP( ceil, 1 )

	_LMAT_DEFINE_REAL_FUNMAP( exp, 1 )
	_LMAT_DEFINE_REAL_FUNMAP( log, 1 )
	_LMAT_DEFINE_REAL_FUNMAP( log10, 1 )
	_LMAT_DEFINE_REAL_FUNMAP( xlogy, 2 )
	_LMAT_DEFINE_REAL_FUNMAP( xlogx, 1 )

	_LMAT_DEFINE_REAL_FUNMAP( sin, 1 )
	_LMAT_DEFINE_REAL_FUNMAP( cos, 1 )
	_LMAT_DEFINE_REAL_FUNMAP( tan, 1 )

	_LMAT_DEFINE_REAL_FUNMAP( asin, 1 )
	_LMAT_DEFINE_REAL_FUNMAP( acos, 1 )
	_LMAT_DEFINE_REAL_FUNMAP( atan, 1 )
	_LMAT_DEFINE_REAL_FUNMAP( atan2, 2 )

	_LMAT_DEFINE_REAL_FUNMAP( sinh, 1 )
	_LMAT_DEFINE_REAL_FUNMAP( cosh, 1 )
	_LMAT_DEFINE_REAL_FUNMAP( tanh, 1 )

	// C++11 real math

#ifdef LMAT_HAS_CXX11_MATH

	_LMAT_DEFINE_REAL_FUNMAP( cbrt, 1 )
	_LMAT_DEFINE_REAL_FUNMAP( hypot, 2 )

	_LMAT_DEFINE_REAL_FUNMAP( round, 1 )
	_LMAT_DEFINE_REAL_FUNMAP( trunc, 1 )

	_LMAT_DEFINE_REAL_FUNMAP( exp2, 1 )
	_LMAT_DEFINE_REAL_FUNMAP( log2, 1 )
	_LMAT_DEFINE_REAL_FUNMAP( expm1, 1 )
	_LMAT_DEFINE_REAL_FUNMAP( log1p, 1 )

	_LMAT_DEFINE_REAL_FUNMAP( asinh, 1 )
	_LMAT_DEFINE_REAL_FUNMAP( acosh, 1 )
	_LMAT_DEFINE_REAL_FUNMAP( atanh, 1 )

	_LMAT_DEFINE_REAL_FUNMAP( erf, 1 )
	_LMAT_DEFINE_REAL_FUNMAP( erfc, 1 )

	_LMAT_DEFINE_REAL_FUNMAP( lgamma, 1 )
	_LMAT_DEFINE_REAL_FUNMAP( tgamma, 1 )

	_LMAT_DEFINE_REAL_FUNMAP( signbit, 1 )
	_LMAT_DEFINE_REAL_FUNMAP( isfinite, 1 )
	_LMAT_DEFINE_REAL_FUNMAP( isinf, 1 )
	_LMAT_DEFINE_REAL_FUNMAP( isnan, 1 )

#endif

	// conditional selection

	template<typename T>
	struct fun_map<ftags::cond_, bool, T, T>
	{
		typedef cond_fun<T> type;
	};

	template<typename T>
	struct fun_map<ftags::cond_, mask_t<T>, T, T>
	{
		typedef cond_fun<T> type;
	};


}

#endif 
