/**
 * @file basic_functors.h
 *
 * Arithmetic and Basic math functors
 * 
 * @author Dahua Lin 
 */

#ifdef _MSC_VER
#pragma once
#endif

#ifndef LIGHTMAT_BASIC_FUNCTORS_H_
#define LIGHTMAT_BASIC_FUNCTORS_H_

#include <light_mat/math/math_base.h>

#include <light_mat/math/fun_tags.h>
#include <light_mat/math/functor_base.h>
#include <light_mat/simd/simd.h>


/************************************************
 *
 *  Macros for defining logical operators
 *
 ************************************************/

#define _LMAT_DEFINE_LOGICAL_FUNCTOR( Name, NA, MExpr, BExpr ) \
	template<typename T> \
	struct Name##_fun { \
		typedef mask_t<T> result_type; \
		LMAT_ENSURE_INLINE \
		result_type operator()(LMAT_REPEAT_ARGS_P##NA(const mask_t<T>& x)) const { return MExpr; } }; \
	template<typename T, typename Kind> \
	struct Name##_fun<simd_pack<T, Kind> > { \
		typedef simd_bpack<T, Kind> result_type; \
		LMAT_ENSURE_INLINE \
		simd_bpack<T, Kind> operator()(LMAT_REPEAT_ARGS_P##NA(const LMAT_SIMD_BPACK_(T, Kind)& x)) const { return MExpr; } }; \
	template<> \
	struct Name##_fun<bool> { \
		typedef bool result_type; \
		LMAT_ENSURE_INLINE \
		bool operator()(LMAT_REPEAT_ARGS_P##NA(bool x)) const { return BExpr; } };

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

#define _LMAT_DEFINE_LOGICAL_FUN(Name, NA, MExpr, BExpr) \
	_LMAT_DEFINE_LOGICAL_FUNCTOR( Name, NA, MExpr, BExpr ) \
	_LMAT_DEFINE_SIMD_SUPPORT( ftags::Name##_, Name##_fun ) \
	_LMAT_DEFINE_LOGICAL_FUNMAP_##NA( Name )



/************************************************
 *
 *  Specific definitions
 *
 ************************************************/

namespace lmat {

	// arithmetics

	_LMAT_DEFINE_GENERIC_MATH_FUN_EX( add, 2, x1 + x2 )
	_LMAT_DEFINE_GENERIC_MATH_FUN_EX( sub, 2, x1 - x2 )
	_LMAT_DEFINE_GENERIC_MATH_FUN_EX( mul, 2, x1 * x2 )
	_LMAT_DEFINE_GENERIC_MATH_FUN_EX( div, 2, x1 / x2 )
	_LMAT_DEFINE_GENERIC_MATH_FUN_EX( neg, 1, -x1 )

	_LMAT_DEFINE_REAL_MATH_FUN( fma, 3 )

	_LMAT_DEFINE_GENERIC_MATH_FUN( max, 2 )
	_LMAT_DEFINE_GENERIC_MATH_FUN( min, 2 )
	_LMAT_DEFINE_GENERIC_MATH_FUN( clamp, 3 )

	// simple power

	_LMAT_DEFINE_GENERIC_MATH_FUN( abs, 1 )
	_LMAT_DEFINE_GENERIC_MATH_FUN( sqr, 1 )

	_LMAT_DEFINE_REAL_MATH_FUN( cube, 1 )
	_LMAT_DEFINE_REAL_MATH_FUN( rcp, 1 )
	_LMAT_DEFINE_REAL_MATH_FUN( sqrt, 1 )
	_LMAT_DEFINE_REAL_MATH_FUN( rsqrt, 1 )

	// comparison

	_LMAT_DEFINE_GENERIC_PRED_FUN_EX( eq, 2, x1 == x2 )
	_LMAT_DEFINE_GENERIC_PRED_FUN_EX( ne, 2, x1 != x2 )
	_LMAT_DEFINE_GENERIC_PRED_FUN_EX( lt, 2, x1 < x2 )
	_LMAT_DEFINE_GENERIC_PRED_FUN_EX( le, 2, x1 <= x2 )
	_LMAT_DEFINE_GENERIC_PRED_FUN_EX( gt, 2, x1 > x2 )
	_LMAT_DEFINE_GENERIC_PRED_FUN_EX( ge, 2, x1 >= x2 )

	// logical

	_LMAT_DEFINE_LOGICAL_FUN( logical_not, 1, ~x1, !x1 )
	_LMAT_DEFINE_LOGICAL_FUN( logical_and, 2, x1 & x2, x1 && x2 )
	_LMAT_DEFINE_LOGICAL_FUN( logical_or,  2, x1 | x2, x1 || x2 )
	_LMAT_DEFINE_LOGICAL_FUN( logical_eq,  2, x1 == x2, x1 == x2 )
	_LMAT_DEFINE_LOGICAL_FUN( logical_ne,  2, x1 != x2, x1 != x2 )

	// rounding

	_LMAT_DEFINE_REAL_MATH_FUN( floor, 1 )
	_LMAT_DEFINE_REAL_MATH_FUN( ceil, 1 )
	_LMAT_DEFINE_REAL_MATH_FUN( round, 1 )
	_LMAT_DEFINE_REAL_MATH_FUN( trunc, 1 )

	// FP classification

	_LMAT_DEFINE_REAL_PRED_FUN( signbit, 1 )
	_LMAT_DEFINE_REAL_PRED_FUN( isinf, 1 )
	_LMAT_DEFINE_REAL_PRED_FUN( isnan, 1 )
	_LMAT_DEFINE_REAL_PRED_FUN( isfinite, 1 )


	// conditional selection

	template<typename T>
	struct cond_fun
	{
		typedef T result_type;

		LMAT_ENSURE_INLINE
		T operator() (bool b, const T& x, const T& y) const
		{
			return math::cond(b, x, y);
		}

		LMAT_ENSURE_INLINE
		T operator() (const mask_t<T>& m, const T& x, const T& y) const
		{
			return math::cond(m, x, y);
		}
	};

	template<typename T, typename Kind>
	struct cond_fun<simd_pack<T, Kind> >
	{
		typedef simd_pack<T, Kind> result_type;
		typedef simd_bpack<T, Kind> mask_type;

		LMAT_ENSURE_INLINE
		result_type operator() (const mask_type& m,
				const result_type& x, const result_type& y) const
		{
			return math::cond(m, x, y);
		}
	};

	_LMAT_DEFINE_SIMD_SUPPORT( ftags::cond_, cond_fun )

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


