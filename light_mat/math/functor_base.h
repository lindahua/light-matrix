/**
 * @file functor_base.h
 *
 * @brief The basis of evaluation functors
 *
 * @author Dahua Lin
 */

#ifdef _MSC_VER
#pragma once
#endif

#ifndef LIGHTMAT_FUNCTOR_BASE_H_
#define LIGHTMAT_FUNCTOR_BASE_H_

#include <light_mat/common/basic_defs.h>
#include <light_mat/simd/simd_base.h>

namespace lmat {

	template<typename Fun, typename Kind>
	struct is_simdizable : public meta::false_ { };

	template<typename Fun, typename Kind>
	struct simdize_map;

	// generic result of predicates

	template<typename T>
	struct pred_result
	{
		typedef mask_t<T> type;
	};

	template<typename T>
	struct pred_result<mask_t<T> >
	{
		typedef mask_t<T> type;
	};

	template<>
	struct pred_result<bool>
	{
		typedef bool type;
	};

	template<typename T, typename Kind>
	struct pred_result<simd_pack<T, Kind> >
	{
		typedef simd_bpack<T, Kind> type;
	};

	template<typename T, typename Kind>
	struct pred_result<simd_bpack<T, Kind> >
	{
		typedef simd_bpack<T, Kind> type;
	};

	// functor map

	template<typename FTag, typename... T>
	struct fun_map;

	template<typename FTag, typename... T>
	struct fun_result
	{
		typedef typename std::result_of<typename fun_map<FTag, T...>::type(T...)>::type type;
	};
}


/************************************************
 *
 *  Macros to support SIMDization
 *
 ************************************************/

#define LMAT_DECL_SIMDIZABLE_ON_REAL(FunT) \
		template<typename Kind> \
		struct is_simdizable<FunT<float>, Kind> : public std::true_type { }; \
		template<typename Kind> \
		struct is_simdizable<FunT<double>, Kind> : public std::true_type { };

#define LMAT_DEF_TRIVIAL_SIMDIZE_MAP(FunT) \
		template<typename Kind> \
		struct simdize_map<FunT<float>, Kind> { \
			typedef FunT<simd_pack<float, Kind> > type; \
			LMAT_ENSURE_INLINE \
			static type get(FunT<float> ) { return type(); } \
		}; \
		template<typename Kind> \
		struct simdize_map<FunT<double>, Kind> { \
			typedef FunT<simd_pack<double, Kind> > type; \
			LMAT_ENSURE_INLINE \
			static type get(FunT<double> ) { return type(); } \
		};

#define LMAT_DEF_SIMD_SUPPORT( FunT ) \
	LMAT_DECL_SIMDIZABLE_ON_REAL( FunT ) \
	LMAT_DEF_TRIVIAL_SIMDIZE_MAP( FunT )


/************************************************
 *
 *  Macros for functor definitions
 *
 ************************************************/

// fun-map

#define LMAT_DEF_GENERIC_FUNMAP( FTag, NA, Functor ) \
	template<typename T> \
	struct fun_map<FTag, LMAT_REPEAT_ARGS_##NA(T)> { \
		typedef Functor<T> type; };

#define LMAT_DEF_REAL_FUNMAP( FTag, NA, Functor ) \
	template<> \
	struct fun_map<FTag, LMAT_REPEAT_ARGS_##NA(float)> { \
		typedef Functor<float> type; }; \
	template<> \
	struct fun_map<FTag, LMAT_REPEAT_ARGS_##NA(double)> { \
		typedef Functor<double> type; };

// Useful macros to define functors

#define LMAT_DEF_GENERIC_MATH_FUNCTOR( FTag, NA, Functor, Expr ) \
	template<typename T> \
	struct Functor { \
		typedef T result_type; \
		LMAT_ENSURE_INLINE \
		result_type operator()(LMAT_REPEAT_ARGS_P##NA(const T& x)) const \
		{ return Expr; } }; \

#define LMAT_DEF_GENERIC_PRED_FUNCTOR( FTag, NA, Functor, Expr ) \
	template<typename T> \
	struct Functor { \
		typedef typename pred_result<T>::type result_type; \
		LMAT_ENSURE_INLINE \
		result_type operator()(LMAT_REPEAT_ARGS_P##NA(const T& x)) const \
		{ return Expr; } };

// Macros to define both functor and fun-map together

#define LMAT_DEF_GENERIC_MATH_FUN( FTag, NA, Functor, Expr ) \
	LMAT_DEF_GENERIC_MATH_FUNCTOR( FTag, NA, Functor, Expr ) \
	LMAT_DEF_GENERIC_FUNMAP( FTag, NA, Functor )

#define LMAT_DEF_GENERIC_PRED_FUN( FTag, NA, Functor, Expr ) \
	LMAT_DEF_GENERIC_PRED_FUNCTOR( FTag, NA, Functor, Expr ) \
	LMAT_DEF_GENERIC_FUNMAP( FTag, NA, Functor )

#define LMAT_DEF_REAL_MATH_FUN( FTag, NA, Functor, Expr ) \
	LMAT_DEF_GENERIC_MATH_FUNCTOR( FTag, NA, Functor, Expr ) \
	LMAT_DEF_REAL_FUNMAP( FTag, NA, Functor )

#define LMAT_DEF_REAL_PRED_FUN( FTag, NA, Functor, Expr ) \
	LMAT_DEF_GENERIC_PRED_FUNCTOR( FTag, NA, Functor, Expr ) \
	LMAT_DEF_REAL_FUNMAP( FTag, NA, Functor )


/************************************************
 *
 *  Internal Macros
 *
 ************************************************/

#define _LMAT_DEFINE_SIMD_SUPPORT(FTag, FunT) \
	LMAT_DEF_TRIVIAL_SIMDIZE_MAP( FunT ) \
	template<typename Kind> \
	struct is_simdizable<FunT<float>, Kind> : public meta::has_simd_support<FTag, float, Kind> { }; \
	template<typename Kind> \
	struct is_simdizable<FunT<double>, Kind> : public meta::has_simd_support<FTag, double, Kind> { };

#define _LMAT_DEFINE_GENERIC_MATH_FUN_EX( Name, NA, Expr ) \
	LMAT_DEF_GENERIC_MATH_FUN( ftags::Name##_, NA, Name##_fun, Expr ) \
	_LMAT_DEFINE_SIMD_SUPPORT( ftags::Name##_, Name##_fun )

#define _LMAT_DEFINE_GENERIC_PRED_FUN_EX( Name, NA, Expr ) \
	LMAT_DEF_GENERIC_PRED_FUN( ftags::Name##_, NA, Name##_fun, Expr ) \
	_LMAT_DEFINE_SIMD_SUPPORT( ftags::Name##_, Name##_fun )

#define _LMAT_DEFINE_GENERIC_MATH_FUN( Name, NA ) \
	LMAT_DEF_GENERIC_MATH_FUN( ftags::Name##_, NA, Name##_fun, math::Name(LMAT_REPEAT_ARGS_P##NA(x)) ) \
	_LMAT_DEFINE_SIMD_SUPPORT( ftags::Name##_, Name##_fun )

#define _LMAT_DEFINE_GENERIC_PRED_FUN( Name, NA ) \
	LMAT_DEF_GENERIC_PRED_FUN( ftags::Name##_, NA, Name##_fun, math::Name(LMAT_REPEAT_ARGS_P##NA(x)) ) \
	_LMAT_DEFINE_SIMD_SUPPORT( ftags::Name##_, Name##_fun )

#define _LMAT_DEFINE_REAL_MATH_FUN( Name, NA ) \
	LMAT_DEF_REAL_MATH_FUN( ftags::Name##_, NA, Name##_fun, math::Name(LMAT_REPEAT_ARGS_P##NA(x)) ) \
	_LMAT_DEFINE_SIMD_SUPPORT( ftags::Name##_, Name##_fun )

#define _LMAT_DEFINE_REAL_PRED_FUN( Name, NA ) \
	LMAT_DEF_REAL_PRED_FUN( ftags::Name##_, NA, Name##_fun, math::Name(LMAT_REPEAT_ARGS_P##NA(x)) ) \
	_LMAT_DEFINE_SIMD_SUPPORT( ftags::Name##_, Name##_fun )


#endif
