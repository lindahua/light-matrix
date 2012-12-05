/*
 * @file functor_base.h
 *
 * The basic definitions for functors
 *
 * @author Dahua Lin
 */

#ifndef _MSC_VER
#pragma once
#endif

#ifndef LIGHTMAT_FUNCTOR_BASE_H_
#define LIGHTMAT_FUNCTOR_BASE_H_

#include <light_mat/common/basic_defs.h>
#include <functional>

namespace lmat
{
	// kernel categories

	struct scalar_ker { };
	struct sse_ker { };
	struct avx_ker { };

	// operation type testing

	template<typename Tag>
	struct is_op_tag
	{
		static const bool value = false;
	};

	// maps

	template<typename Tag, class TList>
	struct op_result;

	template<typename Tag, typename Ker, class TList>
	struct op_fun;

}


/************************************************
 *
 *  Macros for tag declaration
 *
 ************************************************/

// generic numerical operation

#define LMAT_DECLARE_NUMERIC_UNARY_OP(Tag) \
	struct Tag { }; \
	template<> struct is_op_tag<Tag> { static const bool value = true; }; \
	template<typename T> \
	struct op_result<Tag, LMAT_TYPELIST_1(T) > { typedef T type; };

#define LMAT_DECLARE_NUMERIC_BINARY_OP(Tag) \
	struct Tag { }; \
	template<> struct is_op_tag<Tag> { static const bool value = true; }; \
	template<typename T> \
	struct op_result<Tag, LMAT_TYPELIST_2(T, T) > { typedef T type; };

#define LMAT_DECLARE_REAL_UNARY_OP(Op) \
	struct Op { }; \
	template<> struct is_op_tag<Op> { static const bool value = true; }; \
	template<> struct op_result<Op, LMAT_TYPELIST_1(float) > { typedef float type; }; \
	template<> struct op_result<Op, LMAT_TYPELIST_1(double) > { typedef double type; };

#define LMAT_DECLARE_REAL_BINARY_OP(Op) \
	struct Op { }; \
	template<> struct is_op_tag<Op> { static const bool value = true; }; \
	template<> struct op_result<Op, LMAT_TYPELIST_2(float, float) > { typedef float type; }; \
	template<> struct op_result<Op, LMAT_TYPELIST_2(double, double) > { typedef double type; };

// operation on real numbers

#define LMAT_DECLARE_NUMERIC_UNARY_OP_S(Tag, SType) \
	struct Tag { }; \
	template<> struct is_op_tag<Tag> { static const bool value = true; }; \
	template<typename T> \
	struct op_result<Tag, LMAT_TYPELIST_1(T) > { typedef SType type; };

#define LMAT_DECLARE_NUMERIC_BINARY_OP_S(Tag, SType) \
	struct Tag { }; \
	template<> struct is_op_tag<Tag> { static const bool value = true; }; \
	template<typename T> \
	struct op_result<Tag, LMAT_TYPELIST_2(T, T) > { typedef SType type; };

#define LMAT_DECLARE_REAL_UNARY_OP_S(Tag, SType) \
	struct Tag { }; \
	template<> struct is_op_tag<Tag> { static const bool value = true; }; \
	template<> struct op_result<Tag, LMAT_TYPELIST_1(float) > { typedef SType type; }; \
	template<> struct op_result<Tag, LMAT_TYPELIST_2(double) > { typedef bool type; };

#define LMAT_DECLARE_REAL_BINARY_OP_S(Tag, SType) \
	struct Tag { }; \
	template<> struct is_op_tag<Tag> { static const bool value = true; }; \
	template<> struct op_result<Tag, LMAT_TYPELIST_2(float, float) > { typedef SType type; }; \
	template<> struct op_result<Tag, LMAT_TYPELIST_2(double, double) > { typedef SType type; };

// logical operations

#define LMAT_DECLARE_LOGICAL_UNARY_OP(Tag) \
	struct Tag { }; \
	template<> struct is_op_tag<Tag> { static const bool value = true; }; \
	template<> struct op_result<Tag, LMAT_TYPELIST_1(bool) > { typedef bool type; };

#define LMAT_DECLARE_LOGICAL_BINARY_OP(Tag) \
	struct Tag { }; \
	template<> struct is_op_tag<Tag> { static const bool value = true; }; \
	template<> struct op_result<Tag, LMAT_TYPELIST_2(bool, bool) > { typedef bool type; };



/************************************************
 *
 *  Macros for scalar functor definition
 *
 ************************************************/

#define LMAT_DEFINE_UNARY_FUNCTOR_GENERIC( Fname, Expr ) \
	template<typename T> struct Fname##_fun { \
		typedef typename op_result<Fname##_t, LMAT_TYPELIST_1(T) >::type result_type; \
		LMAT_ENSURE_INLINE Fname##_fun() { } \
		LMAT_ENSURE_INLINE Fname##_fun( Fname##_t ) { } \
		LMAT_ENSURE_INLINE result_type operator() (const T& x) const { return Expr; } \
	};

#define LMAT_DEFINE_BINARY_FUNCTOR_GENERIC( Fname, Expr ) \
	template<typename T> struct Fname##_fun { \
		typedef typename op_result<Fname##_t, LMAT_TYPELIST_2(T, T) >::type result_type; \
		LMAT_ENSURE_INLINE Fname##_fun() { } \
		LMAT_ENSURE_INLINE Fname##_fun( Fname##_t ) { } \
		LMAT_ENSURE_INLINE result_type operator() (const T& x, const T& y) const { return Expr; } \
	};


// for numeric operations

#define LMAT_DEFINE_NUMERIC_UNARY_FUNCTOR( Fname, Expr ) \
	LMAT_DEFINE_UNARY_FUNCTOR_GENERIC( Fname, Expr ) \
	template<typename T> \
	struct op_fun<Fname##_t, scalar_ker, LMAT_TYPELIST_1(T) > { typedef Fname##_fun<T> type; };

#define LMAT_DEFINE_NUMERIC_BINARY_FUNCTOR( Fname, Expr ) \
	LMAT_DEFINE_BINARY_FUNCTOR_GENERIC( Fname, Expr ) \
	template<typename T> \
	struct op_fun<Fname##_t, scalar_ker, LMAT_TYPELIST_2(T, T) > { typedef Fname##_fun<T> type; };

// for real-number operations

#define LMAT_DEFINE_REAL_UNARY_FUNCTOR( Fname, Expr ) \
	LMAT_DEFINE_UNARY_FUNCTOR_GENERIC( Fname, Expr ) \
	template<> struct op_fun<Fname##_t, scalar_ker, LMAT_TYPELIST_1(float) > { typedef Fname##_fun<float> type; }; \
	template<> struct op_fun<Fname##_t, scalar_ker, LMAT_TYPELIST_1(double) > { typedef Fname##_fun<double> type; };

#define LMAT_DEFINE_REAL_BINARY_FUNCTOR( Fname, Expr ) \
	LMAT_DEFINE_BINARY_FUNCTOR_GENERIC( Fname, Expr ) \
	template<> struct op_fun<Fname##_t, scalar_ker, LMAT_TYPELIST_2(float, float) > { typedef Fname##_fun<float> type; }; \
	template<> struct op_fun<Fname##_t, scalar_ker, LMAT_TYPELIST_2(double, double) > { typedef Fname##_fun<double> type; };

// for logical operations

#define LMAT_DEFINE_LOGICAL_UNARY_FUNCTOR( Fname, Expr ) \
	struct Fname##_fun { \
		typedef bool result_type; \
		LMAT_ENSURE_INLINE Fname##_fun() { } \
		LMAT_ENSURE_INLINE Fname##_fun( Fname##_t ) { } \
		LMAT_ENSURE_INLINE result_type operator() (const bool& x) const { return Expr; } \
	}; \
	template<> struct op_fun<Fname##_t, scalar_ker, LMAT_TYPELIST_1(bool) > { typedef Fname##_fun type; };

#define LMAT_DEFINE_LOGICAL_BINARY_FUNCTOR( Fname, Expr ) \
	struct Fname##_fun { \
		typedef bool result_type; \
		LMAT_ENSURE_INLINE Fname##_fun() { } \
		LMAT_ENSURE_INLINE Fname##_fun( Fname##_t ) { } \
		LMAT_ENSURE_INLINE result_type operator() (const bool& x, const bool& y) const { return Expr; } \
	}; \
	template<> struct op_fun<Fname##_t, scalar_ker, LMAT_TYPELIST_2(bool, bool) > { typedef Fname##_fun type; };


#endif 




