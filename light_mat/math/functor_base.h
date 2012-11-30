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

	struct any_kernel_t { };
	struct scalar_kernel_t { };
	struct simd_kernel_t { };

	// operation type testing

	template<typename Op>
	struct is_unary_op
	{
		static const bool value = false;
	};

	template<typename Op>
	struct is_binary_op
	{
		static const bool value = false;
	};

	template<typename Op>
	struct is_ternary_op
	{
		static const bool value = false;
	};

	// maps

	template<typename Op, typename T>
	struct unary_op_result;

	template<typename Op, typename T1, typename T2>
	struct binary_op_result;

	template<typename Op, typename T1, typename T2, typename T3>
	struct ternary_op_result;


	template<typename Op, typename Ker, typename T>
	struct unary_op_fun;

	template<typename Op, typename Ker, typename T1, typename T2>
	struct binary_op_fun;

	template<typename Op, typename Ker, typename T1, typename T2, typename T3>
	struct ternary_op_fun;

}

/************************************************
 *
 *  Macros for tag declaration
 *
 ************************************************/

// generic numerical operation

#define LMAT_DECLARE_NUMERIC_UNARY_OP(Tag) \
	struct Tag { }; \
	template<> struct is_unary_op<Tag> { static const bool value = true; }; \
	template<typename T> \
	struct unary_op_result<Tag, T> { typedef T type; };

#define LMAT_DECLARE_NUMERIC_BINARY_OP(Tag) \
	struct Tag { }; \
	template<> struct is_binary_op<Tag> { static const bool value = true; }; \
	template<typename T> \
	struct binary_op_result<Tag, T, T> { typedef T type; };

#define LMAT_DECLARE_REAL_UNARY_OP(Op) \
	struct Op { }; \
	template<> struct is_unary_op<Op> { static const bool value = true; }; \
	template<> struct unary_op_result<Op, float> { typedef float type; }; \
	template<> struct unary_op_result<Op, double> { typedef double type; };

#define LMAT_DECLARE_REAL_BINARY_OP(Op) \
	struct Op { }; \
	template<> struct is_binary_op<Op> { static const bool value = true; }; \
	template<> struct binary_op_result<Op, float, float> { typedef float type; }; \
	template<> struct binary_op_result<Op, double, double> { typedef double type; };

// operation on real numbers

#define LMAT_DECLARE_NUMERIC_UNARY_OP_S(Tag, SType) \
	struct Tag { }; \
	template<> struct is_unary_op<Tag> { static const bool value = true; }; \
	template<typename T> \
	struct unary_op_result<Tag, T> { typedef SType type; };

#define LMAT_DECLARE_NUMERIC_BINARY_OP_S(Tag, SType) \
	struct Tag { }; \
	template<> struct is_binary_op<Tag> { static const bool value = true; }; \
	template<typename T> \
	struct binary_op_result<Tag, T, T> { typedef SType type; };

#define LMAT_DECLARE_REAL_UNARY_OP_S(Tag, SType) \
	struct Tag { }; \
	template<> struct is_unary_op<Tag> { static const bool value = true; }; \
	template<> struct unary_op_result<Tag, float> { typedef SType type; }; \
	template<> struct unary_op_result<Tag, double> { typedef bool type; };

#define LMAT_DECLARE_REAL_BINARY_OP_S(Tag, SType) \
	struct Tag { }; \
	template<> struct is_binary_op<Tag> { static const bool value = true; }; \
	template<> struct binary_op_result<Tag, float, float> { typedef SType type; }; \
	template<> struct binary_op_result<Tag, double, double> { typedef SType type; };

// logical operations

#define LMAT_DECLARE_LOGICAL_UNARY_OP(Tag) \
	struct Tag { }; \
	template<> struct is_unary_op<Tag> { static const bool value = true; }; \
	template<> struct unary_op_result<Tag, bool> { typedef bool type; };

#define LMAT_DECLARE_LOGICAL_BINARY_OP(Tag) \
	struct Tag { }; \
	template<> struct is_binary_op<Tag> { static const bool value = true; }; \
	template<> struct binary_op_result<Tag, bool, bool> { typedef bool type; };



/************************************************
 *
 *  Macros for scalar functor definition
 *
 ************************************************/

#define LMAT_DEFINE_UNARY_FUNCTOR_GENERIC( Fname, Expr ) \
	template<typename T> struct Fname##_fun { \
		typedef typename unary_op_result<Fname##_t, T>::type result_type; \
		LMAT_ENSURE_INLINE Fname##_fun() { } \
		LMAT_ENSURE_INLINE Fname##_fun( Fname##_t ) { } \
		LMAT_ENSURE_INLINE result_type operator() (const T& x) const { return Expr; } \
	};

#define LMAT_DEFINE_BINARY_FUNCTOR_GENERIC( Fname, Expr ) \
	template<typename T> struct Fname##_fun { \
		typedef typename binary_op_result<Fname##_t, T, T>::type result_type; \
		LMAT_ENSURE_INLINE Fname##_fun() { } \
		LMAT_ENSURE_INLINE Fname##_fun( Fname##_t ) { } \
		LMAT_ENSURE_INLINE T operator() (const T& x, const T& y) const { return Expr; } \
	};


// for numeric operations

#define LMAT_DEFINE_NUMERIC_UNARY_FUNCTOR( Fname, Expr ) \
	LMAT_DEFINE_UNARY_FUNCTOR_GENERIC( Fname, Expr ) \
	template<typename T> \
	struct unary_op_fun<Fname##_t, scalar_kernel_t, T> { typedef Fname##_fun<T> type; };

#define LMAT_DEFINE_NUMERIC_BINARY_FUNCTOR( Fname, Expr ) \
	LMAT_DEFINE_BINARY_FUNCTOR_GENERIC( Fname, Expr ) \
	template<typename T> \
	struct binary_op_fun<Fname##_t, scalar_kernel_t, T, T> { typedef Fname##_fun<T> type; };

// for real-number operations

#define LMAT_DEFINE_REAL_UNARY_FUNCTOR( Fname, Expr ) \
	LMAT_DEFINE_UNARY_FUNCTOR_GENERIC( Fname, Expr ) \
	template<> struct unary_op_fun<Fname##_t, scalar_kernel_t, float> { typedef Fname##_fun<float> type; }; \
	template<> struct unary_op_fun<Fname##_t, scalar_kernel_t, double> { typedef Fname##_fun<double> type; };

#define LMAT_DEFINE_REAL_BINARY_FUNCTOR( Fname, Expr ) \
	LMAT_DEFINE_BINARY_FUNCTOR_GENERIC( Fname, Expr ) \
	template<> struct binary_op_fun<Fname##_t, scalar_kernel_t, float, float> { typedef Fname##_fun<float> type; }; \
	template<> struct binary_op_fun<Fname##_t, scalar_kernel_t, double, double> { typedef Fname##_fun<double> type; };

// for logical operations

#define LMAT_DEFINE_LOGICAL_UNARY_FUNCTOR( Fname, Expr ) \
	struct Fname##_fun { \
		typedef bool result_type; \
		LMAT_ENSURE_INLINE Fname##_fun() { } \
		LMAT_ENSURE_INLINE Fname##_fun( Fname##_t ) { } \
		LMAT_ENSURE_INLINE result_type operator() (const bool& x) const { return Expr; } \
	}; \
	struct binary_op_fun<Op, scalar_kernel_t, bool> { typedef Fname##_fun type; };

#define LMAT_DEFINE_LOGICAL_BINARY_FUNCTOR( Fname, Expr ) \
	struct Fname##_fun { \
		typedef bool result_type; \
		LMAT_ENSURE_INLINE Fname##_fun() { } \
		LMAT_ENSURE_INLINE Fname##_fun( Fname##_t ) { } \
		LMAT_ENSURE_INLINE result_type operator() (const bool& x, const bool& y) const { return Expr; } \
	}; \
	struct binary_op_fun<Op, scalar_kernel_t, bool> { typedef Fname##_fun type; };


#endif 




