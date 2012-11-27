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
 *  Useful macros
 *
 ************************************************/

// for numeric operations

#define LMAT_DEFINE_NUMERIC_UNARY_OP(Op) \
	struct Op { }; \
	template<> struct is_unary_op<Op> { static const bool value = true; }; \
	template<typename T> \
	struct unary_op_result<Op, T> { typedef T type; };

#define LMAT_DEFINE_NUMERIC_BINARY_OP(Op) \
	struct Op { }; \
	template<> struct is_binary_op<Op> { static const bool value = true; }; \
	template<typename T> \
	struct binary_op_result<Op, T, T> { typedef T type; };

#define LMAT_DEFINE_NUMERIC_UNARY_FUNMAP(Op, Ker, TFun) \
	template<typename T> \
	struct unary_op_fun<Op, scalar_kernel_t, T> { typedef TFun<T> type; };

#define LMAT_DEFINE_NUMERIC_BINARY_FUNMAP(Op, Ker, TFun) \
	template<typename T> \
	struct binary_op_fun<Op, scalar_kernel_t, T, T> { typedef TFun<T> type; };


// for real-number operations

#define LMAT_DEFINE_REAL_UNARY_OP(Op) \
	struct Op { }; \
	template<> struct is_unary_op<Op> { static const bool value = true; }; \
	template<> struct unary_op_result<Op, float> { typedef float type; }; \
	template<> struct unary_op_result<Op, double> { typedef double type; };

#define LMAT_DEFINE_REAL_BINARY_OP(Op) \
	struct Op { }; \
	template<> struct is_binary_op<Op> { static const bool value = true; }; \
	template<> struct binary_op_result<Op, float, float> { typedef float type; }; \
	template<> struct binary_op_result<Op, double, double> { typedef double type; };

#define LMAT_DEFINE_REAL_UNARY_FUNMAP(Op, Ker, TFun) \
	template<> struct unary_op_fun<Op, scalar_kernel_t, float> { typedef TFun<float> type; }; \
	template<> struct unary_op_fun<Op, scalar_kernel_t, double> { typedef TFun<double> type; };

#define LMAT_DEFINE_REAL_BINARY_FUNMAP(Op, Ker, TFun) \
	template<> struct binary_op_fun<Op, scalar_kernel_t, float, float> { typedef TFun<float> type; }; \
	template<> struct binary_op_fun<Op, scalar_kernel_t, double, double> { typedef TFun<double> type; };

// for general predicates

#define LMAT_DEFINE_UNARY_PRED_OP(Op) \
	struct Op { }; \
	template<> struct is_unary_op<Op> { static const bool value = true; }; \
	template<typename T> \
	struct unary_op_result<Op, T> { typedef bool type; };

#define LMAT_DEFINE_BINARY_PRED_OP(Op) \
	struct Op { }; \
	template<> struct is_binary_op<Op> { static const bool value = true; }; \
	template<typename T> \
	struct binary_op_result<Op, T, T> { typedef bool type; };

#define LMAT_DEFINE_UNARY_PRED_FUNMAP(Op, Ker, TFun) \
	template<typename T> \
	struct unary_op_fun<Op, scalar_kernel_t, T> { typedef TFun<T> type; };

#define LMAT_DEFINE_BINARY_PRED_FUNMAP(Op, Ker, TFun) \
	template<typename T> \
	struct binary_op_fun<Op, scalar_kernel_t, T, T> { typedef TFun<T> type; };


// for logical operations

#define LMAT_DEFINE_LOGICAL_UNARY_OP(Op) \
	struct Op { }; \
	template<> struct is_unary_op<Op> { static const bool value = true; }; \
	template<> struct unary_op_result<Op, bool> { typedef bool type; };

#define LMAT_DEFINE_LOGICAL_BINARY_OP(Op) \
	struct Op { }; \
	template<> struct is_binary_op<Op> { static const bool value = true; }; \
	template<> struct binary_op_result<Op, bool, bool> { typedef bool type; };

#define LMAT_DEFINE_LOGICAL_UNARY_FUNMAP(Op, Ker, TFun) \
	template<> struct unary_op_fun<Op, scalar_kernel_t, bool> { typedef TFun type; };

#define LMAT_DEFINE_LOGICAL_BINARY_FUNMAP(Op, Ker, TFun) \
	template<> struct binary_op_fun<Op, scalar_kernel_t, bool, bool> { typedef TFun type; };



#endif 




