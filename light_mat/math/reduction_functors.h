/**
 * @file reduction_functors.h
 *
 * Basic reduction functors
 *
 * @author Dahua Lin
 */

#ifdef _MSC_VER
#pragma once
#endif

#ifndef LIGHTMAT_REDUCTION_FUNCTORS_H_
#define LIGHTMAT_REDUCTION_FUNCTORS_H_

#include <light_mat/common/arg_check.h>
#include <light_mat/math/functor_base.h>
#include <light_mat/math/math_base.h>

namespace lmat
{
	template<typename Op>
	struct is_unary_reduc_op
	{
		static const bool value = false;
	};

	template<typename Op>
	struct is_binary_reduc_op
	{
		static const bool value = false;
	};

	template<typename Op, typename T>
	struct unary_reduc_result;

	template<typename Op, typename T1, typename T2>
	struct binary_reduc_result;

	template<typename Op, typename Ker, typename T>
	struct unary_reduc_fun;

	template<typename Op, typename Ker, typename T1, typename T2>
	struct binary_reduc_fun;

}

// useful macros

#define LMAT_DEFINE_UNARY_NUMERIC_REDUC_OP(Op) \
	struct Op { }; \
	template<> struct is_unary_reduc_op<Op> { static const bool value = true; }; \
	template<typename T> \
	struct unary_reduc_result<Op, T> { typedef T type; };

#define LMAT_DEFINE_BINARY_NUMERIC_REDUC_OP(Op) \
	struct Op { }; \
	template<> struct is_binary_reduc_op<Op> { static const bool value = true; }; \
	template<typename T> \
	struct binary_reduc_result<Op, T, T> { typedef T type; };

#define LMAT_DEFINE_UNARY_REAL_REDUC_OP(Op) \
	struct Op { }; \
	template<> struct is_unary_reduc_op<Op> { static const bool value = true; }; \
	struct unary_reduc_result<Op, float> { typedef float type; }; \
	struct unary_reduc_result<Op, double> { typedef double type; };

#define LMAT_DEFINE_BINARY_REAL_REDUC_OP(Op) \
	struct Op { }; \
	template<> struct is_binary_reduc_op<Op> { static const bool value = true; }; \
	struct binary_reduc_result<Op, float, float> { typedef float type; }; \
	struct binary_reduc_result<Op, double, double> { typedef double type; };

#define LMAT_DEFINE_UNARY_REDUCTOR(Name, EmptyVal, InitExpr, CombExpr ) \
	template<typename T> \
	struct Name##_fun { \
		LMAT_ENSURE_INLINE \
		T empty_value() const { return EmptyVal; } \
		LMAT_ENSURE_INLINE \
		T transform(const T& x) const { return InitExpr; } \
		LMAT_ENSURE_INLINE \
		T combine(const T& x, const T& y) const { return CombExpr; } \
	}; \
	template<typename T> \
	struct unary_reduc_fun<Name##_t, scalar_kernel_t, T> { \
		typedef Name##_fun<T> type; \
	};

#define LMAT_DEFINE_BINARY_REDUCTOR(Name, EmptyVal, InitExpr, CombExpr ) \
	template<typename T> \
	struct Name##_fun { \
		LMAT_ENSURE_INLINE \
		T empty_value() const { return EmptyVal; } \
		LMAT_ENSURE_INLINE \
		T transform(const T& x, const T& y) const { return InitExpr; } \
		LMAT_ENSURE_INLINE \
		T combine(const T& x, const T& y) const { return CombExpr; } \
	}; \
	template<typename T> \
	struct binary_reduc_fun<Name##_t, scalar_kernel_t, T, T> { \
		typedef Name##_fun<T> type; \
	};


namespace lmat
{
	// operator tags

	LMAT_DEFINE_UNARY_NUMERIC_REDUC_OP( sum_t )
	LMAT_DEFINE_UNARY_NUMERIC_REDUC_OP( maximum_t )
	LMAT_DEFINE_UNARY_NUMERIC_REDUC_OP( minimum_t )

	LMAT_DEFINE_UNARY_REAL_REDUC_OP( L1norm_t )
	LMAT_DEFINE_UNARY_REAL_REDUC_OP( sqL2norm_t )
	LMAT_DEFINE_UNARY_REAL_REDUC_OP( Linfnorm_t )

	LMAT_DEFINE_UNARY_REAL_REDUC_OP( logsum_t )
	LMAT_DEFINE_UNARY_REAL_REDUC_OP( entropy_t )

	LMAT_DEFINE_BINARY_REAL_REDUC_OP( dot_t )

	template<typename T>
	inline T no_empty_value(const char *msg)
	{
		throw invalid_operation(msg);
	};


	/********************************************
	 *
	 *  reduction functors
	 *
	 ********************************************/

	// sum, maximum, and minimum

	LMAT_DEFINE_UNARY_REDUCTOR( sum, T(0), x, x + y )

	LMAT_DEFINE_UNARY_REDUCTOR( maximum,
			no_empty_value<T>("maximum is not allowed for empty array"), x, (math::max)(x, y) )

	LMAT_DEFINE_UNARY_REDUCTOR( minimum,
			no_empty_value<T>("minimum is not allowed for empty array"), x, (math::min)(x, y) )

	// L-norms

	LMAT_DEFINE_UNARY_REDUCTOR( L1norm, T(0), math::abs(x), x + y )
	LMAT_DEFINE_UNARY_REDUCTOR( sqL2norm, T(0), math::sqr(x), x + y )
	LMAT_DEFINE_UNARY_REDUCTOR( Linfnorm, T(0), math::abs(x), (math::max)(x, y) )

	// logsum & entropy

	namespace math
	{
		LMAT_ENSURE_INLINE inline float xlogx(float x) { return x > 0 ? x * math::log(x) : 0.f; }
		LMAT_ENSURE_INLINE inline double xlogx(double x) { return x > 0 ? x * math::log(x) : 0.0; }
	}

	LMAT_DEFINE_UNARY_REDUCTOR( logsum, T(0), math::log(x), x + y )
	LMAT_DEFINE_UNARY_REDUCTOR( entropy, T(0), - math::xlogx(x), x + y )

	// dot

	LMAT_DEFINE_BINARY_REDUCTOR( dot, T(0), x * y, x + y )

}

#endif





