/**
 * @file reductors.h
 *
 * Devices to support reduction
 *
 * @author Dahua Lin
 */

#ifdef _MSC_VER
#pragma once
#endif

#ifndef LIGHTMAT_REDUCTORS_H_
#define LIGHTMAT_REDUCTORS_H_

#include <light_mat/common/arg_check.h>
#include <light_mat/math/functor_base.h>
#include <light_mat/math/math_base.h>

#include <light_mat/math/arith_functors.h>
#include <light_mat/math/emath_functors.h>

#include <limits>


#define LMAT_DECLARE_REDUCTOR( Name ) \
	struct Name##_t { }; \
	template<> struct is_reductor_tag< Name##_t > { static const bool value = true; };

#define LMAT_DEFINE_UNARY_REDUCTOR( Name, RT, TF, CF, PF, EmptyVal ) \
	template<typename T> \
	struct reductor_traits<Name##_t, LMAT_TYPELIST_1(T) > { \
		typedef RT result_type; \
		typedef TF term_fun_tag; \
		typedef CF combine_fun_tag; \
		typedef PF post_fun_tag; \
		static result_type empty() { return EmptyVal; } };

#define LMAT_DEFINE_BINARY_REDUCTOR( Name, RT, TF, CF, PF, EmptyVal ) \
	template<typename T> \
	struct reductor_traits<Name##_t, LMAT_TYPELIST_2(T, T) > { \
		typedef RT result_type; \
		typedef TF term_fun_tag; \
		typedef CF combine_fun_tag; \
		typedef PF post_fun_tag; \
		static result_type empty() { return EmptyVal; } };


namespace lmat
{

	template<typename Tag>
	struct is_reductor_tag
	{
		static const bool value = false;
	};


	template<typename Tag, class QLst>
	struct reductor_traits;

	struct divide_by_dim { };

	/********************************************
	 *
	 *  declarations
	 *
	 ********************************************/

	LMAT_DECLARE_REDUCTOR( sum )
	LMAT_DECLARE_REDUCTOR( maximum )
	LMAT_DECLARE_REDUCTOR( minimum )
	LMAT_DECLARE_REDUCTOR( mean )

	LMAT_DECLARE_REDUCTOR( L1norm )
	LMAT_DECLARE_REDUCTOR( sqL2norm )
	LMAT_DECLARE_REDUCTOR( L2norm )
	LMAT_DECLARE_REDUCTOR( Linfnorm )

	LMAT_DECLARE_REDUCTOR( logsum )
	LMAT_DECLARE_REDUCTOR( entropy )

	LMAT_DECLARE_REDUCTOR( dot )


	/********************************************
	 *
	 *  definitions
	 *
	 ********************************************/

	namespace internal
	{
		template<typename T, bool HasInf> struct _maximum_empty;
		template<typename T, bool HasInf> struct _minimum_empty;

		template<typename T>
		struct _maximum_empty<T, false>
		{
			LMAT_ENSURE_INLINE
			static T get() { return std::numeric_limits<T>::min(); }
		};

		template<typename T>
		struct _maximum_empty<T, true>
		{
			LMAT_ENSURE_INLINE
			static T get() { return - std::numeric_limits<T>::infinity(); }
		};

		template<typename T>
		struct _minimum_empty<T, false>
		{
			LMAT_ENSURE_INLINE
			static T get() { return std::numeric_limits<T>::max(); }
		};

		template<typename T>
		struct _minimum_empty<T, true>
		{
			LMAT_ENSURE_INLINE
			static T get() { return std::numeric_limits<T>::infinity(); }
		};

		template<typename T>
		struct maximum_empty
		{
			LMAT_ENSURE_INLINE
			static T get()
			{
				return _maximum_empty<T, is_floating_point<T>::value>::get();
			}
		};

		template<typename T>
		struct minimum_empty
		{
			LMAT_ENSURE_INLINE
			static T get()
			{
				return _minimum_empty<T, is_floating_point<T>::value>::get();
			}
		};

	}


	/********************************************
	 *
	 *  definitions
	 *
	 ********************************************/

	LMAT_DEFINE_UNARY_REDUCTOR( sum, 		T, id_t, add_t, id_t, T(0))
	LMAT_DEFINE_UNARY_REDUCTOR( maximum, 	T, id_t, max_t, id_t, internal::maximum_empty<T>::get())
	LMAT_DEFINE_UNARY_REDUCTOR( minimum, 	T, id_t, min_t, id_t, internal::minimum_empty<T>::get())
	LMAT_DEFINE_UNARY_REDUCTOR( mean, 		T, id_t, add_t, divide_by_dim, std::numeric_limits<T>::quiet_NaN())

	LMAT_DEFINE_UNARY_REDUCTOR( L1norm, 	T, abs_t, add_t, id_t, T(0))
	LMAT_DEFINE_UNARY_REDUCTOR( sqL2norm, 	T, sqr_t, add_t, id_t, T(0))
	LMAT_DEFINE_UNARY_REDUCTOR( L2norm, 	T, sqr_t, add_t, sqrt_t, T(0))
	LMAT_DEFINE_UNARY_REDUCTOR( Linfnorm, 	T, abs_t, max_t, id_t, T(0))

	LMAT_DEFINE_UNARY_REDUCTOR( logsum, 	T, log_t, add_t, id_t, T(0))
	LMAT_DEFINE_UNARY_REDUCTOR( entropy, 	T, xlogx_t, add_t, negate_t, T(0))

	LMAT_DEFINE_BINARY_REDUCTOR( dot, 		T, multiply_t, add_t, id_t, T(0) )

}

#endif





