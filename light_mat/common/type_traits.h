/**
 * @file type_traits.h
 *
 * @brief Import of TR1 type traits.
 *
 * @author Dahua Lin
 */

#ifdef _MSC_VER
#pragma once
#endif

#ifndef LIGHTMAT_TYPE_TRAITS_H_
#define LIGHTMAT_TYPE_TRAITS_H_

#include <light_mat/config/config.h>

#ifdef LMAT_USE_C11_STDLIB
#include <type_traits>
#else
#include <tr1/type_traits>
#endif

namespace lmat
{
	using LMAT_TR1::is_void;
	using LMAT_TR1::is_integral;
	using LMAT_TR1::is_floating_point;
	using LMAT_TR1::is_array;
	using LMAT_TR1::is_pointer;
	using LMAT_TR1::is_reference;
	using LMAT_TR1::is_member_function_pointer;
	using LMAT_TR1::is_member_object_pointer;
	using LMAT_TR1::is_enum;
	using LMAT_TR1::is_union;
	using LMAT_TR1::is_class;
	using LMAT_TR1::is_function;

	using LMAT_TR1::is_arithmetic;
	using LMAT_TR1::is_fundamental;
	using LMAT_TR1::is_object;
	// using LMAT_TR1::is_scalar;
	using LMAT_TR1::is_compound;
	using LMAT_TR1::is_member_pointer;

	using LMAT_TR1::is_signed;
	using LMAT_TR1::is_unsigned;

	using LMAT_TR1::is_const;
	using LMAT_TR1::is_volatile;
	using LMAT_TR1::is_pod;
	// using LMAT_TR1::is_empty;
	using LMAT_TR1::is_polymorphic;
	using LMAT_TR1::is_abstract;

	using LMAT_TR1::is_same;
	using LMAT_TR1::is_base_of;

	using LMAT_TR1::remove_const;
	using LMAT_TR1::remove_volatile;
	using LMAT_TR1::remove_cv;
	using LMAT_TR1::remove_pointer;
	using LMAT_TR1::remove_reference;

	using LMAT_TR1::add_const;
	using LMAT_TR1::add_volatile;
	using LMAT_TR1::add_cv;

	using LMAT_TR1::alignment_of;


	template<typename S, typename T>
	struct is_implicitly_convertible_int
	{
		static const bool value =
				is_integral<S>::value && is_integral<T>::value &&
				(
					( is_signed<S>::value == is_signed<T>::value && sizeof(S) <= sizeof(T) ) ||
					( is_signed<T>::value && sizeof(S) < sizeof(T))
				);
	};


	template<typename S, typename T>
	struct is_implicitly_convertible_real
	{
		static const bool value =
				(is_same<S, float>::value && (is_same<T, float>::value || is_same<T, double>::value))
				||
				(is_same<S, double>::value && (is_same<T, double>::value));
	};


	/**
	 * Indicates that implicitly converting
	 * value_type from S to D is ok.
	 */
	template<typename S_, typename T_>
	struct is_implicitly_convertible
	{
		typedef typename remove_cv<S_>::type S;
		typedef typename remove_cv<T_>::type T;

		static const bool value =
			is_same<S, T>::value
			||
			is_implicitly_convertible_int<S, T>::value
			||
			is_implicitly_convertible_real<S, T>::value
			||
			(is_integral<S>::value && is_floating_point<T>::value &&
					sizeof(S) <= sizeof(T)
			)
			||
			(is_same<S, bool>::value && is_integral<T>::value);
	};

}

#endif


