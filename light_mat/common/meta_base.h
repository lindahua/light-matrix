/**
 * @file meta_base.h
 *
 * @brief Basic facilities for meta-programming.
 *
 * @author Dahua Lin
 */

#ifdef _MSC_VER
#pragma once
#endif

#ifndef LIGHTMAT_META_BASE_H_
#define LIGHTMAT_META_BASE_H_

#include <light_mat/common/prim_types.h>

namespace lmat {  namespace meta {


	/********************************************
	 *
	 *  Compile-time evaluation
	 *
	 ********************************************/

	// type function

	template<typename T>
	struct id_
	{
		typedef T type;
	};

	// logical operations

	template<typename CondT>
	struct not_
	{
		static const bool value = !CondT::value;
	};

	template<typename CondT1, typename CondT2>
	struct and_
	{
		static const bool value = CondT1::value && CondT2::value;
	};

	template<typename CondT1, typename CondT2>
	struct or_
	{
		static const bool value = CondT1::value || CondT2::value;
	};

	template<typename CondT1, typename CondT2>
	struct xor_
	{
		static const bool value =
				(CondT1::value && !CondT2::value) ||
				(CondT2::value && !CondT1::value);
	};


	// arithmetic calculations on integer

	template<typename ValT1, typename ValT2>
	struct add_
	{
		static const int value = ValT1::value + ValT2::value;
	};

	template<typename ValT1, typename ValT2>
	struct sub_
	{
		static const int value = ValT1::value - ValT2::value;
	};

	template<typename ValT1, typename ValT2>
	struct mul_
	{
		static const int value = ValT1::value * ValT2::value;
	};

	template<typename ValT1, typename ValT2>
	struct max_
	{
		static const int value = (ValT1::value > ValT2::value ? ValT1::value : ValT2::value);
	};

	template<typename ValT1, typename ValT2>
	struct min_
	{
		static const int value = (ValT1::value < ValT2::value ? ValT1::value : ValT2::value);
	};


	// comparison operations

	template<typename ValT1, typename ValT2>
	struct eq_
	{
		static const bool value = ValT1::value == ValT2::value;
	};

	template<typename ValT1, typename ValT2>
	struct ne_
	{
		static const bool value = ValT1::value != ValT2::value;
	};

	template<typename ValT1, typename ValT2>
	struct gt_
	{
		static const bool value = ValT1::value > ValT2::value;
	};

	template<typename ValT1, typename ValT2>
	struct ge_
	{
		static const bool value = ValT1::value >= ValT2::value;
	};

	template<typename ValT1, typename ValT2>
	struct lt_
	{
		static const bool value = ValT1::value < ValT2::value;
	};

	template<typename ValT1, typename ValT2>
	struct le_
	{
		static const bool value = ValT1::value <= ValT2::value;
	};


	/********************************************
	 *
	 *  Conditional statements
	 *
	 ********************************************/

	// meta if statement

	template<bool Cond, typename Ttrue, typename Tfalse>
	struct if_c;

	template<typename Ttrue, typename Tfalse>
	struct if_c<true, Ttrue, Tfalse>
	{
		typedef Ttrue type;
	};

	template<typename Ttrue, typename Tfalse>
	struct if_c<false, Ttrue, Tfalse>
	{
		typedef Tfalse type;
	};

	template<typename CondT, typename Ttrue, typename Tfalse>
	struct if_
	{
		typedef typename if_c<CondT::value, Ttrue, Tfalse>::type type;
	};

	// enable_if

	template<bool Cond, typename T>
	struct enable_if_c;

	template<typename T>
	struct enable_if_c<true, T>
	{
		typedef T type;
	};

	template<typename T>
	struct enable_if_c<false, T> { };

	template<typename CondT, typename T>
	struct enable_if
	{
		typedef typename enable_if_c<CondT::value, T>::type type;
	};

	template<bool Cond, int Val>
	struct enable_int_if_c;

	template<int Val>
	struct enable_int_if_c<true, Val>
	{
		static const int value = Val;
	};

	template<int Val>
	struct enable_int_if_c<false, Val> { };

	template<typename CondT, int Val>
	struct enable_int_if
	{
		static const int value = enable_int_if_c<CondT::value, Val>::value;
	};


	/********************************************
	 *
	 *  compile-time list
	 *
	 ********************************************/

	template<typename... T>
	struct type_list { };

	template<class List> struct len_;

	template<typename T0, typename... T>
	struct len_<type_list<T0, T...> >
	{
		static const int value = 1 + len_<type_list<T...> >::value;
	};

	template<>
	struct len_<type_list<> >
	{
		static const int value = 0;
	};

	template<class List, int I> struct get_;

	template<typename T0, typename... T>
	struct get_<type_list<T0, T...>, 0>
	{
		typedef T0 type;
	};

	template<typename T0, typename... T, int I>
	struct get_<type_list<T0, T...>, I>
	{
		typedef typename get_<type_list<T...>, I-1>::type type;
	};


	template<template<typename U> class Fun, class List> struct map_;

	template<template<typename U> class Fun, typename... T>
	struct map_<Fun, type_list<T...> >
	{
		typedef type_list<typename Fun<T>::type ...> type;
	};


} } // lmat::meta


#endif


