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
	 *  Meta functions
	 *
	 ********************************************/

	// basic meta types

	template<typename T>
	struct type_
	{
		typedef T type;
	};

	typedef std::true_type true_;
	typedef std::false_type false_;

	template<bool B>
	struct bool_ : public std::integral_constant<bool, B> { };

	template<int I>
	struct int_ : public std::integral_constant<int, I> { };

	template<typename W>
	struct id_
	{
		typedef typename W::type type;
	};

	// integer arithmetics

	template<typename W1, typename W2>
	struct add_ : public int_<W1::value + W2::value> { };

	template<typename W1, typename W2>
	struct sub_ : public int_<W1::value - W2::value> { };

	template<typename W1, typename W2>
	struct mul_ : public int_<W1::value * W2::value> { };

	template<typename W1, typename W2>
	struct max_ : public int_<(W1::value > W2::value) ? W1::value : W2::value> { };

	template<typename W1, typename W2>
	struct min_ : public int_<(W1::value < W2::value) ? W1::value : W2::value> { };

	// logical operations

	// specific treatment to achieve short-circuit effect for and_ | or_

	namespace internal
	{
		template<bool L, typename R> struct _and_;

		template<typename R>
		struct _and_<true, R> : public R { };

		template<typename R>
		struct _and_<false, R> : public false_ { };

		template<bool L, typename R> struct _or_;

		template<typename R>
		struct _or_<true, R> : public true_ { };

		template<typename R>
		struct _or_<false, R> : public R { };
	}


	template<typename W>
	struct not_ : public bool_<!W::value> { };

	template<typename W1, typename W2>
	struct and_ : public internal::_and_<W1::value, W2> { };

	template<typename W1, typename W2>
	struct or_ : public internal::_or_<W1::value, W2> { };


	// comparison operations

	template<typename W1, typename W2>
	struct eq_ : public bool_<W1::value == W2::value> { };

	template<typename W1, typename W2>
	struct ne_ : public bool_<W1::value != W2::value> { };

	template<typename W1, typename W2>
	struct ge_ : public bool_<(W1::value >= W2::value)> { };

	template<typename W1, typename W2>
	struct gt_ : public bool_<(W1::value > W2::value)> { };

	template<typename W1, typename W2>
	struct le_ : public bool_<(W1::value <= W2::value)> { };

	template<typename W1, typename W2>
	struct lt_ : public bool_<(W1::value < W2::value)> { };


	/********************************************
	 *
	 *  Conditional statements
	 *
	 ********************************************/

	template<typename Cond, typename T_true, typename T_false>
	struct if_ : public std::conditional<Cond::value, T_true, T_false> { };

	template<typename Cond, typename T=void>
	struct enable_if_ : public std::enable_if<Cond::value, T> { };

	typedef true_ otherwise_;

	template<
		typename C1, typename T1,
		typename C2=otherwise_, typename T2=void,
		typename C3=otherwise_, typename T3=void,
		typename C4=otherwise_, typename T4=void,
		typename C5=otherwise_, typename T5=void,
		typename C6=otherwise_, typename T6=void>
	struct select_
	{
		typedef
			typename std::conditional<C1::value, T1,
			typename std::conditional<C2::value, T2,
			typename std::conditional<C3::value, T3,
			typename std::conditional<C4::value, T4,
			typename std::conditional<C5::value, T5,
			typename std::conditional<C6::value, T6, void
			>::type >::type >::type >::type >::type >::type
			type;
	};

	template<typename... T> struct common_;

	template<typename T>
	struct common_<T>
	{
		typedef T type;
	};

	template<typename T>
	struct common_<T, T>
	{
		typedef T type;
	};

	template<typename T1, typename T2>
	struct common_<T1, T2>
	{
		static_assert(std::is_same<T1, T2>::value,
			"T1 and T2 must be the same type");
	};

	template<typename T0, typename... T>
	struct common_<T0, T...>
	{
		typedef typename
				common_<T0, typename common_<T...>::type >::type type;
	};


	/********************************************
	 *
	 *  folding
	 *
	 ********************************************/

	template<template<typename U1, typename U2> class BOp, typename... W> struct fold_;

	template<template<typename U1, typename U2> class BOp, typename W>
	struct fold_<BOp, W> : public W { };

	template<template<typename U1, typename U2> class BOp, typename W1, typename W2>
	struct fold_<BOp, W1, W2> : public BOp<W1, W2> { };

	template<template<typename U1, typename U2> class BOp, typename W1, typename W2, typename... WR>
	struct fold_<BOp, W1, W2, WR...> : public BOp<BOp<W1, W2>, fold_<BOp, WR...> > { };


	template<typename... W> struct all_ : public fold_<and_, W...> { };

	template<typename... W> struct any_ : public fold_<or_, W...> { };

	template<typename... W> struct sum_ : public fold_<add_, W...> { };

	template<typename... W> struct prod_ : public fold_<mul_, W...> { };

	template<typename... W> struct maximum_ : public fold_<max_, W...> { };

	template<typename... W> struct minimum_ : public fold_<min_, W...> { };

} } // lmat::meta


namespace lmat
{
	// import of some common meta types to lmat namespace

	using meta::type_;

	using meta::int_;

	using meta::bool_;
	using meta::true_;
	using meta::false_;

}


#endif


