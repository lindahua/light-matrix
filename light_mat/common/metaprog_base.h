/**
 * @file metaprog_base.h
 *
 * @brief Basic facilities for meta-programming.
 *
 * @author Dahua Lin
 */

#ifdef _MSC_VER
#pragma once
#endif

#ifndef LIGHTMAT_METAPROG_BASE_H_
#define LIGHTMAT_METAPROG_BASE_H_

#include <light_mat/config/config.h>
#include <light_mat/common/type_traits.h>

namespace lmat
{
	struct true_t
	{
		static const bool value = true;
	};

	struct false_t
	{
		static const bool value = false;
	};

	template<typename T>
	struct type{ };

	template<int I>
	struct fix_int
	{
		static const int value = I;
	};


	// compile-time logical operations

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


	// common type

	template<typename T1, typename T2>
	struct common_type
	{
		typedef typename enable_if<is_same<T1, T2>, T1>::type type;
	};

	template<int N1, int N2>
	struct common_int
	{
		static const int value = enable_int_if_c<N1 == N2, N1>::value;
	};


}

#endif


