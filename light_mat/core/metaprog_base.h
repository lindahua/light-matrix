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

namespace lmat
{
	struct nil_type { };

	struct true_t
	{
		static const bool value = true;
	};

	struct false_t
	{
		static const bool value = false;
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

	// and, or

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

}

#endif
