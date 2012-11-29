/**
 * @file tie.h
 *
 * The devices to tie multiple output arguments
 * 
 * @author Dahua Lin 
 */

#ifdef _MSC_VER
#pragma once
#endif

#ifndef LIGHTMAT_TIE_H_
#define LIGHTMAT_TIE_H_

#include <light_mat/common/prim_types.h>

namespace lmat
{
	template<typename T1, typename T2, typename T3=nil_t> struct tie_t;

	template<typename T1, typename T2, typename T3>
	struct tie_t
	{
		T1& arg1;
		T2& arg2;
		T3& arg3;

		LMAT_ENSURE_INLINE
		tie_t(T1& a1, T2& a2, T3& a3)
		: arg1(a1), arg2(a2), arg3(a3) { }
	};

	template<typename T1, typename T2>
	struct tie_t<T1, T2, nil_t>
	{
		T1& arg1;
		T2& arg2;

		LMAT_ENSURE_INLINE
		tie_t(T1& a1, T2& a2)
		: arg1(a1), arg2(a2) { }
	};

}

#endif 
