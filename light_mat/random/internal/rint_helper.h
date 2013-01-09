/**
 * @file rint_helper.h
 *
 * Helpers for processing random integer
 * 
 * @author Dahua Lin 
 */

#ifdef _MSC_VER
#pragma once
#endif

#ifndef LIGHTMAT_RINT_HELPER_H_
#define LIGHTMAT_RINT_HELPER_H_

#include <light_mat/common/prim_types.h>

namespace lmat { namespace internal {

	template<typename TI, bool IsSigned>
	struct _rint_helper;

	template<typename TI>
	struct _rint_helper<TI, false>
	{
		LMAT_ENSURE_INLINE
		static bool is_nonneg(TI x) { return true; }
	};

	template<typename TI>
	struct _rint_helper<TI, true>
	{
		LMAT_ENSURE_INLINE
		static bool is_nonneg(TI x) { return x >= 0; }
	};

	template<typename TI>
	LMAT_ENSURE_INLINE
	inline bool is_nonneg_int(TI x)
	{
		return _rint_helper<TI, std::is_signed<TI>::value>::is_nonneg(x);
	}

} }

#endif 
