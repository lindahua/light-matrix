/**
 * @file mat_reduce.h
 *
 * Reduction on matrices
 * 
 * @author Dahua Lin 
 */

#ifdef _MSC_VER
#pragma once
#endif

#ifndef LIGHTMAT_MAT_REDUCE_H_
#define LIGHTMAT_MAT_REDUCE_H_

#include "internal/mat_reduce_internal.h"

namespace lmat
{

	template<typename T, class Mat>
	LMAT_ENSURE_INLINE
	inline T sum(const IRegularMatrix<Mat, T>& mat)
	{
		static_assert(meta::is_continuous<Mat>::value,
				"mat should be a compile-time-continuous matrix.");

		unsigned int len = (unsigned int)mat.nelems();
		T r;

		if (len > 0)
		{
			r = internal::sum_impl((unsigned int)mat.nelems(),
					contvec_reader<T, atags::simd<T, default_simd_kind> >(mat.ptr_data()));
		}
		else
		{
			r = T(0);
		}

		return r;
	}



}

#endif 
