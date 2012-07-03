/*
 * @file bool_functors.h
 *
 * Logical mask functors
 *
 * @author Dahua Lin
 */

#ifndef _MSC_VER
#pragma once
#endif

#ifndef LIGHTMAT_MASK_FUNCTORS_H_
#define LIGHTMAT_MASK_FUNCTORS_H_

#include <light_mat/math/functor_base.h>

namespace lmat
{
	template<typename T>
	struct mask_not_op : public unary_mask_ewise_functor<T>
	{
		LMAT_ENSURE_INLINE mask_t<T> operator() (const mask_t<T>& x) const
		{
			return ~x;
		}
	};

	template<typename T>
	struct mask_and_op : public binary_mask_ewise_functor<T>
	{
		LMAT_ENSURE_INLINE mask_t<T> operator() (const mask_t<T>& x, const mask_t<T>& y) const
		{
			return x & y;
		}
	};

	template<typename T>
	struct mask_or_op : public binary_mask_ewise_functor<T>
	{
		LMAT_ENSURE_INLINE mask_t<T> operator() (const mask_t<T>& x, const mask_t<T>& y) const
		{
			return x | y;
		}
	};

	template<typename T>
	struct mask_xor_op : public binary_mask_ewise_functor<T>
	{
		LMAT_ENSURE_INLINE mask_t<T> operator() (const mask_t<T>& x, const mask_t<T>& y) const
		{
			return x ^ y;
		}
	};


	// declaration as ewise functors

	LMAT_DECLARE_AS_UNARY_EWISE_TFUNCTOR ( mask_not_op, true )
	LMAT_DECLARE_AS_BINARY_EWISE_TFUNCTOR( mask_and_op, true )
	LMAT_DECLARE_AS_BINARY_EWISE_TFUNCTOR( mask_or_op,  true )
	LMAT_DECLARE_AS_BINARY_EWISE_TFUNCTOR( mask_xor_op, true )
}

#endif 
