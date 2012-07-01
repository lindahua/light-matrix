/*
 * @file bool_functors.h
 *
 * Boolean functors
 *
 * @author Dahua Lin
 */

#ifndef _MSC_VER
#pragma once
#endif

#ifndef LIGHTMAT_LOGICAL_FUNCTORS_H_
#define LIGHTMAT_LOGICAL_FUNCTORS_H_

#include <light_mat/math/functor_base.h>

namespace lmat
{
	struct not_op : public unary_bool_ewise_functor
	{
		LMAT_ENSURE_INLINE bool operator() (bool x) const
		{
			return !x;
		}
	};

	struct and_op : public binary_bool_ewise_functor
	{
		LMAT_ENSURE_INLINE bool operator() (bool x, bool y) const
		{
			return x && y;
		}
	};

	struct or_op : public binary_bool_ewise_functor
	{
		LMAT_ENSURE_INLINE bool operator() (bool x, bool y) const
		{
			return x || y;
		}
	};


	// declaration as ewise functors

	LMAT_DECLARE_AS_UNARY_EWISE_FUNCTOR( not_op, false )
	LMAT_DECLARE_AS_BINARY_EWISE_FUNCTOR( and_op, false )
	LMAT_DECLARE_AS_BINARY_EWISE_FUNCTOR( or_op, false )
}

#endif 
