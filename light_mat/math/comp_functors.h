/*
 * @file comp_functors.h
 *
 * Element-wise comparison functors
 *
 * @author Dahua Lin
 */

#ifndef _MSC_VER
#pragma once
#endif

#ifndef LIGHTMAT_COMP_FUNCTORS_H_
#define LIGHTMAT_COMP_FUNCTORS_H_

#include <light_mat/math/functor_base.h>

namespace lmat
{
	template<typename T>
	struct eq_op : public binary_predicate_ewise_functor<T>
	{
		LMAT_ENSURE_INLINE bool operator() (const T& x, const T& y) const
		{
			return x == y;
		}
	};

	template<typename T>
	struct ne_op : public binary_predicate_ewise_functor<T>
	{
		LMAT_ENSURE_INLINE bool operator() (const T& x, const T& y) const
		{
			return x != y;
		}
	};

	template<typename T>
	struct lt_op : public binary_predicate_ewise_functor<T>
	{
		LMAT_ENSURE_INLINE bool operator() (const T& x, const T& y) const
		{
			return x < y;
		}
	};

	template<typename T>
	struct le_op : public binary_predicate_ewise_functor<T>
	{
		LMAT_ENSURE_INLINE bool operator() (const T& x, const T& y) const
		{
			return x <= y;
		}
	};


	template<typename T>
	struct gt_op : public binary_predicate_ewise_functor<T>
	{
		LMAT_ENSURE_INLINE bool operator() (const T& x, const T& y) const
		{
			return x > y;
		}
	};

	template<typename T>
	struct ge_op : public binary_predicate_ewise_functor<T>
	{
		LMAT_ENSURE_INLINE bool operator() (const T& x, const T& y) const
		{
			return x >= y;
		}
	};


	// declaration as ewise functors

	LMAT_DECLARE_AS_BINARY_EWISE_TFUNCTOR( eq_op, false )
	LMAT_DECLARE_AS_BINARY_EWISE_TFUNCTOR( ne_op, false )
	LMAT_DECLARE_AS_BINARY_EWISE_TFUNCTOR( ge_op, false )
	LMAT_DECLARE_AS_BINARY_EWISE_TFUNCTOR( gt_op, false )
	LMAT_DECLARE_AS_BINARY_EWISE_TFUNCTOR( le_op, false )
	LMAT_DECLARE_AS_BINARY_EWISE_TFUNCTOR( lt_op, false )

}

#endif 
