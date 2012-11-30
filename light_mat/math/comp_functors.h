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
	// declarations

	LMAT_DECLARE_NUMERIC_BINARY_OP_S( eq_t, bool )
	LMAT_DECLARE_NUMERIC_BINARY_OP_S( ne_t, bool )

	LMAT_DECLARE_NUMERIC_BINARY_OP_S( gt_t, bool )
	LMAT_DECLARE_NUMERIC_BINARY_OP_S( ge_t, bool )

	LMAT_DECLARE_NUMERIC_BINARY_OP_S( lt_t, bool )
	LMAT_DECLARE_NUMERIC_BINARY_OP_S( le_t, bool )

	/********************************************
	 *
	 *  scalar functors
	 *
	 ********************************************/

	LMAT_DEFINE_NUMERIC_BINARY_FUNCTOR( eq, x == y )
	LMAT_DEFINE_NUMERIC_BINARY_FUNCTOR( ne, x != y )
	LMAT_DEFINE_NUMERIC_BINARY_FUNCTOR( gt, x > y )
	LMAT_DEFINE_NUMERIC_BINARY_FUNCTOR( ge, x >= y )
	LMAT_DEFINE_NUMERIC_BINARY_FUNCTOR( lt, x < y )
	LMAT_DEFINE_NUMERIC_BINARY_FUNCTOR( le, x <= y )

}

#endif 
