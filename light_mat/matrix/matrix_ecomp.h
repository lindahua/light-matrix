/*
 * @file matrix_ecomp.h
 *
 * Element-wise comparison operations
 *
 * @author Dahua Lin
 */

#ifndef _MSC_VER
#pragma once
#endif

#ifndef LIGHTMAT_MATRIX_ECOMP_H_
#define LIGHTMAT_MATRIX_ECOMP_H_

#include <light_mat/matrix/matrix_ewise_expr.h>
#include <light_mat/math/comp_functors.h>

namespace lmat
{
	/********************************************
	 *
	 *  Expression Type mapping
	 *
	 ********************************************/

	LMAT_DECLARE_BINARY_TYPE_MAP_EX( eq, eq_op )
	LMAT_DECLARE_BINARY_TYPE_MAP_EX( ne, ne_op )
	LMAT_DECLARE_BINARY_TYPE_MAP_EX( lt, lt_op )
	LMAT_DECLARE_BINARY_TYPE_MAP_EX( le, le_op )
	LMAT_DECLARE_BINARY_TYPE_MAP_EX( gt, gt_op )
	LMAT_DECLARE_BINARY_TYPE_MAP_EX( ge, ge_op )

	/********************************************
	 *
	 *  Specific Expressions
	 *
	 ********************************************/

	LMAT_DEFINE_BINARY_MATFUNCTION_EX( eq, operator ==, eq_op )
	LMAT_DEFINE_BINARY_MATFUNCTION_EX( ne, operator !=, ne_op )
	LMAT_DEFINE_BINARY_MATFUNCTION_EX( lt, operator <, lt_op )
	LMAT_DEFINE_BINARY_MATFUNCTION_EX( le, operator <=, le_op )
	LMAT_DEFINE_BINARY_MATFUNCTION_EX( gt, operator >, gt_op )
	LMAT_DEFINE_BINARY_MATFUNCTION_EX( ge, operator >=, ge_op )
}

#endif 
