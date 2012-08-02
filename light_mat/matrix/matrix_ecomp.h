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
	 *  Specific Expressions
	 *
	 ********************************************/

	LMAT_DEFINE_BINARY_MATFUNCTION_EX( operator ==, eq_op )
	LMAT_DEFINE_BINARY_MATFUNCTION_EX( operator !=, ne_op )
	LMAT_DEFINE_BINARY_MATFUNCTION_EX( operator <, lt_op )
	LMAT_DEFINE_BINARY_MATFUNCTION_EX( operator <=, le_op )
	LMAT_DEFINE_BINARY_MATFUNCTION_EX( operator >, gt_op )
	LMAT_DEFINE_BINARY_MATFUNCTION_EX( operator >=, ge_op )
}

#endif 
