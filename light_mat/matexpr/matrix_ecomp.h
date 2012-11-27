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

#include <light_mat/matexpr/matrix_ewise_eval.h>
#include <light_mat/math/comp_functors.h>

namespace lmat
{
	LMAT_DEFINE_BINARY_MATFUNCTION( operator ==, eq_t )
	LMAT_DEFINE_BINARY_MATFUNCTION( operator !=, ne_t )
	LMAT_DEFINE_BINARY_MATFUNCTION( operator >=, ge_t )
	LMAT_DEFINE_BINARY_MATFUNCTION( operator >,  gt_t )
	LMAT_DEFINE_BINARY_MATFUNCTION( operator <=, le_t )
	LMAT_DEFINE_BINARY_MATFUNCTION( operator <,  lt_t )
}

#endif 
