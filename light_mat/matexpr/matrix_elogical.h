/*
 * @file matrix_elogical.h
 *
 * Element-wise logical operations
 *
 * @author Dahua Lin
 */

#ifndef _MSC_VER
#pragma once
#endif

#ifndef LIGHTMAT_MATRIX_ELOGICAL_H_
#define LIGHTMAT_MATRIX_ELOGICAL_H_

#include <light_mat/matexpr/matrix_ewise_eval.h>
#include <light_mat/math/logical_functors.h>


namespace lmat
{
	LMAT_DEFINE_UNARY_MATFUNCTION_S ( operator ~, not_t, bool )
	LMAT_DEFINE_BINARY_MATFUNCTION_S( operator &, and_t, bool )
	LMAT_DEFINE_BINARY_MATFUNCTION_S( operator |, or_t,  bool )
	LMAT_DEFINE_BINARY_MATFUNCTION_S( operator ^, xor_t, bool )
}

#endif 
