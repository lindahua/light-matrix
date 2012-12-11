/**
 * @file mateval_fwd.h
 *
 * Forward declarations for matrix evaluation module
 * 
 * @author Dahua Lin 
 */

#ifdef _MSC_VER
#pragma once
#endif

#ifndef LIGHTMAT_MATEVAL_FWD_H_
#define LIGHTMAT_MATEVAL_FWD_H_

#include <light_mat/matrix/matrix_concepts.h>

namespace lmat
{
	// Accessor interfaces

	template<class Derived, typename T> class IVecAccessor;

	template<class Derived, typename T> class IMultiColAccessor;

	// specific accessors

	template<typename T, typename Kind> class const_accessor;
	template<typename T, typename Kind> class continuous_linear_accessor;
	template<typename T, typename Kind> class step_linear_accessor;
	template<typename T, typename Kind> class multi_col_accessor;
	template<typename T, typename Kind> class multi_stepcol_accessor;

	// accessor map

	template<class Mat, typename Kind> struct vec_accessor_map;
	template<class Mat, typename Kind> struct multicol_accessor_map;

}

#endif 
