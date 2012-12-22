/**
 * @file matexpr_fwd.h
 *
 * Forward declaration of matrix expressions
 * 
 * @author Dahua Lin 
 */

#ifdef _MSC_VER
#pragma once
#endif

#ifndef LIGHTMAT_MATEXPR_FWD_H_
#define LIGHTMAT_MATEXPR_FWD_H_

#include <light_mat/matrix/matrix_classes.h>

namespace lmat
{
	template<typename Tag, typename... Args> class map_expr;

	template<typename Arg> class transpose_expr;
}


#endif 
