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

namespace lmat
{

	template<typename... QArgs> class map_expr;

	template<typename QArg> class transpose_expr;

}


#endif 
