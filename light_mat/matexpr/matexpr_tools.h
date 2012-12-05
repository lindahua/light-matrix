/**
 * @file matrix_expr_tools.h
 *
 * @brief Tools for expression construction
 *
 * @author Dahua Lin
 */

#ifndef LIGHTMAT_MATEXPR_TOOLS_H_
#define LIGHTMAT_MATEXPR_TOOLS_H_

#include <light_mat/common/expr_base.h>
#include <light_mat/matrix/matrix_properties.h>

namespace lmat
{

	/********************************************
	 *
	 *  common dimensions
	 *
	 ********************************************/

	template<class QArg1>
	LMAT_ENSURE_INLINE
	index_t common_nrows(
			const tied_forwarder< LMAT_TYPELIST_1(QArg1) >& tfwd )
	{
		return tfwd.arg1_fwd.arg.nrows();
	}


	template<class QArg1, class QArg2>
	LMAT_ENSURE_INLINE
	index_t common_nrows(
			const tied_forwarder< LMAT_TYPELIST_2(QArg1, QArg2) >& tfwd )
	{
		return common_nrows( tfwd.arg1_fwd.arg, tfwd.arg2_fwd.arg );
	}

	template<class QArg1, class QArg2, class QArg3>
	LMAT_ENSURE_INLINE
	index_t common_nrows(
			const tied_forwarder< LMAT_TYPELIST_3(QArg1, QArg2, QArg3) >& tfwd )
	{
		return common_nrows( tfwd.arg1_fwd.arg, tfwd.arg2_fwd.arg, tfwd.arg3_fwd.arg );
	}

	template<class QArg1, class QArg2, class QArg3, class QArg4>
	LMAT_ENSURE_INLINE
	index_t common_nrows(
			const tied_forwarder< LMAT_TYPELIST_4(QArg1, QArg2, QArg3, QArg4) >& tfwd )
	{
		return common_nrows( tfwd.arg1_fwd.arg, tfwd.arg2_fwd.arg, tfwd.arg3_fwd.arg, tfwd.arg4_fwd.arg );
	}


	template<class QArg1>
	LMAT_ENSURE_INLINE
	index_t common_ncols(
			const tied_forwarder< LMAT_TYPELIST_1(QArg1) >& tfwd )
	{
		return tfwd.arg1_fwd.arg.ncolumns();
	}


	template<class QArg1, class QArg2>
	LMAT_ENSURE_INLINE
	index_t common_ncols(
			const tied_forwarder< LMAT_TYPELIST_2(QArg1, QArg2) >& tfwd )
	{
		return common_ncols( tfwd.arg1_fwd.arg, tfwd.arg2_fwd.arg );
	}

	template<class QArg1, class QArg2, class QArg3>
	LMAT_ENSURE_INLINE
	index_t common_ncols(
			const tied_forwarder< LMAT_TYPELIST_3(QArg1, QArg2, QArg3) >& tfwd )
	{
		return common_ncols( tfwd.arg1_fwd.arg, tfwd.arg2_fwd.arg, tfwd.arg3_fwd.arg );
	}

	template<class QArg1, class QArg2, class QArg3, class QArg4>
	LMAT_ENSURE_INLINE
	index_t common_ncols(
			const tied_forwarder< LMAT_TYPELIST_4(QArg1, QArg2, QArg3, QArg4) >& tfwd )
	{
		return common_ncols( tfwd.arg1_fwd.arg, tfwd.arg2_fwd.arg, tfwd.arg3_fwd.arg, tfwd.arg4_fwd.arg );
	}

}

#endif /* MATRIX_EXPR_TOOLS_H_ */
