/**
 * @file mat_allany.h
 *
 * @brief all and any reduction on matrices
 *
 * @author Dahua Lin
 */

#ifndef LIGHTMAT_MAT_ALLANY_H_
#define LIGHTMAT_MAT_ALLANY_H_

#include "internal/mat_allany_internal.h"
#include <light_mat/mateval/mat_pred.h>

namespace lmat
{

	/********************************************
	 *
	 *  full reduction
	 *
	 ********************************************/

	template<typename T, class Mat>
	inline bool all(const IEWiseMatrix<Mat, mask_t<T> >& mat, bool val=true)
	{
		typedef typename preferred_macc_policy<Mat>::type policy_t;

		if (val)
			return internal::all_(mat.shape(), mat, policy_t());
		else
			return !internal::any_(mat.shape(), mat, policy_t());
	}


	template<class Mat>
	inline bool all(const IEWiseMatrix<Mat, bool>& mat, bool val=true)
	{
		typedef typename preferred_macc_policy<Mat>::type policy_t;

		if (val)
			return internal::all_(mat.shape(), mat, policy_t());
		else
			return !internal::any_(mat.shape(), mat, policy_t());
	}


	template<typename T, class Mat>
	inline bool any(const IEWiseMatrix<Mat, mask_t<T> >& mat, bool val=true)
	{
		typedef typename preferred_macc_policy<Mat>::type policy_t;

		if (val)
			return internal::any_(mat.shape(), mat, policy_t());
		else
			return !internal::all_(mat.shape(), mat, policy_t());
	}


	template<class Mat>
	inline bool any(const IEWiseMatrix<Mat, bool>& mat, bool val=true)
	{
		typedef typename preferred_macc_policy<Mat>::type policy_t;

		if (val)
			return internal::any_(mat.shape(), mat, policy_t());
		else
			return !internal::all_(mat.shape(), mat, policy_t());
	}


}

#endif
