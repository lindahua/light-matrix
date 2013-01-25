/**
 * @file mat_allany.h
 *
 * @brief all and any reduction on matrices
 *
 * @author Dahua Lin
 */

#ifdef _MSC_VER
#pragma once
#endif

#ifndef LIGHTMAT_MAT_ALLANY_H_
#define LIGHTMAT_MAT_ALLANY_H_

#include "internal/mat_allany_internal.h"
#include <light_mat/matexpr/mat_pred.h>

namespace lmat
{

	/********************************************
	 *
	 *  full reduction
	 *
	 ********************************************/

	namespace internal
	{
		template<class Mat>
		struct allany_policy
		{
			typedef typename preferred_macc_policy<
					typename meta::shape<Mat>::type,
					copy_kernel<typename meta::value_type_of<Mat>::type>,
					Mat>::type type;
		};
	}


	template<typename T, class Mat>
	inline bool all(const IEWiseMatrix<Mat, mask_t<T> >& mat, bool val=true)
	{
		typedef typename internal::allany_policy<Mat>::type policy_t;

		if (val)
			return internal::all_(mat.shape(), type_<T>(), mat, policy_t());
		else
			return !internal::any_(mat.shape(), type_<T>(), mat, policy_t());
	}

	template<class Mat>
	inline bool all(const IEWiseMatrix<Mat, bool>& mat, bool val=true)
	{
		typedef typename internal::allany_policy<Mat>::type policy_t;

		if (val)
			return internal::all_(mat.shape(), type_<bool>(), mat, policy_t());
		else
			return !internal::any_(mat.shape(), type_<bool>(), mat, policy_t());
	}

	template<typename T, class Mat>
	inline bool any(const IEWiseMatrix<Mat, mask_t<T> >& mat, bool val=true)
	{
		typedef typename internal::allany_policy<Mat>::type policy_t;

		if (val)
			return internal::any_(mat.shape(), type_<T>(), mat, policy_t());
		else
			return !internal::all_(mat.shape(), type_<T>(), mat, policy_t());
	}

	template<class Mat>
	inline bool any(const IEWiseMatrix<Mat, bool>& mat, bool val=true)
	{
		typedef typename internal::allany_policy<Mat>::type policy_t;

		if (val)
			return internal::any_(mat.shape(), type_<bool>(), mat, policy_t());
		else
			return !internal::all_(mat.shape(), type_<bool>(), mat, policy_t());
	}


	/********************************************
	 *
	 *  colwise reduction
	 *
	 ********************************************/

	template<typename T, class Mat, class DMat>
	inline void colwise_all(const IEWiseMatrix<Mat, mask_t<T> >& mat,
			IRegularMatrix<DMat, bool>& dmat, bool val=true)
	{
		typedef default_simd_kind kind;
		const bool use_simd = supports_simd<Mat, kind>::value;
		typedef typename std::conditional<use_simd, simd_<kind>, scalar_>::type U;

		LMAT_CHECK_DIMS( dmat.nelems() == mat.ncolumns() )
		internal::colwise_all_(mat.shape(), type_<T>(), mat.derived(), dmat.derived(), val, U());
	}

	template<class Mat, class DMat>
	inline void colwise_all(const IEWiseMatrix<Mat, bool>& mat,
			IRegularMatrix<DMat, bool>& dmat, bool val=true)
	{
		LMAT_CHECK_DIMS( dmat.nelems() == mat.ncolumns() )
		internal::colwise_all_(mat.shape(), type_<bool>(), mat.derived(), dmat.derived(), val, scalar_());
	}

	template<typename T, class Mat, class DMat>
	inline void colwise_any(const IEWiseMatrix<Mat, mask_t<T> >& mat,
			IRegularMatrix<DMat, bool>& dmat, bool val=true)
	{
		typedef default_simd_kind kind;
		const bool use_simd = supports_simd<Mat, kind>::value;
		typedef typename std::conditional<use_simd, simd_<kind>, scalar_>::type U;

		LMAT_CHECK_DIMS( dmat.nelems() == mat.ncolumns() )
		internal::colwise_any_(mat.shape(), type_<T>(), mat.derived(), dmat.derived(), val, U());
	}

	template<class Mat, class DMat>
	inline void colwise_any(const IEWiseMatrix<Mat, bool>& mat,
			IRegularMatrix<DMat, bool>& dmat, bool val=true)
	{

		LMAT_CHECK_DIMS( dmat.nelems() == mat.ncolumns() )
		internal::colwise_any_(mat.shape(), type_<bool>(), mat.derived(), dmat.derived(), val, scalar_());
	}
}

#endif
