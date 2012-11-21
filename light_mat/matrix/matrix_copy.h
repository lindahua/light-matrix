/*
 * @file matrix_copy.h
 *
 * Functions to copy or cast matrices
 *
 * @author Dahua Lin
 */

#ifndef _MSC_VER
#pragma once
#endif

#ifndef LIGHTMAT_MATRIX_COPY_H_
#define LIGHTMAT_MATRIX_COPY_H_


#include "bits/matrix_copy_internal.h"

namespace lmat
{
	template<typename T, class RMat>
	LMAT_ENSURE_INLINE
	void copy(const T *ps, IDenseMatrix<RMat, T>& dst)
	{
		typedef typename detail::mat_copier_p_map<RMat>::type copier_t;
		copier_t::copy(ps, dst.derived());
	}

	template<typename T, class LMat>
	LMAT_ENSURE_INLINE
	void copy(const IDenseMatrix<LMat, T>& src, T* pd)
	{
		typedef typename detail::mat_copier_p_map<LMat>::type copier_t;
		copier_t::copy(src.derived(), pd);
	}


	template<typename T, class LMat, class RMat>
	inline
	void copy(const IDenseMatrix<LMat, T>& src, IDenseMatrix<RMat, T>& dst)
	{
		LMAT_CHECK_SAME_SHAPE(src, dst)

		typedef typename detail::mat_copier_map<LMat, RMat>::type copier_t;
		copier_t::copy(src.derived(), dst.derived());
	}


	template<typename T, class DMat>
	LMAT_ENSURE_INLINE
	DMat& operator << (IDenseMatrix<DMat, T>& dmat, const T *src)
	{
		copy(src, dmat.derived());
		return dmat.derived();
	}


	template<class SMat, class DMat>
	struct matrix_copy_scheme
	{
		const SMat& smat;
		DMat& dmat;

		LMAT_ENSURE_INLINE
		matrix_copy_scheme(const SMat& sm, DMat& dm)
		: smat(sm), dmat(dm) { }

		LMAT_ENSURE_INLINE
		void evaluate()
		{
			typedef typename detail::mat_copier_map<SMat, DMat>::type copier_t;
			copier_t::copy(smat, dmat);
		}
	};


}

#endif 
