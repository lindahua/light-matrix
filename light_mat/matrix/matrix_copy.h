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

#include "internal/matrix_copy_internal.h"

namespace lmat
{
	template<typename T, class RMat>
	LMAT_ENSURE_INLINE
	inline void copy(const T *ps, IRegularMatrix<RMat, T>& dst)
	{
		internal::copy(ps, dst.derived(), internal::get_copy_scheme(dst));
	}

	template<typename T, class LMat>
	LMAT_ENSURE_INLINE
	inline void copy(const IRegularMatrix<LMat, T>& src, T* pd)
	{
		internal::copy(src.derived(), pd, internal::get_copy_scheme(src));
	}


	template<typename T, class LMat, class RMat>
	LMAT_ENSURE_INLINE
	inline void copy(const IRegularMatrix<LMat, T>& src, IRegularMatrix<RMat, T>& dst)
	{
		internal::copy(src.derived(), dst.derived(), internal::get_copy_scheme(src, dst));
	}

	template<typename T, class DMat>
	LMAT_ENSURE_INLINE
	inline DMat& operator << (IRegularMatrix<DMat, T>& dmat, const T *src)
	{
		copy(src, dmat.derived());
		return dmat.derived();
	}

	template<typename T, class SExpr, class DMat>
	LMAT_ENSURE_INLINE
	inline typename meta::enable_if< meta::is_mat_assignable<SExpr, DMat>, void>::type
	evaluate(const IRegularMatrix<SExpr, T>& sexpr, IRegularMatrix<DMat, T>& dmat)
	{
		copy(sexpr.derived(), dmat.derived());
	}

}

#endif 
