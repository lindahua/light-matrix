/**
 * @file mat_enorms.h
 *
 * @brief L-norms of matrices
 *
 * @author Dahua Lin
 */

#ifdef _MSC_VER
#pragma once
#endif

#ifndef LIGHTMAT_MAT_ENORMS_H_
#define LIGHTMAT_MAT_ENORMS_H_

#include <light_mat/mateval/mat_reduce.h>
#include <light_mat/matexpr/mat_arith.h>

namespace lmat
{

	namespace norms
	{
		struct L1_ { };
		struct L2_ { };
		struct Linf_ { };
	}

	// full reduction

	template<typename T, class Mat>
	LMAT_ENSURE_INLINE
	inline T norm(const IEWiseMatrix<Mat, T>& mat, norms::L1_)
	{
		return asum(mat.derived());
	}

	template<typename T, class Mat>
	LMAT_ENSURE_INLINE
	inline T norm(const IEWiseMatrix<Mat, T>& mat, norms::L2_)
	{
		return math::sqrt(sqsum(mat));
	}

	template<typename T, class Mat>
	LMAT_ENSURE_INLINE
	inline T norm(const IEWiseMatrix<Mat, T>& mat, norms::Linf_)
	{
		return amax(mat);
	}


	// colwise reduction

	template<typename T, class Mat, class DMat>
	LMAT_ENSURE_INLINE
	inline void colwise_norm(const IEWiseMatrix<Mat, T>& mat, IRegularMatrix<DMat, T>& dmat, norms::L1_)
	{
		colwise_asum(mat, dmat);
	}

	template<typename T, class Mat, class DMat>
	LMAT_ENSURE_INLINE
	inline void colwise_norm(const IEWiseMatrix<Mat, T>& mat, IRegularMatrix<DMat, T>& dmat, norms::L2_)
	{
		colwise_sqsum(mat, dmat);
		dmat.derived() = sqrt(dmat);
	}

	template<typename T, class Mat, class DMat>
	LMAT_ENSURE_INLINE
	inline void colwise_norm(const IEWiseMatrix<Mat, T>& mat, IRegularMatrix<DMat, T>& dmat, norms::Linf_)
	{
		colwise_amax(mat, dmat);
	}


	// rowwise reduction

	template<typename T, class Mat, class DMat>
	LMAT_ENSURE_INLINE
	inline void rowwise_norm(const IEWiseMatrix<Mat, T>& mat, IRegularMatrix<DMat, T>& dmat, norms::L1_)
	{
		rowwise_asum(mat, dmat);
	}

	template<typename T, class Mat, class DMat>
	LMAT_ENSURE_INLINE
	inline void rowwise_norm(const IEWiseMatrix<Mat, T>& mat, IRegularMatrix<DMat, T>& dmat, norms::L2_)
	{
		rowwise_sqsum(mat, dmat);
		dmat.derived() = sqrt(dmat);
	}

	template<typename T, class Mat, class DMat>
	LMAT_ENSURE_INLINE
	inline void rowwise_norm(const IEWiseMatrix<Mat, T>& mat, IRegularMatrix<DMat, T>& dmat, norms::Linf_)
	{
		rowwise_amax(mat, dmat);
	}


}

#endif /* MAT_ENORMS_H_ */
