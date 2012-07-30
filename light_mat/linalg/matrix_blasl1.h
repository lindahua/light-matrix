/**
 * @file matrix_blasl1.h
 *
 * BLAS Level 1 on matrices
 *
 * @author Dahua Lin
 */

#ifdef _MSC_VER
#pragma once
#endif

#ifndef LIGHTMAT_MATRIX_BLASL1_H_
#define LIGHTMAT_MATRIX_BLASL1_H_

#include "bits/matrix_blasl1_internal.h"

namespace lmat { namespace blas {

	template<class Mat>
	LMAT_ENSURE_INLINE
	float asum(const IMatrixXpr<Mat, float>& x)
	{
		typedef typename lmat::detail::unary_blasl1_internal_map<Mat>::type intern_t;
		return intern_t::asum(x);
	}

	template<class Mat>
	LMAT_ENSURE_INLINE
	double asum(const IMatrixXpr<Mat, double>& x)
	{
		typedef typename lmat::detail::unary_blasl1_internal_map<Mat>::type intern_t;
		return intern_t::asum(x);
	}


	template<class XMat, class YMat>
	LMAT_ENSURE_INLINE
	void axpy(const float& a, const IMatrixXpr<XMat, float>& x, IDenseMatrix<YMat, float>& y)
	{
		typedef typename lmat::detail::binary_blasl1_internal_map<XMat, YMat>::type intern_t;
		intern_t::axpy(a, x, y);
	}

	template<class XMat, class YMat>
	LMAT_ENSURE_INLINE
	void axpy(const double& a, const IMatrixXpr<XMat, double>& x, IDenseMatrix<YMat, double>& y)
	{
		typedef typename lmat::detail::binary_blasl1_internal_map<XMat, YMat>::type intern_t;
		intern_t::axpy(a, x, y);
	}


	template<class XMat, class YMat>
	LMAT_ENSURE_INLINE
	float dot(const IMatrixXpr<XMat, float>& x, const IMatrixXpr<YMat, float>& y)
	{
		typedef typename lmat::detail::binary_blasl1_internal_map<XMat, YMat>::type intern_t;
		return intern_t::dot(x, y);
	}

	template<class XMat, class YMat>
	LMAT_ENSURE_INLINE
	double dot(const IMatrixXpr<XMat, double>& x, const IMatrixXpr<YMat, double>& y)
	{
		typedef typename lmat::detail::binary_blasl1_internal_map<XMat, YMat>::type intern_t;
		return intern_t::dot(x, y);
	}


	template<class Mat>
	LMAT_ENSURE_INLINE
	float nrm2(const IMatrixXpr<Mat, float>& x)
	{
		typedef typename lmat::detail::unary_blasl1_internal_map<Mat>::type intern_t;
		return intern_t::nrm2(x);
	}

	template<class Mat>
	LMAT_ENSURE_INLINE
	double nrm2(const IMatrixXpr<Mat, double>& x)
	{
		typedef typename lmat::detail::unary_blasl1_internal_map<Mat>::type intern_t;
		return intern_t::nrm2(x);
	}


	template<class XMat, class YMat>
	LMAT_ENSURE_INLINE
	void rot(IDenseMatrix<XMat, float>& x, IDenseMatrix<YMat, float>& y, const float& c, const float& s)
	{
		typedef typename lmat::detail::binary_blasl1_internal_map<XMat, YMat>::type intern_t;
		intern_t::rot(x, y, c, s);
	}

	template<class XMat, class YMat>
	LMAT_ENSURE_INLINE
	void rot(IDenseMatrix<XMat, double>& x, IDenseMatrix<YMat, double>& y, const double& c, const double& s)
	{
		typedef typename lmat::detail::binary_blasl1_internal_map<XMat, YMat>::type intern_t;
		intern_t::rot(x, y, c, s);
	}

} }

#endif /* MATRIX_BLASL1_H_ */
