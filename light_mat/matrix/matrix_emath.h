/**
 * @file matrix_emath.h
 *
 * Element-wise math functions on matrix
 *
 * @author Dahua Lin
 */

#ifdef _MSC_VER
#pragma once
#endif

#ifndef LIGHTMAT_MATRIX_EMATH_H_
#define LIGHTMAT_MATRIX_EMATH_H_

#include <light_mat/matrix/matrix_ewise_eval.h>
#include <light_mat/math/emath_functors.h>

namespace lmat
{

	template<typename T, class LMat, class RMat>
	LMAT_ENSURE_INLINE
	inline binary_ewise_expr<max_fun<T>, LMat, RMat>
	(max)(const IMatrixXpr<LMat, T>& A, const IMatrixXpr<RMat, T>& B)
	{
		return ewise_expr(max_fun<T>(), A, B);
	}

	template<typename T, class LMat>
	LMAT_ENSURE_INLINE
	inline binary_fix2_ewise_expr<max_fun<T>, LMat>
	(max)(const IMatrixXpr<LMat, T>& A, const T& b)
	{
		typedef typename const_mat_same_form<LMat>::type cst_t;
		return ewise_expr(max_fun<T>(), A, b);
	}

	template<typename T, class RMat>
	LMAT_ENSURE_INLINE
	inline binary_fix1_ewise_expr<max_fun<T>, RMat>
	(max)(const T& a, const IMatrixXpr<RMat, T>& B)
	{
		typedef typename const_mat_same_form<RMat>::type cst_t;
		return ewise_expr(max_fun<T>(), a, B);
	}


	template<typename T, class LMat, class RMat>
	LMAT_ENSURE_INLINE
	inline binary_ewise_expr<min_fun<T>, LMat, RMat>
	(min)(const IMatrixXpr<LMat, T>& A, const IMatrixXpr<RMat, T>& B)
	{
		return ewise_expr(min_fun<T>(), A, B);
	}

	template<typename T, class LMat>
	LMAT_ENSURE_INLINE
	inline binary_fix2_ewise_expr<min_fun<T>, LMat>
	(min)(const IMatrixXpr<LMat, T>& A, const T& b)
	{
		typedef typename const_mat_same_form<LMat>::type cst_t;
		return ewise_expr(min_fun<T>(), A, b);
	}

	template<typename T, class RMat>
	LMAT_ENSURE_INLINE
	inline binary_fix1_ewise_expr<min_fun<T>, RMat>
	(min)(const T& a, const IMatrixXpr<RMat, T>& B)
	{
		typedef typename const_mat_same_form<RMat>::type cst_t;
		return ewise_expr(min_fun<T>(), a, B);
	}


	template<typename T, class LMat, class RMat>
	LMAT_ENSURE_INLINE
	inline binary_ewise_expr<pow_fun<T>, LMat, RMat>
	pow(const IMatrixXpr<LMat, T>& A, const IMatrixXpr<RMat, T>& B)
	{
		return ewise_expr(pow_fun<T>(), A, B);
	}

	template<typename T, class LMat>
	LMAT_ENSURE_INLINE
	inline binary_fix2_ewise_expr<pow_fun<T>, LMat>
	pow(const IMatrixXpr<LMat, T>& A, const T& b)
	{
		typedef typename const_mat_same_form<LMat>::type cst_t;
		return ewise_expr(pow_fun<T>(), A, b);
	}


	template<typename T, class Mat>
	LMAT_ENSURE_INLINE
	inline unary_ewise_expr<floor_fun<T>, Mat>
	floor(const IMatrixXpr<Mat, T>& A)
	{
		return ewise_expr(floor_fun<T>(), A);
	}

	template<typename T, class Mat>
	LMAT_ENSURE_INLINE
	inline unary_ewise_expr<ceil_fun<T>, Mat>
	ceil(const IMatrixXpr<Mat, T>& A)
	{
		return ewise_expr(ceil_fun<T>(), A);
	}


	template<typename T, class Mat>
	LMAT_ENSURE_INLINE
	inline unary_ewise_expr<exp_fun<T>, Mat>
	exp(const IMatrixXpr<Mat, T>& A)
	{
		return ewise_expr(exp_fun<T>(), A);
	}

	template<typename T, class Mat>
	LMAT_ENSURE_INLINE
	inline unary_ewise_expr<log_fun<T>, Mat>
	log(const IMatrixXpr<Mat, T>& A)
	{
		return ewise_expr(log_fun<T>(), A);
	}

	template<typename T, class Mat>
	LMAT_ENSURE_INLINE
	inline unary_ewise_expr<log10_fun<T>, Mat>
	log10(const IMatrixXpr<Mat, T>& A)
	{
		return ewise_expr(log10_fun<T>(), A);
	}


	template<typename T, class Mat>
	LMAT_ENSURE_INLINE
	inline unary_ewise_expr<sin_fun<T>, Mat>
	sin(const IMatrixXpr<Mat, T>& A)
	{
		return ewise_expr(sin_fun<T>(), A);
	}

	template<typename T, class Mat>
	LMAT_ENSURE_INLINE
	inline unary_ewise_expr<cos_fun<T>, Mat>
	cos(const IMatrixXpr<Mat, T>& A)
	{
		return ewise_expr(cos_fun<T>(), A);
	}

	template<typename T, class Mat>
	LMAT_ENSURE_INLINE
	inline unary_ewise_expr<tan_fun<T>, Mat>
	tan(const IMatrixXpr<Mat, T>& A)
	{
		return ewise_expr(tan_fun<T>(), A);
	}


	template<typename T, class Mat>
	LMAT_ENSURE_INLINE
	inline unary_ewise_expr<asin_fun<T>, Mat>
	asin(const IMatrixXpr<Mat, T>& A)
	{
		return ewise_expr(asin_fun<T>(), A);
	}

	template<typename T, class Mat>
	LMAT_ENSURE_INLINE
	inline unary_ewise_expr<acos_fun<T>, Mat>
	acos(const IMatrixXpr<Mat, T>& A)
	{
		return ewise_expr(acos_fun<T>(), A);
	}

	template<typename T, class Mat>
	LMAT_ENSURE_INLINE
	inline unary_ewise_expr<atan_fun<T>, Mat>
	atan(const IMatrixXpr<Mat, T>& A)
	{
		return ewise_expr(atan_fun<T>(), A);
	}

	template<typename T, class LMat, class RMat>
	LMAT_ENSURE_INLINE
	inline binary_ewise_expr<atan2_fun<T>, LMat, RMat>
	atan2(const IMatrixXpr<LMat, T>& A, const IMatrixXpr<RMat, T>& B)
	{
		return ewise_expr(atan2_fun<T>(), A, B);
	}

	template<typename T, class Mat>
	LMAT_ENSURE_INLINE
	inline unary_ewise_expr<sinh_fun<T>, Mat>
	sinh(const IMatrixXpr<Mat, T>& A)
	{
		return ewise_expr(sinh_fun<T>(), A);
	}

	template<typename T, class Mat>
	LMAT_ENSURE_INLINE
	inline unary_ewise_expr<cosh_fun<T>, Mat>
	cosh(const IMatrixXpr<Mat, T>& A)
	{
		return ewise_expr(cosh_fun<T>(), A);
	}

	template<typename T, class Mat>
	LMAT_ENSURE_INLINE
	inline unary_ewise_expr<tanh_fun<T>, Mat>
	tanh(const IMatrixXpr<Mat, T>& A)
	{
		return ewise_expr(tanh_fun<T>(), A);
	}

#ifdef LMAT_HAS_CXX11_MATH

	template<typename T, class Mat>
	LMAT_ENSURE_INLINE
	inline unary_ewise_expr<round_fun<T>, Mat>
	round(const IMatrixXpr<Mat, T>& A)
	{
		return ewise_expr(round_fun<T>(), A);
	}

	template<typename T, class Mat>
	LMAT_ENSURE_INLINE
	inline unary_ewise_expr<trunc_fun<T>, Mat>
	trunc(const IMatrixXpr<Mat, T>& A)
	{
		return ewise_expr(trunc_fun<T>(), A);
	}

	template<typename T, class Mat>
	LMAT_ENSURE_INLINE
	inline unary_ewise_expr<cbrt_fun<T>, Mat>
	cbrt(const IMatrixXpr<Mat, T>& A)
	{
		return ewise_expr(cbrt_fun<T>(), A);
	}

	template<typename T, class LMat, class RMat>
	LMAT_ENSURE_INLINE
	inline binary_ewise_expr<hypot_fun<T>, LMat, RMat>
	hypot(const IMatrixXpr<LMat, T>& A, const IMatrixXpr<RMat, T>& B)
	{
		return ewise_expr(hypot_fun<T>(), A, B);
	}


	template<typename T, class Mat>
	LMAT_ENSURE_INLINE
	inline unary_ewise_expr<exp2_fun<T>, Mat>
	exp2(const IMatrixXpr<Mat, T>& A)
	{
		return ewise_expr(exp2_fun<T>(), A);
	}

	template<typename T, class Mat>
	LMAT_ENSURE_INLINE
	inline unary_ewise_expr<log2_fun<T>, Mat>
	log2(const IMatrixXpr<Mat, T>& A)
	{
		return ewise_expr(log2_fun<T>(), A);
	}

	template<typename T, class Mat>
	LMAT_ENSURE_INLINE
	inline unary_ewise_expr<expm1_fun<T>, Mat>
	expm1(const IMatrixXpr<Mat, T>& A)
	{
		return ewise_expr(expm1_fun<T>(), A);
	}

	template<typename T, class Mat>
	LMAT_ENSURE_INLINE
	inline unary_ewise_expr<log1p_fun<T>, Mat>
	log1p(const IMatrixXpr<Mat, T>& A)
	{
		return ewise_expr(log1p_fun<T>(), A);
	}


	template<typename T, class Mat>
	LMAT_ENSURE_INLINE
	inline unary_ewise_expr<asinh_fun<T>, Mat>
	asinh(const IMatrixXpr<Mat, T>& A)
	{
		return ewise_expr(asinh_fun<T>(), A);
	}

	template<typename T, class Mat>
	LMAT_ENSURE_INLINE
	inline unary_ewise_expr<acosh_fun<T>, Mat>
	acosh(const IMatrixXpr<Mat, T>& A)
	{
		return ewise_expr(acosh_fun<T>(), A);
	}

	template<typename T, class Mat>
	LMAT_ENSURE_INLINE
	inline unary_ewise_expr<atanh_fun<T>, Mat>
	atanh(const IMatrixXpr<Mat, T>& A)
	{
		return ewise_expr(atanh_fun<T>(), A);
	}

#endif

}

#endif /* MATRIX_EMATH_H_ */


