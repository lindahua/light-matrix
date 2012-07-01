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
	inline typename binary_ewise_expr_map<max_fun<T>, LMat, RMat>::type
	(max)(const IMatrixXpr<LMat, T>& A, const IMatrixXpr<RMat, T>& B)
	{
		return ewise(max_fun<T>(), A.derived(), B.derived());
	}

	template<typename T, class LMat>
	LMAT_ENSURE_INLINE
	inline typename binary_fix2_ewise_expr_map<max_fun<T>, LMat>::type
	(max)(const IMatrixXpr<LMat, T>& A, const T& b)
	{
		return ewise(max_fun<T>(), A.derived(), b);
	}

	template<typename T, class RMat>
	LMAT_ENSURE_INLINE
	inline typename binary_fix1_ewise_expr_map<max_fun<T>, RMat>::type
	(max)(const T& a, const IMatrixXpr<RMat, T>& B)
	{
		return ewise(max_fun<T>(), a, B.derived());
	}


	template<typename T, class LMat, class RMat>
	LMAT_ENSURE_INLINE
	inline typename binary_ewise_expr_map<min_fun<T>, LMat, RMat>::type
	(min)(const IMatrixXpr<LMat, T>& A, const IMatrixXpr<RMat, T>& B)
	{
		return ewise(min_fun<T>(), A.derived(), B.derived());
	}

	template<typename T, class LMat>
	LMAT_ENSURE_INLINE
	inline typename binary_fix2_ewise_expr_map<min_fun<T>, LMat>::type
	(min)(const IMatrixXpr<LMat, T>& A, const T& b)
	{
		return ewise(min_fun<T>(), A.derived(), b);
	}

	template<typename T, class RMat>
	LMAT_ENSURE_INLINE
	inline typename binary_fix1_ewise_expr_map<min_fun<T>, RMat>::type
	(min)(const T& a, const IMatrixXpr<RMat, T>& B)
	{
		return ewise(min_fun<T>(), a, B.derived());
	}


	template<typename T, class LMat, class RMat>
	LMAT_ENSURE_INLINE
	inline typename binary_ewise_expr_map<pow_fun<T>, LMat, RMat>::type
	pow(const IMatrixXpr<LMat, T>& A, const IMatrixXpr<RMat, T>& B)
	{
		return ewise(pow_fun<T>(), A.derived(), B.derived());
	}

	template<typename T, class LMat>
	LMAT_ENSURE_INLINE
	inline typename binary_fix2_ewise_expr_map<pow_fun<T>, LMat>::type
	pow(const IMatrixXpr<LMat, T>& A, const T& b)
	{
		return ewise(pow_fun<T>(), A.derived(), b);
	}


	template<typename T, class Mat>
	LMAT_ENSURE_INLINE
	inline typename unary_ewise_expr_map<floor_fun<T>, Mat>::type
	floor(const IMatrixXpr<Mat, T>& A)
	{
		return ewise(floor_fun<T>(), A.derived());
	}

	template<typename T, class Mat>
	LMAT_ENSURE_INLINE
	inline typename unary_ewise_expr_map<ceil_fun<T>, Mat>::type
	ceil(const IMatrixXpr<Mat, T>& A)
	{
		return ewise(ceil_fun<T>(), A.derived());
	}


	template<typename T, class Mat>
	LMAT_ENSURE_INLINE
	inline typename unary_ewise_expr_map<exp_fun<T>, Mat>::type
	exp(const IMatrixXpr<Mat, T>& A)
	{
		return ewise(exp_fun<T>(), A.derived());
	}

	template<typename T, class Mat>
	LMAT_ENSURE_INLINE
	inline typename unary_ewise_expr_map<log_fun<T>, Mat>::type
	log(const IMatrixXpr<Mat, T>& A)
	{
		return ewise(log_fun<T>(), A.derived());
	}

	template<typename T, class Mat>
	LMAT_ENSURE_INLINE
	inline typename unary_ewise_expr_map<log10_fun<T>, Mat>::type
	log10(const IMatrixXpr<Mat, T>& A)
	{
		return ewise(log10_fun<T>(), A.derived());
	}


	template<typename T, class Mat>
	LMAT_ENSURE_INLINE
	inline typename unary_ewise_expr_map<sin_fun<T>, Mat>::type
	sin(const IMatrixXpr<Mat, T>& A)
	{
		return ewise(sin_fun<T>(), A.derived());
	}

	template<typename T, class Mat>
	LMAT_ENSURE_INLINE
	inline typename unary_ewise_expr_map<cos_fun<T>, Mat>::type
	cos(const IMatrixXpr<Mat, T>& A)
	{
		return ewise(cos_fun<T>(), A.derived());
	}

	template<typename T, class Mat>
	LMAT_ENSURE_INLINE
	inline typename unary_ewise_expr_map<tan_fun<T>, Mat>::type
	tan(const IMatrixXpr<Mat, T>& A)
	{
		return ewise(tan_fun<T>(), A.derived());
	}


	template<typename T, class Mat>
	LMAT_ENSURE_INLINE
	inline typename unary_ewise_expr_map<asin_fun<T>, Mat>::type
	asin(const IMatrixXpr<Mat, T>& A)
	{
		return ewise(asin_fun<T>(), A.derived());
	}

	template<typename T, class Mat>
	LMAT_ENSURE_INLINE
	inline typename unary_ewise_expr_map<acos_fun<T>, Mat>::type
	acos(const IMatrixXpr<Mat, T>& A)
	{
		return ewise(acos_fun<T>(), A.derived());
	}

	template<typename T, class Mat>
	LMAT_ENSURE_INLINE
	inline typename unary_ewise_expr_map<atan_fun<T>, Mat>::type
	atan(const IMatrixXpr<Mat, T>& A)
	{
		return ewise(atan_fun<T>(), A.derived());
	}

	template<typename T, class LMat, class RMat>
	LMAT_ENSURE_INLINE
	inline typename binary_ewise_expr_map<atan2_fun<T>, LMat, RMat>::type
	atan2(const IMatrixXpr<LMat, T>& A, const IMatrixXpr<RMat, T>& B)
	{
		return ewise(atan2_fun<T>(), A.derived(), B.derived());
	}

	template<typename T, class Mat>
	LMAT_ENSURE_INLINE
	inline typename unary_ewise_expr_map<sinh_fun<T>, Mat>::type
	sinh(const IMatrixXpr<Mat, T>& A)
	{
		return ewise(sinh_fun<T>(), A.derived());
	}

	template<typename T, class Mat>
	LMAT_ENSURE_INLINE
	inline typename unary_ewise_expr_map<cosh_fun<T>, Mat>::type
	cosh(const IMatrixXpr<Mat, T>& A)
	{
		return ewise(cosh_fun<T>(), A.derived());
	}

	template<typename T, class Mat>
	LMAT_ENSURE_INLINE
	inline typename unary_ewise_expr_map<tanh_fun<T>, Mat>::type
	tanh(const IMatrixXpr<Mat, T>& A)
	{
		return ewise(tanh_fun<T>(), A.derived());
	}

#ifdef LMAT_HAS_CXX11_MATH

	template<typename T, class Mat>
	LMAT_ENSURE_INLINE
	inline typename unary_ewise_expr_map<round_fun<T>, Mat>::type
	round(const IMatrixXpr<Mat, T>& A)
	{
		return ewise(round_fun<T>(), A.derived());
	}

	template<typename T, class Mat>
	LMAT_ENSURE_INLINE
	inline typename unary_ewise_expr_map<trunc_fun<T>, Mat>::type
	trunc(const IMatrixXpr<Mat, T>& A)
	{
		return ewise(trunc_fun<T>(), A.derived());
	}

	template<typename T, class Mat>
	LMAT_ENSURE_INLINE
	inline typename unary_ewise_expr_map<cbrt_fun<T>, Mat>::type
	cbrt(const IMatrixXpr<Mat, T>& A)
	{
		return ewise(cbrt_fun<T>(), A.derived());
	}

	template<typename T, class LMat, class RMat>
	LMAT_ENSURE_INLINE
	inline typename binary_ewise_expr_map<hypot_fun<T>, LMat, RMat>::type
	hypot(const IMatrixXpr<LMat, T>& A, const IMatrixXpr<RMat, T>& B)
	{
		return ewise(hypot_fun<T>(), A.derived(), B.derived());
	}


	template<typename T, class Mat>
	LMAT_ENSURE_INLINE
	inline typename unary_ewise_expr_map<exp2_fun<T>, Mat>::type
	exp2(const IMatrixXpr<Mat, T>& A)
	{
		return ewise(exp2_fun<T>(), A.derived());
	}

	template<typename T, class Mat>
	LMAT_ENSURE_INLINE
	inline typename unary_ewise_expr_map<log2_fun<T>, Mat>::type
	log2(const IMatrixXpr<Mat, T>& A)
	{
		return ewise(log2_fun<T>(), A.derived());
	}

	template<typename T, class Mat>
	LMAT_ENSURE_INLINE
	inline typename unary_ewise_expr_map<expm1_fun<T>, Mat>::type
	expm1(const IMatrixXpr<Mat, T>& A)
	{
		return ewise(expm1_fun<T>(), A.derived());
	}

	template<typename T, class Mat>
	LMAT_ENSURE_INLINE
	inline typename unary_ewise_expr_map<log1p_fun<T>, Mat>::type
	log1p(const IMatrixXpr<Mat, T>& A)
	{
		return ewise(log1p_fun<T>(), A.derived());
	}


	template<typename T, class Mat>
	LMAT_ENSURE_INLINE
	inline typename unary_ewise_expr_map<asinh_fun<T>, Mat>::type
	asinh(const IMatrixXpr<Mat, T>& A)
	{
		return ewise(asinh_fun<T>(), A.derived());
	}

	template<typename T, class Mat>
	LMAT_ENSURE_INLINE
	inline typename unary_ewise_expr_map<acosh_fun<T>, Mat>::type
	acosh(const IMatrixXpr<Mat, T>& A)
	{
		return ewise(acosh_fun<T>(), A.derived());
	}

	template<typename T, class Mat>
	LMAT_ENSURE_INLINE
	inline typename unary_ewise_expr_map<atanh_fun<T>, Mat>::type
	atanh(const IMatrixXpr<Mat, T>& A)
	{
		return ewise(atanh_fun<T>(), A.derived());
	}

#endif

}

#endif /* MATRIX_EMATH_H_ */


