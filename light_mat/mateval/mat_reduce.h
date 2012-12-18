/**
 * @file mat_reduce.h
 *
 * Reduction on matrices
 * 
 * @author Dahua Lin 
 */

#ifdef _MSC_VER
#pragma once
#endif

#ifndef LIGHTMAT_MAT_REDUCE_H_
#define LIGHTMAT_MAT_REDUCE_H_

#include "internal/mat_reduce_internal.h"


#define LMAT_DEFINE_FULL_REDUCTION_1( Name, Reduc, ScaFun, EmptyVal ) \
	template<typename T, class Mat> \
	LMAT_ENSURE_INLINE \
	inline T Name(const IRegularMatrix<Mat, T>& mat) { \
		typedef default_simd_kind kind; \
		unsigned int len = (unsigned int)reduc_get_length(mat); \
		T r; \
		if (len > 0) { \
			internal::Reduc##_impl_x(len, r, kind(), ScaFun<T>(), \
					contvec_reader<T, atags::simd<T, kind> >(mat.ptr_data())); } \
		else { r = EmptyVal; } \
		return r; }

#define LMAT_DEFINE_FULL_REDUCTION_2( Name, Reduc, ScaFun, EmptyVal ) \
	template<typename T, class Mat1, class Mat2> \
	LMAT_ENSURE_INLINE \
	inline T Name(const IRegularMatrix<Mat1, T>& mat1, const IRegularMatrix<Mat2, T>& mat2) { \
		typedef default_simd_kind kind; \
		unsigned int len = (unsigned int)reduc_get_length(mat1, mat2); \
		T r; \
		if (len > 0) { \
			internal::Reduc##_impl_x(len, r, kind(), ScaFun<T>(), \
					contvec_reader<T, atags::simd<T, kind> >(mat1.ptr_data()), \
					contvec_reader<T, atags::simd<T, kind> >(mat2.ptr_data())); } \
		else { r = EmptyVal; } \
		return r; }



namespace lmat
{

	/********************************************
	 *
	 *  auxiliary function
	 *
	 ********************************************/

	template<typename T, class Mat>
	LMAT_ENSURE_INLINE
	inline index_t reduc_get_length(const IRegularMatrix<Mat, T>& mat)
	{
		static_assert(meta::is_continuous<Mat>::value,
				"mat should be a compile-time-continuous matrix.");

		return mat.nelems();
	}

	template<typename T, class Mat1, class Mat2>
	LMAT_ENSURE_INLINE
	inline index_t reduc_get_length(const IRegularMatrix<Mat1, T>& mat1, const IRegularMatrix<Mat2, T>& mat2)
	{
		static_assert(meta::is_continuous<Mat1>::value,
				"mat1 should be a compile-time-continuous matrix.");

		static_assert(meta::is_continuous<Mat2>::value,
				"mat2 should be a compile-time-continuous matrix.");

		return common_shape(mat1.derived(), mat2.derived()).nelems();
	}


	/********************************************
	 *
	 *  basic reduction function
	 *
	 ********************************************/

	template<typename T, class Mat>
	LMAT_ENSURE_INLINE
	inline T sum(const IRegularMatrix<Mat, T>& mat)
	{
		typedef default_simd_kind kind;

		unsigned int len = (unsigned int)reduc_get_length(mat);
		T r;

		if (len > 0)
		{

			r = internal::sum_impl(len, kind(),
					contvec_reader<T, atags::simd<T, kind> >(mat.ptr_data()));
		}
		else
		{
			r = T(0);
		}

		return r;
	}

	template<typename T, class Mat>
	LMAT_ENSURE_INLINE
	inline T maximum(const IRegularMatrix<Mat, T>& mat)
	{
		typedef default_simd_kind kind;

		unsigned int len = (unsigned int)reduc_get_length(mat);
		T r;

		if (len > 0)
		{
			r = internal::maximum_impl(len, kind(),
					contvec_reader<T, atags::simd<T, default_simd_kind> >(mat.ptr_data()));
		}
		else
		{
			r = - std::numeric_limits<T>::infinity();
		}

		return r;
	}

	template<typename T, class Mat>
	LMAT_ENSURE_INLINE
	inline T minimum(const IRegularMatrix<Mat, T>& mat)
	{
		typedef default_simd_kind kind;

		unsigned int len = (unsigned int)reduc_get_length(mat);
		T r;

		if (len > 0)
		{
			r = internal::minimum_impl(len, kind(),
					contvec_reader<T, atags::simd<T, default_simd_kind> >(mat.ptr_data()));
		}
		else
		{
			r = std::numeric_limits<T>::infinity();
		}

		return r;
	}


	/********************************************
	 *
	 *  derived reduction function
	 *
	 ********************************************/

	template<typename T, class Mat>
	LMAT_ENSURE_INLINE
	inline T mean(const IRegularMatrix<Mat, T>& mat)
	{
		unsigned int len = (unsigned int)reduc_get_length(mat);

		if (len > 0)
		{
			return sum(mat) * (T(1) / T(len));
		}
		else return std::numeric_limits<T>::quiet_NaN();
	}

	LMAT_DEFINE_FULL_REDUCTION_1( L1norm, sum, math::abs_fun, T(0) )
	LMAT_DEFINE_FULL_REDUCTION_1( sqL2norm, sum, math::sqr_fun, T(0) )
	LMAT_DEFINE_FULL_REDUCTION_1( Linfnorm, maximum, math::abs_fun, T(0) )

	template<typename T, class Mat>
	LMAT_ENSURE_INLINE
	inline T L2norm(const IRegularMatrix<Mat, T>& mat)
	{
		return math::sqrt(sqL2norm(mat));
	}

	LMAT_DEFINE_FULL_REDUCTION_2( diff_L1norm, sum, math::diff_abs_fun, T(0) )
	LMAT_DEFINE_FULL_REDUCTION_2( diff_sqL2norm, sum, math::diff_sqr_fun, T(0) )
	LMAT_DEFINE_FULL_REDUCTION_2( diff_Linfnorm, maximum, math::diff_abs_fun, T(0) )

	template<typename T, class Mat1, class Mat2>
	LMAT_ENSURE_INLINE
	inline T diff_L2norm(const IRegularMatrix<Mat1, T>& mat1, const IRegularMatrix<Mat2, T>& mat2)
	{
		return std::sqrt(diff_sqL2norm(mat1, mat2));
	}


}

#endif 
