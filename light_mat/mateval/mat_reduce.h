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


#define LMAT_DEFINE_BASIC_FULL_REDUCTION( Name ) \
	template<typename T, class Mat> \
	LMAT_ENSURE_INLINE \
	inline T Name(const IRegularMatrix<Mat, T>& mat) { \
		typedef default_simd_kind kind; \
		dimension<meta::nelems<Mat>::value> dim = reduc_get_length(mat); \
		T r; \
		if (dim.value() > 0) { \
			r = internal::Name##_(dim, atags::simd<T, kind>(), in_(mat.derived())); } \
		else { \
			r = empty_values<T>::Name(); } \
		return r; }

#define LMAT_DEFINE_BASIC_COLWISE_REDUCTION( Name ) \
	template<typename T, class Mat, class DMat> \
	LMAT_ENSURE_INLINE \
	inline void colwise_##Name(const IRegularMatrix<Mat, T>& mat, IRegularMatrix<DMat, T>& dmat) { \
		typedef default_simd_kind kind; \
		typename meta::shape<Mat>::type shape = reduc_get_shape(mat); \
		LMAT_CHECK_DIMS( dmat.nelems() == shape.ncolumns() ); \
		if (shape.nrows() > 0) { \
			internal::colwise_##Name##_(shape, atags::simd<T, kind>(), dmat.derived(), in_(mat.derived())); } \
		else { \
			fill(dmat.derived(), empty_values<T>::Name()); } }


#define LMAT_DEFINE_FULL_REDUCTION_1( Name, Reduc, ScaFun, EmptyVal ) \
	template<typename T, class Mat> \
	LMAT_ENSURE_INLINE \
	inline T Name(const IRegularMatrix<Mat, T>& mat) { \
		typedef default_simd_kind kind; \
		dimension<meta::nelems<Mat>::value> dim = reduc_get_length(mat); \
		T r; \
		if (dim.value() > 0) { \
			r = internal::Reduc##x_(dim, atags::simd<T, kind>(), ScaFun<T>(), in_(mat.derived())); } \
		else { r = T(0); } \
		return r; }

#define LMAT_DEFINE_FULL_REDUCTION_2( Name, Reduc, ScaFun, EmptyVal ) \
	template<typename T, class Mat1, class Mat2> \
	LMAT_ENSURE_INLINE \
	inline T Name(const IRegularMatrix<Mat1, T>& mat1, const IRegularMatrix<Mat2, T>& mat2) { \
		typedef default_simd_kind kind; \
		dimension<meta::common_nelems<Mat1, Mat2>::value> dim = reduc_get_length(mat1, mat2); \
		T r; \
		if (dim.value() > 0) { \
			r = internal::Reduc##x_(dim, atags::simd<T, kind>(), ScaFun<T>(), \
					in_(mat1.derived()), in_(mat2.derived())); } \
		else { r = T(0); } \
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


	template<typename T, class Mat>
	LMAT_ENSURE_INLINE
	inline typename meta::shape<Mat>::type
	reduc_get_shape(const IRegularMatrix<Mat, T>& mat)
	{
		static_assert(meta::is_percol_continuous<Mat>::value,
				"mat should be a compile-time percol-continuous matrix.");

		return mat.shape();
	}

	template<typename T, class Mat1, class Mat2>
	LMAT_ENSURE_INLINE
	inline typename meta::common_shape<Mat1, Mat2>::type
	reduc_get_shape(const IRegularMatrix<Mat1, T>& mat1, const IRegularMatrix<Mat2, T>& mat2)
	{
		static_assert(meta::is_continuous<Mat1>::value,
				"mat1 should be a compile-time percol-continuous matrix.");

		static_assert(meta::is_continuous<Mat2>::value,
				"mat2 should be a compile-time percol-continuous matrix.");

		return common_shape(mat1, mat2);
	}


	template<typename T>
	struct empty_values
	{
		LMAT_ENSURE_INLINE
		static T sum() { return T(0); }

		LMAT_ENSURE_INLINE
		static T mean() { return std::numeric_limits<T>::quiet_NaN(); }

		LMAT_ENSURE_INLINE
		static T maximum() { return - std::numeric_limits<T>::infinity(); }

		LMAT_ENSURE_INLINE
		static T minimum() { return std::numeric_limits<T>::infinity(); }
	};


	/********************************************
	 *
	 *  basic reduction function
	 *
	 ********************************************/

	LMAT_DEFINE_BASIC_FULL_REDUCTION( sum )
	LMAT_DEFINE_BASIC_FULL_REDUCTION( mean )
	LMAT_DEFINE_BASIC_FULL_REDUCTION( maximum )
	LMAT_DEFINE_BASIC_FULL_REDUCTION( minimum )

	LMAT_DEFINE_BASIC_COLWISE_REDUCTION( sum )
	LMAT_DEFINE_BASIC_COLWISE_REDUCTION( mean )
	LMAT_DEFINE_BASIC_COLWISE_REDUCTION( maximum )
	LMAT_DEFINE_BASIC_COLWISE_REDUCTION( minimum )


	/********************************************
	 *
	 *  derived reduction function
	 *
	 ********************************************/

	LMAT_DEFINE_FULL_REDUCTION_1( asum,  sum,     math::abs_fun, T(0) )
	LMAT_DEFINE_FULL_REDUCTION_1( amean, mean,    math::abs_fun, T(0) )
	LMAT_DEFINE_FULL_REDUCTION_1( amax,  maximum, math::abs_fun, T(0) )
	LMAT_DEFINE_FULL_REDUCTION_1( sqsum, sum,     math::sqr_fun, T(0) )

	LMAT_DEFINE_FULL_REDUCTION_2( diff_asum,  sum,     math::diff_abs_fun, T(0) )
	LMAT_DEFINE_FULL_REDUCTION_2( diff_amean, mean,    math::diff_abs_fun, T(0) )
	LMAT_DEFINE_FULL_REDUCTION_2( diff_amax,  maximum, math::diff_abs_fun, T(0) )
	LMAT_DEFINE_FULL_REDUCTION_2( diff_sqsum, sum,     math::diff_sqr_fun, T(0) )

	LMAT_DEFINE_FULL_REDUCTION_2( dot, sum, math::mul_fun, T(0) )


	/********************************************
	 *
	 *  norms
	 *
	 ********************************************/

	namespace norms
	{
		struct L1_ { };
		struct L2_ { };
		struct Linf_ { };
	}

	template<typename T, class Mat>
	LMAT_ENSURE_INLINE
	inline T norm(const IRegularMatrix<Mat, T>& mat, norms::L1_)
	{
		return asum(mat.derived());
	}

	template<typename T, class Mat>
	LMAT_ENSURE_INLINE
	inline T norm(const IRegularMatrix<Mat, T>& mat, norms::L2_)
	{
		return math::sqrt(sqsum(mat));
	}

	template<typename T, class Mat>
	LMAT_ENSURE_INLINE
	inline T norm(const IRegularMatrix<Mat, T>& mat, norms::Linf_)
	{
		return amax(mat);
	}


	template<typename T, class Mat1, class Mat2>
	LMAT_ENSURE_INLINE
	inline T diff_norm( const IRegularMatrix<Mat1, T>& mat1,
						const IRegularMatrix<Mat2, T>& mat2, norms::L1_)
	{
		return diff_asum(mat1, mat2);
	}

	template<typename T, class Mat1, class Mat2>
	LMAT_ENSURE_INLINE
	inline T diff_norm( const IRegularMatrix<Mat1, T>& mat1,
						const IRegularMatrix<Mat2, T>& mat2, norms::L2_)
	{
		return math::sqrt(diff_sqsum(mat1, mat2));
	}

	template<typename T, class Mat1, class Mat2>
	LMAT_ENSURE_INLINE
	inline T diff_norm( const IRegularMatrix<Mat1, T>& mat1,
						const IRegularMatrix<Mat2, T>& mat2, norms::Linf_)
	{
		return diff_amax(mat1, mat2);
	}
}

#endif 
