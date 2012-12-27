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


/********************************************
 *
 *  macros to define reduction functions
 *
 ********************************************/

// basic reduction

#define LMAT_DEFINE_BASIC_FULL_REDUCTION( Name ) \
	template<typename T, class A> \
	LMAT_ENSURE_INLINE \
	inline T Name(const IEWiseMatrix<A, T>& a) { \
		typedef typename internal::full_reduc_policy<Name##_folder<T>, A>::atag atag; \
		dimension<meta::nelems<A>::value> dim = internal::reduc_get_length(a); \
		return dim.value() > 0 ? \
				fold(Name##_folder<T>(), atag())(dim, in_(a.derived())) : \
				internal::empty_values<T>::Name(); }


#define LMAT_DEFINE_BASIC_COLWISE_REDUCTION( Name ) \
	template<typename T, class Mat, class DMat> \
	LMAT_ENSURE_INLINE \
	inline void colwise_##Name(const IEWiseMatrix<Mat, T>& mat, IRegularMatrix<DMat, T>& dmat) { \
		typedef default_simd_kind kind; \
		typename meta::shape<Mat>::type shape = internal::reduc_get_shape(mat); \
		LMAT_CHECK_DIMS( dmat.nelems() == shape.ncolumns() ); \
		if (shape.nrows() > 0) { \
			internal::colwise_##Name##_(shape, atags::simd<kind>(), dmat.derived(), in_(mat.derived())); } \
		else { \
			fill(dmat.derived(), internal::empty_values<T>::Name()); } }

#define LMAT_DEFINE_BASIC_ROWWISE_REDUCTION( Name ) \
	template<typename T, class Mat, class DMat> \
	void rowwise_##Name(const IEWiseMatrix<Mat, T>& mat, IRegularMatrix<DMat, T>& dmat) { \
		typedef default_simd_kind kind; \
		typename meta::shape<Mat>::type shape = internal::reduc_get_shape(mat); \
		LMAT_CHECK_DIMS( dmat.nelems() == shape.nrows() ); \
		if (shape.ncolumns() > 0) { \
			internal::rowwise_##Name##_(shape, atags::simd<kind>(), dmat.derived(), in_(mat.derived())); } \
		else { fill(dmat, internal::empty_values<T>::Name()); } }


// derived reduction

#define LMAT_DEFINE_FULL_REDUCTION_1( Name, Reduc, TExpr, EmptyVal ) \
	template<typename T, class A> \
	LMAT_ENSURE_INLINE \
	inline T Name(const IEWiseMatrix<A, T>& a) { \
		dimension<meta::nelems<A>::value> dim = internal::reduc_get_length(a); \
		return dim.value() > 0 ? Reduc(TExpr) : EmptyVal; \
	}

#define LMAT_DEFINE_FULL_REDUCTION_2( Name, Reduc, TExpr, EmptyVal ) \
	template<typename T, class A, class B> \
	LMAT_ENSURE_INLINE \
	inline T Name(const IEWiseMatrix<A, T>& a, const IEWiseMatrix<B, T>& b) { \
		dimension<meta::common_nelems<A, B>::value> dim = internal::reduc_get_length(a, b); \
		return dim.value() > 0 ? Reduc(TExpr) : EmptyVal; }


#define LMAT_DEFINE_COLWISE_REDUCTION_1( Name, Reduc, TExpr, EmptyVal ) \
	template<typename T, class A, class DMat> \
	LMAT_ENSURE_INLINE \
	inline void colwise_##Name(const IEWiseMatrix<A, T>& a, IRegularMatrix<DMat, T>& dmat) { \
		typename meta::shape<A>::type shape = internal::reduc_get_shape(a); \
		LMAT_CHECK_DIMS( dmat.nelems() == shape.ncolumns() ); \
		if (shape.nrows() > 0) { \
			colwise_##Reduc(TExpr, dmat); \
		} \
		else { fill(dmat.derived(), EmptyVal); } }

#define LMAT_DEFINE_COLWISE_REDUCTION_2( Name, Reduc, TExpr, EmptyVal ) \
	template<typename T, class A, class B, class DMat> \
	LMAT_ENSURE_INLINE \
	inline void colwise_##Name(const IEWiseMatrix<A, T>& a, const IEWiseMatrix<B, T>& b, \
			IRegularMatrix<DMat, T>& dmat) { \
		typename meta::common_shape<A, B>::type shape = internal::reduc_get_shape(a, b); \
		LMAT_CHECK_DIMS( dmat.nelems() == shape.ncolumns() ); \
		if (shape.nrows() > 0) { \
			colwise_##Reduc(TExpr, dmat); \
		} \
		else { fill(dmat.derived(), EmptyVal); } }


#define LMAT_DEFINE_ROWWISE_REDUCTION_1( Name, Reduc, TExpr, EmptyVal ) \
	template<typename T, class A, class DMat> \
	LMAT_ENSURE_INLINE \
	inline void rowwise_##Name(const IEWiseMatrix<A, T>& a, IRegularMatrix<DMat, T>& dmat) { \
		typename meta::shape<A>::type shape = internal::reduc_get_shape(a); \
		LMAT_CHECK_DIMS( dmat.nelems() == shape.nrows() ); \
		if (shape.ncolumns() > 0) { \
			rowwise_##Reduc(TExpr, dmat); \
		} \
		else { fill(dmat.derived(), EmptyVal); } }

#define LMAT_DEFINE_ROWWISE_REDUCTION_2( Name, Reduc, TExpr, EmptyVal ) \
	template<typename T, class A, class B, class DMat> \
	LMAT_ENSURE_INLINE \
	inline void rowwise_##Name(const IEWiseMatrix<A, T>& a, const IEWiseMatrix<B, T>& b, \
			IRegularMatrix<DMat, T>& dmat) { \
		typename meta::common_shape<A, B>::type shape = internal::reduc_get_shape(a, b); \
		LMAT_CHECK_DIMS( dmat.nelems() == shape.nrows() ); \
		if (shape.ncolumns() > 0) { \
			rowwise_##Reduc(TExpr, dmat); \
		} \
		else { fill(dmat.derived(), EmptyVal); } }



namespace lmat
{
	/********************************************
	 *
	 *  basic reduction function
	 *
	 ********************************************/

	// full reduction

	LMAT_DEFINE_BASIC_FULL_REDUCTION( sum )
	LMAT_DEFINE_BASIC_FULL_REDUCTION( maximum )
	LMAT_DEFINE_BASIC_FULL_REDUCTION( minimum )

	template<typename T, class A>
	LMAT_ENSURE_INLINE
	inline T mean(const IEWiseMatrix<A, T>& a)
	{
		dimension<meta::nelems<A>::value> dim = internal::reduc_get_length(a);
		return dim.value() > 0 ?
				sum(a) / T(dim.value()) :
				internal::empty_values<T>::mean();
	}


	// colwise reduction

	LMAT_DEFINE_BASIC_COLWISE_REDUCTION( sum )
	LMAT_DEFINE_BASIC_COLWISE_REDUCTION( mean )
	LMAT_DEFINE_BASIC_COLWISE_REDUCTION( maximum )
	LMAT_DEFINE_BASIC_COLWISE_REDUCTION( minimum )

	// rowwise reduction

	LMAT_DEFINE_BASIC_ROWWISE_REDUCTION( sum )
	LMAT_DEFINE_BASIC_ROWWISE_REDUCTION( mean )
	LMAT_DEFINE_BASIC_ROWWISE_REDUCTION( maximum )
	LMAT_DEFINE_BASIC_ROWWISE_REDUCTION( minimum )


	/********************************************
	 *
	 *  derived reduction function
	 *
	 ********************************************/

	// full reduction

	LMAT_DEFINE_FULL_REDUCTION_1( asum,  sum,     abs(a), T(0) )
	LMAT_DEFINE_FULL_REDUCTION_1( amean, mean,    abs(a), T(0) )
	LMAT_DEFINE_FULL_REDUCTION_1( amax,  maximum, abs(a), T(0) )
	LMAT_DEFINE_FULL_REDUCTION_1( sqsum, sum,     sqr(a), T(0) )

	LMAT_DEFINE_FULL_REDUCTION_2( diff_asum,  sum,     abs(a - b), T(0) )
	LMAT_DEFINE_FULL_REDUCTION_2( diff_amean, mean,    abs(a - b), T(0) )
	LMAT_DEFINE_FULL_REDUCTION_2( diff_amax,  maximum, abs(a - b), T(0) )
	LMAT_DEFINE_FULL_REDUCTION_2( diff_sqsum, sum,     sqr(a - b), T(0) )

	LMAT_DEFINE_FULL_REDUCTION_2( dot, sum, a * b, T(0) )

	// colwise reduction

	LMAT_DEFINE_COLWISE_REDUCTION_1( asum,  sum,     abs(a), T(0) )
	LMAT_DEFINE_COLWISE_REDUCTION_1( amean, mean,    abs(a), T(0) )
	LMAT_DEFINE_COLWISE_REDUCTION_1( amax,  maximum, abs(a), T(0) )
	LMAT_DEFINE_COLWISE_REDUCTION_1( sqsum, sum,     sqr(a), T(0) )

	LMAT_DEFINE_COLWISE_REDUCTION_2( diff_asum,  sum,     abs(a - b), T(0) )
	LMAT_DEFINE_COLWISE_REDUCTION_2( diff_amean, mean,    abs(a - b), T(0) )
	LMAT_DEFINE_COLWISE_REDUCTION_2( diff_amax,  maximum, abs(a - b), T(0) )
	LMAT_DEFINE_COLWISE_REDUCTION_2( diff_sqsum, sum,     sqr(a - b), T(0) )

	LMAT_DEFINE_COLWISE_REDUCTION_2( dot, sum, a * b, T(0) )

	// rowwise reduction

	LMAT_DEFINE_ROWWISE_REDUCTION_1( asum,  sum,     abs(a), T(0) )
	LMAT_DEFINE_ROWWISE_REDUCTION_1( amean, mean,    abs(a), T(0) )
	LMAT_DEFINE_ROWWISE_REDUCTION_1( amax,  maximum, abs(a), T(0) )
	LMAT_DEFINE_ROWWISE_REDUCTION_1( sqsum, sum,     sqr(a), T(0) )

	LMAT_DEFINE_ROWWISE_REDUCTION_2( diff_asum,  sum,     abs(a - b), T(0) )
	LMAT_DEFINE_ROWWISE_REDUCTION_2( diff_amean, mean,    abs(a - b), T(0) )
	LMAT_DEFINE_ROWWISE_REDUCTION_2( diff_amax,  maximum, abs(a - b), T(0) )
	LMAT_DEFINE_ROWWISE_REDUCTION_2( diff_sqsum, sum,     sqr(a - b), T(0) )

	LMAT_DEFINE_ROWWISE_REDUCTION_2( dot, sum, a * b, T(0) )


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


	template<typename T, class Mat1, class Mat2>
	LMAT_ENSURE_INLINE
	inline T diff_norm( const IEWiseMatrix<Mat1, T>& mat1,
						const IEWiseMatrix<Mat2, T>& mat2, norms::L1_)
	{
		return diff_asum(mat1, mat2);
	}

	template<typename T, class Mat1, class Mat2>
	LMAT_ENSURE_INLINE
	inline T diff_norm( const IEWiseMatrix<Mat1, T>& mat1,
						const IEWiseMatrix<Mat2, T>& mat2, norms::L2_)
	{
		return math::sqrt(diff_sqsum(mat1, mat2));
	}

	template<typename T, class Mat1, class Mat2>
	LMAT_ENSURE_INLINE
	inline T diff_norm( const IEWiseMatrix<Mat1, T>& mat1,
						const IEWiseMatrix<Mat2, T>& mat2, norms::Linf_)
	{
		return diff_amax(mat1, mat2);
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


	template<typename T, class Mat1, class Mat2, class DMat>
	LMAT_ENSURE_INLINE
	inline void colwise_diff_norm(
			const IEWiseMatrix<Mat1, T>& mat1, const IEWiseMatrix<Mat2, T>& mat2,
			IRegularMatrix<DMat, T>& dmat, norms::L1_)
	{
		colwise_asum(mat1.derived() - mat2.derived(), dmat);
	}

	template<typename T, class Mat1, class Mat2, class DMat>
	LMAT_ENSURE_INLINE
	inline void colwise_diff_norm(
			const IEWiseMatrix<Mat1, T>& mat1, const IEWiseMatrix<Mat2, T>& mat2,
			IRegularMatrix<DMat, T>& dmat, norms::L2_)
	{
		colwise_diff_sqsum(mat1, mat2, dmat);
		dmat.derived() = sqrt(dmat);
	}

	template<typename T, class Mat1, class Mat2, class DMat>
	LMAT_ENSURE_INLINE
	inline void colwise_diff_norm(
			const IEWiseMatrix<Mat1, T>& mat1, const IEWiseMatrix<Mat2, T>& mat2,
			IRegularMatrix<DMat, T>& dmat, norms::Linf_)
	{
		colwise_diff_amax(mat1, mat2, dmat);
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

	template<typename T, class Mat1, class Mat2, class DMat>
	LMAT_ENSURE_INLINE
	inline void rowwise_diff_norm(
			const IEWiseMatrix<Mat1, T>& mat1, const IEWiseMatrix<Mat2, T>& mat2,
			IRegularMatrix<DMat, T>& dmat, norms::L1_)
	{
		rowwise_diff_asum(mat1, mat2, dmat);
	}

	template<typename T, class Mat1, class Mat2, class DMat>
	LMAT_ENSURE_INLINE
	inline void rowwise_diff_norm(
			const IEWiseMatrix<Mat1, T>& mat1, const IEWiseMatrix<Mat2, T>& mat2,
			IRegularMatrix<DMat, T>& dmat, norms::L2_)
	{
		rowwise_diff_sqsum(mat1, mat2, dmat);
		dmat.derived() = sqrt(dmat);
	}

	template<typename T, class Mat1, class Mat2, class DMat>
	LMAT_ENSURE_INLINE
	inline void rowwise_diff_norm(
			const IEWiseMatrix<Mat1, T>& mat1, const IEWiseMatrix<Mat2, T>& mat2,
			IRegularMatrix<DMat, T>& dmat, norms::Linf_)
	{
		rowwise_diff_amax(mat1, mat2, dmat);
	}

}

#endif 
