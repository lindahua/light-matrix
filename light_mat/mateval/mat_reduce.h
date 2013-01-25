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
		return a.nelems() > 0 ? \
				fold(Name##_kernel<T>())(a.shape(), in_(a)) : \
				internal::empty_values<T>::Name(); }

#define LMAT_DEFINE_BASIC_COLWISE_REDUCTION( Name ) \
	template<typename T, class A, class DMat> \
	inline void colwise_##Name(const IEWiseMatrix<A, T>& a, IRegularMatrix<DMat, T>& dmat) { \
		auto shape = a.shape(); \
		if (shape.nrows() > 0) { \
			internal::colwise_fold_impl(shape, Name##_kernel<T>(), dmat, a ); } \
		else { fill(dmat, internal::empty_values<T>::Name()); } }

#define LMAT_DEFINE_BASIC_ROWWISE_REDUCTION( Name ) \
	template<typename T, class A, class DMat> \
	inline void rowwise_##Name(const IEWiseMatrix<A, T>& a, IRegularMatrix<DMat, T>& dmat) { \
		auto shape = a.shape(); \
		if (shape.ncolumns() > 0) { \
			internal::rowwise_fold_impl(shape, Name##_kernel<T>(), dmat, a); } \
		else { \
			fill(dmat, internal::empty_values<T>::Name()); } }


// derived reduction

#define LMAT_DEFINE_FULL_REDUCTION_1( Name, Reduc, TExpr, EmptyVal ) \
	template<typename T, class A> \
	LMAT_ENSURE_INLINE \
	inline T Name(const IEWiseMatrix<A, T>& a) { \
		dimension<meta::nelems<A>::value> dim = internal::reduc_get_length(a); \
		return dim.value() > 0 ? Reduc(TExpr) : EmptyVal; }

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

	using internal::colwise_fold_getter;
	using internal::make_colwise_fold_getter;

	LMAT_DEFINE_BASIC_COLWISE_REDUCTION( sum )
	LMAT_DEFINE_BASIC_COLWISE_REDUCTION( maximum )
	LMAT_DEFINE_BASIC_COLWISE_REDUCTION( minimum )

	template<typename T, class A, class DMat>
	inline void colwise_mean(const IEWiseMatrix<A, T>& a, IRegularMatrix<DMat, T>& dmat)
	{
		auto shape = internal::reduc_get_shape(a);
		if (shape.nrows() > 0)
		{
			colwise_sum(a, dmat);
			dmat *= math::rcp((T)shape.nrows());
		}
		else
		{
			fill(dmat, internal::empty_values<T>::mean());
		}
	}

	// rowwise reduction

	LMAT_DEFINE_BASIC_ROWWISE_REDUCTION( sum )
	LMAT_DEFINE_BASIC_ROWWISE_REDUCTION( maximum )
	LMAT_DEFINE_BASIC_ROWWISE_REDUCTION( minimum )

	template<typename T, class A, class DMat>
	inline void rowwise_mean(const IEWiseMatrix<A, T>& a, IRegularMatrix<DMat, T>& dmat)
	{
		auto shape = internal::reduc_get_shape(a);
		if (shape.ncolumns() > 0)
		{
			rowwise_sum(a, dmat);
			dmat *= math::rcp((T)shape.ncolumns());
		}
		else
		{
			fill(dmat, internal::empty_values<T>::mean());
		}
	}


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


}

#endif 
