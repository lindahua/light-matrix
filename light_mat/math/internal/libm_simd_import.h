/**
 * @file libm_simd_import.h
 *
 * Import of AMD LibM SIMD functions
 * 
 * @author Dahua Lin 
 */

#ifdef _MSC_VER
#pragma once
#endif

#ifndef LIGHTMAT_LIBM_SIMD_IMPORT_H_
#define LIGHTMAT_LIBM_SIMD_IMPORT_H_

#include <light_mat/math/simd_ops.h>

#define LMAT_LIBM_SSE_F( name ) amd_vrs4_##name##f
#define LMAT_LIBM_SSE_D( name ) amd_vrd2_##name

#define LMAT_DECLARE_LIBM_SSE_EXTERN1( name ) \
	__m128  LMAT_LIBM_SSE_F(name)( __m128 ); \
	__m128d LMAT_LIBM_SSE_D(name)( __m128d );

#define LMAT_DECLARE_LIBM_SSE_EXTERN2( name ) \
	__m128  LMAT_LIBM_SSE_F(name)( __m128,  __m128  ); \
	__m128d LMAT_LIBM_SSE_D(name)( __m128d, __m128d );


/************************************************
 *
 *  Declaration of external LibM functions
 *
 ************************************************/

extern "C"
{
	// power functions

	LMAT_DECLARE_LIBM_SSE_EXTERN2( pow )
	LMAT_DECLARE_LIBM_SSE_EXTERN1( cbrt )

	// exp & log

	LMAT_DECLARE_LIBM_SSE_EXTERN1( exp )
	LMAT_DECLARE_LIBM_SSE_EXTERN1( log )
	LMAT_DECLARE_LIBM_SSE_EXTERN1( log10 )

	LMAT_DECLARE_LIBM_SSE_EXTERN1( exp2 )
	LMAT_DECLARE_LIBM_SSE_EXTERN1( log2 )
	LMAT_DECLARE_LIBM_SSE_EXTERN1( exp10 )
	LMAT_DECLARE_LIBM_SSE_EXTERN1( expm1 )
	LMAT_DECLARE_LIBM_SSE_EXTERN1( log1p )

	// trigonometry

	LMAT_DECLARE_LIBM_SSE_EXTERN1( sin )
	LMAT_DECLARE_LIBM_SSE_EXTERN1( cos )
	LMAT_DECLARE_LIBM_SSE_EXTERN1( tan )
}


/************************************************
 *
 *  Import as LMAT functions
 *
 ************************************************/

#define LMAT_IMPORT_LIBM_SSE1( Name ) \
	LMAT_ENSURE_INLINE \
	inline sse_f32pk Name( const sse_f32pk& a ) { \
		return LMAT_LIBM_SSE_F(Name)(a); } \
	LMAT_ENSURE_INLINE \
	inline sse_f64pk Name( const sse_f64pk& a ) { \
		return LMAT_LIBM_SSE_D(Name)(a); }

#define LMAT_IMPORT_LIBM_SSE2( Name ) \
	LMAT_ENSURE_INLINE \
	inline sse_f32pk Name( const sse_f32pk& a, const sse_f32pk& b ) { \
		return LMAT_LIBM_SSE_F(Name)(a, b); } \
	LMAT_ENSURE_INLINE \
	inline sse_f64pk Name( const sse_f64pk& a, const sse_f64pk& b ) { \
		return LMAT_LIBM_SSE_D(Name)(a, b); }

namespace lmat { namespace math {

	// power functions

	LMAT_IMPORT_LIBM_SSE2( pow )
	LMAT_IMPORT_LIBM_SSE1( cbrt )

	// exp & log

	LMAT_IMPORT_LIBM_SSE1( exp )
	LMAT_IMPORT_LIBM_SSE1( log )
	LMAT_IMPORT_LIBM_SSE1( log10 )

	LMAT_IMPORT_LIBM_SSE1( exp2 )
	LMAT_IMPORT_LIBM_SSE1( log2 )
	LMAT_IMPORT_LIBM_SSE1( exp10 )
	LMAT_IMPORT_LIBM_SSE1( expm1 )
	LMAT_IMPORT_LIBM_SSE1( log1p )

	// trigonometry

	LMAT_IMPORT_LIBM_SSE1( sin )
	LMAT_IMPORT_LIBM_SSE1( cos )
	LMAT_IMPORT_LIBM_SSE1( tan )

	// xlogy

	LMAT_ENSURE_INLINE
	inline sse_f32pk xlogy(const sse_f32pk& a, const sse_f32pk& b)
	{
		sse_f32pk z = sse_f32pk::zeros();
		return cond(a > z, log(b), z) * a;
	}

	LMAT_ENSURE_INLINE
	inline sse_f64pk xlogy(const sse_f64pk& a, const sse_f64pk& b)
	{
		sse_f64pk z = sse_f64pk::zeros();
		return cond(a > z, log(b), z) * a;
	}

	// xlogx

	LMAT_ENSURE_INLINE
	inline sse_f32pk xlogx(const sse_f32pk& a)
	{
		return xlogy(a, a);
	}

	LMAT_ENSURE_INLINE
	inline sse_f64pk xlogx(const sse_f64pk& a)
	{
		return xlogy(a, a);
	}

} }


/************************************************
 *
 *  Declaration of SIMD support
 *
 ************************************************/

#define _LMAT_DECLARE_LIBM_SIMD_SUPPORT( name ) LMAT_DEFINE_HAS_SSE_SUPPORT( name )

namespace lmat { namespace meta {

	// power functions

	_LMAT_DECLARE_LIBM_SIMD_SUPPORT( pow_ )
	_LMAT_DECLARE_LIBM_SIMD_SUPPORT( cbrt_ )

	// exp & log

	_LMAT_DECLARE_LIBM_SIMD_SUPPORT( exp_ )
	_LMAT_DECLARE_LIBM_SIMD_SUPPORT( log_ )
	_LMAT_DECLARE_LIBM_SIMD_SUPPORT( log10_ )
	_LMAT_DECLARE_LIBM_SIMD_SUPPORT( xlogy_ )
	_LMAT_DECLARE_LIBM_SIMD_SUPPORT( xlogx_ )

	_LMAT_DECLARE_LIBM_SIMD_SUPPORT( exp2_ )
	_LMAT_DECLARE_LIBM_SIMD_SUPPORT( log2_ )
	_LMAT_DECLARE_LIBM_SIMD_SUPPORT( exp10_ )
	_LMAT_DECLARE_LIBM_SIMD_SUPPORT( expm1_ )
	_LMAT_DECLARE_LIBM_SIMD_SUPPORT( log1p_ )

	// trigonometry

	_LMAT_DECLARE_LIBM_SIMD_SUPPORT( sin_ )
	_LMAT_DECLARE_LIBM_SIMD_SUPPORT( cos_ )
	_LMAT_DECLARE_LIBM_SIMD_SUPPORT( tan_ )

} }



#endif 
