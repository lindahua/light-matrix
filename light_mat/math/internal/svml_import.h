/**
 * @file svml_import.h
 *
 * Import of Intel SVML functions
 * 
 * @author Dahua Lin 
 */

#ifdef _MSC_VER
#pragma once
#endif

#ifndef LIGHTMAT_SVML_IMPORT_H_
#define LIGHTMAT_SVML_IMPORT_H_

#include <light_mat/simd/simd_packs.h>

#ifdef LMAT_USE_INTEL_SVML

/************************************************
 *
 *  Declaration of external SVML functions
 *
 ************************************************/

#define LMAT_SVML_SSE_F( name ) __svml_##name##f4
#define LMAT_SVML_SSE_D( name ) __svml_##name##2

#define LMAT_SVML_AVX_F( name ) __svml_##name##f8
#define LMAT_SVML_AVX_D( name ) __svml_##name##4

#define LMAT_DECLARE_SVML_EXTERN1( name ) \
	__m128  LMAT_SVML_SSE_F(name)( __m128 ); \
	__m128d LMAT_SVML_SSE_D(name)( __m128d ); \
	__m256  LMAT_SVML_AVX_F(name)( __m256 ); \
	__m256d LMAT_SVML_AVX_D(name)( __m256d );

#define LMAT_DECLARE_SVML_EXTERN2( name ) \
	__m128  LMAT_SVML_SSE_F(name)( __m128,  __m128  ); \
	__m128d LMAT_SVML_SSE_D(name)( __m128d, __m128d ); \
	__m256  LMAT_SVML_AVX_F(name)( __m256,  __m256  ); \
	__m256d LMAT_SVML_AVX_D(name)( __m256d, __m256d );


extern "C"
{
	// power functions

	LMAT_DECLARE_SVML_EXTERN2( pow )

	LMAT_DECLARE_SVML_EXTERN2( hypot )
	LMAT_DECLARE_SVML_EXTERN1( cbrt )

	// exp & log

	LMAT_DECLARE_SVML_EXTERN1( exp )
	LMAT_DECLARE_SVML_EXTERN1( log )
	LMAT_DECLARE_SVML_EXTERN1( log10 )

	LMAT_DECLARE_SVML_EXTERN1( exp2 )
	LMAT_DECLARE_SVML_EXTERN1( log2 )
	LMAT_DECLARE_SVML_EXTERN1( exp10 )
	LMAT_DECLARE_SVML_EXTERN1( expm1 )
	LMAT_DECLARE_SVML_EXTERN1( log1p )

	// trigonometry

	LMAT_DECLARE_SVML_EXTERN1( sin )
	LMAT_DECLARE_SVML_EXTERN1( cos )
	LMAT_DECLARE_SVML_EXTERN1( tan )

	LMAT_DECLARE_SVML_EXTERN1( asin )
	LMAT_DECLARE_SVML_EXTERN1( acos )
	LMAT_DECLARE_SVML_EXTERN1( atan )
	LMAT_DECLARE_SVML_EXTERN2( atan2 )

	// hyperbolic

	LMAT_DECLARE_SVML_EXTERN1( sinh )
	LMAT_DECLARE_SVML_EXTERN1( cosh )
	LMAT_DECLARE_SVML_EXTERN1( tanh )

	LMAT_DECLARE_SVML_EXTERN1( asinh )
	LMAT_DECLARE_SVML_EXTERN1( acosh )
	LMAT_DECLARE_SVML_EXTERN1( atanh )

	// special functions

	LMAT_DECLARE_SVML_EXTERN1( erf )
	LMAT_DECLARE_SVML_EXTERN1( erfc )

	LMAT_DECLARE_SVML_EXTERN1( erfinv )
	LMAT_DECLARE_SVML_EXTERN1( erfcinv )

	LMAT_DECLARE_SVML_EXTERN1( cdfnorm )
	LMAT_DECLARE_SVML_EXTERN1( cdfnorminv )
}


/************************************************
 *
 *  Import of LMAT functions
 *
 ************************************************/

#define LMAT_IMPORT_SVML1( Name ) \
	LMAT_ENSURE_INLINE \
	inline sse_f32pk Name( const sse_f32pk& a ) { \
		return LMAT_SVML_SSE_F(Name)(a); } \
	LMAT_ENSURE_INLINE \
	inline sse_f64pk Name( const sse_f64pk& a ) { \
		return LMAT_SVML_SSE_D(Name)(a); } \
	LMAT_ENSURE_INLINE \
	inline avx_f32pk Name( const avx_f32pk& a ) { \
		return LMAT_SVML_AVX_F(Name)(a); } \
	LMAT_ENSURE_INLINE \
	inline avx_f64pk Name( const avx_f64pk& a ) { \
		return LMAT_SVML_AVX_D(Name)(a); }

#define LMAT_IMPORT_SVML2( Name ) \
	LMAT_ENSURE_INLINE \
	inline sse_f32pk Name( const sse_f32pk& a, const sse_f32pk& b ) { \
		return LMAT_SVML_SSE_F(Name)(a, b); } \
	LMAT_ENSURE_INLINE \
	inline sse_f64pk Name( const sse_f64pk& a, const sse_f64pk& b ) { \
		return LMAT_SVML_SSE_D(Name)(a, b); } \
	LMAT_ENSURE_INLINE \
	inline avx_f32pk Name( const avx_f32pk& a, const avx_f32pk& b ) { \
		return LMAT_SVML_AVX_F(Name)(a, b); } \
	LMAT_ENSURE_INLINE \
	inline avx_f64pk Name( const avx_f64pk& a, const avx_f64pk& b ) { \
		return LMAT_SVML_AVX_D(Name)(a, b); }


namespace lmat { namespace math {

	// power functions

	LMAT_IMPORT_SVML2( pow )
	LMAT_IMPORT_SVML1( cbrt )
	LMAT_IMPORT_SVML2( hypot )

	// exp & log

	LMAT_IMPORT_SVML1( exp )
	LMAT_IMPORT_SVML1( log )
	LMAT_IMPORT_SVML1( log10 )

	LMAT_IMPORT_SVML1( exp2 )
	LMAT_IMPORT_SVML1( log2 )
	LMAT_IMPORT_SVML1( exp10 )
	LMAT_IMPORT_SVML1( expm1 )
	LMAT_IMPORT_SVML1( log1p )

	// trigonometry

	LMAT_IMPORT_SVML1( sin )
	LMAT_IMPORT_SVML1( cos )
	LMAT_IMPORT_SVML1( tan )

	LMAT_IMPORT_SVML1( asin )
	LMAT_IMPORT_SVML1( acos )
	LMAT_IMPORT_SVML1( atan )
	LMAT_IMPORT_SVML2( atan2 )

	// hyperbolic

	LMAT_IMPORT_SVML1( sinh )
	LMAT_IMPORT_SVML1( cosh )
	LMAT_IMPORT_SVML1( tanh )

	LMAT_IMPORT_SVML1( asinh )
	LMAT_IMPORT_SVML1( acosh )
	LMAT_IMPORT_SVML1( atanh )

	// special functions

	LMAT_IMPORT_SVML1( erf )
	LMAT_IMPORT_SVML1( erfc )

	LMAT_IMPORT_SVML1( erfinv )
	LMAT_IMPORT_SVML1( erfcinv )

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

	LMAT_ENSURE_INLINE
	inline avx_f32pk xlogy(const avx_f32pk& a, const avx_f32pk& b)
	{
		avx_f32pk z = avx_f32pk::zeros();
		return cond(a > z, log(b), z) * a;
	}

	LMAT_ENSURE_INLINE
	inline avx_f64pk xlogy(const avx_f64pk& a, const avx_f64pk& b)
	{
		avx_f64pk z = avx_f64pk::zeros();
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

	LMAT_ENSURE_INLINE
	inline avx_f32pk xlogx(const avx_f32pk& a)
	{
		return xlogy(a, a);
	}

	LMAT_ENSURE_INLINE
	inline avx_f64pk xlogx(const avx_f64pk& a)
	{
		return xlogy(a, a);
	}

	// norminv

	LMAT_ENSURE_INLINE
	inline sse_f32pk norminv(const sse_f32pk& a)
	{
		return LMAT_SVML_SSE_F(cdfnorminv)(a);
	}

	LMAT_ENSURE_INLINE
	inline sse_f64pk norminv(const sse_f64pk& a)
	{
		return LMAT_SVML_SSE_D(cdfnorminv)(a);
	}

	LMAT_ENSURE_INLINE
	inline avx_f32pk norminv(const avx_f32pk& a)
	{
		return LMAT_SVML_AVX_F(cdfnorminv)(a);
	}

	LMAT_ENSURE_INLINE
	inline avx_f64pk norminv(const avx_f64pk& a)
	{
		return LMAT_SVML_AVX_D(cdfnorminv)(a);
	}

} }


/************************************************
 *
 *  Declaration of SIMD support
 *
 ************************************************/

#define _LMAT_DECLARE_SVML_SIMD_SUPPORT( name ) \
	LMAT_DEFINE_HAS_SSE_SUPPORT( name ) \
	LMAT_DEFINE_HAS_AVX_SUPPORT( name )

namespace lmat { namespace meta {

	// power functions

	_LMAT_DECLARE_SVML_SIMD_SUPPORT( pow_ )
	_LMAT_DECLARE_SVML_SIMD_SUPPORT( cbrt_ )
	_LMAT_DECLARE_SVML_SIMD_SUPPORT( hypot_ )

	// exp & log

	_LMAT_DECLARE_SVML_SIMD_SUPPORT( exp_ )
	_LMAT_DECLARE_SVML_SIMD_SUPPORT( log_ )
	_LMAT_DECLARE_SVML_SIMD_SUPPORT( log10_ )
	_LMAT_DECLARE_SVML_SIMD_SUPPORT( xlogy_ )
	_LMAT_DECLARE_SVML_SIMD_SUPPORT( xlogx_ )

	_LMAT_DECLARE_SVML_SIMD_SUPPORT( exp2_ )
	_LMAT_DECLARE_SVML_SIMD_SUPPORT( log2_ )
	_LMAT_DECLARE_SVML_SIMD_SUPPORT( exp10_ )
	_LMAT_DECLARE_SVML_SIMD_SUPPORT( expm1_ )
	_LMAT_DECLARE_SVML_SIMD_SUPPORT( log1p_ )

	// trigonometry

	_LMAT_DECLARE_SVML_SIMD_SUPPORT( sin_ )
	_LMAT_DECLARE_SVML_SIMD_SUPPORT( cos_ )
	_LMAT_DECLARE_SVML_SIMD_SUPPORT( tan_ )

	_LMAT_DECLARE_SVML_SIMD_SUPPORT( asin_ )
	_LMAT_DECLARE_SVML_SIMD_SUPPORT( acos_ )
	_LMAT_DECLARE_SVML_SIMD_SUPPORT( atan_ )
	_LMAT_DECLARE_SVML_SIMD_SUPPORT( atan2_ )

	// hyperbolic

	_LMAT_DECLARE_SVML_SIMD_SUPPORT( sinh_ )
	_LMAT_DECLARE_SVML_SIMD_SUPPORT( cosh_ )
	_LMAT_DECLARE_SVML_SIMD_SUPPORT( tanh_ )

	_LMAT_DECLARE_SVML_SIMD_SUPPORT( asinh_ )
	_LMAT_DECLARE_SVML_SIMD_SUPPORT( acosh_ )
	_LMAT_DECLARE_SVML_SIMD_SUPPORT( atanh_ )

	// special functions

	_LMAT_DECLARE_SVML_SIMD_SUPPORT( erf_ )
	_LMAT_DECLARE_SVML_SIMD_SUPPORT( erfc_ )
	_LMAT_DECLARE_SVML_SIMD_SUPPORT( erfinv_ )
	_LMAT_DECLARE_SVML_SIMD_SUPPORT( erfcinv_ )

	_LMAT_DECLARE_SVML_SIMD_SUPPORT( norminv_ )

} }

#endif // LMAT_USE_INTEL_SVML


#endif 
