/**
 * @file sse_math_emulate.h
 *
 * Emulation of math functions on SSE packs using scalar
 * 
 * @author Dahua Lin 
 */

#ifdef _MSC_VER
#pragma once
#endif

#ifndef LIGHTMAT_SSE_MATH_EMULATE_H_
#define LIGHTMAT_SSE_MATH_EMULATE_H_

#include <light_mat/math/math_base.h>
#include <light_mat/math/sse_packs.h>

#define LMAT_ENABLE_SIMD_EMULATE

#define LMAT_DEFINE_SSE_MATH_EMULATE_1( Name ) \
	LMAT_ENSURE_INLINE \
	inline sse_f32pk Name##_emulate(const sse_f32pk& a) { \
		LMAT_ALIGN_SSE float r_[4]; \
		for (unsigned int i = 0; i < 4; ++i) { \
			r_[i] = lmat::math::Name(a[i]); \
		} \
		sse_f32pk r; \
		r.load_a(r_); \
		return r; } \
	LMAT_ENSURE_INLINE \
	inline sse_f64pk Name##_emulate(const sse_f64pk& a) { \
		LMAT_ALIGN_SSE double r_[2]; \
		r_[0] = lmat::math::Name(a[0]); \
		r_[1] = lmat::math::Name(a[1]); \
		sse_f64pk r; \
		r.load_a(r_); \
		return r; }

#define LMAT_DEFINE_SSE_MATH_EMULATE_2( Name ) \
	LMAT_ENSURE_INLINE \
	inline sse_f32pk Name##_emulate(const sse_f32pk& a, const sse_f32pk& b) { \
		LMAT_ALIGN_SSE float r_[4]; \
		for (unsigned int i = 0; i < 4; ++i) { \
			r_[i] = lmat::math::Name(a[i], b[i]); \
		} \
		sse_f32pk r; \
		r.load_a(r_); \
		return r; } \
	LMAT_ENSURE_INLINE \
	inline sse_f64pk Name##_emulate(const sse_f64pk& a, const sse_f64pk& b) { \
		LMAT_ALIGN_SSE double r_[2]; \
		r_[0] = lmat::math::Name(a[0], b[0]); \
		r_[1] = lmat::math::Name(a[1], b[1]); \
		sse_f64pk r; \
		r.load_a(r_); \
		return r; }


#define LMAT_ACTIVATE_SSE_MATH_EMULATE_1( Name ) \
	LMAT_ENSURE_INLINE \
	inline sse_f32pk Name( const sse_f32pk& a ) { \
		return internal::Name##_emulate(a); } \
	LMAT_ENSURE_INLINE \
	inline sse_f64pk Name( const sse_f64pk& a ) { \
		return internal::Name##_emulate(a); }


#define LMAT_ACTIVATE_SSE_MATH_EMULATE_2( Name ) \
	LMAT_ENSURE_INLINE \
	inline sse_f32pk Name( const sse_f32pk& a, const sse_f32pk& b ) { \
		return internal::Name##_emulate(a, b); } \
	LMAT_ENSURE_INLINE \
	inline sse_f64pk Name( const sse_f64pk& a, const sse_f64pk& b ) { \
		return internal::Name##_emulate(a, b); }


namespace lmat { namespace math { namespace internal {

	// power, exp, and log

	LMAT_DEFINE_SSE_MATH_EMULATE_2( pow )

	LMAT_DEFINE_SSE_MATH_EMULATE_1( exp )
	LMAT_DEFINE_SSE_MATH_EMULATE_1( log )
	LMAT_DEFINE_SSE_MATH_EMULATE_1( log10 )

	// trigonometry

	LMAT_DEFINE_SSE_MATH_EMULATE_1( sin )
	LMAT_DEFINE_SSE_MATH_EMULATE_1( cos )
	LMAT_DEFINE_SSE_MATH_EMULATE_1( tan )

	LMAT_DEFINE_SSE_MATH_EMULATE_1( asin )
	LMAT_DEFINE_SSE_MATH_EMULATE_1( acos )
	LMAT_DEFINE_SSE_MATH_EMULATE_1( atan )
	LMAT_DEFINE_SSE_MATH_EMULATE_2( atan2 )

	// hyperbolic

	LMAT_DEFINE_SSE_MATH_EMULATE_1( sinh )
	LMAT_DEFINE_SSE_MATH_EMULATE_1( cosh )
	LMAT_DEFINE_SSE_MATH_EMULATE_1( tanh )

#ifdef LMAT_HAS_CXX11_MATH

	// hypot & cbrt

	LMAT_DEFINE_SSE_MATH_EMULATE_2( hypot )
	LMAT_DEFINE_SSE_MATH_EMULATE_1( cbrt )

	// extended exp & log

	LMAT_DEFINE_SSE_MATH_EMULATE_1( exp2 )
	LMAT_DEFINE_SSE_MATH_EMULATE_1( log2 )
	LMAT_DEFINE_SSE_MATH_EMULATE_1( expm1 )
	LMAT_DEFINE_SSE_MATH_EMULATE_1( log1p )

	// inverse hyperbolic

	LMAT_DEFINE_SSE_MATH_EMULATE_1( asinh )
	LMAT_DEFINE_SSE_MATH_EMULATE_1( acosh )
	LMAT_DEFINE_SSE_MATH_EMULATE_1( atanh )

	// error functions

	LMAT_DEFINE_SSE_MATH_EMULATE_1( erf )
	LMAT_DEFINE_SSE_MATH_EMULATE_1( erfc )

	// gamma

	LMAT_DEFINE_SSE_MATH_EMULATE_1( lgamma )
	LMAT_DEFINE_SSE_MATH_EMULATE_1( tgamma )

#endif

} } }

#endif 
