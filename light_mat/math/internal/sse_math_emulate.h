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

#define LMAT_DEFINE_SSE_MATH_EMULATE_1( Name ) \
	LMAT_ENSURE_INLINE \
	inline sse_f32pk Name##_emulate(const sse_f32pk& a) { \
		LMAT_ALIGN_SSE float r_[4]; \
		for (int i = 0; i < 4; ++i) { \
			r_[i] = lmat::math::Name(a.e[i]); \
		} \
		sse_f32pk r; \
		r.load_a(r_); \
		return r; } \
	LMAT_ENSURE_INLINE \
	inline sse_f64pk Name##_emulate(const sse_f64pk& a) { \
		LMAT_ALIGN_SSE double r_[2]; \
		r_[0] = lmat::math::Name(a.e[0]); \
		r_[1] = lmat::math::Name(a.e[1]); \
		sse_f64pk r; \
		r.load_a(r_); \
		return r; }

#define LMAT_DEFINE_SSE_MATH_EMULATE_2( Name ) \
	LMAT_ENSURE_INLINE \
	inline sse_f32pk Name##_emulate(const sse_f32pk& a, const sse_f32pk& b) { \
		LMAT_ALIGN_SSE float r_[4]; \
		for (int i = 0; i < 4; ++i) { \
			r_[i] = lmat::math::Name(a.e[i], b.e[i]); \
		} \
		sse_f32pk r; \
		r.load_a(r_); \
		return r; } \
	LMAT_ENSURE_INLINE \
	inline sse_f64pk Name##_emulate(const sse_f64pk& a, const sse_f64pk& b) { \
		LMAT_ALIGN_SSE double r_[2]; \
		r_[0] = lmat::math::Name(a.e[0], b.e[0]); \
		r_[1] = lmat::math::Name(a.e[1], b.e[1]); \
		sse_f64pk r; \
		r.load_a(r_); \
		return r; }



namespace lmat { namespace math { namespace internal {

	// exp & log

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


} } }

#endif 
