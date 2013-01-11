/**
 * @file avx_math_emulate.h
 *
 * @brief Emulation of math functions on AVX packs using scalar
 *
 * @author Dahua Lin
 */

#ifdef _MSC_VER
#pragma once
#endif

#ifndef LIGHTMAT_AVX_MATH_EMULATE_H_
#define LIGHTMAT_AVX_MATH_EMULATE_H_

#include <light_mat/math/math_base.h>
#include <light_mat/math/avx_packs.h>

#define LMAT_ENABLE_SIMD_EMULATE

#define LMAT_DEFINE_AVX_MATH_EMULATE_1( Name ) \
	LMAT_ENSURE_INLINE \
	inline avx_f32pk Name##_emulate(const avx_f32pk& a) { \
		LMAT_ALIGN_AVX float r_[8]; \
		for (unsigned int i = 0; i < 8; ++i) { \
			r_[i] = lmat::math::Name(a[i]); \
		} \
		avx_f32pk r; \
		r.load_a(r_); \
		return r; } \
	LMAT_ENSURE_INLINE \
	inline avx_f64pk Name##_emulate(const avx_f64pk& a) { \
		LMAT_ALIGN_AVX double r_[4]; \
		for (unsigned int i = 0; i < 4; ++i) { \
			r_[i] = lmat::math::Name(a[i]); \
		} \
		avx_f64pk r; \
		r.load_a(r_); \
		return r; }

#define LMAT_DEFINE_AVX_MATH_EMULATE_2( Name ) \
	LMAT_ENSURE_INLINE \
	inline avx_f32pk Name##_emulate(const avx_f32pk& a, const avx_f32pk& b) { \
		LMAT_ALIGN_AVX float r_[8]; \
		for (unsigned int i = 0; i < 8; ++i) { \
			r_[i] = lmat::math::Name(a[i], b[i]); \
		} \
		avx_f32pk r; \
		r.load_a(r_); \
		return r; } \
	LMAT_ENSURE_INLINE \
	inline avx_f64pk Name##_emulate(const avx_f64pk& a, const avx_f64pk& b) { \
		LMAT_ALIGN_AVX double r_[4]; \
		for (unsigned int i = 0; i < 4; ++i) { \
			r_[i] = lmat::math::Name(a[i], b[i]); \
		} \
		avx_f64pk r; \
		r.load_a(r_); \
		return r; }


#define LMAT_ACTIVATE_AVX_MATH_EMULATE_1( Name ) \
	LMAT_ENSURE_INLINE \
	inline avx_f32pk Name( const avx_f32pk& a ) { \
		return internal::Name##_emulate(a); } \
	LMAT_ENSURE_INLINE \
	inline avx_f64pk Name( const avx_f64pk& a ) { \
		return internal::Name##_emulate(a); }


#define LMAT_ACTIVATE_AVX_MATH_EMULATE_2( Name ) \
	LMAT_ENSURE_INLINE \
	inline avx_f32pk Name( const avx_f32pk& a, const avx_f32pk& b ) { \
		return internal::Name##_emulate(a, b); } \
	LMAT_ENSURE_INLINE \
	inline avx_f64pk Name( const avx_f64pk& a, const avx_f64pk& b ) { \
		return internal::Name##_emulate(a, b); }


namespace lmat { namespace math { namespace internal {

	// power, exp, and log

	LMAT_DEFINE_AVX_MATH_EMULATE_2( pow )

	LMAT_DEFINE_AVX_MATH_EMULATE_1( exp )
	LMAT_DEFINE_AVX_MATH_EMULATE_1( log )
	LMAT_DEFINE_AVX_MATH_EMULATE_1( log10 )

	LMAT_DEFINE_AVX_MATH_EMULATE_2( xlogy )

	// trigonometry

	LMAT_DEFINE_AVX_MATH_EMULATE_1( sin )
	LMAT_DEFINE_AVX_MATH_EMULATE_1( cos )
	LMAT_DEFINE_AVX_MATH_EMULATE_1( tan )

	LMAT_DEFINE_AVX_MATH_EMULATE_1( asin )
	LMAT_DEFINE_AVX_MATH_EMULATE_1( acos )
	LMAT_DEFINE_AVX_MATH_EMULATE_1( atan )
	LMAT_DEFINE_AVX_MATH_EMULATE_2( atan2 )

	// hyperbolic

	LMAT_DEFINE_AVX_MATH_EMULATE_1( sinh )
	LMAT_DEFINE_AVX_MATH_EMULATE_1( cosh )
	LMAT_DEFINE_AVX_MATH_EMULATE_1( tanh )

#ifdef LMAT_HAS_CXX11_MATH

	// hypot & cbrt

	LMAT_DEFINE_AVX_MATH_EMULATE_2( hypot )
	LMAT_DEFINE_AVX_MATH_EMULATE_1( cbrt )

	// extended exp & log

	LMAT_DEFINE_AVX_MATH_EMULATE_1( exp2 )
	LMAT_DEFINE_AVX_MATH_EMULATE_1( log2 )
	LMAT_DEFINE_AVX_MATH_EMULATE_1( expm1 )
	LMAT_DEFINE_AVX_MATH_EMULATE_1( log1p )

	// inverse hyperbolic

	LMAT_DEFINE_AVX_MATH_EMULATE_1( asinh )
	LMAT_DEFINE_AVX_MATH_EMULATE_1( acosh )
	LMAT_DEFINE_AVX_MATH_EMULATE_1( atanh )

#endif

} } }


#endif /* AVX_MATH_EMULATE_H_ */
