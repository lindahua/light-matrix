/**
 * @file simd_math.h
 *
 * SIMD math functions
 * 
 * @author Dahua Lin 
 */

#ifdef _MSC_VER
#pragma once
#endif

#ifndef LIGHTMAT_SIMD_MATH_H_
#define LIGHTMAT_SIMD_MATH_H_

#include <light_mat/math/math.h>
#include <light_mat/math/simd_ops.h>

#if defined(LMAT_USE_INTEL_SVML) && defined(LMAT_USE_AMD_LIBM)
#error SVML and LIBM cannot be used simultaneously.
#endif

#ifdef LMAT_USE_INTEL_SVML
#include "internal/svml_import.h"
#elif LMAT_USE_AMD_LIBM
#endif

#endif 
