/**
 * @file simd_ops.h
 *
 * @brief Overall header that includes relevant SIMD ops headers
 *
 * @author Dahua Lin
 */

#ifdef _MSC_VER
#pragma once
#endif

#ifndef LIGHTMAT_SIMD_OPS_H_
#define LIGHTMAT_SIMD_OPS_H_

#include <light_mat/math/simd_arch.h>
#include <light_mat/math/sse_ops.h>
#ifdef LMAT_HAS_AVX
#include <light_mat/math/avx_ops.h>
#endif

#endif
