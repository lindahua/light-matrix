/**
 * @file simd_packs.h
 *
 * @brief Overall headers to include all SIMD pack classes
 *
 * @author Dahua Lin
 */

#ifdef _MSC_VER
#pragma once
#endif

#ifndef LIGHTMAT_SIMD_PACKS_H_
#define LIGHTMAT_SIMD_PACKS_H_

#include <light_mat/simd/simd_base.h>

#include <light_mat/simd/sse_packs.h>
#include <light_mat/simd/sse_bpacks.h>
#include <light_mat/simd/sse_reduce.h>

#ifdef LMAT_HAS_AVX
#include <light_mat/simd/avx_packs.h>
#include <light_mat/simd/avx_bpacks.h>
#include <light_mat/simd/avx_reduce.h>
#endif

#endif /* SIMD_PACKS_H_ */
