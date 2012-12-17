/**
 * @file simd_packs.h
 *
 * @brief Overall headers to include all SIMD pack classes
 *
 * @author Dahua Lin
 */

#ifndef LIGHTMAT_SIMD_PACKS_H_
#define LIGHTMAT_SIMD_PACKS_H_

#include <light_mat/math/simd_base.h>

#include <light_mat/math/sse_packs.h>
#include <light_mat/math/sse_bpacks.h>
#include <light_mat/math/sse_reduce.h>

#ifdef LMAT_HAS_AVX
#include <light_mat/math/avx_packs.h>
#include <light_mat/math/avx_bpacks.h>
#include <light_mat/math/avx_reduce.h>
#endif

#endif /* SIMD_PACKS_H_ */
