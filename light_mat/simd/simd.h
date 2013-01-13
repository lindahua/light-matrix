/**
 * @file simd.h
 *
 * @brief The overall header for SIMD
 *
 * @author Dahua Lin
 */

#ifdef _MSC_VER
#pragma once
#endif

#ifndef LIGHTMAT_SIMD_H_
#define LIGHTMAT_SIMD_H_

#include <light_mat/simd/simd_base.h>

#include <light_mat/simd/sse.h>

#ifdef LMAT_HAS_AVX
#include <light_mat/simd/avx.h>
#endif

#endif /* SIMD_H_ */
