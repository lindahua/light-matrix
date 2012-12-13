/**
 * @file simd_base.h
 *
 * @brief The basis for SIMD
 *
 * @author Dahua Lin
 */

#ifndef LIGHTMAT_SIMD_BASE_H_
#define LIGHTMAT_SIMD_BASE_H_

#include <light_mat/math/simd_arch.h>
#include <light_mat/common/basic_defs.h>

// system headers for SIMD intrinsics

#if LMAT_SIMD > 7   	// AVX 2
#ifdef __GNUC__
#include <x86intrin.h>
#else
#include <immintrin.h>
#endif // __GNUC__
#elif LMAT_SIMD == 7
#include <immintrin.h> 	// AVX
#elif LMAT_SIMD == 6
#include <nmmintrin.h> 	// SSE4.2
#elif LMAT_SIMD == 5
#include <smmintrin.h> 	// SSE4.1
#elif LMAT_SIMD == 4
#include <tmmintrin.h>	// SSSE3
#elif LMAT_SIMD == 3
#include <pmmintrin.h> 	// SSE3
#elif LMAT_SIMD == 2
#include <emmintrin.h> 	// SSE2
#elif LMAT_SIMD == 1
#include <xmmintrin.h>	// SSE
#endif

namespace lmat { namespace math {

	// SIMD kind

	struct sse_t { };
	struct avx_t { };

#if (LMAT_SIMD < 7)
	typedef sse_t default_simd_kind;
#else
	typedef avx_t default_simd_kind;
#endif

	// forward declaration of classes


	template<typename T, typename Kind> struct simd_traits;

	template<typename T, typename Kind> struct simd_pack;

} }


// Useful macros

#define LMAT_ALIGN_SSE LMAT_ALIGN(16)
#define LMAT_ALIGN_AVX LMAT_ALIGN(32)

#define LMAT_DEFINE_SIMD_TRAITS( Kind, ScalarT, Wid, Bytes ) \
	template<> struct simd_traits<ScalarT, Kind> { \
		typedef ScalarT scalar_type; \
		static const unsigned int pack_width = Wid; \
		static const unsigned int pack_bytes = Bytes; \
		static const unsigned int pack_nbits = Bytes * 8; \
	};


#define LMAT_DEFINE_FOR_SIMD_PACK( Kind, ScalarT, Wid ) \
		typedef Kind simd_kind; \
		typedef ScalarT scalar_type; \
		static const unsigned int pack_width = Wid;


#endif /* SIMD_BASE_H_ */
