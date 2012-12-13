/****************************  instrset.h   **********************************
| Author:        Agner Fog
| Date created:  2012-05-30
| Last modified: 2012-08-01
| Version:       1.02 Beta
| Project:       vector classes
| Description:
| Header file for various common tasks to vector class library:
| > selects the supported instruction set
| > defines integer types
| > defines compiler version macro if GCC compiler used
| > defines template class to represent compile-time integer constant
| > defines template for compile-time error messages
|
| (c) Copyright 2012 GNU General Public License http://www.gnu.org/licenses
\*****************************************************************************/

#ifndef LIGHTMAT_VC_INSTRSET_H
#define LIGHTMAT_VC_INSTRSET_H

#include <stdlib.h> // define abs(int)

// Detect 64 bit mode
#if (defined(_M_AMD64) || defined(_M_X64) || defined(__amd64) ) && ! defined(__x86_64__)
#define __x86_64__ 1  // There are many different macros for this, decide on only one
#endif

// Find instruction set from compiler macros if INSTRSET not defined
// Note: Not all compilers define these macros automatically
#ifndef INSTRSET
#if defined ( __AVX2__ )
#define INSTRSET 8
#elif defined ( __AVX__ )
#define INSTRSET 7
#elif defined ( __SSE4_2__ )
#define INSTRSET 6
#elif defined ( __SSE4_1__ )
#define INSTRSET 5
#elif defined ( __SSSE3__ )
#define INSTRSET 4
#elif defined ( __SSE3__ )
#define INSTRSET 3
#elif defined ( __SSE2__ ) || defined ( __x86_64__ )
#define INSTRSET 2
#elif defined ( __SSE__ )
#define INSTRSET 1
#elif defined ( _M_IX86_FP )           // Defined in MS compiler. 1: SSE, 2: SSE2
#define INSTRSET _M_IX86_FP
#else 
#define INSTRSET 0
#endif // instruction set defines
#endif // INSTRSET

// Include the appropriate header file for intrinsic functions
#if INSTRSET > 7                       // AVX2 and later
#ifdef __GNUC__
#include <x86intrin.h>                 // x86intrin.h includes header files for whatever instruction 
                                       // sets are specified on the compiler command line, such as:
                                       // xopintrin.h, fma4intrin.h
#else
#include <immintrin.h>                 // MS version of immintrin.h covers AVX, AVX2 and FMA3
#endif // __GNUC__
#elif INSTRSET == 7
#include <immintrin.h>                 // AVX
#elif INSTRSET == 6
#include <nmmintrin.h>                 // SSE4.2
#elif INSTRSET == 5
#include <smmintrin.h>                 // SSE4.1
#elif INSTRSET == 4
#include <tmmintrin.h>                 // SSSE3
#elif INSTRSET == 3
#include <pmmintrin.h>                 // SSE3
#elif INSTRSET == 2
#include <emmintrin.h>                 // SSE2
#elif INSTRSET == 1
#include <xmmintrin.h>                 // SSE
#endif // INSTRSET

#if INSTRSET < 2   // SSE2 required
#error Please compile for the SSE2 instruction set or higher
#endif

// AMD  instruction sets
#if defined (__XOP__) || defined (__FMA4__)
#ifdef __GNUC__
#include <x86intrin.h>                 // AMD XOP (Gnu)
#else
#include <ammintrin.h>                 // AMD XOP (Microsoft)
#endif //  __GNUC__
#elif defined (__SSE4A__)              // AMD SSE4A
#include <ammintrin.h>
#endif // __XOP__ 

// FMA3 instruction set
#if defined (__FMA__)
#include <fmaintrin.h> 
#endif // __FMA__ 

// FMA4 instruction set
#if defined (__FMA4__)
#include <fma4intrin.h> // must have both x86intrin.h and fma4intrin.h, don't know why
#endif // __FMA4__ 


#ifdef _MSC_VER                                  // Microsoft compiler or compatible Intel compiler
#include <intrin.h>                              // define _BitScanReverse(int), __cpuid(int[4],int), _xgetbv(int)
#endif // _MSC_VER

// GCC version
#if defined(__GNUC__) && ! defined (GCC_VERSION)
#define GCC_VERSION  ((__GNUC__) * 10000 + (__GNUC_MINOR__) * 100 + (__GNUC_PATCHLEVEL__))
#endif

namespace lmat { namespace vc {

	// functions in instrset_detect.cpp
	int  instrset_detect(void);                      // tells which instruction sets are supported
	bool hasFMA3(void);                              // true if FMA3 instructions supported
	bool hasFMA4(void);                              // true if FMA4 instructions supported
	bool hasXOP (void);                              // true if XOP  instructions supported

	// Template class to represent compile-time integer constant
	template <int32_t  n> class Const_int_t  {};     // represent compile-time signed integer constant
	template <uint32_t n> class Const_uint_t {};     // represent compile-time unsigned integer constant

	// Template for compile-time error messages
	template <bool> class Static_error_check {
	public:  Static_error_check(){};
	};
	template <> class Static_error_check<false> {    // generate compile-time error if false
	private: Static_error_check(){};
	};

} }

#define const_int(n)  (Const_int_t <n>())        // n must be compile-time integer constant
#define const_uint(n) (Const_uint_t<n>())        // n must be compile-time unsigned integer constant

#endif 
