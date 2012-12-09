/**
 * @file platform_config.h
 *
 * The platform-specific configuration for BCSLib
 *
 * @author Dahua Lin
 */

#ifdef _MSC_VER
#pragma once
#endif

#ifndef LIGHTMAT_PLATFORM_CONFIG_H
#define LIGHTMAT_PLATFORM_CONFIG_H

#define LIGHTMAT_MSVC 0x01
#define LIGHTMAT_GCC 0x02
#define LIGHTMAT_CLANG 0x03

#define LIGHTMAT_WIN32 0x11
#define LIGHTMAT_POSIX 0x12

#if (defined(_WIN32) || defined(_WIN64)) && defined(_MSC_VER)
	#if _MSC_VER < 1600
		#error Microsoft Visual C++ of version lower than MSVC 2010 is not supported.
	#endif
	#define LIGHTMAT_COMPILER LIGHTMAT_MSVC

	#define LIGHTMAT_PLATFORM LIGHTMAT_WIN32

	#define LMAT_USE_C11_STDLIB
	#define LMAT_USE_STATIC_ASSERT
	#define LMAT_HAS_NULLPTR
	#define LMAT_HAS_DECLTYPE

#elif (defined(__GNUC__))

	#if (defined(__clang__))
		#if ((__clang_major__ < 2) || (__clang_major__ == 2 && __clang_minor__ < 9))
			#error CLANG of version lower than 2.9.0 is not supported
		#endif
		#define LIGHTMAT_COMPILER LIGHTMAT_CLANG

		#define LMAT_USE_C11_STDLIB
		#define LMAT_USE_STATIC_ASSERT
		#define LMAT_HAS_C99_MATH
		#define LMAT_HAS_DECLTYPE

		#if (__clang_major__ >= 3)
			#define LMAT_HAS_NULLPTR
		#endif

	#else
		#if ((__GNUC__ < 4) || (__GNUC__ == 4 && __GNUC_MINOR__ < 3))
			#error GCC of version lower than 4.3.0 is not supported
		#endif
		#define LIGHTMAT_COMPILER LIGHTMAT_GCC

		#define LMAT_HAS_C99_MATH

		#if (!(defined(__GXX_EXPERIMENTAL_CXX0X__)))
			#error Light-Matrix requires support of C++11 standard.
		#endif
	#endif

	#define LIGHTMAT_PLATFORM LIGHTMAT_POSIX

#else
	#error Light-Matrix can only be used with Microsoft Visual C++, GCC (G++), or clang (clang++).
#endif

#endif

