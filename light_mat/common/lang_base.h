/**
 * @file lang_base.h
 *
 * @brief Basic language support facilities.
 *
 * @author Dahua Lin
 */

#ifdef _MSC_VER
#pragma once
#endif

#ifndef LIGHTMAT_LANG_BASE_H_
#define LIGHTMAT_LANG_BASE_H_

#include <light_mat/config/config.h>

// useful macros

#if LIGHTMAT_PLATFORM == LIGHTMAT_POSIX

	#define LMAT_ALIGN(a) __attribute__((aligned(a)))
	#define LMAT_ENSURE_INLINE __attribute__((always_inline))

#elif LIGHTMAT_PLATFORM == LIGHTMAT_WIN32

	#define LMAT_ALIGN(a) __declspec(align(a))
	#define LMAT_ENSURE_INLINE __forceinline

#endif

#ifdef LMAT_HAS_NULLPTR
#define LMAT_NULL nullptr
#else
#define LMAT_NULL NULL
#endif

#define LMAT_CRTP_REF \
		LMAT_ENSURE_INLINE const Derived& derived() const { return *(static_cast<const Derived*>(this)); } \
		LMAT_ENSURE_INLINE Derived& derived() { return *(static_cast<Derived*>(this)); }

namespace lmat
{
	/**
	 * @brief The base class to ensure derived classes to be non-copyable.
	 */
	class noncopyable
	{
	protected:
		noncopyable() { }
		~noncopyable() { }

	private:
		noncopyable(const noncopyable& );
		noncopyable& operator= (const noncopyable& );
	};

}

#endif
