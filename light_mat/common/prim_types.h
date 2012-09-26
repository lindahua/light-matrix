/**
 * @file prim_types.h
 *
 * @brief Declaration of primitive types
 *
 * @author Dahua Lin
 */

#ifdef _MSC_VER
#pragma once
#endif

#ifndef LIGHTMAT_PRIM_TYPES_H_
#define LIGHTMAT_PRIM_TYPES_H_

#include <light_mat/config/config.h>

#include <cstddef>

#ifdef LMAT_USE_C11_STDLIB
#include <cstdint>
#else
#include <tr1/cstdint>
#endif


// Useful macros

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
	struct nil_t { };

	// primitive types

	using LMAT_TR1::int8_t;
	using LMAT_TR1::int16_t;
	using LMAT_TR1::int32_t;
	using LMAT_TR1::int64_t;

	using LMAT_TR1::uint8_t;
	using LMAT_TR1::uint16_t;
	using LMAT_TR1::uint32_t;
	using LMAT_TR1::uint64_t;

	using std::ptrdiff_t;
	using std::size_t;

#if (LMAT_INDEX_SIZE == 4)
	typedef int32_t index_t;
#elif (LMAT_INDEX_SIZE == 8)
	typedef int64_t index_t;
#else
#error Invalid value for LMAT_INDEX_SIZE
#endif

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
