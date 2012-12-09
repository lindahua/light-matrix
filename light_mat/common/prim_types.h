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
#include <cstdint>
#include <type_traits>

// Useful macros

#if LIGHTMAT_PLATFORM == LIGHTMAT_POSIX

	#define LMAT_ALIGN(a) __attribute__((aligned(a)))
	#define LMAT_ENSURE_INLINE __attribute__((always_inline))

#elif LIGHTMAT_PLATFORM == LIGHTMAT_WIN32

	#define LMAT_ALIGN(a) __declspec(align(a))
	#define LMAT_ENSURE_INLINE __forceinline

#endif


#define LMAT_CRTP_REF \
		LMAT_ENSURE_INLINE const Derived& derived() const { return *(static_cast<const Derived*>(this)); } \
		LMAT_ENSURE_INLINE Derived& derived() { return *(static_cast<Derived*>(this)); }

namespace lmat
{
	struct nil_t { };

	template<int I>
	struct fix_int : std::integral_constant<int, I> { };

	using std::true_type;
	using std::false_type;

	template<typename T> struct type{ };

	// primitive types

	using std::int8_t;
	using std::int16_t;
	using std::int32_t;
	using std::int64_t;

	using std::uint8_t;
	using std::uint16_t;
	using std::uint32_t;
	using std::uint64_t;

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
