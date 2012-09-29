/*
 * @file mask_type.h
 *
 * The mask types
 *
 * @author Dahua Lin
 */

#ifndef _MSC_VER
#pragma once
#endif

#ifndef LIGHTMAT_MASK_TYPE_H_
#define LIGHTMAT_MASK_TYPE_H_

#include <light_mat/common/prim_types.h>
#include <light_mat/common/type_traits.h>

namespace lmat
{

	template<typename T>
	class mask_t
	{
#ifdef LMAT_USE_STATIC_ASSERT
		static_assert(is_integral<T>::value, "T must be an integral type.");
#endif

	private:
		T value;

		LMAT_ENSURE_INLINE
		mask_t(T x) : value(x) { };

	public:
		LMAT_ENSURE_INLINE
		mask_t() { };

		LMAT_ENSURE_INLINE
		mask_t(bool b)
		: value(b ? ~T(0) : T(0)) { };

		LMAT_ENSURE_INLINE
		operator bool() const
		{
			return static_cast<bool>(value);
		}

		LMAT_ENSURE_INLINE
		mask_t operator ~ () const
		{
			return ~value;
		}

		LMAT_ENSURE_INLINE
		mask_t operator & (const mask_t& r) const
		{
			return value & r.value;
		}

		LMAT_ENSURE_INLINE
		mask_t operator | (const mask_t& r) const
		{
			return value | r.value;
		}

		LMAT_ENSURE_INLINE
		mask_t operator ^ (const mask_t& r) const
		{
			return value ^ r.value;
		}
	};


	template<>
	class mask_t<float>
	{
	private:
		union
		{
			float value;
			uint32_t vint;
		};

	private:
		LMAT_ENSURE_INLINE
		mask_t(uint32_t vi) : vint(vi) { }

	public:
		LMAT_ENSURE_INLINE
		mask_t() { };

		LMAT_ENSURE_INLINE
		mask_t(bool b)
		: vint(b ? ~uint32_t(0) : uint32_t(0)) { };

		LMAT_ENSURE_INLINE
		operator bool() const
		{
			return static_cast<bool>(vint);
		}

		LMAT_ENSURE_INLINE
		mask_t operator ~ () const
		{
			return ~vint;
		}

		LMAT_ENSURE_INLINE
		mask_t operator & (const mask_t& r) const
		{
			return vint & r.vint;
		}

		LMAT_ENSURE_INLINE
		mask_t operator | (const mask_t& r) const
		{
			return vint | r.vint;
		}

		LMAT_ENSURE_INLINE
		mask_t operator ^ (const mask_t& r) const
		{
			return vint ^ r.vint;
		}
	};


	template<>
	class mask_t<double>
	{
	private:
		union
		{
			double value;
			uint64_t vint;
		};

	private:
		LMAT_ENSURE_INLINE
		mask_t(uint64_t vi) : vint(vi) { }

	public:
		LMAT_ENSURE_INLINE
		mask_t() { };

		LMAT_ENSURE_INLINE
		mask_t(bool b)
		: vint(b ? ~uint64_t(0) : uint64_t(0)) { };

		LMAT_ENSURE_INLINE
		operator bool() const
		{
			return static_cast<bool>(vint);
		}

		LMAT_ENSURE_INLINE
		mask_t operator ~ () const
		{
			return ~vint;
		}

		LMAT_ENSURE_INLINE
		mask_t operator & (const mask_t& r) const
		{
			return vint & r.vint;
		}

		LMAT_ENSURE_INLINE
		mask_t operator | (const mask_t& r) const
		{
			return vint | r.vint;
		}

		LMAT_ENSURE_INLINE
		mask_t operator ^ (const mask_t& r) const
		{
			return vint ^ r.vint;
		}
	};


	// Typedefs

	typedef mask_t<double> mask_f64_t;
	typedef mask_t<float>  mask_f32_t;
	typedef mask_t<uint8_t> mask_u8_t;
	typedef mask_t<int8_t>  mask_i8_t;
	typedef mask_t<uint16_t> mask_u16_t;
	typedef mask_t<int16_t>  mask_i16_t;
	typedef mask_t<uint32_t> mask_u32_t;
	typedef mask_t<int32_t>  mask_i32_t;
	typedef mask_t<uint64_t> mask_u64_t;
	typedef mask_t<int64_t>  mask_i64_t;

}

#endif 

