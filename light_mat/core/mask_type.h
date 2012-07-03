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

#include <light_mat/core/prim_types.h>
#include <light_mat/core/lang_base.h>
#include <light_mat/core/type_traits.h>

namespace lmat
{

	template<typename T>
	class mask_t
	{
		static_assert(is_integral<T>::value, "T must be an integral type.");

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


	template<typename T>
	struct is_compatible_type<mask_t<T>, bool>
	{
		static const bool value = true;
	};

	template<typename T>
	struct is_compatible_type<bool, mask_t<T> >
	{
		static const bool value = true;
	};
}

#endif 

