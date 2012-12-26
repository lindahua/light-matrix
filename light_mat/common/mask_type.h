/**
 * @file mask_type.h
 *
 * @brief The definition of mask types
 *
 * @author Dahua Lin
 */

#ifndef LIGHTMAT_MASK_TYPE_H_
#define LIGHTMAT_MASK_TYPE_H_

#include <light_mat/common/prim_types.h>

namespace lmat
{
	// mask type

	template<typename T>
	struct mask_t
	{
		bool bvalue;

		LMAT_ENSURE_INLINE
		mask_t() : bvalue(false) { }

		LMAT_ENSURE_INLINE
		mask_t(bool b)
		{
			bvalue = b;
		}

		LMAT_ENSURE_INLINE
		operator bool() const
		{
			return bvalue;
		}

		LMAT_ENSURE_INLINE
		mask_t operator ! () const
		{
			return !bvalue;
		}

		LMAT_ENSURE_INLINE
		mask_t operator ~ () const
		{
			return !bvalue;
		}

		LMAT_ENSURE_INLINE
		mask_t operator & (const mask_t& r) const
		{
			return bvalue & r.bvalue;
		}

		LMAT_ENSURE_INLINE
		mask_t operator & (const bool& r) const
		{
			return bvalue & r;
		}

		LMAT_ENSURE_INLINE
		mask_t& operator &= (const mask_t& r)
		{
			bvalue &= r.bvalue;
			return *this;
		}

		LMAT_ENSURE_INLINE
		mask_t& operator &= (const bool& r)
		{
			bvalue &= r;
			return *this;
		}

		LMAT_ENSURE_INLINE
		mask_t operator | (const mask_t& r) const
		{
			return bvalue | r.bvalue;
		}

		LMAT_ENSURE_INLINE
		mask_t operator | (const bool& r) const
		{
			return bvalue | r;
		}

		LMAT_ENSURE_INLINE
		mask_t& operator |= (const mask_t& r)
		{
			bvalue |= r.bvalue;
			return *this;
		}

		LMAT_ENSURE_INLINE
		mask_t& operator |= (const bool& r)
		{
			bvalue |= r;
			return *this;
		}

		LMAT_ENSURE_INLINE
		mask_t operator == (const mask_t& r) const
		{
			return bvalue == r.bvalue;
		}

		LMAT_ENSURE_INLINE
		mask_t operator == (const bool& r) const
		{
			return bvalue == r;
		}

		LMAT_ENSURE_INLINE
		mask_t operator != (const mask_t& r) const
		{
			return bvalue != r.bvalue;
		}

		LMAT_ENSURE_INLINE
		mask_t operator != (const bool& r) const
		{
			return bvalue != r;
		}
	};


}

#endif /* MASK_TYPE_H_ */
