/**
 * @file memalloc.h
 *
 * @brief Memory allocation devices
 *
 * @author Dahua Lin
 */

#ifdef _MSC_VER
#pragma once
#endif

#ifndef LIGHTMAT_MEMALLOC_H_
#define LIGHTMAT_MEMALLOC_H_

#include <light_mat/common/basic_defs.h>
#include "internal/align_alloc.h"

#include <limits>

namespace lmat
{

	/********************************************
	 *
	 *  Aligned memory allocation
	 *
	 ********************************************/

    template<typename T, unsigned int Align=LMAT_DEFAULT_ALIGNMENT>
    class aligned_allocator
    {
    public:
    	typedef T value_type;
    	typedef T* pointer;
    	typedef T& reference;
    	typedef const T* const_pointer;
    	typedef const T& const_reference;
    	typedef size_t size_type;
    	typedef ptrdiff_t difference_type;

    	template<typename TOther>
    	struct rebind
    	{
    		typedef aligned_allocator<TOther> other;
    	};

    public:
    	LMAT_ENSURE_INLINE
    	aligned_allocator() { }

    	template<typename U>
    	LMAT_ENSURE_INLINE
    	aligned_allocator(const aligned_allocator<U>& r) { }

    	LMAT_ENSURE_INLINE
    	unsigned int alignment() const
    	{
    		return Align;
    	}

    	LMAT_ENSURE_INLINE
    	pointer address( reference x ) const
    	{
    		return &x;
    	}

    	LMAT_ENSURE_INLINE
    	const_pointer address( const_reference x ) const
    	{
    		return &x;
    	}

    	LMAT_ENSURE_INLINE
    	size_type max_size() const
    	{
    		return std::numeric_limits<size_type>::max() / sizeof(value_type);
    	}

    	LMAT_ENSURE_INLINE
    	pointer allocate(size_type n, const void* hint=0)
    	{
    		return (pointer)internal::aligned_allocate(n * sizeof(value_type), Align);
    	}

    	LMAT_ENSURE_INLINE
    	void deallocate(pointer p, size_type)
    	{
    		internal::aligned_release(p);
    	}

    	LMAT_ENSURE_INLINE
    	void construct (pointer p, const_reference val)
    	{
    		new (p) value_type(val);
    	}

    	LMAT_ENSURE_INLINE
    	void destroy (pointer p)
    	{
    		p->~value_type();
    	}

    }; // end class aligned_allocator

}


#endif /* MEMALLOC_H_ */
