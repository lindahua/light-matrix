/**
 * @file memory.h
 *
 * @brief Functions for memory operations and allocation
 *
 * @author Dahua Lin
 */

#ifdef _MSC_VER
#pragma once
#endif

#ifndef LIGHTMAT_MEMORY_H_
#define LIGHTMAT_MEMORY_H_

#include <light_mat/common/basic_defs.h>
#include <cstring>
#include <new>
#include <limits>

#if LIGHTMAT_PLATFORM == LIGHTMAT_POSIX
#include <stdlib.h>
#elif LIGHTMAT_PLATFORM == LIGHTMAT_WIN32
#include <malloc.h>
#endif

namespace lmat
{
	/********************************************
	 *
	 *  Basic memory operations
	 *
	 ********************************************/

	template<typename T>
	LMAT_ENSURE_INLINE inline size_t nbytes(index_t n)
	{
		return static_cast<size_t>(n) * sizeof(T);
	}

	// copy

	template<typename T>
	LMAT_ENSURE_INLINE
	inline void copy_vec(const index_t n, const T *a, T *b)
	{
		std::memcpy(b, a, nbytes<T>(n));
	}

	template<typename T>
	LMAT_ENSURE_INLINE
	inline void copy_vec(const index_t n, const T* a, T* b, const index_t bstep)
	{
		for (index_t i = 0; i < n; ++i) b[i * bstep] = a[i];
	}

	template<typename T>
	LMAT_ENSURE_INLINE
	inline void copy_vec(const index_t n, const T* a, const index_t astep, T* b)
	{
		for (index_t i = 0; i < n; ++i) b[i] = a[i * astep];
	}

	template<typename T>
	LMAT_ENSURE_INLINE
	inline void copy_vec(const index_t n, const T* a, const index_t astep, T* b, const index_t bstep)
	{
		for (index_t i = 0; i < n; ++i) b[i * bstep] = a[i * astep];
	}

	// zero & fill

	template<typename T>
	LMAT_ENSURE_INLINE
	inline void set_zero_value(T& x) { x = T(0); }

	template<typename T>
	LMAT_ENSURE_INLINE
	inline void zero_vec(const index_t n, T *dst)
	{
		std::memset(dst, 0, nbytes<T>(n));
	}

	template<typename T>
	LMAT_ENSURE_INLINE
	inline void zero_vec(const index_t n, T *dst, const index_t step)
	{
		for (index_t i = 0; i < n; ++i) set_zero_value(dst[i * step]);
	}

	template<typename T>
	LMAT_ENSURE_INLINE
	inline void fill_vec(const T& v, const index_t n, T *dst)
	{
		for (index_t i = 0; i < n; ++i) dst[i] = v;
	}

	template<typename T>
	LMAT_ENSURE_INLINE
	inline void fill_vec(const T& v, const index_t n, T *dst, const index_t step)
	{
		for (index_t i = 0; i < n; ++i) dst[i * step] = v;
	}


	/********************************************
	 *
	 *  memory setter
	 *
	 ********************************************/

	template<class Derived>
	class IMemorySetter
	{
	public:
		LMAT_CRTP_REF
	};


	class zero_t : public IMemorySetter<zero_t>
	{
	};

	template<typename T>
	class copy_t : public IMemorySetter<copy_t<T> >
	{
	public:
		LMAT_ENSURE_INLINE copy_t(const T* src)
		: m_src(src)
		{ }

		LMAT_ENSURE_INLINE const T *source() const
		{
			return m_src;
		}

	private:
		const T* m_src;
	};


	template<typename T>
	class fill_t : public IMemorySetter<fill_t<T> >
	{
	public:
		LMAT_ENSURE_INLINE fill_t(const T& v)
		: m_val(v)
		{ }

		LMAT_ENSURE_INLINE T value() const
		{
			return m_val;
		}

	private:
		T m_val;
	};

	LMAT_ENSURE_INLINE
	inline zero_t zero() { return zero_t(); }

	template<typename T>
	LMAT_ENSURE_INLINE
	inline copy_t<T> copy_from(const T* src) { return copy_t<T>(src); }

	template<typename T>
	LMAT_ENSURE_INLINE
	inline fill_t<T> fill(const T& v) { return fill_t<T>(v); }


	template<typename T>
	LMAT_ENSURE_INLINE
	inline void apply(zero_t, index_t n, T *p)
	{
		zero_vec(n, p);
	}

	template<typename T>
	LMAT_ENSURE_INLINE
	inline void apply(const fill_t<T>& op, index_t n, T *p)
	{
		fill_vec(op.value(), n, p);
	}

	template<typename T>
	LMAT_ENSURE_INLINE
	inline void apply(const copy_t<T>& op, index_t n, T *p)
	{
		copy_vec(n, op.source(), p);
	}


	/********************************************
	 *
	 *  Aligned memory allocation
	 *
	 ********************************************/

#if LIGHTMAT_PLATFORM == LIGHTMAT_POSIX

	LMAT_ENSURE_INLINE
	inline void* aligned_allocate(size_t nbytes, unsigned int alignment)
	{
		char* p = 0;
		if (::posix_memalign((void**)(&p), alignment, nbytes) != 0)
		{
			throw std::bad_alloc();
		}
		return p;
	}

	LMAT_ENSURE_INLINE
	inline void aligned_release(void *p)
	{
		::free(p);
	}

#elif LIGHTMAT_PLATFORM == LIGHTMAT_WIN32

	LMAT_ENSURE_INLINE
	inline void* aligned_allocate(size_t nbytes, unsigned int alignment)
	{
		void* p = ::_aligned_malloc(nbytes, alignment));
		if (!p)
		{
			throw std::bad_alloc();
		}
		return p;
	}

	LMAT_ENSURE_INLINE
	inline void aligned_release(void *p)
	{
		::_aligned_free(p);
	}

#endif


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
    		return (pointer)aligned_allocate(n * sizeof(value_type), Align);
    	}

    	LMAT_ENSURE_INLINE
    	void deallocate(pointer p, size_type)
    	{
    		aligned_release(p);
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

#endif /* MEM_OP_H_ */
