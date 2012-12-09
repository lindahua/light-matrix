/**
 * @file block.h
 *
 * Block classes
 *
 * @author Dahua Lin
 */

#ifdef _MSC_VER
#pragma once
#endif

#ifndef LIGHTMAT_BLOCK_H_
#define LIGHTMAT_BLOCK_H_

#include <light_mat/common/memory.h>
#include <light_mat/common/memalloc.h>
#include <algorithm>

namespace lmat
{

    /********************************************
     *
     *   Dynamic block
     *
     ********************************************/

	template<typename T, typename Allocator=aligned_allocator<T> >
	class dblock
	{
	public:
		typedef T value_type;
		typedef Allocator allocator_type;

		typedef typename allocator_type::size_type size_type;
		typedef typename allocator_type::difference_type difference_type;

		typedef typename allocator_type::pointer pointer;
		typedef typename allocator_type::reference reference;
		typedef typename allocator_type::const_pointer const_pointer;
		typedef typename allocator_type::const_reference const_reference;

	public:
		LMAT_ENSURE_INLINE
		explicit dblock(const allocator_type& allocator = allocator_type())
		: m_allocator(allocator)
		, m_len(0)
		, m_ptr(nullptr)
		{
		}

		LMAT_ENSURE_INLINE
		explicit dblock(index_t len, const allocator_type& allocator = allocator_type())
		: m_allocator(allocator)
		, m_len(len)
		, m_ptr(alloc(len))
		{
		}

		template<class Setter>
		LMAT_ENSURE_INLINE
		explicit dblock(index_t len, const IMemorySetter<Setter>& setter,
				const allocator_type& allocator = allocator_type())
		: m_allocator(allocator)
		, m_len(len)
		, m_ptr(alloc(len))
		{
			apply(setter.derived(), len, m_ptr);
		}

		LMAT_ENSURE_INLINE
		dblock(const dblock& s, const allocator_type& allocator = allocator_type())
		: m_allocator(allocator)
		, m_len(s.nelems())
		, m_ptr(alloc(m_len))
		{
			copy_vec(m_len, s.ptr_data(), m_ptr);
		}

		LMAT_ENSURE_INLINE
		dblock(dblock&& s)
		: m_allocator(allocator_type())
		, m_len(0)
		, m_ptr(nullptr)
		{
			swap(s);
		}

		LMAT_ENSURE_INLINE
		~dblock()
		{
			dealloc(m_ptr);
			m_ptr = nullptr;  // make sure that it alarms when access after destructed
		}

		LMAT_ENSURE_INLINE
		void swap(dblock& r)
		{
			using std::swap;

			swap(m_allocator, r.m_allocator);
			swap(m_len, r.m_len);
			swap(m_ptr, r.m_ptr);
		}

		LMAT_ENSURE_INLINE
		dblock& operator = (const dblock& r)
		{
			if (this != &r)
			{
				if (nelems() == r.nelems())  // no need to re-allocate memory
				{
					copy_vec(nelems(), r.ptr_data(), this->ptr_data());
				}
				else
				{
					dblock tmp(r);
					swap(tmp);
				}
			}
			return *this;
		}


		LMAT_ENSURE_INLINE
		dblock& operator = (dblock&& r)
		{
			if (this != &r)
			{
				dealloc(m_ptr);
				m_ptr = nullptr;
				m_len = 0;

				swap(r);
			}
			return *this;
		}

		template<class Setter>
		LMAT_ENSURE_INLINE
		dblock& operator = (const IMemorySetter<Setter>& setter)
		{
			apply(setter.derived(), m_len, m_ptr);
			return *this;
		}


		LMAT_ENSURE_INLINE
		const allocator_type& get_allocator() const
		{
			return m_allocator;
		}

		LMAT_ENSURE_INLINE
		void resize(index_t n)
		{
			if (n != this->nelems())
			{
				reset_size(n);
			}
		}

	public:
		LMAT_ENSURE_INLINE size_type size() const
		{
			return static_cast<size_type>(m_len);
		}

		LMAT_ENSURE_INLINE index_t nelems() const
		{
			return m_len;
		}

		LMAT_ENSURE_INLINE const_pointer ptr_data() const
		{
			return m_ptr;
		}

		LMAT_ENSURE_INLINE pointer ptr_data()
		{
			return m_ptr;
		}

		LMAT_ENSURE_INLINE
		const_reference operator[] (const index_t i) const
		{
			return m_ptr[i];
		}

		LMAT_ENSURE_INLINE
		reference operator[] (const index_t i)
		{
			return m_ptr[i];
		}


	private:
		LMAT_ENSURE_INLINE
		pointer alloc(index_t n)
		{
			return n > 0 ? m_allocator.allocate(size_type(n)) : nullptr;
		}

		LMAT_ENSURE_INLINE
		void dealloc(pointer p)
		{
			if (m_ptr) m_allocator.deallocate(p, (size_type)m_len);
		}

		LMAT_ENSURE_INLINE
		void reset_size(index_t n)
		{
			pointer p = n > 0 ? m_allocator.allocate(size_type(n)) : nullptr;
			dealloc(m_ptr);
			m_ptr = p;
			m_len = n;
		}

	private:
		allocator_type m_allocator;
		index_t m_len;
		pointer m_ptr;

	}; // end class dblock


	template<typename T, class Allocator>
	LMAT_ENSURE_INLINE
	inline void swap(dblock<T, Allocator>& a, dblock<T, Allocator>& b)
	{
		a.swap(b);
	}


    /********************************************
     *
     *   Static Block
     *
     ********************************************/

	template<typename T, int N>
	class sblock
	{
	public:
		typedef T value_type;

		typedef size_t size_type;
		typedef ptrdiff_t difference_type;

		typedef T* pointer;
		typedef T& reference;
		typedef const T* const_pointer;
		typedef const T& const_reference;

		typedef index_t index_type;

	public:
		LMAT_ENSURE_INLINE
		explicit sblock()
		{
		}

		template<class Setter>
		LMAT_ENSURE_INLINE
		explicit sblock(const IMemorySetter<Setter>& setter)
		{
			apply(setter.derived(), N, m_arr);
		}


		LMAT_ENSURE_INLINE
		sblock(const sblock& src)
		{
			copy_vec(N, src.m_arr, m_arr);
		}

		LMAT_ENSURE_INLINE
		sblock& operator = (const sblock& src)
		{
			if (this != &src)
			{
				copy_vec(N, src.m_arr, m_arr);
			}
			return *this;
		}

		template<class Setter>
		LMAT_ENSURE_INLINE
		sblock& operator = (const IMemorySetter<Setter>& setter)
		{
			apply(setter.derived(), N, m_arr);
		}

		LMAT_ENSURE_INLINE
		void swap(sblock& r)
		{
			T tmp[N];
			copy_vec(N, m_arr, tmp);
			copy_vec(N, r.m_arr, m_arr);
			copy_vec(N, tmp, r.m_arr);
		}

	public:
		LMAT_ENSURE_INLINE size_type size() const
		{
			return static_cast<size_t>(N);
		}

		LMAT_ENSURE_INLINE index_type nelems() const
		{
			return N;
		}

		LMAT_ENSURE_INLINE const_pointer ptr_data() const
		{
			return m_arr;
		}

		LMAT_ENSURE_INLINE pointer ptr_data()
		{
			return m_arr;
		}

		LMAT_ENSURE_INLINE
		const_reference operator[] (index_type i) const
		{
			return m_arr[i];
		}

		LMAT_ENSURE_INLINE
		reference operator[] (index_type i)
		{
			return m_arr[i];
		}

	private:
		T m_arr[N];

	}; // end class sblock

	template<typename T, int N>
	LMAT_ENSURE_INLINE
	inline void swap(sblock<T, N>& a, sblock<T, N>& b)
	{
		a.swap(b);
	}

}

#endif /* BLOCK_H_ */
