/**
 * @file block.h
 *
 * Memory block classes.
 *
 * @author Dahua Lin
 */

#ifdef _MSC_VER
#pragma once
#endif

#ifndef LIGHTMAT_BLOCK_H_
#define LIGHTMAT_BLOCK_H_

#include <light_mat/core/mem_op.h>
#include <light_mat/core/mem_alloc.h>

#include <algorithm>

namespace lmat
{
	/**
	 * @brief The interface for memory blocks.
	 */
	template<class Derived, typename T>
	class IBlock
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

		LMAT_CRTP_REF

	public:
		LMAT_ENSURE_INLINE size_type size() const
		{
			return (size_type)(derived().nelems());
		}

		LMAT_ENSURE_INLINE index_type nelems() const
		{
			return derived().nelems();
		}

		LMAT_ENSURE_INLINE const_pointer ptr_begin() const
		{
			return derived().ptr_begin();
		}

		LMAT_ENSURE_INLINE pointer ptr_begin()
		{
			return derived().ptr_begin();
		}

		LMAT_ENSURE_INLINE const_pointer ptr_end() const
		{
			return derived().ptr_end();
		}

		LMAT_ENSURE_INLINE pointer ptr_end()
		{
			return derived().ptr_end();
		}

		LMAT_ENSURE_INLINE const_reference operator[] (index_type i) const
		{
			return derived().elem(i);
		}

		LMAT_ENSURE_INLINE reference operator[] (index_type i)
		{
			return derived().elem(i);
		}

	}; // end class IBlock


	// Operations on blocks

	template<typename T, class Derived>
	LMAT_ENSURE_INLINE
	inline void copy(IBlock<Derived, T>& blk, const T* src)
	{
		copy_mem(blk.nelems(), src, blk.ptr_begin());
	}

	template<typename T, class Derived>
	LMAT_ENSURE_INLINE
	void copy(const IBlock<Derived, T>& blk, T* dst)
	{
		copy_mem(blk.nelems(), blk.ptr_begin(), dst);
	}

	template<typename T, class Derived>
	LMAT_ENSURE_INLINE
	void fill(IBlock<Derived, T>& blk, const T& v)
	{
		fill_mem(blk.nelems(), blk.ptr_begin(), v);
	}

	template<typename T, class Derived>
	LMAT_ENSURE_INLINE
	void zero(IBlock<Derived, T>& blk)
	{
		zero_mem(blk.nelems(), blk.ptr_begin());
	}

	template<typename T, class LDerived, class RDerived>
	LMAT_ENSURE_INLINE
	inline bool is_equal(const IBlock<LDerived, T>& B1, const IBlock<RDerived, T>& B2)
	{
		return B1.nelems() == B2.nelems() && mem_equal(B1.nelems(), B1.ptr_begin(), B2.ptr_begin());
	}

	template<typename T, class Derived>
	LMAT_ENSURE_INLINE
	inline bool all_equal(const IBlock<Derived, T>& B, const T& v)
	{
		return mem_equal(B.nelems(), B.ptr_begin(), v);
	}

	template<typename T, class Derived>
	LMAT_ENSURE_INLINE
	inline bool all_equal(const IBlock<Derived, T>& B, const T* r)
	{
		return mem_equal(B.nelems(), B.ptr_begin(), r);
	}


    /********************************************
     *
     *   block
     *
     ********************************************/

	template<typename T, typename Allocator=aligned_allocator<T> >
	class block : public IBlock<block<T, Allocator>,  T>
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

		typedef index_t index_type;

	public:
		LMAT_ENSURE_INLINE
		explicit block(const allocator_type& allocator = allocator_type())
		: m_allocator(allocator)
		, m_len(0)
		, m_ptr(LMAT_NULL)
		{
		}

		LMAT_ENSURE_INLINE
		explicit block(index_type len, const allocator_type& allocator = allocator_type())
		: m_allocator(allocator)
		, m_len(len)
		, m_ptr(alloc(len))
		{
		}

		LMAT_ENSURE_INLINE
		block(index_type len, const value_type& v, const allocator_type& allocator = allocator_type())
		: m_allocator(allocator)
		, m_len(len)
		, m_ptr(alloc(len))
		{
			fill_mem(len, m_ptr, v);
		}

		LMAT_ENSURE_INLINE
		block(index_type len, const_pointer src, const allocator_type& allocator = allocator_type())
		: m_allocator(allocator)
		, m_len(len)
		, m_ptr(alloc(len))
		{
			copy_mem(m_len, src, m_ptr);
		}

		LMAT_ENSURE_INLINE
		block(const block& s, const allocator_type& allocator = allocator_type())
		: m_allocator(allocator)
		, m_len(s.nelems())
		, m_ptr(alloc(m_len))
		{
			copy_mem(m_len, s.ptr_begin(), m_ptr);
		}

		template<class OtherDerived>
		LMAT_ENSURE_INLINE
		block(const IBlock<OtherDerived, T>& s, const allocator_type& allocator = allocator_type())
		: m_allocator(allocator)
		, m_len(s.nelems())
		, m_ptr(alloc(m_len))
		{
			copy_mem(m_len, s.ptr_begin(), m_ptr);
		}

		LMAT_ENSURE_INLINE
		~block()
		{
			dealloc(m_ptr);
			m_ptr = LMAT_NULL;  // make sure that it alarms when access after destructed
		}

		LMAT_ENSURE_INLINE
		void swap(block& r)
		{
			using std::swap;

			swap(m_allocator, r.m_allocator);
			swap(m_len, r.m_len);
			swap(m_ptr, r.m_ptr);
		}

		LMAT_ENSURE_INLINE
		block& operator = (const block& r)
		{
			if (this != &r) assign(r);
			return *this;
		}

		template<class OtherDerived>
		LMAT_ENSURE_INLINE
		block& operator = (const IBlock<OtherDerived, T>& r)
		{
			assign(r);
			return *this;
		}

		LMAT_ENSURE_INLINE
		const allocator_type& get_allocator() const
		{
			return m_allocator;
		}

		LMAT_ENSURE_INLINE
		void resize(index_type n)
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

		LMAT_ENSURE_INLINE index_type nelems() const
		{
			return m_len;
		}

		LMAT_ENSURE_INLINE const_pointer ptr_begin() const
		{
			return m_ptr;
		}

		LMAT_ENSURE_INLINE pointer ptr_begin()
		{
			return m_ptr;
		}

		LMAT_ENSURE_INLINE const_pointer ptr_end() const
		{
			return m_ptr + m_len;
		}

		LMAT_ENSURE_INLINE pointer ptr_end()
		{
			return m_ptr + m_len;
		}

		LMAT_ENSURE_INLINE const_reference operator[] (index_type i) const
		{
			return m_ptr[i];
		}

		LMAT_ENSURE_INLINE reference operator[] (index_type i)
		{
			return m_ptr[i];
		}


	private:
		LMAT_ENSURE_INLINE
		pointer alloc(index_type n)
		{
			return n > 0 ? m_allocator.allocate(size_type(n)) : LMAT_NULL;
		}

		LMAT_ENSURE_INLINE
		void dealloc(pointer p)
		{
			if (m_ptr) m_allocator.deallocate(p, (size_type)m_len);
		}

		LMAT_ENSURE_INLINE
		void reset_size(index_type n)
		{
			pointer p = n > 0 ? m_allocator.allocate(size_type(n)) : LMAT_NULL;
			dealloc(m_ptr);
			m_ptr = p;
			m_len = n;
		}

		template<class OtherDerived>
		LMAT_ENSURE_INLINE
		void assign(const IBlock<OtherDerived, T>& r)
		{
			if (nelems() == r.nelems())  // no need to re-allocate memory
			{
				copy_mem(nelems(), r.ptr_begin(), this->ptr_begin());
			}
			else
			{
				block tmp(r);
				swap(tmp);
			}
		}

	private:
		allocator_type m_allocator;
		index_type m_len;
		pointer m_ptr;

	}; // end class block


	template<typename T, class Allocator>
	LMAT_ENSURE_INLINE
	void swap(block<T, Allocator>& a, block<T, Allocator>& b)
	{
		a.swap(b);
	}


    /********************************************
     *
     *   Scoped Block
     *
     ********************************************/

	template<typename T, typename Allocator=aligned_allocator<T> >
	class scoped_block : public IBlock<scoped_block<T, Allocator>, T>, private noncopyable
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

		typedef index_t index_type;

	public:
		LMAT_ENSURE_INLINE
		explicit scoped_block(index_type len, const allocator_type& allocator = allocator_type())
		: m_allocator(allocator)
		, m_len(len)
		, m_ptr(alloc(len))
		{
		}

		LMAT_ENSURE_INLINE
		scoped_block(index_type len, const value_type& v, const allocator_type& allocator = allocator_type())
		: m_allocator(allocator)
		, m_len(len)
		, m_ptr(alloc(len))
		{
			fill_mem(len, m_ptr, v);
		}

		LMAT_ENSURE_INLINE
		scoped_block(index_type len, const_pointer src, const allocator_type& allocator = allocator_type())
		: m_allocator(allocator)
		, m_len(len)
		, m_ptr(alloc(len))
		{
			copy_mem(len, src, m_ptr);
		}

		LMAT_ENSURE_INLINE
		~scoped_block()
		{
			dealloc(m_ptr);
		}

		LMAT_ENSURE_INLINE
		const allocator_type& get_allocator() const
		{
			return m_allocator;
		}

	public:
		LMAT_ENSURE_INLINE size_type size() const
		{
			return static_cast<size_type>(m_len);
		}

		LMAT_ENSURE_INLINE index_type nelems() const
		{
			return m_len;
		}

		LMAT_ENSURE_INLINE const_pointer ptr_begin() const
		{
			return m_ptr;
		}

		LMAT_ENSURE_INLINE pointer ptr_begin()
		{
			return m_ptr;
		}

		LMAT_ENSURE_INLINE const_pointer ptr_end() const
		{
			return m_ptr + m_len;
		}

		LMAT_ENSURE_INLINE pointer ptr_end()
		{
			return m_ptr + m_len;
		}

		LMAT_ENSURE_INLINE
		const_reference operator[] (index_type i) const
		{
			return m_ptr[i];
		}

		LMAT_ENSURE_INLINE
		reference operator[] (index_type i)
		{
			return m_ptr[i];
		}

	private:
		LMAT_ENSURE_INLINE
		pointer alloc(index_type n)
		{
			return n > 0 ? m_allocator.allocate(size_type(n)) : LMAT_NULL;
		}

		LMAT_ENSURE_INLINE
		void dealloc(pointer p)
		{
			if (m_ptr) m_allocator.deallocate(p, (size_type)m_len);
		}

	private:
		allocator_type m_allocator;
		index_type m_len;
		pointer m_ptr;

	}; // end class scoped_block



    /********************************************
     *
     *   Static Block
     *
     ********************************************/

	template<typename T, int N>
	class static_block : public IBlock<static_block<T, N>, T>
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
		explicit static_block()
		{
		}

		LMAT_ENSURE_INLINE
		explicit static_block(const value_type& v)
		{
			fill_mem(N, m_arr, v);
		}

		LMAT_ENSURE_INLINE
		explicit static_block(const_pointer src)
		{
			copy_mem(N, src, m_arr);
		}

		LMAT_ENSURE_INLINE
		static_block(const static_block& src)
		{
			copy_mem(N, src.m_arr, m_arr);
		}

		template<class OtherDerived>
		LMAT_ENSURE_INLINE
		static_block(const IBlock<OtherDerived, T>& src)
		{
			copy_mem(N, src.ptr_begin(), m_arr);
		}

		LMAT_ENSURE_INLINE
		static_block& operator = (const static_block& src)
		{
			if (this != &src)
			{
				copy_mem(N, src.m_arr, m_arr);
			}
			return *this;
		}

		template<class OtherDerived>
		LMAT_ENSURE_INLINE
		static_block& operator = (const IBlock<OtherDerived, T>& src)
		{
			copy_mem(N, src.ptr_begin(), m_arr);
			return *this;
		}

		LMAT_ENSURE_INLINE
		void swap(static_block& r)
		{
			T tmp[N];
			copy_mem(N, m_arr, tmp);
			copy_mem(N, r.m_arr, m_arr);
			copy_mem(N, tmp, r.m_arr);
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

		LMAT_ENSURE_INLINE const_pointer ptr_begin() const
		{
			return m_arr;
		}

		LMAT_ENSURE_INLINE pointer ptr_begin()
		{
			return m_arr;
		}

		LMAT_ENSURE_INLINE const_pointer ptr_end() const
		{
			return m_arr + N;
		}

		LMAT_ENSURE_INLINE pointer ptr_end()
		{
			return m_arr + N;
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

	}; // end class static_block

}

#endif /* BLOCK_H_ */
