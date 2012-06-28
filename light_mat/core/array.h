/**
 * @file array.h
 *
 * Array classes
 *
 * @author Dahua Lin
 */

#ifdef _MSC_VER
#pragma once
#endif

#ifndef LIGHTMAT_ARRAY_H_
#define LIGHTMAT_ARRAY_H_

#include <light_mat/core/mem_op.h>
#include <light_mat/core/mem_alloc.h>

#include <algorithm>

namespace lmat
{
	/**
	 * @brief The interface for arrays.
	 */
	template<class Derived, typename T>
	class IArray
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

		LMAT_ENSURE_INLINE const_reference operator[] (const index_type i) const
		{
			return derived().elem(i);
		}

		LMAT_ENSURE_INLINE reference operator[] (const index_type i)
		{
			return derived().elem(i);
		}

	}; // end class IBlock


	// Operations on arrays

	template<typename T, class Derived>
	LMAT_ENSURE_INLINE
	inline void copy(const T* src, IArray<Derived, T>& a)
	{
		copy_mem(a.nelems(), src, a.ptr_begin());
	}

	template<typename T, class Derived>
	LMAT_ENSURE_INLINE
	void copy(const IArray<Derived, T>& a, T* dst)
	{
		copy_mem(a.nelems(), a.ptr_begin(), dst);
	}

	template<typename T, class Derived>
	LMAT_ENSURE_INLINE
	void fill(IArray<Derived, T>& a, const T& v)
	{
		fill_mem(a.nelems(), a.ptr_begin(), v);
	}

	template<typename T, class Derived>
	LMAT_ENSURE_INLINE
	void zero(IArray<Derived, T>& a)
	{
		zero_mem(a.nelems(), a.ptr_begin());
	}

	template<typename T, class LDerived, class RDerived>
	LMAT_ENSURE_INLINE
	inline bool operator == (const IArray<LDerived, T>& a, const IArray<RDerived, T>& b)
	{
		return a.nelems() == b.nelems() && mem_equal(a.nelems(), a.ptr_begin(), b.ptr_begin());
	}

	template<typename T, class LDerived, class RDerived>
	LMAT_ENSURE_INLINE
	inline bool operator != (const IArray<LDerived, T>& a, const IArray<RDerived, T>& b)
	{
		return !(a == b);
	}

	template<typename T, class Derived>
	LMAT_ENSURE_INLINE
	inline bool elems_equal(const IArray<Derived, T>& B, const T& v)
	{
		return mem_equal(B.nelems(), B.ptr_begin(), v);
	}

	template<typename T, class Derived>
	LMAT_ENSURE_INLINE
	inline bool elems_equal(const IArray<Derived, T>& B, const T* r)
	{
		return mem_equal(B.nelems(), B.ptr_begin(), r);
	}


    /********************************************
     *
     *   Dynamic array
     *
     ********************************************/

	template<typename T, typename Allocator=aligned_allocator<T> >
	class darray : public IArray<darray<T, Allocator>,  T>
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
		explicit darray(const allocator_type& allocator = allocator_type())
		: m_allocator(allocator)
		, m_len(0)
		, m_ptr(LMAT_NULL)
		{
		}

		LMAT_ENSURE_INLINE
		explicit darray(index_type len, const allocator_type& allocator = allocator_type())
		: m_allocator(allocator)
		, m_len(len)
		, m_ptr(alloc(len))
		{
		}

		LMAT_ENSURE_INLINE
		darray(index_type len, const value_type& v, const allocator_type& allocator = allocator_type())
		: m_allocator(allocator)
		, m_len(len)
		, m_ptr(alloc(len))
		{
			fill_mem(len, m_ptr, v);
		}

		LMAT_ENSURE_INLINE
		darray(index_type len, const_pointer src, const allocator_type& allocator = allocator_type())
		: m_allocator(allocator)
		, m_len(len)
		, m_ptr(alloc(len))
		{
			copy_mem(m_len, src, m_ptr);
		}

		LMAT_ENSURE_INLINE
		darray(const darray& s, const allocator_type& allocator = allocator_type())
		: m_allocator(allocator)
		, m_len(s.nelems())
		, m_ptr(alloc(m_len))
		{
			copy_mem(m_len, s.ptr_begin(), m_ptr);
		}

		template<class OtherDerived>
		LMAT_ENSURE_INLINE
		darray(const IArray<OtherDerived, T>& s, const allocator_type& allocator = allocator_type())
		: m_allocator(allocator)
		, m_len(s.nelems())
		, m_ptr(alloc(m_len))
		{
			copy_mem(m_len, s.ptr_begin(), m_ptr);
		}

		LMAT_ENSURE_INLINE
		~darray()
		{
			dealloc(m_ptr);
			m_ptr = LMAT_NULL;  // make sure that it alarms when access after destructed
		}

		LMAT_ENSURE_INLINE
		void swap(darray& r)
		{
			using std::swap;

			swap(m_allocator, r.m_allocator);
			swap(m_len, r.m_len);
			swap(m_ptr, r.m_ptr);
		}

		LMAT_ENSURE_INLINE
		darray& operator = (const darray& r)
		{
			if (this != &r) assign(r);
			return *this;
		}

		template<class OtherDerived>
		LMAT_ENSURE_INLINE
		darray& operator = (const IArray<OtherDerived, T>& r)
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

		LMAT_ENSURE_INLINE
		const_reference operator[] (const index_type i) const
		{
			return m_ptr[i];
		}

		LMAT_ENSURE_INLINE
		reference operator[] (const index_type i)
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
		void assign(const IArray<OtherDerived, T>& r)
		{
			if (nelems() == r.nelems())  // no need to re-allocate memory
			{
				copy_mem(nelems(), r.ptr_begin(), this->ptr_begin());
			}
			else
			{
				darray tmp(r);
				swap(tmp);
			}
		}

	private:
		allocator_type m_allocator;
		index_type m_len;
		pointer m_ptr;

	}; // end class darray


	template<typename T, class Allocator>
	LMAT_ENSURE_INLINE
	void swap(darray<T, Allocator>& a, darray<T, Allocator>& b)
	{
		a.swap(b);
	}


    /********************************************
     *
     *   Scoped Array
     *
     ********************************************/

	template<typename T, typename Allocator=aligned_allocator<T> >
	class scoped_array : public IArray<scoped_array<T, Allocator>, T>, private noncopyable
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
		explicit scoped_array(index_type len, const allocator_type& allocator = allocator_type())
		: m_allocator(allocator)
		, m_len(len)
		, m_ptr(alloc(len))
		{
		}

		LMAT_ENSURE_INLINE
		scoped_array(index_type len, const value_type& v, const allocator_type& allocator = allocator_type())
		: m_allocator(allocator)
		, m_len(len)
		, m_ptr(alloc(len))
		{
			fill_mem(len, m_ptr, v);
		}

		LMAT_ENSURE_INLINE
		scoped_array(index_type len, const_pointer src, const allocator_type& allocator = allocator_type())
		: m_allocator(allocator)
		, m_len(len)
		, m_ptr(alloc(len))
		{
			copy_mem(len, src, m_ptr);
		}

		LMAT_ENSURE_INLINE
		~scoped_array()
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
		const_reference operator[] (const index_type i) const
		{
			return m_ptr[i];
		}

		LMAT_ENSURE_INLINE
		reference operator[] (const index_type i)
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

	}; // end class scoped_array



    /********************************************
     *
     *   Static Array
     *
     ********************************************/

	template<typename T, int N>
	class sarray : public IArray<sarray<T, N>, T>
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
		explicit sarray()
		{
		}

		LMAT_ENSURE_INLINE
		explicit sarray(const value_type& v)
		{
			fill_mem(N, m_arr, v);
		}

		LMAT_ENSURE_INLINE
		explicit sarray(const_pointer src)
		{
			copy_mem(N, src, m_arr);
		}

		LMAT_ENSURE_INLINE
		sarray(const sarray& src)
		{
			copy_mem(N, src.m_arr, m_arr);
		}

		template<class OtherDerived>
		LMAT_ENSURE_INLINE
		sarray(const IArray<OtherDerived, T>& src)
		{
			copy_mem(N, src.ptr_begin(), m_arr);
		}

		LMAT_ENSURE_INLINE
		sarray& operator = (const sarray& src)
		{
			if (this != &src)
			{
				copy_mem(N, src.m_arr, m_arr);
			}
			return *this;
		}

		template<class OtherDerived>
		LMAT_ENSURE_INLINE
		sarray& operator = (const IArray<OtherDerived, T>& src)
		{
			copy_mem(N, src.ptr_begin(), m_arr);
			return *this;
		}

		LMAT_ENSURE_INLINE
		void swap(sarray& r)
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

	}; // end class sarray

	template<typename T, int N>
	LMAT_ENSURE_INLINE
	void swap(sarray<T, N>& a, sarray<T, N>& b)
	{
		a.swap(b);
	}

}

#endif /* BLOCK_H_ */
