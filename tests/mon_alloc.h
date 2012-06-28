/**
 * @file mon_alloc.h
 *
 * @brief Monitored memory allocation
 *
 * @author Dahua Lin
 */

#ifdef _MSC_VER
#pragma once
#endif

#ifndef MON_ALLOC_H_
#define MON_ALLOC_H_

#include <light_mat/core/basic_defs.h>
#include <map>

namespace lmat { namespace test {

	class memory_allocation_monitor : private noncopyable
	{
	public:
		memory_allocation_monitor()
		: m_pending_bytes(0)
		{
		}

		~memory_allocation_monitor()
		{
		}

		bool has_pending() const
		{
			return !m_ptrmap.empty();
		}

		size_t num_pending_sections() const
		{
			return m_ptrmap.size();
		}

		size_t num_pending_bytes() const
		{
			return m_pending_bytes;
		}

		bool verify(void *p, size_t nbytes) const
		{
			char *pc = static_cast<char*>(p);
			map_type::const_iterator it = m_ptrmap.find(pc);

			if (it != m_ptrmap.end())
			{
				const std::pair<char*, size_t>& e = *it;
				return e.second == nbytes;
			}
			else return false;
		}

		void* request(size_t nbytes)
		{
			char *p = new char[nbytes];

			m_ptrmap.insert(std::make_pair(p, nbytes));
			m_pending_bytes += nbytes;

			return p;
		}

		void release(void *p, size_t nbytes)
		{
			char *pc = static_cast<char*>(p);

			map_type::iterator it = m_ptrmap.find(pc);
			if (it != m_ptrmap.end())
			{
				const std::pair<char*, size_t>& e = *it;
				if (e.second == nbytes)
				{
					m_ptrmap.erase(it);
					m_pending_bytes -= nbytes;
				}
				else
				{
					throw std::runtime_error("The size to be released does not match the record.");
				}
			}
			else
			{
				throw std::runtime_error("captured invalid pointer to be released.");
			}

			delete[] static_cast<char*>(p);
		}

	private:
		size_t m_pending_bytes;

		typedef std::map<char*, size_t> map_type;
		map_type m_ptrmap;

	}; // end class memory_allocation_monitor



	extern memory_allocation_monitor global_memory_allocation_monitor;

    template<typename T>
    class monitored_allocator
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
    		typedef monitored_allocator<TOther> other;
    	};

    public:
    	monitored_allocator()
    	{
    	}

    	template<typename U>
    	monitored_allocator(const monitored_allocator<U>& r)
    	{
    	}

    	pointer address( reference x ) const
    	{
    		return &x;
    	}

    	const_pointer address( const_reference x ) const
    	{
    		return &x;
    	}

    	size_type max_size() const
    	{
    		return std::numeric_limits<size_type>::max() / sizeof(value_type);
    	}

    	pointer allocate(size_type n, const void* hint=0)
    	{
    		return static_cast<pointer>(
    				global_memory_allocation_monitor.request(n * sizeof(value_type)));
    	}

    	void deallocate(pointer p, size_type n)
    	{
    		global_memory_allocation_monitor.release(p, n * sizeof(value_type));
    	}

    	void construct (pointer p, const_reference val)
    	{
    		new ((void*)p) value_type(val);
    	}

    	void destroy (pointer p)
    	{
    		p->~value_type();
    	}

    }; // end class monitored_allocator



} }

#endif /* MON_ALLOC_H_ */
