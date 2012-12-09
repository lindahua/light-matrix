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

namespace lmat
{
	/********************************************
	 *
	 *  Iterator to access interval memory
	 *
	 ********************************************/

	template<typename T> // T can be a const type
	class step_ptr_t
	{
	public:
		LMAT_ENSURE_INLINE
		step_ptr_t(T *p, const index_t& step)
		: _p(p), _step(step)
		{ }

		LMAT_ENSURE_INLINE
		index_t step() const
		{
			return _step;
		}

		LMAT_ENSURE_INLINE
		T& operator[] (index_t i) const
		{
			return _p[i * _step];
		}

		LMAT_ENSURE_INLINE
		T& operator * () const
		{
			return *_p;
		}

		LMAT_ENSURE_INLINE
		T* operator -> () const
		{
			return _p;
		}

		LMAT_ENSURE_INLINE
		step_ptr_t& operator++ ()
		{
			_p += _step;
			return *this;
		}

		LMAT_ENSURE_INLINE
		step_ptr_t operator-- ()
		{
			_p -= _step;
			return *this;
		}

		LMAT_ENSURE_INLINE
		step_ptr_t operator++ (int)
		{
			step_ptr_t old(*this);
			_p += _step;
			return old;
		}

		LMAT_ENSURE_INLINE
		step_ptr_t& operator-- (int)
		{
			step_ptr_t old(*this);
			_p -= _step;
			return old;
		}

		LMAT_ENSURE_INLINE
		step_ptr_t& operator += (index_t n)
		{
			_p += _step * n;
			return *this;
		}

		LMAT_ENSURE_INLINE
		step_ptr_t& operator -= (index_t n)
		{
			_p -= _step * n;
			return *this;
		}

		LMAT_ENSURE_INLINE
		step_ptr_t operator + (index_t n) const
		{
			return step_ptr_t(_p + _step * n, _step);
		}

		LMAT_ENSURE_INLINE
		step_ptr_t operator - (index_t n) const
		{
			return step_ptr_t(_p - _step * n, _step);
		}

		LMAT_ENSURE_INLINE
		index_t operator - (const step_ptr_t& r) const
		{
			return (r._p - _p) / _step;
		}

		LMAT_ENSURE_INLINE
		bool operator == (const step_ptr_t& r) const
		{
			return _p == r._p;
		}

		LMAT_ENSURE_INLINE
		bool operator != (const step_ptr_t& r) const
		{
			return _p != r._p;
		}

	private:
		T * _p;
		const index_t _step;
	};


	template<typename T>
	LMAT_ENSURE_INLINE
	inline step_ptr_t<T> step_ptr(T* p, index_t step)
	{
		return step_ptr_t<T>(p, step);
	}


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
	inline void copy_vec(index_t n, const T *a, T *b)
	{
		std::memcpy(b, a, nbytes<T>(n));
	}

	template<typename SrcIter, typename DstIter>
	LMAT_ENSURE_INLINE
	inline void copy_vec(index_t n, SrcIter a, DstIter b)
	{
		for (index_t i = 0; i < n; ++i) b[i] = a[i];
	}

	// zero & fill

	template<typename T>
	LMAT_ENSURE_INLINE
	inline void set_zero_value(T& x) { x = T(0); }

	template<typename T>
	LMAT_ENSURE_INLINE
	inline void zero_vec(index_t n, T *dst)
	{
		std::memset(dst, 0, nbytes<T>(n));
	}

	template<typename DstIter>
	LMAT_ENSURE_INLINE
	inline void zero_vec(index_t n, DstIter dst)
	{
		for (index_t i = 0; i < n; ++i) set_zero_value(dst[i]);
	}

	template<typename DstIter, typename T>
	LMAT_ENSURE_INLINE
	inline void fill_vec(const index_t n, DstIter dst, const T& v)
	{
		for (index_t i = 0; i < n; ++i) dst[i] = v;
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
		fill_vec(n, p, op.value());
	}

	template<typename T>
	LMAT_ENSURE_INLINE
	inline void apply(const copy_t<T>& op, index_t n, T *p)
	{
		copy_vec(n, op.source(), p);
	}


}

#endif /* MEM_OP_H_ */
