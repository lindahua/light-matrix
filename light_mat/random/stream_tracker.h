/**
 * @file stream_tracker.h
 *
 * @brief A stream tracker to support the implementation of cache-based streams
 *
 * @author Dahua Lin
 */

#ifdef _MSC_VER
#pragma once
#endif

#ifndef LIGHTMAT_STREAM_TRACKER_H_
#define LIGHTMAT_STREAM_TRACKER_H_

#include <light_mat/common/basic_defs.h>

namespace lmat { namespace random {

	template<typename T> // T is the unit type
	class stream_tracker
	{
		static const unsigned int U = (unsigned int)sizeof(T);

	public:
		LMAT_ENSURE_INLINE
		explicit stream_tracker(size_t len)
		: m_len(len), m_i(len) { }

		LMAT_ENSURE_INLINE
		size_t length() const  // in terms of #units
		{
			return m_len;
		}

		LMAT_ENSURE_INLINE
		size_t offset() const // in terms of units
		{
			return m_i;
		}

		LMAT_ENSURE_INLINE
		size_t remain() const // in terms of units
		{
			return m_i < m_len ? m_len - m_i : 0;
		}

		LMAT_ENSURE_INLINE
		bool is_end() const
		{
			return m_i >= m_len;
		}

		LMAT_ENSURE_INLINE
		void set_offset(size_t i)
		{
			m_i = i;
		}

		LMAT_ENSURE_INLINE
		void forward(size_t n) // n units
		{
			m_i += n;
		}

		LMAT_ENSURE_INLINE
		void set_end()  // call this after the stream is exhausted, but next state is not ready
		{
			m_i = m_len;
		}

		LMAT_ENSURE_INLINE  // call this after the state is refresh
		void rewind()
		{
			m_i = 0;
		}

		LMAT_ENSURE_INLINE
		void to_boundary(bdtags::dbl)
		{
			if (m_i & 1) ++m_i;
		}

		LMAT_ENSURE_INLINE
		void to_boundary(bdtags::quad)
		{
			if (m_i & 3) m_i = ((m_i >> 2) + 1) << 2;
		}

		LMAT_ENSURE_INLINE
		void to_boundary(bdtags::oct)
		{
			if (m_i & 7) m_i = ((m_i >> 3) + 1) << 3;
		}

		LMAT_ENSURE_INLINE
		void to_boundary(bdtags::hex)
		{
			if (m_i & 15) m_i = ((m_i >> 4) + 1) << 4;
		}

		LMAT_ENSURE_INLINE
		void forward_bytes(size_t bytes)  // pre-condition: bytes <= remain() * U
		{
			size_t n = int_div<U>::quo(bytes);
			if (int_div<U>::rem(bytes)) ++n;

			forward(n);
		}

	private:
		const size_t m_len;
		size_t m_i;
	};

} }

#endif
