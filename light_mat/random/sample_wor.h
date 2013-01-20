/**
 * @file sample_wor.h
 *
 * @brief Sampling without replacement
 *
 * @author Dahua Lin
 */

#ifdef _MSC_VER
#pragma once
#endif

#ifndef LIGHTMAT_SAMPLE_WOR_H_
#define LIGHTMAT_SAMPLE_WOR_H_

#include <light_mat/matrix/dense_matrix.h>
#include <light_mat/random/uniform_int_distr.h>
#include <unordered_set>

namespace lmat { namespace random {

	template<typename TI=uint32_t, class RStream=default_rand_stream>
	class rand_shuffle_enumerator
	{
	public:
		typedef TI value_type;

		LMAT_ENSURE_INLINE
		rand_shuffle_enumerator(RStream& rstream, TI n)
		: m_rstream(rstream), m_seq((index_t)n), m_i(0)
		{
			reset();
		}

		void reset()
		{
			const index_t n = m_seq.nelems();
			for (index_t i = 0; i < n; ++i)
			{
				m_seq[i] = static_cast<TI>(i);
			}
			m_i = 0;
		}

		LMAT_ENSURE_INLINE
		index_t length() const
		{
			return m_seq.nelems();
		}

		LMAT_ENSURE_INLINE
		index_t remain() const
		{
			return m_seq.nelems() - m_i;
		}

		LMAT_ENSURE_INLINE
		bool is_end() const
		{
			return m_i == m_seq.nelems();
		}

		LMAT_ENSURE_INLINE
		TI next()
		{
			std_uniform_int_distr<index_t> d(remain() - 1);
			index_t i = m_i + d(m_rstream);

			TI t = m_seq[i];
			m_seq[i] = m_seq[m_i];
			m_seq[m_i++] = t;

			return t;
		}

	private:
		RStream& m_rstream;
		dense_col<TI> m_seq;
		index_t m_i;
	};


	template<typename TI=uint32_t, class RStream=default_rand_stream>
	class past_avoid_rand_enumerator
	{
	public:
		typedef TI value_type;

		LMAT_ENSURE_INLINE
		past_avoid_rand_enumerator(RStream& rstream, TI n)
		: m_rstream(rstream), m_past((size_t)n)
		, m_len((index_t)n), m_i(0), m_distr((TI)(n - 1))
		{
			reset();
		}

		void reset()
		{
			m_past.clear();
			m_i = 0;
		}

		LMAT_ENSURE_INLINE
		index_t length() const
		{
			return m_len;
		}

		LMAT_ENSURE_INLINE
		index_t remain() const
		{
			return m_len - m_i;
		}

		LMAT_ENSURE_INLINE
		bool is_end() const
		{
			return m_i == m_len;
		}

		TI next()
		{
			TI x(0);

			for(;;)
			{
				x = m_distr(m_rstream);

				if (m_past.find(x) == m_past.end())
				{
					m_past.insert(x);
					break;
				}
			}

			++ m_i;
			return x;
		}

	private:
		RStream& m_rstream;
		std::unordered_set<TI> m_past;
		index_t m_len;
		index_t m_i;
		std_uniform_int_distr<TI> m_distr;
	};



} }

#endif
