/**
 * @file discrete_distr.h
 *
 * @brief Discrete distribution PRNG
 *
 * @author Dahua Lin
 */

#ifdef _MSC_VER
#pragma once
#endif

#ifndef LIGHTMAT_DISCRETE_DISTR_H_
#define LIGHTMAT_DISCRETE_DISTR_H_

#include <light_mat/random/uniform_real_distr.h>
#include <light_mat/matrix/dense_matrix.h>
#include <algorithm>


namespace lmat { namespace random {

	/********************************************
	 *
	 *  convenient functions
	 *
	 ********************************************/

	template<typename TI, typename TReal>
	LMAT_ENSURE_INLINE
	inline TI dd_draw(TI n, const TReal* weights, const TReal& tw, const TReal& u)  // u ~ [0, 1)
	{
		TReal v = tw * u;
		TI i(0);
		TI max_i = n - 1;

		TReal cw = weights[0];

		while(v > cw && i < max_i)
			cw += weights[++i];

		return i;
	}

	template<typename TI, typename TReal, class RStream>
	LMAT_ENSURE_INLINE
	inline TI dd_draw(TI n, const TReal* weights, const TReal& tw, RStream& rs)
	{
		std_uniform_real_distr<double> ud;
		return dd_draw(n, weights, tw, ud(rs));
	}


	/********************************************
	 *
	 *  internal implementation
	 *
	 ********************************************/

	namespace internal
	{
		template<typename TI, typename Method>
		struct discrete_distr_impl;

		template<typename TI>
		struct discrete_distr_impl<TI, naive_>
		{
		public:
			template<typename InputIter>
			LMAT_ENSURE_INLINE
			explicit discrete_distr_impl(InputIter first, InputIter last)
			: m_weights(), m_total(0.0), m_n(0)
			{
				// scan
				for (InputIter it = first; it != last; ++it, ++m_n)
					m_total += double(*it);

				m_inv_total = 1.0 / m_total;

				// copy
				m_weights.require_size((index_t)m_n);
				std::copy_n(first, (size_t)m_n, m_weights.ptr_data());
			}

			template<class RStream>
			LMAT_ENSURE_INLINE
			TI operator() (RStream& rs) const
			{
				return dd_draw(m_n, m_weights.ptr_data(), m_total, rs);
			}

			LMAT_ENSURE_INLINE
			TI n() const
			{
				return m_n;
			}

			LMAT_ENSURE_INLINE
			double p(TI x) const
			{
				return m_weights[x] * m_inv_total;
			}

		private:
			dense_col<double> m_weights;
			double m_total;
			double m_inv_total;
			TI m_n;
		};
	}

	template<typename TI, typename Method>
	class discrete_distr  // Uniform over [0, n)
	{
		typedef internal::discrete_distr_impl<TI, Method> impl_t;

	public:
		typedef TI result_type;

		template<typename InputIter>
		LMAT_ENSURE_INLINE
		explicit discrete_distr(InputIter first, InputIter last)
		: m_impl(first, last) { }

		LMAT_ENSURE_INLINE
		TI n() const
		{
			return m_impl.n();
		}

		LMAT_ENSURE_INLINE
		double p(TI x) const
		{
			return m_impl.p(x);
		}

		template<class RStream>
		LMAT_ENSURE_INLINE
		TI operator() (RStream& rs) const
		{
			return m_impl(rs);
		}

	private:
		impl_t m_impl;
	};


} }

#endif
