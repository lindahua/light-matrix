/**
 * @file gamma_distr_internal.h
 *
 * @brief Internal implementation of PRNG for gamma distribution
 *
 * @author Dahua Lin
 */

#ifdef _MSC_VER
#pragma once
#endif

#ifndef LIGHTMAT_GAMMA_DISTR_INTERNAL_H_
#define LIGHTMAT_GAMMA_DISTR_INTERNAL_H_

#include <light_mat/random/uniform_real_distr.h>
#include <light_mat/random/exponential_distr.h>
#include <light_mat/math/math.h>


namespace lmat { namespace random { namespace internal {

	template<typename T, typename Method>
	struct std_gamma_distr_impl;

	/********************************************
	 *
	 *  Basic implementation
	 *  (borrowed from libc++)
	 *
	 ********************************************/

	template<typename T>
	struct _basic_gamma_vgen_params
	{
		LMAT_ENSURE_INLINE
		_basic_gamma_vgen_params(const T& a_)
		{
			a = a_;
			ra = math::rcp(a);
			b = a - T(1);
			rb = math::rcp(b);
			c = T(3) * a - T(0.75);
		}

		T a, ra;
		T b, rb;
		T c;
	};

	template<typename T, class RStream>
	T basic_gamma_variate_gen(RStream& rs, const _basic_gamma_vgen_params<T>& p)
	{
		std_exponential_distr<T> eg;

	    T x;
	    const T a = p.a;
	    const T ra = p.ra;
	    const T b = p.b;
	    const T rb = p.rb;
	    const T c = p.c;

	    if (a == 1)
	    {
	        x = eg(rs);
	    }
	    else if (a > 1)
	    {
	        while (true)
	        {
	            const T u = rand_real<T>::c0o1(rs);
	            const T v = rand_real<T>::c0o1(rs);
	            const T w = u * (1 - u);
	            if (w != 0)
	            {
	                const T y = math::sqrt(c / w) * (u - T(0.5));
	                x = b + y;
	                if (x >= 0)
	                {
	                    const T z = T(64) * w * w * w * v * v;
	                    if (z <= T(1) - T(2) * y * y / x)
	                        break;
	                    if (math::log(z) <= 2 * (b * math::log(x * rb) - y))
	                        break;
	                }
	            }
	        }
	    }
	    else  // a < 1
	    {
	        while (true)
	        {
	            const T u = rand_real<T>::c0o1(rs);
	            const T es = eg(rs);
	            if (u <= 1 - a)
	            {
	                x = math::pow(u, ra);
	                if (x <= es)
	                    break;
	            }
	            else
	            {
	                const T e = - math::log((1-u) * ra);
	                x = math::pow(1 - a + a * e, ra);
	                if (x <= e + es)
	                    break;
	            }
	        }
	    }

	    return x;
	}


	template<typename T>
	struct std_gamma_distr_impl<T, basic_>
	{
		_basic_gamma_vgen_params<T> params;

		LMAT_ENSURE_INLINE
		std_gamma_distr_impl(const T& a)
		: params(a) { }

		LMAT_ENSURE_INLINE
		T alpha() const
		{
			return params.a;
		}

		template<class RStream>
		LMAT_ENSURE_INLINE
		T operator() (RStream& rs) const
		{
			return basic_gamma_variate_gen(rs, params);
		}
	};


} } }

#endif
