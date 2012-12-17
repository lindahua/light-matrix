/**
 * @file mat_reduce_internal.h
 *
 * Internal implementation of matrix reduction
 * 
 * @author Dahua Lin 
 */

#ifdef _MSC_VER
#pragma once
#endif

#ifndef LIGHTMAT_MAT_REDUCE_INTERNAL_H_
#define LIGHTMAT_MAT_REDUCE_INTERNAL_H_

#include <light_mat/mateval/vec_accessors.h>
#include <light_mat/math/math_functors.h>
#include <light_mat/math/simd.h>


namespace lmat { namespace internal {

	template<class Reader>
	inline typename Reader::scalar_type
	sum_impl(unsigned int len, const Reader& rd)
	{
		typedef typename Reader::scalar_type T;
		typedef typename Reader::pack_type pack_t;

		const unsigned int pw = pack_t::pack_width;
		unsigned int npacks = int_div<pw>::quo(len);
		T r;
		unsigned int i;

		if (npacks)
		{
			pack_t a0, a1, a2, a3;

			unsigned int m4 = npacks >> 2;
			if (m4)
			{
				const unsigned int w4 = pw << 2;
				unsigned int l4 = m4 * w4;

				a0 = rd.pack(0);
				a1 = rd.pack(pw);
				a2 = rd.pack(pw * 2);
				a3 = rd.pack(pw * 3);

				for (i = w4; i < l4; i += w4)
				{
					a0 += rd.pack(i);
					a1 += rd.pack(i + pw);
					a2 += rd.pack(i + pw * 2);
					a3 += rd.pack(i + pw * 3);
				}

				a0 += a2;
				a1 += a3;

				if (npacks & 2)
				{
					a0 += rd.pack(i);
					a1 += rd.pack(i + pw);
					i += pw * 2;
				}

				a0 += a1;

				if (npacks & 1)
				{
					a0 += rd.pack(i);
					i += pw;
				}
			}
			else if (npacks >> 1)
			{
				a0 = rd.pack(0) + rd.pack(pw);
				i = pw << 1;

				if (npacks & 1)
				{
					a0 += rd.pack(i);
					i += pw;
				}
			}
			else
			{
				a0 = rd.pack(0);
				i = pw;
			}

			r = sum(a0);
		}
		else
		{
			r = T(0);
			i = 0;
		}

		for (; i < len; ++i) r += rd.scalar(i);

		return r;
	}



} }

#endif 
