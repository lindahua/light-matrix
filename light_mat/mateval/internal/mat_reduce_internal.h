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

	template<typename T>
	struct _sum_rfun
	{
		template<class Pack>
		LMAT_ENSURE_INLINE
		T reduce(const Pack& p) const
		{
			return math::sum(p);
		}

		template<class A>
		LMAT_ENSURE_INLINE
		void fold(A& a, const A& b) const
		{
			a += b;
		}
	};

	template<typename T>
	struct _maximum_rfun
	{
		template<class Pack>
		LMAT_ENSURE_INLINE
		T reduce(const Pack& p) const
		{
			return math::maximum(p);
		}

		template<class A>
		LMAT_ENSURE_INLINE
		void fold(A& a, const A& b) const
		{
			a = math::max(a, b);
		}
	};

	template<typename T>
	struct _minimum_rfun
	{
		template<class Pack>
		LMAT_ENSURE_INLINE
		T reduce(const Pack& p) const
		{
			return math::minimum(p);
		}

		template<class A>
		LMAT_ENSURE_INLINE
		void fold(A& a, const A& b) const
		{
			a = math::min(a, b);
		}
	};


	template<class RFun, typename RT, typename Kind, class Reader>
	inline void fold_impl(unsigned int len, RFun rfun, RT& r_out, Kind kind, const Reader& rd)
	{
		typedef math::simd_pack<RT, Kind> pack_t;

		const unsigned int pw = pack_t::pack_width;
		unsigned int npacks = int_div<pw>::quo(len);
		unsigned int i;
		RT r;

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
					rfun.fold(a0, rd.pack(i));
					rfun.fold(a1, rd.pack(i + pw));
					rfun.fold(a2, rd.pack(i + pw * 2));
					rfun.fold(a3, rd.pack(i + pw * 3));
				}

				rfun.fold(a0, a2);
				rfun.fold(a1, a3);

				if (npacks & 2)
				{
					rfun.fold(a0, rd.pack(i));
					rfun.fold(a1, rd.pack(i + pw));
					i += pw * 2;
				}

				rfun.fold(a0, a1);

				if (npacks & 1)
				{
					rfun.fold(a0, rd.pack(i));
					i += pw;
				}
			}
			else if (npacks >> 1)
			{
				a0 = rd.pack(0);
				rfun.fold(a0, rd.pack(pw));
				i = pw << 1;

				if (npacks & 1)
				{
					rfun.fold(a0, rd.pack(i));
					i += pw;
				}
			}
			else
			{
				rfun.fold(a0, rd.pack(0));
				i = pw;
			}

			r = rfun.reduce(a0);
		}
		else
		{
			r = rd.scalar(0);
			i = 1;
		}

		for (; i < len; ++i) rfun.fold(r, rd.scalar(i));
		r_out = r;
	}


	template<class RFun, typename RT, typename Kind, typename TFun, typename... Reader>
	inline void foldx_impl(unsigned int len, RFun rfun, RT& r_out, Kind, TFun tfun, const Reader&... rds)
	{
		typedef math::simd_pack<RT, Kind> pack_t;

		const unsigned int pw = pack_t::pack_width;
		unsigned int npacks = int_div<pw>::quo(len);
		unsigned int i;
		RT r;

		if (npacks)
		{
			pack_t a0, a1, a2, a3;
			pack_t t0, t1, t2, t3;

			unsigned int m4 = npacks >> 2;
			if (m4)
			{
				const unsigned int w4 = pw << 2;
				unsigned int l4 = m4 * w4;

				a0 = tfun(rds.pack(0)...);
				a1 = tfun(rds.pack(pw)...);
				a2 = tfun(rds.pack(pw * 2)...);
				a3 = tfun(rds.pack(pw * 3)...);

				for (i = w4; i < l4; i += w4)
				{
					t0 = tfun(rds.pack(i)...);
					t1 = tfun(rds.pack(i + pw)...);
					t2 = tfun(rds.pack(i + pw * 2)...);
					t3 = tfun(rds.pack(i + pw * 3)...);

					rfun.fold(a0, t0);
					rfun.fold(a1, t1);
					rfun.fold(a2, t2);
					rfun.fold(a3, t3);
				}

				rfun.fold(a0, a2);
				rfun.fold(a1, a3);

				if (npacks & 2)
				{
					t0 = tfun(rds.pack(i)...);
					t1 = tfun(rds.pack(i + pw)...);

					rfun.fold(a0, t0);
					rfun.fold(a1, t1);
					i += pw * 2;
				}

				rfun.fold(a0, a1);

				if (npacks & 1)
				{
					rfun.fold(a0, tfun(rds.pack(i)...));
					i += pw;
				}
			}
			else if (npacks >> 1)
			{
				a0 = tfun(rds.pack(0)...);
				rfun.fold(a0, tfun(rds.pack(pw)...));
				i = pw << 1;

				if (npacks & 1)
				{
					rfun.fold(a0, tfun(rds.pack(i)...));
					i += pw;
				}
			}
			else
			{
				rfun.fold(a0, tfun(rds.pack(0)...));
				i = pw;
			}

			r = rfun.reduce(a0);
		}
		else
		{
			r = tfun(rds.scalar(0)...);
			i = 1;
		}

		for (; i < len; ++i) rfun.fold(r, tfun(rds.scalar(i)...));
		r_out = r;
	}


	// Basic reduction function

	template<typename Kind, class Reader>
	inline typename Reader::scalar_type
	sum_impl(unsigned int len, Kind, const Reader& rd)
	{
		typedef typename Reader::scalar_type T;

		T r;
		fold_impl(len, _sum_rfun<T>(), r, Kind(), rd);
		return r;
	}

	template<typename Kind, class Reader>
	inline typename Reader::scalar_type
	maximum_impl(unsigned int len, Kind, const Reader& rd)
	{
		typedef typename Reader::scalar_type T;

		T r;
		fold_impl(len, _maximum_rfun<T>(), r, Kind(), rd);
		return r;
	}

	template<typename Kind, class Reader>
	inline typename Reader::scalar_type
	minimum_impl(unsigned int len, Kind, const Reader& rd)
	{
		typedef typename Reader::scalar_type T;

		T r;
		fold_impl(len, _minimum_rfun<T>(), r, Kind(), rd);
		return r;
	}

	// Extended reduction function

	template<typename TFun, typename RT, typename Kind, typename... Reader>
	inline void sum_impl_x(unsigned int len, RT& r, Kind, const TFun& tfun, const Reader&... rds)
	{
		foldx_impl(len, _sum_rfun<RT>(), r, Kind(), tfun, rds...);
	}

	template<typename TFun, typename RT, typename Kind, typename... Reader>
	inline void maximum_impl_x(unsigned int len, RT& r, Kind, const TFun& tfun, const Reader&... rds)
	{
		foldx_impl(len, _maximum_rfun<RT>(), r, Kind(), tfun, rds...);
	}

	template<typename TFun, typename RT, typename Kind, typename... Reader>
	inline void minimum_impl_x(unsigned int len, RT& r, Kind, const TFun& tfun, const Reader&... rds)
	{
		foldx_impl(len, _minimum_rfun<RT>(), r, Kind(), tfun, rds...);
	}


} }

#endif 
