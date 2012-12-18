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

#include <light_mat/mateval/multicol_accessors.h>
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


	/********************************************
	 *
	 *  core implementation
	 *
	 ********************************************/

	template<int N, typename T, typename SKind, class RFun, class Reader>
	inline T fold_impl(const dimension<N>& dim, atags::simd<T, SKind>,
			RFun rfun, const Reader& rd)
	{
		typedef math::simd_pack<T, SKind> pack_t;

		const index_t pw = (index_t)pack_t::pack_width;

		const index_t len = dim.value();
		const index_t npacks = int_div<pack_t::pack_width>::quo(len);
		index_t i;
		T r;

		if (npacks)
		{
			pack_t a0, a1, a2, a3;

			const index_t m4 = npacks >> 2;
			if (m4)
			{
				const index_t w4 = pw << 2;
				const index_t l4 = m4 * w4;

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
				a0 = rd.pack(0);
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
		return r;
	}


	template<int N, typename T, typename SKind, class RFun, typename TFun, typename... Reader>
	inline T foldx_impl(const dimension<N>& dim, atags::simd<T, SKind>,
			RFun rfun, const TFun& tfun, const Reader&... rds)
	{
		typedef math::simd_pack<T, SKind> pack_t;
		const index_t pw = (index_t)pack_t::pack_width;

		const index_t len = dim.value();
		const index_t npacks = int_div<pack_t::pack_width>::quo(len);
		index_t i;
		T r;

		if (npacks)
		{
			pack_t a0, a1, a2, a3;
			pack_t t0, t1, t2, t3;

			const index_t m4 = npacks >> 2;
			if (m4)
			{
				const index_t w4 = pw << 2;
				const index_t l4 = m4 * w4;

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
				a0 = tfun(rds.pack(0)...);
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
		return r;
	}



	/********************************************
	 *
	 *  full reduction
	 *
	 ********************************************/

	template<int N, typename T, typename Kind, class Wrap>
	inline T sum_(const dimension<N>& dim, atags::simd<T, Kind> u, const Wrap& wrap)
	{
		return fold_impl(dim, u, _sum_rfun<T>(), make_vec_accessor(u, wrap));
	}

	template<int N, typename T, typename Kind, class Wrap>
	inline T mean_(const dimension<N>& dim, atags::simd<T, Kind> u, const Wrap& wrap)
	{
		T r = sum_(dim, u, wrap);
		return r / T(dim.value());
	}

	template<int N, typename T, typename Kind, class Wrap>
	inline T maximum_(const dimension<N>& dim, atags::simd<T, Kind> u, const Wrap& wrap)
	{
		return fold_impl(dim, u, _maximum_rfun<T>(), make_vec_accessor(u, wrap));
	}

	template<int N, typename T, typename Kind, class Wrap>
	inline T minimum_(const dimension<N>& dim, atags::simd<T, Kind> u, const Wrap& wrap)
	{
		return fold_impl(dim, u, _minimum_rfun<T>(), make_vec_accessor(u, wrap));
	}


	template<int N, typename T, typename Kind, typename TFun, typename... Wrap>
	inline T sumx_(const dimension<N>& dim, atags::simd<T, Kind> u, const TFun& tfun, const Wrap&... wraps)
	{
		return foldx_impl(dim, u, _sum_rfun<T>(), tfun, make_vec_accessor(u, wraps)...);
	}

	template<int N, typename T, typename Kind, typename TFun, typename... Wrap>
	inline T meanx_(const dimension<N>& dim, atags::simd<T, Kind> u, const TFun& tfun, const Wrap&... wraps)
	{
		T r = sumx_(dim, u, tfun, wraps...);
		return r / T(dim.value());
	}

	template<int N, typename T, typename Kind, typename TFun, typename... Wrap>
	inline T maximumx_(const dimension<N>& dim, atags::simd<T, Kind> u, const TFun& tfun, const Wrap&... wraps)
	{
		return foldx_impl(dim, u, _maximum_rfun<T>(), tfun, make_vec_accessor(u, wraps)...);
	}

	template<int N, typename T, typename Kind, typename TFun, typename... Wrap>
	inline T minimumx_(const dimension<N>& dim, atags::simd<T, Kind> u, const TFun& tfun, const Wrap&... wraps)
	{
		return foldx_impl(dim, u, _minimum_rfun<T>(), tfun, make_vec_accessor(u, wraps)...);
	}

	/********************************************
	 *
	 *  colwise reduction
	 *
	 ********************************************/

	template<int M, int N, typename U, class RFun, class DMat, class MultiColReader>
	inline void colwise_fold_impl(const matrix_shape<M, N>& shape, U u,
			RFun rfun, DMat& dmat, const MultiColReader& rd)
	{
		dimension<M> col_dim(shape.nrows());
		const index_t n = shape.ncolumns();

		for (index_t j = 0; j < n; ++j)
		{
			dmat[j] = fold_impl(col_dim, u, rfun, rd.col(j));
		}
	}


	template<int M, int N, typename T, typename Kind, class DMat, class Wrap>
	inline void colwise_sum_(const matrix_shape<M, N>& shape, atags::simd<T, Kind> u, DMat& dmat, const Wrap& wrap)
	{
		colwise_fold_impl(shape, u, _sum_rfun<T>(), dmat, make_multicol_accessor(u, wrap));
	}


} }

#endif 
