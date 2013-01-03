/**
 * @file mat_allany_internal.h
 *
 * @brief Internal implementation of all/any reduction
 *
 * @author Dahua Lin
 */

#ifdef _MSC_VER
#pragma once
#endif

#ifndef LIGHTMAT_MAT_ALLANY_INTERNAL_H_
#define LIGHTMAT_MAT_ALLANY_INTERNAL_H_

#include <light_mat/mateval/macc_policy.h>
#include <light_mat/math/sse_reduce.h>
#ifdef LMAT_HAS_AVX
#include <light_mat/math/avx_reduce.h>
#endif

namespace lmat { namespace internal {


	template<int N, typename T, class Reader>
	LMAT_ENSURE_INLINE
	inline bool all_impl(const dimension<N>& dim, type_<T>, atags::scalar, const Reader& rd, index_t i=0)
	{
		const index_t len = dim.value();

		for (; i < len; ++i)
		{
			if (!rd.scalar(i)) return false;
		}

		return true;
	}


	template<int N, typename T, typename SKind, class Reader>
	LMAT_ENSURE_INLINE
	inline bool all_impl(const dimension<N>& dim, type_<T>, atags::simd<SKind>, const Reader& rd)
	{
		typedef typename math::simd_bpack<T, SKind> pack_t;

		const index_t pw = (index_t)pack_t::pack_width;
		const index_t len = dim.value();
		const index_t npacks = (index_t)int_div<pack_t::pack_width>::quo((size_t)len);
		index_t i = 0;

		if (npacks)
		{
			const index_t l = npacks * pw;
			for (i = 0; i < l; i += pw)
			{
				if (math::any_false(rd.pack(i))) return false;
			}
		}

		return all_impl(dim, type_<T>(), atags::scalar(), rd, i);
	}


	template<int N, typename T, class Reader>
	LMAT_ENSURE_INLINE
	inline bool any_impl(const dimension<N>& dim, type_<T>, atags::scalar, const Reader& rd, index_t i=0)
	{
		const index_t len = dim.value();

		for (; i < len; ++i)
		{
			if (rd.scalar(i)) return true;
		}

		return false;
	}


	template<int N, typename T, typename SKind, class Reader>
	LMAT_ENSURE_INLINE
	inline bool any_impl(const dimension<N>& dim, type_<T>, atags::simd<SKind>, const Reader& rd)
	{
		typedef typename math::simd_bpack<T, SKind> pack_t;

		const index_t pw = (index_t)pack_t::pack_width;
		const index_t len = dim.value();
		const index_t npacks = (index_t)int_div<pack_t::pack_width>::quo((size_t)len);
		index_t i = 0;

		if (npacks)
		{
			const index_t l = npacks * pw;
			for (i = 0; i < l; i += pw)
			{
				if (math::any_true(rd.pack(i))) return true;
			}
		}

		return any_impl(dim, type_<T>(), atags::scalar(), rd, i);
	}


	template<int M, int N, typename T, typename VT, class Mat, typename U>
	inline bool all_(const matrix_shape<M, N>& shape, type_<T>, const IEWiseMatrix<Mat, VT>& mat, linear_macc<U>)
	{
		dimension<M * N> dim(shape.nelems());
		return all_impl(dim, type_<T>(), U(), make_vec_accessor(U(), in_(mat.derived())));
	}

	template<int M, int N, typename T, typename VT, class Mat, typename U>
	inline bool all_(const matrix_shape<M, N>& shape, type_<T>, const IEWiseMatrix<Mat, VT>& mat, percol_macc<U>)
	{
		if (shape.nelems() > 0)
		{
			dimension<M> col_dim(shape.nrows());
			auto rd = make_multicol_accessor(U(), in_(mat.derived()));

			const index_t n = shape.ncolumns();
			for (index_t j = 0; j < n; ++j)
			{
				if (!all_impl(col_dim, type_<T>(), U(), rd.col(j)))
					return false;
			}
		}

		return true;
	}


	template<int M, int N, typename T, typename VT, class Mat, typename U>
	inline bool any_(const matrix_shape<M, N>& shape, type_<T>, const IEWiseMatrix<Mat, VT>& mat, linear_macc<U>)
	{
		dimension<M * N> dim(shape.nelems());
		return any_impl(dim, type_<T>(), U(), make_vec_accessor(U(), in_(mat.derived())));
	}

	template<int M, int N, typename T, typename VT, class Mat, typename U>
	inline bool any_(const matrix_shape<M, N>& shape, type_<T>, const IEWiseMatrix<Mat, VT>& mat, percol_macc<U>)
	{
		if (shape.nelems() > 0)
		{
			dimension<M> col_dim(shape.nrows());
			auto rd = make_multicol_accessor(U(), in_(mat.derived()));

			const index_t n = shape.ncolumns();
			for (index_t j = 0; j < n; ++j)
			{
				if (any_impl(col_dim, type_<T>(), U(), rd.col(j)))
					return true;
			}
		}

		return false;
	}


	template<int M, int N, typename T, class A, class D, typename U>
	inline void colwise_all_(const matrix_shape<M, N>& shape, type_<T>, const A& a, D& d, bool expect_val, U)
	{
		const index_t n = a.ncolumns();
		auto rd = make_multicol_accessor(U(), in_(a));
		dimension<M> col_dim(shape.nrows());

		if (expect_val)
		{
			for (index_t j = 0; j < n; ++j)
			{
				d[j] = all_impl(col_dim, type_<T>(), U(), rd.col(j));
			}
		}
		else
		{
			for (index_t j = 0; j < n; ++j)
			{
				d[j] = !any_impl(col_dim, type_<T>(), U(), rd.col(j));
			}
		}
	}

	template<int M, int N, typename T, class A, class D, typename U>
	inline void colwise_any_(const matrix_shape<M, N>& shape, type_<T>, const A& a, D& d, bool expect_val, U)
	{
		const index_t n = a.ncolumns();
		auto rd = make_multicol_accessor(U(), in_(a));
		dimension<M> col_dim(shape.nrows());

		if (expect_val)
		{
			for (index_t j = 0; j < n; ++j)
			{
				d[j] = any_impl(col_dim, type_<T>(), U(), rd.col(j));
			}
		}
		else
		{
			for (index_t j = 0; j < n; ++j)
			{
				d[j] = !all_impl(col_dim, type_<T>(), U(), rd.col(j));
			}
		}
	}


} }

#endif /* MAT_ALLANY_INTERNAL_H_ */
