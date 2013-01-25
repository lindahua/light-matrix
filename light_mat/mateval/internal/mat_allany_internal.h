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

#include <light_mat/mateval/common_kernels.h>
#include <light_mat/mateval/macc_policy.h>
#include <light_mat/simd/simd.h>

namespace lmat { namespace internal {


	template<index_t N, typename T, class Reader>
	LMAT_ENSURE_INLINE
	inline bool all_impl(const dimension<N>& dim, type_<T>, scalar_, const Reader& rd, index_t i=0)
	{
		const index_t len = dim.value();

		for (; i < len; ++i)
		{
			if (!rd.scalar(i)) return false;
		}

		return true;
	}


	template<index_t N, typename T, typename SKind, class Reader>
	LMAT_ENSURE_INLINE
	inline bool all_impl(const dimension<N>& dim, type_<T>, simd_<SKind>, const Reader& rd)
	{
		typedef simd_bpack<T, SKind> pack_t;

		const index_t pw = (index_t)pack_t::pack_width;
		const index_t len = dim.value();
		const index_t npacks = (index_t)int_div<pack_t::pack_width>::quo((size_t)len);
		index_t i = 0;

		if (npacks)
		{
			const index_t l = npacks * pw;
			for (i = 0; i < l; i += pw)
			{
				if (any_false(rd.pack(i))) return false;
			}
		}

		return all_impl(dim, type_<T>(), scalar_(), rd, i);
	}


	template<index_t N, typename T, class Reader>
	LMAT_ENSURE_INLINE
	inline bool any_impl(const dimension<N>& dim, type_<T>, scalar_, const Reader& rd, index_t i=0)
	{
		const index_t len = dim.value();

		for (; i < len; ++i)
		{
			if (rd.scalar(i)) return true;
		}

		return false;
	}


	template<index_t N, typename T, typename SKind, class Reader>
	LMAT_ENSURE_INLINE
	inline bool any_impl(const dimension<N>& dim, type_<T>, simd_<SKind>, const Reader& rd)
	{
		typedef simd_bpack<T, SKind> pack_t;

		const index_t pw = (index_t)pack_t::pack_width;
		const index_t len = dim.value();
		const index_t npacks = (index_t)int_div<pack_t::pack_width>::quo((size_t)len);
		index_t i = 0;

		if (npacks)
		{
			const index_t l = npacks * pw;
			for (i = 0; i < l; i += pw)
			{
				if (any_true(rd.pack(i))) return true;
			}
		}

		return any_impl(dim, type_<T>(), scalar_(), rd, i);
	}


	template<index_t M, index_t N, typename T, typename VT, class Mat, typename U>
	inline bool all_(const matrix_shape<M, N>& shape, type_<T>, const IEWiseMatrix<Mat, VT>& mat, macc_<linear_, U>)
	{
		dimension<M * N> dim(shape.nelems());
		return all_impl(dim, type_<T>(), U(), make_vec_accessor(U(), in_(mat.derived())));
	}

	template<index_t M, index_t N, typename T, typename VT, class Mat, typename U>
	inline bool all_(const matrix_shape<M, N>& shape, type_<T>, const IEWiseMatrix<Mat, VT>& mat, macc_<percol_, U>)
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


	template<index_t M, index_t N, typename T, typename VT, class Mat, typename U>
	inline bool any_(const matrix_shape<M, N>& shape, type_<T>, const IEWiseMatrix<Mat, VT>& mat, macc_<linear_, U>)
	{
		dimension<M * N> dim(shape.nelems());
		return any_impl(dim, type_<T>(), U(), make_vec_accessor(U(), in_(mat.derived())));
	}

	template<index_t M, index_t N, typename T, typename VT, class Mat, typename U>
	inline bool any_(const matrix_shape<M, N>& shape, type_<T>, const IEWiseMatrix<Mat, VT>& mat, macc_<percol_, U>)
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


	template<index_t M, index_t N, typename T, class A, class D, typename U>
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

	template<index_t M, index_t N, typename T, class A, class D, typename U>
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
