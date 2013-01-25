/*
 * @file matrix_copy.h
 *
 * Functions to copy or cast matrices
 *
 * @author Dahua Lin
 */

#ifdef _MSC_VER
#pragma once
#endif

#ifndef LIGHTMAT_MATRIX_COPY_H_
#define LIGHTMAT_MATRIX_COPY_H_

#include "internal/matrix_copy_internal.h"

namespace lmat
{
	template<typename T, class RMat>
	LMAT_ENSURE_INLINE
	inline void copy(const T *ps, IRegularMatrix<RMat, T>& dst)
	{
		internal::copy(ps, dst.derived(), internal::get_copy_scheme(dst));
	}

	template<typename T, class LMat>
	LMAT_ENSURE_INLINE
	inline void copy(const IRegularMatrix<LMat, T>& src, T* pd)
	{
		internal::copy(src.derived(), pd, internal::get_copy_scheme(src));
	}


	template<typename T, class LMat, class RMat>
	LMAT_ENSURE_INLINE
	inline void copy(const IRegularMatrix<LMat, T>& src, IRegularMatrix<RMat, T>& dst)
	{
		internal::copy(src.derived(), dst.derived(), internal::get_copy_scheme(src, dst));
	}

	template<typename T, class DMat>
	LMAT_ENSURE_INLINE
	inline DMat& operator << (IRegularMatrix<DMat, T>& dmat, const T *src)
	{
		copy(src, dmat.derived());
		return dmat.derived();
	}

	template<typename T, class SExpr, class DMat>
	LMAT_ENSURE_INLINE
	inline void evaluate(const IRegularMatrix<SExpr, T>& sexpr, IRegularMatrix<DMat, T>& dmat)
	{
		copy(sexpr.derived(), dmat.derived());
	}


	/********************************************
	 *
	 *  copy part of the matrix
	 *
	 ********************************************/


	// copy all elements with j - i >= k

	template<typename T, class SMat, class DMat>
	inline void copy_triu(const IRegularMatrix<SMat, T>& smat, IRegularMatrix<DMat, T>& dmat, index_t k=0)
	{
		static_assert(meta::is_percol_contiguous<SMat>::value, "smat should be percol contiguous");
		static_assert(meta::is_percol_contiguous<DMat>::value, "dmat should be percol contiguous");

		index_t m = smat.nrows();
		index_t n = smat.ncolumns();

		index_t dm = dmat.nrows();
		index_t dn = dmat.ncolumns();

		LMAT_CHECK_DIMS( (dm >= m || dm >= n-k) && dn == n )

		const T *ps = smat.ptr_data();
		T *pd = dmat.ptr_data();

		index_t ss = smat.col_stride();
		index_t ds = dmat.col_stride();

		index_t j0 = k >= 0 ? k : 0;
		if (j0 > 0)
		{
			ps += j0 * ss;
			pd += j0 * ds;
		}

		if (n <= m + k)
		{
			for (index_t j = j0; j < n; ++j, ps += ss, pd += ds)
				copy_vec(j-k+1, ps, pd);
		}
		else
		{
			for (index_t j = j0; j < m+k; ++j, ps += ss, pd += ds)
				copy_vec(j-k+1, ps, pd);

			for (index_t j = m+k; j < n; ++j, ps += ss, pd += ds)
				copy_vec(m, ps, pd);
		}
	}


	// copy all elements with j - i <= k

	template<typename T, class SMat, class DMat>
	inline void copy_tril(const IRegularMatrix<SMat, T>& smat, IRegularMatrix<DMat, T>& dmat, index_t k=0)
	{
		static_assert(meta::is_percol_contiguous<SMat>::value, "smat should be percol contiguous");
		static_assert(meta::is_percol_contiguous<DMat>::value, "dmat should be percol contiguous");

		index_t m = smat.nrows();
		index_t n = smat.ncolumns();

		index_t dm = dmat.nrows();
		index_t dn = dmat.ncolumns();

		LMAT_CHECK_DIMS( dm == m && (dn >= n || dn >= m + k) )

		const T *ps = smat.ptr_data();
		T *pd = dmat.ptr_data();

		index_t ss = smat.col_stride();
		index_t ds = dmat.col_stride();

		for (index_t j = 0; j <= k; ++j, ps += ss, pd += ds)
			copy_vec(m, ps, pd);

		index_t j0 = k >= -1 ? k + 1 : 0;
		for (index_t j = j0; j <= m+k; ++j, ps += ss, pd += ds)
		{
			index_t di = j-k;
			copy_vec(m-di, ps+di, pd+di);
		}
	}


}

#endif 
