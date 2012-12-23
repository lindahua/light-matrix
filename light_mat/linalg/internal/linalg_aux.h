/**
 * @file linalg_aux.h
 *
 * Auxiliary function to support linear algebra routines
 * 
 * @author Dahua Lin 
 */

#ifdef _MSC_VER
#pragma once
#endif

#ifndef LIGHTMAT_LINALG_AUX_H_
#define LIGHTMAT_LINALG_AUX_H_

#include <light_mat/linalg/linalg_fwd.h>

namespace lmat { namespace internal {

	template<typename T, class Mat>
	LMAT_ENSURE_INLINE
	inline index_t get_vector_intv(const IRegularMatrix<Mat, T>& mat)
	{
		index_t intv;

		if (meta::is_continuous<Mat>::value)
		{
			intv = 1;
		}
		else if (is_column(mat))
		{
			intv = mat.row_stride();
		}
		else if (is_row(mat))
		{
			intv = mat.col_stride();
		}
		else
		{
			throw invalid_argument("The input is not a proper vector.");
		}

		return intv;
	}


	template<typename T, class Mat>
	LMAT_ENSURE_INLINE
	inline index_t scale_columns(index_t m, index_t n, index_t lda, double *A, double *s)
	{
		for (index_t j = 0; j < n; ++j)
		{
			double *a = A + j * lda;
			for (index_t i = 0; i < m; ++i) a[i] *= s[i];
		}
	}


	template<typename T, class SMat, class DMat>
	inline void get_triu(index_t n, const IRegularMatrix<SMat, T>& a, IRegularMatrix<DMat, T>& b, bool uni=false)
	{
		const T *pa = a.ptr_data();
		T *pb = b.ptr_data();
		index_t as = a.col_stride();
		index_t bs = b.col_stride();

		if (uni)
		{
			pb[0] = T(1);
			pa += as;
			pb += bs;

			for (index_t j = 1; j < n; ++j, pa += as, pb += bs)
			{
				copy_vec(j, pa, pb);
				pb[j] = T(1);
			}
		}
		else
		{
			pb[0] = pa[0];
			pa += as;
			pb += bs;

			for (index_t j = 1; j < n; ++j, pa += as, pb += bs)
			{
				copy_vec(j+1, pa, pb);
			}
		}
	}

	template<typename T, class SMat, class DMat>
	inline void get_tril(index_t n, const IRegularMatrix<SMat, T>& a, IRegularMatrix<DMat, T>& b, bool uni=false)
	{
		const T *pa = a.ptr_data();
		T *pb = b.ptr_data();
		index_t as = a.col_stride();
		index_t bs = b.col_stride();

		if (uni)
		{
			for (index_t j = 0; j < n-1; ++j, pa += as, pb += bs)
			{
				index_t l = j+1;

				pb[j] = T(1);
				copy_vec(n-l, pa+l, pb+l);
			}

			pb[n-1] = T(1);
		}
		else
		{
			for (index_t j = 0; j < n-1; ++j, pa += as, pb += bs)
			{
				copy_vec(n-j, pa+j, pb+j);
			}

			pb[n-1] = pa[n-1];
		}
	}


	template<typename T, class Mat>
	inline void complete_sym(index_t n, IRegularMatrix<Mat, T>& a, char uplo)
	{
		T *pa = a.ptr_data();
		index_t as = a.col_stride();

		if (uplo == 'L' || uplo == 'l')
		{
			pa += as;
			for (index_t j = 1; j < n; ++j, pa += as)
			{
				for (index_t i = 0; i < j; ++i) pa[i] = a.elem(j, i);
			}
		}
		else
		{
			for (index_t j = 0; j < n-1; ++j, pa += as)
			{
				for (index_t i = j+1; i < n; ++i) pa[i] = a.elem(j, i);
			}
		}
	}

} }

#endif 
