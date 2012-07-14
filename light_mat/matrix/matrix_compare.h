/**
 * @file matrix_compare.h
 *
 * Functions for matrix comparison
 *
 * @author Dahua Lin
 */

#ifdef _MSC_VER
#pragma once
#endif

#ifndef LIGHTMAT_MATRIX_COMPARE_H_
#define LIGHTMAT_MATRIX_COMPARE_H_

#include <light_mat/matrix/matrix_properties.h>
#include "bits/matrix_compare_internal.h"
#include <cmath>

namespace lmat
{
	template<typename T, class LMat, class RMat>
	LMAT_ENSURE_INLINE
	bool is_equal(const IDenseMatrix<LMat, T>& a, const IDenseMatrix<RMat, T>& b)
	{
		const int M = binary_ct_rows<LMat, RMat>::value;
		const int N = binary_ct_cols<LMat, RMat>::value;
		typedef typename detail::mat_comparer<T, M, N>::type comparer_t;

		return has_same_size(a, b) &&
				comparer_t::is_equal(a.nrows(), a.ncolumns(),
				a.ptr_data(), a.lead_dim(),
				b.ptr_data(), b.lead_dim());
	}

	template<typename T>
	LMAT_ENSURE_INLINE
	bool is_approx_scalar(const T& a, const T& b, const T& tol)
	{
		return std::fabs(a - b) <= tol;
	}

	template<typename T, class LMat, class RMat>
	inline
	bool is_approx(const IDenseMatrix<LMat, T>& a, const IDenseMatrix<RMat, T>& b, const T& tol)
	{
		if (has_same_size(a, b))
		{
			const index_t m = a.nrows();
			const index_t n = a.ncolumns();

			if (n == 1)
			{
				const T *ca = a.ptr_data();
				const T *cb = b.ptr_data();

				for (index_t i = 0; i < m; ++i)
				{
					if (!is_approx_scalar(ca[i], cb[i], tol))
						return false;
				}
			}
			else
			{
				for (index_t j = 0; j < n; ++j)
				{
					const T *ca = a.ptr_col(j);
					const T *cb = b.ptr_col(j);

					for (index_t i = 0; i < m; ++i)
					{
						if (!is_approx_scalar(ca[i], cb[i], tol))
							return false;
					}
				}
			}

			return true;
		}
		else return false;
	}

}

#endif
