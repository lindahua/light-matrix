/**
 * @file matrix_iter.h
 *
 * @brief Matrx iterators
 *
 * @author Dahua Lin
 */

#ifdef _MSC_VER
#pragma once
#endif

#ifndef LIGHTMAT_MATRIX_ITER_H_
#define LIGHTMAT_MATRIX_ITER_H_

#include "internal/matrix_iter_internal.h"

namespace lmat
{
	template<class Mat>
	struct matrix_iter
	{
		typedef internal::matrix_iter_helper<Mat,
				typename internal::iter_tag<Mat>::type> iter_helper;

		typedef internal::matrix_coliter_helper<Mat,
				typename internal::coliter_tag<Mat>::type> coliter_helper;

		typedef typename iter_helper::const_iterator const_iterator;
		typedef typename iter_helper::iterator iterator;

		typedef typename coliter_helper::const_iterator col_const_iterator;
		typedef typename coliter_helper::iterator col_iterator;

		LMAT_ENSURE_INLINE
		static const_iterator begin(const Mat& a)
		{
			return iter_helper::begin(a);
		}

		LMAT_ENSURE_INLINE
		static const_iterator end(const Mat& a)
		{
			return iter_helper::end(a);
		}

		LMAT_ENSURE_INLINE
		static iterator begin(Mat& a)
		{
			return iter_helper::begin(a);
		}

		LMAT_ENSURE_INLINE
		static iterator end(Mat& a)
		{
			return iter_helper::end(a);
		}

		LMAT_ENSURE_INLINE
		static col_const_iterator col_begin(const Mat& a, index_t j)
		{
			return coliter_helper::begin(a, j);
		}

		LMAT_ENSURE_INLINE
		static col_const_iterator col_end(const Mat& a, index_t j)
		{
			return coliter_helper::end(a, j);
		}

		LMAT_ENSURE_INLINE
		static col_iterator col_begin(Mat& a, index_t j)
		{
			return coliter_helper::begin(a, j);
		}

		LMAT_ENSURE_INLINE
		static col_iterator col_end(Mat& a, index_t j)
		{
			return coliter_helper::end(a, j);
		}

	};

}

#endif
