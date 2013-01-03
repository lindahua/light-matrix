/**
 * @file matrix_colviews_internal.h
 *
 * @brief Internal implementation of matrix colviews
 *
 * @author Dahua Lin
 */

#ifdef _MSC_VER
#pragma once
#endif

#ifndef LIGHTMAT_MATRIX_COLVIEWS_INTERNAL_H_
#define LIGHTMAT_MATRIX_COLVIEWS_INTERNAL_H_

#include <light_mat/matrix/matrix_meta.h>

namespace lmat {  namespace internal {

	template<class Mat, class Range, bool IsPerColCont, bool IsReadOnly> struct colview_helper;

	// whole

	template<class Mat>
	struct colview_helper<Mat, whole, true, true>
	{
		typedef typename matrix_traits<Mat>::value_type value_type;

		static const int ctrows = meta::nrows<Mat>::value;
		typedef cref_matrix<value_type, ctrows, 1> type;

		LMAT_ENSURE_INLINE
		static type get(const Mat& mat, const index_t j, whole)
		{
			return type(mat.ptr_col(j), mat.nrows(), 1);
		}
	};

	template<class Mat>
	struct colview_helper<Mat, whole, true, false>
	{
		typedef typename matrix_traits<Mat>::value_type value_type;

		static const int ctrows = meta::nrows<Mat>::value;
		typedef ref_matrix<value_type, ctrows, 1> type;

		LMAT_ENSURE_INLINE
		static type get(Mat& mat, const index_t j, whole)
		{
			return type(mat.ptr_col(j), mat.nrows(), 1);
		}
	};

	template<class Mat>
	struct colview_helper<Mat, whole, false, true>
	{
		typedef typename matrix_traits<Mat>::value_type value_type;

		static const int ctrows = meta::nrows<Mat>::value;
		typedef cref_grid<value_type, ctrows, 1> type;

		LMAT_ENSURE_INLINE
		static type get(const Mat& mat, const index_t j, whole)
		{
			return type(mat.ptr_col(j), mat.nrows(), 1, mat.row_stride(), 0);
		}
	};

	template<class Mat>
	struct colview_helper<Mat, whole, false, false>
	{
		typedef typename matrix_traits<Mat>::value_type value_type;

		static const int ctrows = meta::nrows<Mat>::value;
		typedef ref_grid<value_type, ctrows, 1> type;

		LMAT_ENSURE_INLINE
		static type get(Mat& mat, const index_t j, whole)
		{
			return type(mat.ptr_col(j), mat.nrows(), 1, mat.row_stride(), 0);
		}
	};


	// range

	template<class Mat>
	struct colview_helper<Mat, range, true, true>
	{
		typedef typename matrix_traits<Mat>::value_type value_type;

		typedef cref_matrix<value_type, 0, 1> type;

		LMAT_ENSURE_INLINE
		static type get(const Mat& mat, const index_t j, const range& rgn)
		{
			return type(mat.ptr_col(j) + rgn.begin_index(), rgn.num(), 1);
		}
	};

	template<class Mat>
	struct colview_helper<Mat, range, true, false>
	{
		typedef typename matrix_traits<Mat>::value_type value_type;

		typedef ref_matrix<value_type, 0, 1> type;

		LMAT_ENSURE_INLINE
		static type get(Mat& mat, const index_t j, const range& rgn)
		{
			return type(mat.ptr_col(j) + rgn.begin_index(), rgn.num(), 1);
		}
	};

	template<class Mat>
	struct colview_helper<Mat, range, false, true>
	{
		typedef typename matrix_traits<Mat>::value_type value_type;

		typedef cref_grid<value_type, 0, 1> type;

		LMAT_ENSURE_INLINE
		static type get(const Mat& mat, const index_t j, const range& rgn)
		{
			const index_t rs = mat.row_stride();
			return type(mat.ptr_col(j) + rgn.begin_index() * rs, rgn.num(), 1, rs, 0);
		}
	};

	template<class Mat>
	struct colview_helper<Mat, range, false, false>
	{
		typedef typename matrix_traits<Mat>::value_type value_type;

		typedef ref_grid<value_type, 0, 1> type;

		LMAT_ENSURE_INLINE
		static type get(Mat& mat, const index_t j, const range& rgn)
		{
			const index_t rs = mat.row_stride();
			return type(mat.ptr_col(j) + rgn.begin_index() * rs, rgn.num(), 1, rs, 0);
		}
	};


	// step_range

	template<class Mat>
	struct colview_helper<Mat, step_range, true, true>
	{
		typedef typename matrix_traits<Mat>::value_type value_type;

		typedef cref_grid<value_type, 0, 1> type;

		LMAT_ENSURE_INLINE
		static type get(const Mat& mat, const index_t j, const step_range& rgn)
		{
			return type(mat.ptr_col(j) + rgn.begin_index(), rgn.num(), 1, rgn.step(), 0);
		}
	};

	template<class Mat>
	struct colview_helper<Mat, step_range, true, false>
	{
		typedef typename matrix_traits<Mat>::value_type value_type;

		typedef ref_grid<value_type, 0, 1> type;

		LMAT_ENSURE_INLINE
		static type get(Mat& mat, const index_t j, const step_range& rgn)
		{
			return type(mat.ptr_col(j) + rgn.begin_index(), rgn.num(), 1, rgn.step(), 0);
		}
	};

	template<class Mat>
	struct colview_helper<Mat, step_range, false, true>
	{
		typedef typename matrix_traits<Mat>::value_type value_type;

		typedef cref_grid<value_type, 0, 1> type;

		LMAT_ENSURE_INLINE
		static type get(const Mat& mat, const index_t j, const step_range& rgn)
		{
			const index_t rs = mat.row_stride();
			return type(mat.ptr_col(j) + rgn.begin_index() * rs, rgn.num(), 1, rgn.step() * rs, 0);
		}
	};

	template<class Mat>
	struct colview_helper<Mat, step_range, false, false>
	{
		typedef typename matrix_traits<Mat>::value_type value_type;

		typedef ref_grid<value_type, 0, 1> type;

		LMAT_ENSURE_INLINE
		static type get(Mat& mat, const index_t j, const step_range& rgn)
		{
			const index_t rs = mat.row_stride();
			return type(mat.ptr_col(j) + rgn.begin_index() * rs, rgn.num(), 1, rgn.step() * rs, 0);
		}
	};


} }

#endif /* MATRIX_SUBVIEWS_INTERNAL_H_ */
