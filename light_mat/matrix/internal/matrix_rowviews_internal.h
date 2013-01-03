/**
 * @file matrix_rowviews_internal.h
 *
 * @brief Internal implementation of matrix rowviews
 *
 * @author Dahua Lin
 */

#ifdef _MSC_VER
#pragma once
#endif

#ifndef LIGHTMAT_MATRIX_ROWVIEWS_INTERNAL_H_
#define LIGHTMAT_MATRIX_ROWVIEWS_INTERNAL_H_

#include <light_mat/matrix/matrix_meta.h>

namespace lmat {  namespace internal {

	template<class Mat, class Range, bool IsPerRowCont, bool IsReadOnly> struct rowview_helper;

	// whole

	template<class Mat>
	struct rowview_helper<Mat, whole, true, true>
	{
		typedef typename matrix_traits<Mat>::value_type value_type;

		static const int ctcols = meta::ncols<Mat>::value;
		typedef cref_matrix<value_type, 1, ctcols> type;

		LMAT_ENSURE_INLINE
		static type get(const Mat& mat, const index_t i, whole)
		{
			return type(mat.ptr_row(i), 1, mat.ncolumns());
		}
	};

	template<class Mat>
	struct rowview_helper<Mat, whole, true, false>
	{
		typedef typename matrix_traits<Mat>::value_type value_type;

		static const int ctcols = meta::ncols<Mat>::value;
		typedef ref_matrix<value_type, 1, ctcols> type;

		LMAT_ENSURE_INLINE
		static type get(Mat& mat, const index_t i, whole)
		{
			return type(mat.ptr_row(i), 1, mat.ncolumns());
		}
	};

	template<class Mat>
	struct rowview_helper<Mat, whole, false, true>
	{
		typedef typename matrix_traits<Mat>::value_type value_type;

		static const int ctcols = meta::ncols<Mat>::value;
		typedef cref_block<value_type, 1, ctcols> type;

		LMAT_ENSURE_INLINE
		static type get(const Mat& mat, const index_t i, whole)
		{
			return type(mat.ptr_row(i), 1, mat.ncolumns(), mat.col_stride());
		}
	};

	template<class Mat>
	struct rowview_helper<Mat, whole, false, false>
	{
		typedef typename matrix_traits<Mat>::value_type value_type;

		static const int ctcols = meta::ncols<Mat>::value;
		typedef ref_block<value_type, 1, ctcols> type;

		LMAT_ENSURE_INLINE
		static type get(Mat& mat, const index_t i, whole)
		{
			return type(mat.ptr_row(i), 1, mat.ncolumns(), mat.col_stride());
		}
	};


	// range

	template<class Mat>
	struct rowview_helper<Mat, range, true, true>
	{
		typedef typename matrix_traits<Mat>::value_type value_type;

		static const int ctcols = meta::ncols<Mat>::value;
		typedef cref_matrix<value_type, 1, 0> type;

		LMAT_ENSURE_INLINE
		static type get(const Mat& mat, const index_t i, const range& rgn)
		{
			return type(mat.ptr_row(i) + rgn.begin_index(), 1, rgn.num());
		}
	};

	template<class Mat>
	struct rowview_helper<Mat, range, true, false>
	{
		typedef typename matrix_traits<Mat>::value_type value_type;

		static const int ctcols = meta::ncols<Mat>::value;
		typedef ref_matrix<value_type, 1, 0> type;

		LMAT_ENSURE_INLINE
		static type get(Mat& mat, const index_t i, const range& rgn)
		{
			return type(mat.ptr_row(i) + rgn.begin_index(), 1, rgn.num());
		}
	};

	template<class Mat>
	struct rowview_helper<Mat, range, false, true>
	{
		typedef typename matrix_traits<Mat>::value_type value_type;

		static const int ctcols = meta::ncols<Mat>::value;
		typedef cref_block<value_type, 1, 0> type;

		LMAT_ENSURE_INLINE
		static type get(const Mat& mat, const index_t i, const range& rgn)
		{
			const index_t cs = mat.col_stride();
			return type(mat.ptr_row(i) + rgn.begin_index() * cs, 1, rgn.num(), cs);
		}
	};

	template<class Mat>
	struct rowview_helper<Mat, range, false, false>
	{
		typedef typename matrix_traits<Mat>::value_type value_type;

		static const int ctcols = meta::ncols<Mat>::value;
		typedef ref_block<value_type, 1, 0> type;

		LMAT_ENSURE_INLINE
		static type get(Mat& mat, const index_t i, const range& rgn)
		{
			const index_t cs = mat.col_stride();
			return type(mat.ptr_row(i) + rgn.begin_index() * cs, 1, rgn.num(), cs);
		}
	};


	// step_range

	template<class Mat>
	struct rowview_helper<Mat, step_range, true, true>
	{
		typedef typename matrix_traits<Mat>::value_type value_type;

		static const int ctcols = meta::ncols<Mat>::value;
		typedef cref_block<value_type, 1, 0> type;

		LMAT_ENSURE_INLINE
		static type get(const Mat& mat, const index_t i, const step_range& rgn)
		{
			return type(mat.ptr_row(i) + rgn.begin_index(), 1, rgn.num(), rgn.step());
		}
	};

	template<class Mat>
	struct rowview_helper<Mat, step_range, true, false>
	{
		typedef typename matrix_traits<Mat>::value_type value_type;

		static const int ctcols = meta::ncols<Mat>::value;
		typedef ref_block<value_type, 1, 0> type;

		LMAT_ENSURE_INLINE
		static type get(Mat& mat, const index_t i, const step_range& rgn)
		{
			return type(mat.ptr_row(i) + rgn.begin_index(), 1, rgn.num(), rgn.step());
		}
	};

	template<class Mat>
	struct rowview_helper<Mat, step_range, false, true>
	{
		typedef typename matrix_traits<Mat>::value_type value_type;

		static const int ctcols = meta::ncols<Mat>::value;
		typedef cref_block<value_type, 1, 0> type;

		LMAT_ENSURE_INLINE
		static type get(const Mat& mat, const index_t i, const step_range& rgn)
		{
			const index_t cs = mat.col_stride();
			return type(mat.ptr_row(i) + rgn.begin_index() * cs, 1, rgn.num(), rgn.step() * cs);
		}
	};

	template<class Mat>
	struct rowview_helper<Mat, step_range, false, false>
	{
		typedef typename matrix_traits<Mat>::value_type value_type;

		static const int ctcols = meta::ncols<Mat>::value;
		typedef ref_block<value_type, 1, 0> type;

		LMAT_ENSURE_INLINE
		static type get(Mat& mat, const index_t i, const step_range& rgn)
		{
			const index_t cs = mat.col_stride();
			return type(mat.ptr_row(i) + rgn.begin_index() * cs, 1, rgn.num(), rgn.step() * cs);
		}
	};



} }


#endif /* MATRIX_ROWVIEWS_INTERNAL_H_ */
