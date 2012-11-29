/**
 * @file matrix_matviews_internal.h
 *
 * @brief Internal implementation of matrix mat-views
 *
 * @author Dahua Lin
 */

#ifndef LIGHTMAT_MATRIX_MATVIEWS_INTERNAL_H_
#define LIGHTMAT_MATRIX_MATVIEWS_INTERNAL_H_

#include <light_mat/matrix/matrix_meta.h>

namespace lmat {  namespace internal {

	// ContLevel
	// 0 : non continuous at each column
	// 1 : continuous per-column
	// 2 : continuous as whole

	template<class Mat, class ColRange, class RowRange, int ContLevel, bool IsReadOnly> struct matview_helper;

	template<class Mat>
	struct matview_cont_level
	{
		static const int value =
				ct_is_continuous<Mat>::value ? 2 : (ct_is_percol_continuous<Mat>::value ? 1 : 0);
	};


	/********************************************
	 *
	 *  continuous matrix
	 *
	 ********************************************/

	// whole x whole

	template<class Mat>
	struct matview_helper<Mat, whole, whole, 2, true>
	{
		typedef typename matrix_traits<Mat>::value_type value_type;

		static const int M = ct_rows<Mat>::value;
		static const int N = ct_cols<Mat>::value;
		typedef cref_matrix<value_type, M, N> type;

		LMAT_ENSURE_INLINE
		static type get(const Mat& mat, whole, whole)
		{
			return type(mat.ptr_data(), mat.nrows(), mat.ncolumns());
		}
	};

	template<class Mat>
	struct matview_helper<Mat, whole, whole, 2, false>
	{
		typedef typename matrix_traits<Mat>::value_type value_type;

		static const int M = ct_rows<Mat>::value;
		static const int N = ct_cols<Mat>::value;
		typedef ref_matrix<value_type, M, N> type;

		LMAT_ENSURE_INLINE
		static type get(Mat& mat, whole, whole)
		{
			return type(mat.ptr_data(), mat.nrows(), mat.ncolumns());
		}
	};

	// whole x range

	template<class Mat>
	struct matview_helper<Mat, whole, range, 2, true>
	{
		typedef typename matrix_traits<Mat>::value_type value_type;

		static const int M = ct_rows<Mat>::value;
		typedef cref_matrix<value_type, M, 0> type;

		LMAT_ENSURE_INLINE
		static type get(const Mat& mat, whole, const range& c)
		{
			return type(mat.ptr_col(c.begin_index()), mat.nrows(), c.num());
		}
	};

	template<class Mat>
	struct matview_helper<Mat, whole, range, 2, false>
	{
		typedef typename matrix_traits<Mat>::value_type value_type;

		static const int M = ct_rows<Mat>::value;
		typedef ref_matrix<value_type, M, 0> type;

		LMAT_ENSURE_INLINE
		static type get(Mat& mat, whole, const range& c)
		{
			return type(mat.ptr_col(c.begin_index()), mat.nrows(), c.num());
		}
	};


	// whole x step

	template<class Mat>
	struct matview_helper<Mat, whole, step_range, 2, true>
	{
		typedef typename matrix_traits<Mat>::value_type value_type;

		static const int M = ct_rows<Mat>::value;
		typedef cref_block<value_type, M, 0> type;

		LMAT_ENSURE_INLINE
		static type get(const Mat& mat, whole, const step_range& c)
		{
			return type(mat.ptr_col(c.begin_index()), mat.nrows(), c.num(),
					mat.col_stride() * c.step());
		}
	};

	template<class Mat>
	struct matview_helper<Mat, whole, step_range, 2, false>
	{
		typedef typename matrix_traits<Mat>::value_type value_type;

		static const int M = ct_rows<Mat>::value;
		typedef ref_block<value_type, M, 0> type;

		LMAT_ENSURE_INLINE
		static type get(Mat& mat, whole, const step_range& c)
		{
			return type(mat.ptr_col(c.begin_index()), mat.nrows(), c.num(),
					mat.col_stride() * c.step());
		}
	};


	// range x whole

	template<class Mat, int L>
	struct matview_helper<Mat, range, whole, L, true>
	{
		typedef typename matrix_traits<Mat>::value_type value_type;

		static const int N = ct_cols<Mat>::value;
		typedef cref_block<value_type, 0, N> type;

		LMAT_ENSURE_INLINE
		static type get(const Mat& mat, const range& r, whole)
		{
			return type(mat.ptr_row(r.begin_index()), r.num(), mat.ncolumns(),
					mat.col_stride());
		}
	};

	template<class Mat, int L>
	struct matview_helper<Mat, range, whole, L, false>
	{
		typedef typename matrix_traits<Mat>::value_type value_type;

		static const int N = ct_cols<Mat>::value;
		typedef ref_block<value_type, 0, N> type;

		LMAT_ENSURE_INLINE
		static type get(Mat& mat, const range& r, whole)
		{
			return type(mat.ptr_row(r.begin_index()), r.num(), mat.ncolumns(),
					mat.col_stride());
		}
	};


	// range x range

	template<class Mat, int L>
	struct matview_helper<Mat, range, range, L, true>
	{
		typedef typename matrix_traits<Mat>::value_type value_type;

		typedef cref_block<value_type, 0, 0> type;

		LMAT_ENSURE_INLINE
		static type get(const Mat& mat, const range& r, const range& c)
		{
			const index_t ri = r.begin_index();
			const index_t ci = c.begin_index();

			return type(&(mat(ri, ci)), r.num(), c.num(), mat.col_stride());
		}
	};

	template<class Mat, int L>
	struct matview_helper<Mat, range, range, L, false>
	{
		typedef typename matrix_traits<Mat>::value_type value_type;

		typedef ref_block<value_type, 0, 0> type;

		LMAT_ENSURE_INLINE
		static type get(Mat& mat, const range& r, const range& c)
		{
			const index_t ri = r.begin_index();
			const index_t ci = c.begin_index();

			return type(&(mat(ri, ci)), r.num(), c.num(), mat.col_stride());
		}
	};


	// range x step

	template<class Mat, int L>
	struct matview_helper<Mat, range, step_range, L, true>
	{
		typedef typename matrix_traits<Mat>::value_type value_type;

		typedef cref_block<value_type, 0, 0> type;

		LMAT_ENSURE_INLINE
		static type get(const Mat& mat, const range& r, const step_range& c)
		{
			const index_t ri = r.begin_index();
			const index_t ci = c.begin_index();

			return type(&(mat(ri, ci)), r.num(), c.num(), mat.col_stride() * c.step());
		}
	};

	template<class Mat, int L>
	struct matview_helper<Mat, range, step_range, L, false>
	{
		typedef typename matrix_traits<Mat>::value_type value_type;

		typedef ref_block<value_type, 0, 0> type;

		LMAT_ENSURE_INLINE
		static type get(Mat& mat, const range& r, const step_range& c)
		{
			const index_t ri = r.begin_index();
			const index_t ci = c.begin_index();

			return type(&(mat(ri, ci)), r.num(), c.num(), mat.col_stride() * c.step());
		}
	};


	// step x whole

	template<class Mat, int L>
	struct matview_helper<Mat, step_range, whole, L, true>
	{
		typedef typename matrix_traits<Mat>::value_type value_type;

		static const int N = ct_cols<Mat>::value;
		typedef cref_grid<value_type, 0, N> type;

		LMAT_ENSURE_INLINE
		static type get(const Mat& mat, const step_range& r, whole)
		{
			return type(mat.ptr_row(r.begin_index()), r.num(), mat.ncolumns(),
					mat.row_stride() * r.step(), mat.col_stride());
		}
	};

	template<class Mat, int L>
	struct matview_helper<Mat, step_range, whole, L, false>
	{
		typedef typename matrix_traits<Mat>::value_type value_type;

		static const int N = ct_cols<Mat>::value;
		typedef ref_grid<value_type, 0, N> type;

		LMAT_ENSURE_INLINE
		static type get(Mat& mat, const step_range& r, whole)
		{
			return type(mat.ptr_row(r.begin_index()), r.num(), mat.ncolumns(),
					mat.row_stride() * r.step(), mat.col_stride());
		}
	};


	// step x range

	template<class Mat, int L>
	struct matview_helper<Mat, step_range, range, L, true>
	{
		typedef typename matrix_traits<Mat>::value_type value_type;

		typedef cref_grid<value_type, 0, 0> type;

		LMAT_ENSURE_INLINE
		static type get(const Mat& mat, const step_range& r, const range& c)
		{
			const index_t ri = r.begin_index();
			const index_t ci = c.begin_index();

			return type(&mat(ri, ci), r.num(), c.num(),
					mat.row_stride() * r.step(), mat.col_stride());
		}
	};

	template<class Mat, int L>
	struct matview_helper<Mat, step_range, range, L, false>
	{
		typedef typename matrix_traits<Mat>::value_type value_type;

		typedef ref_grid<value_type, 0, 0> type;

		LMAT_ENSURE_INLINE
		static type get(Mat& mat, const step_range& r, const range& c)
		{
			const index_t ri = r.begin_index();
			const index_t ci = c.begin_index();

			return type(&mat(ri, ci), r.num(), c.num(),
					mat.row_stride() * r.step(), mat.col_stride());
		}
	};


	// step x step

	template<class Mat, int L>
	struct matview_helper<Mat, step_range, step_range, L, true>
	{
		typedef typename matrix_traits<Mat>::value_type value_type;

		typedef cref_grid<value_type, 0, 0> type;

		LMAT_ENSURE_INLINE
		static type get(const Mat& mat, const step_range& r, const step_range& c)
		{
			const index_t ri = r.begin_index();
			const index_t ci = c.begin_index();

			return type(&mat(ri, ci), r.num(), c.num(),
					mat.row_stride() * r.step(), mat.col_stride() * c.step());
		}
	};

	template<class Mat, int L>
	struct matview_helper<Mat, step_range, step_range, L, false>
	{
		typedef typename matrix_traits<Mat>::value_type value_type;

		typedef ref_grid<value_type, 0, 0> type;

		LMAT_ENSURE_INLINE
		static type get(Mat& mat, const step_range& r, const step_range& c)
		{
			const index_t ri = r.begin_index();
			const index_t ci = c.begin_index();

			return type(&mat(ri, ci), r.num(), c.num(),
					mat.row_stride() * r.step(), mat.col_stride() * c.step());
		}
	};




	/********************************************
	 *
	 *  per-col continuous matrix (specialize)
	 *
	 ********************************************/

	// whole x whole

	template<class Mat>
	struct matview_helper<Mat, whole, whole, 1, true>
	{
		typedef typename matrix_traits<Mat>::value_type value_type;

		static const int M = ct_rows<Mat>::value;
		static const int N = ct_cols<Mat>::value;
		typedef cref_block<value_type, M, N> type;

		LMAT_ENSURE_INLINE
		static type get(const Mat& mat, whole, whole)
		{
			return type(mat.ptr_data(), mat.nrows(), mat.ncolumns(), mat.col_stride());
		}
	};

	template<class Mat>
	struct matview_helper<Mat, whole, whole, 1, false>
	{
		typedef typename matrix_traits<Mat>::value_type value_type;

		static const int M = ct_rows<Mat>::value;
		static const int N = ct_cols<Mat>::value;
		typedef ref_block<value_type, M, N> type;

		LMAT_ENSURE_INLINE
		static type get(Mat& mat, whole, whole)
		{
			return type(mat.ptr_data(), mat.nrows(), mat.ncolumns(), mat.col_stride());
		}
	};


	// whole x range

	template<class Mat>
	struct matview_helper<Mat, whole, range, 1, true>
	{
		typedef typename matrix_traits<Mat>::value_type value_type;

		static const int M = ct_rows<Mat>::value;
		typedef cref_block<value_type, M, 0> type;

		LMAT_ENSURE_INLINE
		static type get(const Mat& mat, whole, const range& c)
		{
			return type(mat.ptr_col(c.begin_index()), mat.nrows(), c.num(), mat.col_stride());
		}
	};

	template<class Mat>
	struct matview_helper<Mat, whole, range, 1, false>
	{
		typedef typename matrix_traits<Mat>::value_type value_type;

		static const int M = ct_rows<Mat>::value;
		typedef ref_block<value_type, M, 0> type;

		LMAT_ENSURE_INLINE
		static type get(Mat& mat, whole, const range& c)
		{
			return type(mat.ptr_col(c.begin_index()), mat.nrows(), c.num(), mat.col_stride());
		}
	};


	// whole x step

	template<class Mat>
	struct matview_helper<Mat, whole, step_range, 1, true>
	{
		typedef typename matrix_traits<Mat>::value_type value_type;

		static const int M = ct_rows<Mat>::value;
		typedef cref_block<value_type, M, 0> type;

		LMAT_ENSURE_INLINE
		static type get(const Mat& mat, whole, const step_range& c)
		{
			return type(mat.ptr_col(c.begin_index()), mat.nrows(), c.num(),
					mat.col_stride() * c.step());
		}
	};

	template<class Mat>
	struct matview_helper<Mat, whole, step_range, 1, false>
	{
		typedef typename matrix_traits<Mat>::value_type value_type;

		static const int M = ct_rows<Mat>::value;
		typedef ref_block<value_type, M, 0> type;

		LMAT_ENSURE_INLINE
		static type get(Mat& mat, whole, const step_range& c)
		{
			return type(mat.ptr_col(c.begin_index()), mat.nrows(), c.num(),
					mat.col_stride() * c.step());
		}
	};




	/********************************************
	 *
	 *  non-continuous matrix
	 *
	 ********************************************/

	// whole x whole

	template<class Mat>
	struct matview_helper<Mat, whole, whole, 0, true>
	{
		typedef typename matrix_traits<Mat>::value_type value_type;

		static const int M = ct_rows<Mat>::value;
		static const int N = ct_cols<Mat>::value;
		typedef cref_grid<value_type, M, N> type;

		LMAT_ENSURE_INLINE
		static type get(const Mat& mat, whole, whole)
		{
			return type(mat.ptr_data(), mat.nrows(), mat.ncolumns(),
					mat.row_stride(), mat.col_stride());
		}
	};

	template<class Mat>
	struct matview_helper<Mat, whole, whole, 0, false>
	{
		typedef typename matrix_traits<Mat>::value_type value_type;

		static const int M = ct_rows<Mat>::value;
		static const int N = ct_cols<Mat>::value;
		typedef ref_grid<value_type, M, N> type;

		LMAT_ENSURE_INLINE
		static type get(Mat& mat, whole, whole)
		{
			return type(mat.ptr_data(), mat.nrows(), mat.ncolumns(),
					mat.row_stride(), mat.col_stride());
		}
	};


	// whole x range

	template<class Mat>
	struct matview_helper<Mat, whole, range, 0, true>
	{
		typedef typename matrix_traits<Mat>::value_type value_type;

		static const int M = ct_rows<Mat>::value;
		static const int N = ct_cols<Mat>::value;
		typedef cref_grid<value_type, M, 0> type;

		LMAT_ENSURE_INLINE
		static type get(const Mat& mat, whole, const range& c)
		{
			return type(mat.ptr_col(c.begin_index()), mat.nrows(), c.num(),
					mat.row_stride(), mat.col_stride());
		}
	};

	template<class Mat>
	struct matview_helper<Mat, whole, range, 0, false>
	{
		typedef typename matrix_traits<Mat>::value_type value_type;

		static const int M = ct_rows<Mat>::value;
		static const int N = ct_cols<Mat>::value;
		typedef ref_grid<value_type, M, 0> type;

		LMAT_ENSURE_INLINE
		static type get(Mat& mat, whole, const range& c)
		{
			return type(mat.ptr_col(c.begin_index()), mat.nrows(), c.num(),
					mat.row_stride(), mat.col_stride());
		}
	};


	// whole x step_range

	template<class Mat>
	struct matview_helper<Mat, whole, step_range, 0, true>
	{
		typedef typename matrix_traits<Mat>::value_type value_type;

		static const int M = ct_rows<Mat>::value;
		typedef cref_grid<value_type, M, 0> type;

		LMAT_ENSURE_INLINE
		static type get(const Mat& mat, whole, const step_range& c)
		{
			return type(mat.ptr_col(c.begin_index()), mat.nrows(), c.num(),
					mat.row_stride(), mat.col_stride() * c.step());
		}
	};

	template<class Mat>
	struct matview_helper<Mat, whole, step_range, 0, false>
	{
		typedef typename matrix_traits<Mat>::value_type value_type;

		static const int M = ct_rows<Mat>::value;
		typedef ref_grid<value_type, M, 0> type;

		LMAT_ENSURE_INLINE
		static type get(Mat& mat, whole, const step_range& c)
		{
			return type(mat.ptr_col(c.begin_index()), mat.nrows(), c.num(),
					mat.row_stride(), mat.col_stride() * c.step());
		}
	};


	// range x whole

	template<class Mat>
	struct matview_helper<Mat, range, whole, 0, true>
	{
		typedef typename matrix_traits<Mat>::value_type value_type;

		static const int N = ct_cols<Mat>::value;
		typedef cref_grid<value_type, 0, N> type;

		LMAT_ENSURE_INLINE
		static type get(const Mat& mat, const range& r, whole)
		{
			return type(mat.ptr_row(r.begin_index()), r.num(), mat.ncolumns(),
					mat.row_stride(), mat.col_stride());
		}
	};

	template<class Mat>
	struct matview_helper<Mat, range, whole, 0, false>
	{
		typedef typename matrix_traits<Mat>::value_type value_type;

		static const int N = ct_cols<Mat>::value;
		typedef ref_grid<value_type, 0, N> type;

		LMAT_ENSURE_INLINE
		static type get(Mat& mat, const range& r, whole)
		{
			return type(mat.ptr_row(r.begin_index()), r.num(), mat.ncolumns(),
					mat.row_stride(), mat.col_stride());
		}
	};


	// range x range

	template<class Mat>
	struct matview_helper<Mat, range, range, 0, true>
	{
		typedef typename matrix_traits<Mat>::value_type value_type;

		typedef cref_grid<value_type, 0, 0> type;

		LMAT_ENSURE_INLINE
		static type get(const Mat& mat, const range& r, const range& c)
		{
			const index_t ri = r.begin_index();
			const index_t ci = c.begin_index();

			return type(&mat(ri, ci), r.num(), c.num(),
					mat.row_stride(), mat.col_stride());
		}
	};

	template<class Mat>
	struct matview_helper<Mat, range, range, 0, false>
	{
		typedef typename matrix_traits<Mat>::value_type value_type;

		typedef ref_grid<value_type, 0, 0> type;

		LMAT_ENSURE_INLINE
		static type get(Mat& mat, const range& r, const range& c)
		{
			const index_t ri = r.begin_index();
			const index_t ci = c.begin_index();

			return type(&mat(ri, ci), r.num(), c.num(),
					mat.row_stride(), mat.col_stride());
		}
	};


	// range x step_range

	template<class Mat>
	struct matview_helper<Mat, range, step_range, 0, true>
	{
		typedef typename matrix_traits<Mat>::value_type value_type;

		typedef cref_grid<value_type, 0, 0> type;

		LMAT_ENSURE_INLINE
		static type get(const Mat& mat, const range& r, const step_range& c)
		{
			const index_t ri = r.begin_index();
			const index_t ci = c.begin_index();

			return type(&mat(ri, ci), r.num(), c.num(),
					mat.row_stride(), mat.col_stride() * c.step());
		}
	};

	template<class Mat>
	struct matview_helper<Mat, range, step_range, 0, false>
	{
		typedef typename matrix_traits<Mat>::value_type value_type;

		typedef ref_grid<value_type, 0, 0> type;

		LMAT_ENSURE_INLINE
		static type get(Mat& mat, const range& r, const step_range& c)
		{
			const index_t ri = r.begin_index();
			const index_t ci = c.begin_index();

			return type(&mat(ri, ci), r.num(), c.num(),
					mat.row_stride(), mat.col_stride() * c.step());
		}
	};

} }

#endif /* MATRIX_MATVIEWS_INTERNAL_H_ */
