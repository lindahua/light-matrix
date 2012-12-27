/**
 * @file matrix_asvec_internal.h
 *
 * @brief Internal implementation of as-vector views
 *
 * @author Dahua Lin
 */

#ifndef LIGHTMAT_MATRIX_ASVEC_INTERNAL_H_
#define LIGHTMAT_MATRIX_ASVEC_INTERNAL_H_

#include <light_mat/matrix/matrix_meta.h>

namespace lmat { namespace internal {

	template<class Mat>
	struct as_vec_indicator
	{
		static const int value =
				meta::is_continuous<Mat>::value ? 1 :
				(meta::is_col<Mat>::value ? 2 :
				(meta::is_row<Mat>::value ? 3 : 0));
	};


	/********************************************
	 *
	 *   as column
	 *
	 ********************************************/

	template<class Mat, int I, bool IsReadOnly>
	struct _as_col_map;

	// continuous

	template<class Mat>
	struct _as_col_map<Mat, 1, true>
	{
		typedef typename matrix_traits<Mat>::value_type T;
		typedef cref_matrix<T, meta::nelems<Mat>::value, 1> type;

		LMAT_ENSURE_INLINE static type get(const Mat& mat)
		{
			return type(mat.ptr_data(), mat.nelems(), 1);
		}
	};

	template<class Mat>
	struct _as_col_map<Mat, 1, false>
	{
		typedef typename matrix_traits<Mat>::value_type T;
		typedef ref_matrix<T, meta::nelems<Mat>::value, 1> type;

		LMAT_ENSURE_INLINE static type get(Mat& mat)
		{
			return type(mat.ptr_data(), mat.nelems(), 1);
		}
	};

	// generic column

	template<class Mat>
	struct _as_col_map<Mat, 2, true>
	{
		typedef typename matrix_traits<Mat>::value_type T;
		typedef cref_grid<T, meta::nrows<Mat>::value, 1> type;

		LMAT_ENSURE_INLINE static type get(const Mat& mat)
		{
			return type(mat.ptr_data(), mat.nrows(), 1, mat.row_stride(), 0);
		}
	};

	template<class Mat>
	struct _as_col_map<Mat, 2, false>
	{
		typedef typename matrix_traits<Mat>::value_type T;
		typedef ref_grid<T, meta::nrows<Mat>::value, 1> type;

		LMAT_ENSURE_INLINE static type get(Mat& mat)
		{
			return type(mat.ptr_data(), mat.nrows(), 1, mat.row_stride(), 0);
		}
	};

	// generic row

	template<class Mat>
	struct _as_col_map<Mat, 3, true>
	{
		typedef typename matrix_traits<Mat>::value_type T;
		typedef cref_grid<T, meta::ncols<Mat>::value, 1> type;

		LMAT_ENSURE_INLINE static type get(const Mat& mat)
		{
			return type(mat.ptr_data(), mat.ncolumns(), 1, mat.col_stride(), 0);
		}
	};

	template<class Mat>
	struct _as_col_map<Mat, 3, false>
	{
		typedef typename matrix_traits<Mat>::value_type T;
		typedef ref_grid<T, meta::ncols<Mat>::value, 1> type;

		LMAT_ENSURE_INLINE static type get(Mat& mat)
		{
			return type(mat.ptr_data(), mat.ncolumns(), 1, mat.col_stride(), 0);
		}
	};

	// wrapper map

	template<class Mat>
	struct as_col_map
	{
		static const int I = as_vec_indicator<Mat>::value;
		static const bool is_readonly = meta::is_readonly<Mat>::value;

		typedef _as_col_map<Mat, I, true> chelper_t;
		typedef _as_col_map<Mat, I, is_readonly> helper_t;

		typedef typename chelper_t::type const_type;
		typedef typename helper_t::type type;

		typedef typename
				meta::if_<meta::is_readonly<Mat>,
					const_type,
					dense_mutable_view<type>
				>::type _type;

		LMAT_ENSURE_INLINE
		static const_type get(const Mat& mat)
		{
			return chelper_t::get(mat);
		}

		LMAT_ENSURE_INLINE
		static _type get(Mat& mat)
		{
			return helper_t::get(mat);
		}
	};


	/********************************************
	 *
	 *   as column
	 *
	 ********************************************/

	template<class Mat, int I, bool IsReadOnly>
	struct _as_row_map;

	// continuous

	template<class Mat>
	struct _as_row_map<Mat, 1, true>
	{
		typedef typename matrix_traits<Mat>::value_type T;
		typedef cref_matrix<T, 1, meta::nelems<Mat>::value> type;

		LMAT_ENSURE_INLINE static type get(const Mat& mat)
		{
			return type(mat.ptr_data(), 1, mat.nelems());
		}
	};

	template<class Mat>
	struct _as_row_map<Mat, 1, false>
	{
		typedef typename matrix_traits<Mat>::value_type T;
		typedef ref_matrix<T, 1, meta::nelems<Mat>::value> type;

		LMAT_ENSURE_INLINE static type get(Mat& mat)
		{
			return type(mat.ptr_data(), 1, mat.nelems());
		}
	};

	// generic column

	template<class Mat>
	struct _as_row_map<Mat, 2, true>
	{
		typedef typename matrix_traits<Mat>::value_type T;
		typedef cref_block<T, 1, meta::nrows<Mat>::value> type;

		LMAT_ENSURE_INLINE static type get(const Mat& mat)
		{
			return type(mat.ptr_data(), 1, mat.nrows(), mat.row_stride());
		}
	};

	template<class Mat>
	struct _as_row_map<Mat, 2, false>
	{
		typedef typename matrix_traits<Mat>::value_type T;
		typedef ref_block<T, 1, meta::nrows<Mat>::value> type;

		LMAT_ENSURE_INLINE static type get(Mat& mat)
		{
			return type(mat.ptr_data(), 1, mat.nrows(), mat.row_stride());
		}
	};

	// generic row

	template<class Mat>
	struct _as_row_map<Mat, 3, true>
	{
		typedef typename matrix_traits<Mat>::value_type T;
		typedef cref_block<T, 1, meta::ncols<Mat>::value> type;

		LMAT_ENSURE_INLINE static type get(const Mat& mat)
		{
			return type(mat.ptr_data(), 1, mat.ncolumns(), mat.col_stride());
		}
	};

	template<class Mat>
	struct _as_row_map<Mat, 3, false>
	{
		typedef typename matrix_traits<Mat>::value_type T;
		typedef ref_block<T, 1, meta::ncols<Mat>::value> type;

		LMAT_ENSURE_INLINE static type get(Mat& mat)
		{
			return type(mat.ptr_data(), 1, mat.ncolumns(), mat.col_stride());
		}
	};


	// wrapper map

	template<class Mat>
	struct as_row_map
	{
		static const int I = as_vec_indicator<Mat>::value;
		static const bool is_readonly = meta::is_readonly<Mat>::value;

		typedef _as_row_map<Mat, I, true> chelper_t;
		typedef _as_row_map<Mat, I, is_readonly> helper_t;

		typedef typename chelper_t::type const_type;
		typedef typename helper_t::type type;

		typedef typename
				meta::if_<meta::is_readonly<Mat>,
					const_type,
					dense_mutable_view<type>
				>::type _type;

		LMAT_ENSURE_INLINE
		static const_type get(const Mat& mat)
		{
			return chelper_t::get(mat);
		}

		LMAT_ENSURE_INLINE
		static _type get(Mat& mat)
		{
			return helper_t::get(mat);
		}
	};

} }

#endif /* MATRIX_ASVEC_INTERNAL_H_ */
