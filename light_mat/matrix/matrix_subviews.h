/**
 * @file matrix_subviews.h
 *
 * Facilities to create matrix sub-views
 *
 * @author Dahua Lin
 */

#ifdef _MSC_VER
#pragma once
#endif

#ifndef LIGHTMAT_MATRIX_SUBVIEWS_H_
#define LIGHTMAT_MATRIX_SUBVIEWS_H_

#include "internal/matrix_colviews_internal.h"
#include "internal/matrix_rowviews_internal.h"
#include "internal/matrix_matviews_internal.h"

namespace lmat
{
	// column views

	template<class Mat, class Rgn>
	struct colview_map
	{
		static const bool is_percol_cont = meta::is_percol_continuous<Mat>::value;
		static const bool is_readonly = meta::is_readonly<Mat>::value;

		typedef internal::colview_helper<Mat, Rgn, is_percol_cont, true> chelper_t;
		typedef internal::colview_helper<Mat, Rgn, is_percol_cont, is_readonly> helper_t;

		typedef typename chelper_t::type const_type;
		typedef typename helper_t::type type;

		typedef typename
				meta::if_<meta::is_readonly<Mat>,
					const_type,
					dense_mutable_view<type>
				>::type _type;

		LMAT_ENSURE_INLINE
		static const_type get(const Mat& mat, const index_t j, const Rgn& rgn)
		{
			return chelper_t::get(mat, j, rgn);
		}

		LMAT_ENSURE_INLINE
		static _type get(Mat& mat, const index_t j, const Rgn& rgn)
		{
			return helper_t::get(mat, j, rgn);
		}
	};


	// row views

	template<class Mat, class Rgn>
	struct rowview_map
	{
		static const bool is_perrow_cont = meta::is_continuous<Mat>::value && meta::is_row<Mat>::value;
		static const bool is_readonly = meta::is_readonly<Mat>::value;

		typedef internal::rowview_helper<Mat, Rgn, is_perrow_cont, true> chelper_t;
		typedef internal::rowview_helper<Mat, Rgn, is_perrow_cont, is_readonly> helper_t;

		typedef typename chelper_t::type const_type;
		typedef typename helper_t::type type;

		typedef typename
				meta::if_<meta::is_readonly<Mat>,
					const_type,
					dense_mutable_view<type>
				>::type _type;

		LMAT_ENSURE_INLINE
		static const_type get(const Mat& mat, const index_t i, const Rgn& rgn)
		{
			return chelper_t::get(mat, i, rgn);
		}

		LMAT_ENSURE_INLINE
		static _type get(Mat& mat, const index_t i, const Rgn& rgn)
		{
			return helper_t::get(mat, i, rgn);
		}
	};


	// vec views

	struct NO_VECVIEWS_FOR_NON_COMPILE_TIME_VECTORS { };

	template<class Mat, class Rgn>
	struct vecview_map
	{
		typedef typename
				meta::if_<meta::is_col<Mat>,
					colview_map<Mat, Rgn>,
					typename
					meta::if_<meta::is_row<Mat>,
						rowview_map<Mat, Rgn>,
						NO_VECVIEWS_FOR_NON_COMPILE_TIME_VECTORS
					>::type
				>::type intern_t;

		typedef typename intern_t::const_type const_type;
		typedef typename intern_t::type type;
		typedef typename intern_t::_type _type;

		LMAT_ENSURE_INLINE
		static const_type get(const Mat& mat, const Rgn& rgn)
		{
			return intern_t::get(mat, 0, rgn);
		}

		LMAT_ENSURE_INLINE
		static _type get(Mat& mat, const Rgn& rgn)
		{
			return intern_t::get(mat, 0, rgn);
		}
	};


	// diagonal views

	template<class Mat>
	struct diagview_map
	{
		typedef typename matrix_traits<Mat>::value_type value_type;

		static const int M = meta::nrows<Mat>::value;
		static const int N = meta::ncols<Mat>::value;
		static const int L = M < N ? M : N;

		typedef cref_grid<value_type, L, 1> const_type;

		typedef typename
				meta::if_<meta::is_readonly<Mat>,
					cref_grid<value_type, L, 1>,
					 ref_grid<value_type, L, 1>
				>::type type;

		typedef typename
				meta::if_<meta::is_readonly<Mat>,
					const_type,
					dense_mutable_view<type>
				>::type _type;

		LMAT_ENSURE_INLINE
		static const_type get(const Mat& mat)
		{
			const index_t m = mat.nrows();
			const index_t n = mat.ncolumns();
			const index_t len = m < n ? m : n;
			const index_t step = mat.row_stride() + mat.col_stride();

			return const_type(mat.ptr_data(), len, 1, step, 0);
		}

		LMAT_ENSURE_INLINE
		static _type get(Mat& mat)
		{
			const index_t m = mat.nrows();
			const index_t n = mat.ncolumns();
			const index_t len = m < n ? m : n;
			const index_t step = mat.row_stride() + mat.col_stride();

			return type(mat.ptr_data(), len, 1, step, 0);
		}
	};



	// mat views

	template<class Mat, class RowRange, class ColRange>
	struct matview_map
	{
		static const int cont_level = internal::matview_cont_level<Mat>::value;
		static const bool is_readonly = meta::is_readonly<Mat>::value;

		typedef internal::matview_helper<Mat, RowRange, ColRange, cont_level, true> chelper_t;
		typedef internal::matview_helper<Mat, RowRange, ColRange, cont_level, is_readonly> helper_t;

		typedef typename chelper_t::type const_type;
		typedef typename helper_t::type type;

		typedef typename
				meta::if_<meta::is_readonly<Mat>,
					const_type,
					dense_mutable_view<type>
				>::type _type;

		LMAT_ENSURE_INLINE
		static const_type get(const Mat& mat, const RowRange &r, const ColRange& c)
		{
			return chelper_t::get(mat, r, c);
		}

		LMAT_ENSURE_INLINE
		static _type get(Mat& mat, const RowRange &r, const ColRange& c)
		{
			return helper_t::get(mat, r, c);
		}
	};


}

#endif /* MATRIX_SUBVIEWS_H_ */
