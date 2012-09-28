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

#include <light_mat/matrix/matrix_meta.h>

namespace lmat
{
	/********************************************
	 *
	 *  column views
	 *
	 ********************************************/

	template<class Mat>
	struct colview_map<Mat, whole>
	{
		typedef typename matrix_traits<Mat>::value_type value_type;

		static const int ctrows = ct_rows<Mat>::value;

		typedef cref_matrix<value_type, ctrows, 1> const_type;
		typedef  ref_matrix<value_type, ctrows, 1> non_const_type;

		typedef typename
				if_<is_readonly_mat<Mat>,
					const_type,
					non_const_type
				>::type _type;

		typedef typename
				if_<is_readonly_mat<Mat>,
					const_type,
					dense_mutable_view<non_const_type>
				>::type type;

		LMAT_ENSURE_INLINE
		static const_type get(const Mat& mat, const index_t j, whole)
		{
			return const_type(mat.ptr_col(j), mat.nrows(), 1);
		}

		LMAT_ENSURE_INLINE
		static type get(Mat& mat, const index_t j, whole)
		{
			return _type(mat.ptr_col(j), mat.nrows(), 1);
		}
	};


	template<class Mat>
	struct colview_map<Mat, range>
	{
		typedef typename matrix_traits<Mat>::value_type value_type;

		typedef cref_matrix<value_type, 0, 1> const_type;
		typedef  ref_matrix<value_type, 0, 1> non_const_type;

		typedef typename
				if_<is_readonly_mat<Mat>,
					const_type,
					non_const_type
				>::type _type;

		typedef typename
				if_<is_readonly_mat<Mat>,
					const_type,
					dense_mutable_view<non_const_type>
				>::type type;

		LMAT_ENSURE_INLINE
		static const_type get(const Mat& mat, const index_t j, const range &rg)
		{
			return const_type(mat.ptr_col(j) + rg.begin_index(), rg.num(), 1);
		}

		LMAT_ENSURE_INLINE
		static type get(Mat& mat, const index_t j, const range& rg)
		{
			return _type(mat.ptr_col(j) + rg.begin_index(), rg.num(), 1);
		}
	};



	/********************************************
	 *
	 *  row views
	 *
	 ********************************************/

	template<class Mat>
	struct rowview_map<Mat, whole>
	{
		typedef typename matrix_traits<Mat>::value_type value_type;

		static const int ctcols = ct_cols<Mat>::value;

		typedef cref_matrix_ex<value_type, 1, ctcols> const_type;
		typedef  ref_matrix_ex<value_type, 1, ctcols> non_const_type;

		typedef typename
				if_<is_readonly_mat<Mat>,
					const_type,
					non_const_type
				>::type _type;

		typedef typename
				if_<is_readonly_mat<Mat>,
					const_type,
					dense_mutable_view<non_const_type>
				>::type type;

		LMAT_ENSURE_INLINE
		static const_type get(const Mat& mat, const index_t i, whole)
		{
			return const_type(mat.ptr_data() + i, 1, mat.ncolumns(), mat.lead_dim());
		}

		LMAT_ENSURE_INLINE
		static type get(Mat& mat, const index_t i, whole)
		{
			return _type(mat.ptr_data() + i, 1, mat.ncolumns(), mat.lead_dim());
		}

	};


	template<class Mat>
	struct rowview_map<Mat, range>
	{
		typedef typename matrix_traits<Mat>::value_type value_type;

		typedef cref_matrix_ex<value_type, 1, 0> const_type;
		typedef  ref_matrix_ex<value_type, 1, 0> non_const_type;

		typedef typename
				if_<is_readonly_mat<Mat>,
					const_type,
					non_const_type
				>::type _type;

		typedef typename
				if_<is_readonly_mat<Mat>,
					const_type,
					dense_mutable_view<non_const_type>
				>::type type;

		LMAT_ENSURE_INLINE
		static const_type get(const Mat& mat, const index_t i, const range& rg)
		{
			return const_type(mat.ptr_col(rg.begin_index()) + i, 1, rg.num(), mat.lead_dim());
		}

		LMAT_ENSURE_INLINE
		static type get(Mat& mat, const index_t i, const range& rg)
		{
			return _type(mat.ptr_col(rg.begin_index()) + i, 1, rg.num(), mat.lead_dim());
		}

	};



	/********************************************
	 *
	 *  subviews
	 *
	 ********************************************/

	namespace detail
	{
		template<class Mat, int CTCols, bool IsCont> struct multicol_helper;

		template<class Mat, int CTCols>
		struct multicol_helper<Mat, CTCols, true>
		{
			typedef typename matrix_traits<Mat>::value_type value_type;
			static const int ctrows = ct_rows<Mat>::value;

			typedef cref_matrix<value_type, ctrows, CTCols> const_type;
			typedef  ref_matrix<value_type, ctrows, CTCols> non_const_type;

			typedef typename
					if_<is_readonly_mat<Mat>,
						const_type,
						non_const_type
					>::type _type;

			typedef typename
					if_<is_readonly_mat<Mat>,
						const_type,
						dense_mutable_view<non_const_type>
					>::type type;

			LMAT_ENSURE_INLINE
			static const_type get(const Mat& mat, const index_t j, const index_t n)
			{
				return const_type(mat.ptr_col(j), mat.nrows(), n);
			}

			LMAT_ENSURE_INLINE
			static type get(Mat& mat, const index_t j, const index_t n)
			{
				return _type(mat.ptr_col(j), mat.nrows(), n);
			}
		};


		template<class Mat, int CTCols>
		struct multicol_helper<Mat, CTCols, false>
		{
			typedef typename matrix_traits<Mat>::value_type value_type;
			static const int ctrows = ct_rows<Mat>::value;

			typedef cref_matrix_ex<value_type, ctrows, CTCols> const_type;
			typedef  ref_matrix_ex<value_type, ctrows, CTCols> non_const_type;

			typedef typename
					if_<is_readonly_mat<Mat>,
						const_type,
						non_const_type
					>::type _type;

			typedef typename
					if_<is_readonly_mat<Mat>,
						const_type,
						dense_mutable_view<non_const_type>
					>::type type;

			LMAT_ENSURE_INLINE
			static const_type get(const Mat& mat, const index_t j, const index_t n)
			{
				return const_type(mat.ptr_col(j), mat.nrows(), n, mat.lead_dim());
			}

			LMAT_ENSURE_INLINE
			static type get(Mat& mat, const index_t j, const index_t n)
			{
				return _type(mat.ptr_col(j), mat.nrows(), n, mat.lead_dim());
			}
		};
	}

	template<class Mat>
	struct matview_map<Mat, whole, whole>
	{
		static const bool is_continuous = ct_has_continuous_layout<Mat>::value;
		typedef detail::multicol_helper<Mat, ct_cols<Mat>::value, is_continuous> helper_t;

		typedef typename helper_t::const_type const_type;
		typedef typename helper_t::type type;

		LMAT_ENSURE_INLINE
		static const_type get(const Mat& mat, whole, whole)
		{
			return helper_t::get(mat, 0, mat.ncolumns());
		}

		LMAT_ENSURE_INLINE
		static type get(Mat& mat, whole, whole)
		{
			return helper_t::get(mat, 0, mat.ncolumns());
		}
	};


	template<class Mat>
	struct matview_map<Mat, whole, range>
	{
		static const bool is_continuous = ct_has_continuous_layout<Mat>::value;
		typedef detail::multicol_helper<Mat, 0, is_continuous> helper_t;

		typedef typename helper_t::const_type const_type;
		typedef typename helper_t::type type;

		LMAT_ENSURE_INLINE
		static const_type get(const Mat& mat, whole, const range& rg)
		{
			return helper_t::get(mat, rg.begin_index(), rg.num());
		}

		LMAT_ENSURE_INLINE
		static type get(Mat& mat, whole, const range& rg)
		{
			return helper_t::get(mat, rg.begin_index(), rg.num());
		}
	};


	template<class Mat>
	struct matview_map<Mat, range, whole>
	{
		typedef typename matrix_traits<Mat>::value_type value_type;
		static const int ctcols = ct_cols<Mat>::value;

		typedef cref_matrix_ex<value_type, 0, ctcols> const_type;
		typedef  ref_matrix_ex<value_type, 0, ctcols> non_const_type;

		typedef typename
				if_<is_readonly_mat<Mat>,
					const_type,
					non_const_type
				>::type _type;

		typedef typename
				if_<is_readonly_mat<Mat>,
					const_type,
					dense_mutable_view<non_const_type>
				>::type type;

		LMAT_ENSURE_INLINE
		static const_type get(const Mat& mat, const range& rg, whole)
		{
			return const_type(mat.ptr_data() + rg.begin_index(), rg.num(),
					mat.ncolumns(), mat.lead_dim());
		}

		LMAT_ENSURE_INLINE
		static type get(Mat& mat, const range& rg, whole)
		{
			return _type(mat.ptr_data() + rg.begin_index(), rg.num(),
					mat.ncolumns(), mat.lead_dim());
		}
	};


	template<class Mat>
	struct matview_map<Mat, range, range>
	{
		typedef typename matrix_traits<Mat>::value_type value_type;

		typedef cref_matrix_ex<value_type, 0, 0> const_type;
		typedef  ref_matrix_ex<value_type, 0, 0> non_const_type;

		typedef typename
				if_<is_readonly_mat<Mat>,
					const_type,
					non_const_type
				>::type _type;

		typedef typename
				if_<is_readonly_mat<Mat>,
					const_type,
					dense_mutable_view<non_const_type>
				>::type type;

		LMAT_ENSURE_INLINE
		static const_type get(const Mat& mat, const range& rrg, const range& crg)
		{
			return const_type(
					mat.ptr_col(crg.begin_index()) + rrg.begin_index(),
					rrg.num(), crg.num(), mat.lead_dim());
		}

		LMAT_ENSURE_INLINE
		static type get(Mat& mat, const range& rrg, const range& crg)
		{
			return _type(
					mat.ptr_col(crg.begin_index()) + rrg.begin_index(),
					rrg.num(), crg.num(), mat.lead_dim());
		}
	};

}

#endif /* MATRIX_SUBVIEWS_H_ */
