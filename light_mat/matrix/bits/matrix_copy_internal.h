/*
 * @file matrix_copy_internal.h
 *
 * Internal implementation of matrix copy
 *
 * @author Dahua Lin
 */

#ifndef _MSC_VER
#pragma once
#endif

#ifndef LIGHTMAT_MATRIX_COPY_INTERNAL_H_
#define LIGHTMAT_MATRIX_COPY_INTERNAL_H_

#include <light_mat/matrix/matrix_properties.h>

namespace lmat { namespace detail {


	/********************************************
	 *
	 *  Copy between matrix and pointer
	 *
	 ********************************************/

	template<typename T, class Mat>
	struct scalar_copier_p
	{
		LMAT_ENSURE_INLINE
		static void copy(const Mat& mat, T *dst)
		{
			*dst = *mat.ptr_data();
		}

		LMAT_ENSURE_INLINE
		static void copy(const T *src, Mat& mat)
		{
			*mat.ptr_data() = *src;
		}
	};

	template<typename T, class Mat>
	struct cc_copier_p
	{
		LMAT_ENSURE_INLINE
		static void copy(const Mat& mat, T *dst)
		{
			copy_vec(mat.nelems(), mat.ptr_data(), dst);
		}

		LMAT_ENSURE_INLINE
		static void copy(const T *src, Mat& mat)
		{
			copy_vec(mat.nelems(), src, mat.ptr_data());
		}
	};


	template<typename T, class Mat>
	struct column_copier_p
	{
		LMAT_ENSURE_INLINE
		static void copy(const Mat& mat, T *dst)
		{
			if (mat.row_stride() == 1)
			{
				copy_vec(mat.nrows(), mat.ptr_data(), dst);
			}
			else
			{
				copy_vec(mat.nrows(), mat.ptr_data(), mat.row_stride(), dst);
			}
		}

		LMAT_ENSURE_INLINE
		static void copy(const T *src, Mat& mat)
		{
			if (mat.row_stride() == 1)
			{
				copy_vec(mat.nrows(), src, mat.ptr_data());
			}
			else
			{
				copy_vec(mat.nrows(), src, mat.ptr_data(), mat.row_stride());
			}
		}
	};


	template<typename T, class Mat>
	struct row_copier_p
	{
		LMAT_ENSURE_INLINE
		static void copy(const Mat& mat, T *dst)
		{
			if (mat.col_stride() == 1)
			{
				copy_vec(mat.ncolumns(), mat.ptr_data(), dst);
			}
			else
			{
				copy_vec(mat.ncolumns(), mat.ptr_data(), mat.col_stride(), dst);
			}
		}

		LMAT_ENSURE_INLINE
		static void copy(const T *src, Mat& mat)
		{
			if (mat.col_stride() == 1)
			{
				copy_vec(mat.ncolumns(), src, mat.ptr_data());
			}
			else
			{
				copy_vec(mat.ncolumns(), src, mat.ptr_data(), mat.col_stride());
			}
		}
	};


	template<typename T, class Mat>
	struct genmat_copier_p
	{
		static void copy(const Mat& mat, T *dst)
		{
			const index_t m = mat.nrows();
			const index_t n = mat.ncolumns();

			if (n == 1)
			{
				column_copier_p<T, Mat>::copy(mat, dst);
			}
			else if (m == 1)
			{
				row_copier_p<T, Mat>::copy(mat, dst);
			}
			else
			{
				const index_t rs = mat.row_stride();

				if (rs == 1)
				{
					for (index_t j = 0; j < n; ++j, dst += m)
					{
						copy_vec(m, mat.ptr_col(j), dst);
					}
				}
				else
				{
					for (index_t j = 0; j < n; ++j, dst += m)
					{
						copy_vec(m, mat.ptr_col(j), rs, dst);
					}
				}
			}
		}

		static void copy(const T *src, Mat& mat)
		{
			const index_t m = mat.nrows();
			const index_t n = mat.ncolumns();

			if (n == 1)
			{
				column_copier_p<T, Mat>::copy(src, mat);
			}
			else if (m == 1)
			{
				row_copier_p<T, Mat>::copy(src, mat);
			}
			else
			{
				const index_t rs = mat.row_stride();

				if (rs == 1)
				{
					for (index_t j = 0; j < n; ++j, src += m)
					{
						copy_vec(m, src, mat.ptr_col(j));
					}
				}
				else
				{
					for (index_t j = 0; j < n; ++j, src += m)
					{
						copy_vec(m, src, mat.ptr_col(j), rs);
					}
				}
			}
		}
	};


	template<class Mat>
	struct mat_copier_p_map
	{
		typedef typename matrix_traits<Mat>::value_type T;
		static const int M = ct_rows<Mat>::value;
		static const int N = ct_cols<Mat>::value;

		typedef typename
				if_c<M == 1 && N == 1,
					scalar_copier_p<T, Mat>,
					typename
					if_c<ct_is_continuous<Mat>::value,
						cc_copier_p<T, Mat>,
						typename
						if_c<N == 1,
							column_copier_p<T, Mat>,
							typename
							if_c<M == 1,
								row_copier_p<T, Mat>,
								genmat_copier_p<T, Mat>
							>::type
						>::type
					>::type
				>::type type;
	};



	/********************************************
	 *
	 *  Copy between two matrices
	 *
	 ********************************************/

	template<typename T, class SMat, class DMat>
	struct scalar_copier
	{
		LMAT_ENSURE_INLINE
		static void copy(const SMat& s, DMat& t)
		{
			*t.ptr_data() = *s.ptr_data();
		}
	};


	template<typename T, class SMat, class DMat>
	struct cc_copier
	{
		LMAT_ENSURE_INLINE
		static void copy(const SMat& s, DMat& t)
		{
			copy_vec(s.nelems(), s.ptr_data(), t.ptr_data());
		}
	};

	template<typename T, class SMat, class DMat>
	struct column_copier
	{
		LMAT_ENSURE_INLINE
		static void copy(const SMat& s, DMat& t)
		{
			const index_t step_s = s.row_stride();
			const index_t step_t = t.row_stride();

			if (step_s == 1)
			{
				if (step_t == 1)
				{
					copy_vec(s.nrows(), s.ptr_data(), t.ptr_data());
				}
				else
				{
					copy_vec(s.nrows(), s.ptr_data(), t.ptr_data(), step_t);
				}
			}
			else
			{
				if (step_t == 1)
				{
					copy_vec(s.nrows(), s.ptr_data(), step_s, t.ptr_data());
				}
				else
				{
					copy_vec(s.nrows(), s.ptr_data(), step_s, t.ptr_data(), step_t);
				}
			}


		}
	};


	template<typename T, class SMat, class DMat>
	struct row_copier
	{
		LMAT_ENSURE_INLINE
		static void copy(const SMat& s, DMat& t)
		{
			const index_t step_s = s.col_stride();
			const index_t step_t = t.col_stride();

			if (step_s == 1)
			{
				if (step_t == 1)
				{
					copy_vec(s.ncolumns(), s.ptr_data(), t.ptr_data());
				}
				else
				{
					copy_vec(s.ncolumns(), s.ptr_data(), t.ptr_data(), step_t);
				}
			}
			else
			{
				if (step_t == 1)
				{
					copy_vec(s.ncolumns(), s.ptr_data(), step_s, t.ptr_data());
				}
				else
				{
					copy_vec(s.ncolumns(), s.ptr_data(), step_s, t.ptr_data(), step_t);
				}
			}
		}
	};

	template<typename T, class SMat, class DMat>
	struct genmat_copier
	{
		static void copy(const SMat& s, DMat& t)
		{
			const index_t m = s.nrows();
			const index_t n = s.ncolumns();

			if (n == 1)
			{
				column_copier<T, SMat, DMat>::copy(s, t);
			}
			else if (m == 1)
			{
				row_copier<T, SMat, DMat>::copy(s, t);
			}
			else
			{
				const index_t step_s = s.row_stride();
				const index_t step_t = t.row_stride();

				if (step_s == 1)
				{
					if (step_t == 1)
					{
						if (s.col_stride() == m && t.col_stride() == m)
						{
							copy_vec(s.nelems(), s.ptr_data(), t.ptr_data());
						}
						else
						{
							for (index_t j = 0; j < n; ++j)
								copy_vec(m, s.ptr_col(j), t.ptr_col(j));
						}
					}
					else
					{
						for (index_t j = 0; j < n; ++j)
							copy_vec(m, s.ptr_col(j), t.ptr_col(j), step_t);
					}
				}
				else
				{
					if (step_t == 1)
					{
						for (index_t j = 0; j < n; ++j)
							copy_vec(m, s.ptr_col(j), step_s, t.ptr_col(j));
					}
					else
					{
						for (index_t j = 0; j < n; ++j)
							copy_vec(m, s.ptr_col(j), step_s, t.ptr_col(j), step_t);
					}
				}
			}
		}
	};


	template<class SMat, class DMat>
	struct mat_copier_map
	{
		typedef typename binary_value_type<SMat, DMat>::type T;

		static const int M = binary_ct_rows<SMat, DMat>::value;
		static const int N = binary_ct_cols<SMat, DMat>::value;

		typedef typename
				if_c<M == 1 && N == 1,
					scalar_copier<T, SMat, DMat>,
					typename
					if_c<ct_is_continuous<SMat>::value && ct_is_continuous<DMat>::value,
						cc_copier<T, SMat, DMat>,
						typename
						if_c<N == 1,
							column_copier<T, SMat, DMat>,
							typename
							if_c<M == 1,
								row_copier<T, SMat, DMat>,
								genmat_copier<T, SMat, DMat>
							>::type
						>::type
					>::type
				>::type type;
	};




} }

#endif 
