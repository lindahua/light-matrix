/**
 * @file matrix_fill_internal.h
 *
 * Internal implementation of matrix filling
 *
 * @author Dahua Lin
 */

#ifdef _MSC_VER
#pragma once
#endif

#ifndef LIGHTMAT_MATRIX_FILL_INTERNAL_H_
#define LIGHTMAT_MATRIX_FILL_INTERNAL_H_

#include <light_mat/matrix/matrix_properties.h>

namespace lmat { namespace detail {

	template<typename T, class DMat>
	struct scalar_filler
	{
		LMAT_ENSURE_INLINE
		static void zero(DMat& dmat)
		{
			set_zero_value(*dmat.ptr_data());
		}

		LMAT_ENSURE_INLINE
		static void fill(const T& v, DMat& dmat)
		{
			*dmat.ptr_data() = v;
		}
	};


	template<typename T, class DMat>
	struct column_filler
	{
		LMAT_ENSURE_INLINE
		static void zero(DMat& dmat)
		{
			if (dmat.row_stride() == 1)
			{
				zero_vec(dmat.nrows(), dmat.ptr_data());
			}
			else
			{
				zero_vec(dmat.nrows(), dmat.ptr_data(), dmat.row_stride());
			}
		}

		LMAT_ENSURE_INLINE
		static void fill(const T& v, DMat& dmat)
		{
			if (dmat.row_stride() == 1)
			{
				fill_vec(v, dmat.nrows(), dmat.ptr_data());
			}
			else
			{
				fill_vec(v, dmat.nrows(), dmat.ptr_data(), dmat.row_stride());
			}
		}
	};


	template<typename T, class DMat>
	struct row_filler
	{
		LMAT_ENSURE_INLINE
		static void zero(DMat& dmat)
		{
			if (dmat.col_stride() == 1)
			{
				zero_vec(dmat.ncolumns(), dmat.ptr_data());
			}
			else
			{
				zero_vec(dmat.ncolumns(), dmat.ptr_data(), dmat.col_stride());
			}
		}

		LMAT_ENSURE_INLINE
		static void fill(const T& v, DMat& dmat)
		{
			if (dmat.col_stride() == 1)
			{
				fill_vec(v, dmat.ncolumns(), dmat.ptr_data());
			}
			else
			{
				fill_vec(v, dmat.ncolumns(), dmat.ptr_data(), dmat.col_stride());
			}
		}
	};


	template<typename T, class DMat>
	struct mat_filler
	{
		LMAT_ENSURE_INLINE
		static void zero(DMat& dmat)
		{
			const index_t m = dmat.nrows();
			const index_t n = dmat.ncolumns();

			if (n == 1)
			{
				column_filler<T, DMat>::zero(dmat);
			}
			else if (m == 1)
			{
				row_filler<T, DMat>::zero(dmat);
			}
			else
			{
				const index_t rs = dmat.row_stride();

				if (rs == 1)
				{
					for (index_t j = 0; j < n; ++j)
					{
						zero_vec(m, dmat.ptr_col(j));
					}
				}
				else
				{
					for (index_t j = 0; j < n; ++j)
					{
						zero_vec(m, dmat.ptr_col(j), rs);
					}
				}
			}
		}

		LMAT_ENSURE_INLINE
		static void fill(const T& v, DMat& dmat)
		{
			const index_t m = dmat.nrows();
			const index_t n = dmat.ncolumns();

			if (n == 1)
			{
				column_filler<T, DMat>::fill(v, dmat);
			}
			else if (m == 1)
			{
				row_filler<T, DMat>::fill(v, dmat);
			}
			else
			{
				const index_t rs = dmat.row_stride();

				if (rs == 1)
				{
					for (index_t j = 0; j < n; ++j)
					{
						fill_vec(v, m, dmat.ptr_col(j));
					}
				}
				else
				{
					for (index_t j = 0; j < n; ++j)
					{
						fill_vec(v, m, dmat.ptr_col(j), rs);
					}
				}
			}
		}
	};


	template<class DMat>
	struct mat_filler_map
	{
		typedef typename matrix_traits<DMat>::value_type T;
		static const int M = ct_rows<DMat>::value;
		static const int N = ct_cols<DMat>::value;

		typedef typename
				if_c<N == 1,
					typename
					if_c<M == 1,
						scalar_filler<T, DMat>,
						column_filler<T, DMat>
					>::type,
					typename
					if_c<M == 1,
						row_filler<T, DMat>,
						mat_filler<T, DMat>
					>::type
				>::type type;
	};



} }

#endif /* MATRIX_FILL_INTERNAL_H_ */


