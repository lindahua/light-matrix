/**
 * @file matrix_compare_internal.h
 *
 * Internal implementation of matrix comparison
 *
 * @author Dahua Lin
 */

#ifdef _MSC_VER
#pragma once
#endif

#ifndef LIGHTMAT_MATRIX_COMPARE_INTERNAL_H_
#define LIGHTMAT_MATRIX_COMPARE_INTERNAL_H_

#include <light_mat/matrix/matrix_properties.h>
#include <light_mat/common/vec_algs.h>

namespace lmat { namespace detail {

	/********************************************
	 *
	 *  Matrix equality
	 *
	 ********************************************/

	template<typename T, class LMat, class RMat>
	struct scalar_comparer
	{
		LMAT_ENSURE_INLINE
		static bool all_equal(const LMat& a, const RMat& b)
		{
			return *a.ptr_data() == *b.ptr_data();
		}
	};


	template<typename T, class LMat, class RMat>
	struct cc_comparer
	{
		LMAT_ENSURE_INLINE
		static bool all_equal(const LMat& a, const RMat& b)
		{
			return vec_equal(a.nelems(), a.ptr_data(), b.ptr_data());
		}
	};

	template<typename T, class LMat, class RMat>
	struct column_comparer
	{
		LMAT_ENSURE_INLINE
		static bool all_equal(const LMat& a, const RMat& b)
		{
			const index_t step_s = a.row_stride();
			const index_t step_t = b.row_stride();

			if (step_s == 1)
			{
				if (step_t == 1)
				{
					return vec_equal(a.nrows(), a.ptr_data(), b.ptr_data());
				}
				else
				{
					return vec_equal(a.nrows(), a.ptr_data(), b.ptr_data(), step_t);
				}
			}
			else
			{
				if (step_t == 1)
				{
					return vec_equal(a.nrows(), a.ptr_data(), step_s, b.ptr_data());
				}
				else
				{
					return vec_equal(a.nrows(), a.ptr_data(), step_s, b.ptr_data(), step_t);
				}
			}
		}
	};


	template<typename T, class LMat, class RMat>
	struct row_comparer
	{
		LMAT_ENSURE_INLINE
		static bool all_equal(const LMat& a, const RMat& b)
		{
			const index_t step_s = a.col_stride();
			const index_t step_t = b.col_stride();

			if (step_s == 1)
			{
				if (step_t == 1)
				{
					return vec_equal(a.ncolumns(), a.ptr_data(), b.ptr_data());
				}
				else
				{
					return vec_equal(a.ncolumns(), a.ptr_data(), b.ptr_data(), step_t);
				}
			}
			else
			{
				if (step_t == 1)
				{
					return vec_equal(a.ncolumns(), a.ptr_data(), step_s, b.ptr_data());
				}
				else
				{
					return vec_equal(a.ncolumns(), a.ptr_data(), step_s, b.ptr_data(), step_t);
				}
			}
		}
	};

	template<typename T, class LMat, class RMat>
	struct genmat_comparer
	{
		static bool all_equal(const LMat& a, const RMat& b)
		{
			const index_t m = a.nrows();
			const index_t n = a.ncolumns();

			if (n == 1)
			{
				return column_comparer<T, LMat, RMat>::all_equal(a, b);
			}
			else if (m == 1)
			{
				return row_comparer<T, LMat, RMat>::all_equal(a, b);
			}
			else
			{
				const index_t step_s = a.row_stride();
				const index_t step_t = b.row_stride();

				if (step_s == 1)
				{
					if (step_t == 1)
					{
						if (a.col_stride() == m && b.col_stride() == m)
						{
							return vec_equal(a.nelems(), a.ptr_data(), b.ptr_data());
						}
						else
						{
							for (index_t j = 0; j < n; ++j)
							{
								if (!vec_equal(m, a.ptr_col(j), b.ptr_col(j)))
									return false;
							}
							return true;
						}
					}
					else
					{
						for (index_t j = 0; j < n; ++j)
						{
							if (!vec_equal(m, a.ptr_col(j), b.ptr_col(j), step_t))
								return false;
						}
						return true;
					}
				}
				else
				{
					if (step_t == 1)
					{
						for (index_t j = 0; j < n; ++j)
						{
							if (!vec_equal(m, a.ptr_col(j), step_s, b.ptr_col(j)))
								return false;
						}
						return true;
					}
					else
					{
						for (index_t j = 0; j < n; ++j)
						{
							if (!vec_equal(m, a.ptr_col(j), step_s, b.ptr_col(j), step_t))
								return false;
						}
						return true;
					}
				}
			}
		}
	};


	template<class LMat, class RMat>
	struct mat_comparer_map
	{
		typedef typename binary_value_type<LMat, RMat>::type T;

		static const int M = binary_ct_rows<LMat, RMat>::value;
		static const int N = binary_ct_cols<LMat, RMat>::value;

		typedef typename
				if_c<M == 1 && N == 1,
					scalar_comparer<T, LMat, RMat>,
					typename
					if_c<ct_is_continuous<LMat>::value && ct_is_continuous<RMat>::value,
						cc_comparer<T, LMat, RMat>,
						typename
						if_c<N == 1,
							column_comparer<T, LMat, RMat>,
							typename
							if_c<M == 1,
								row_comparer<T, LMat, RMat>,
								genmat_comparer<T, LMat, RMat>
							>::type
						>::type
					>::type
				>::type type;
	};


	/********************************************
	 *
	 *  Matrix approximation
	 *
	 ********************************************/

	template<typename T, class LMat, class RMat>
	struct scalar_approx_comparer
	{
		LMAT_ENSURE_INLINE
		static bool all_approx(const LMat& a, const RMat& b, const T& tol)
		{
			return is_approx(*a.ptr_data(), *b.ptr_data(), tol);
		}
	};


	template<typename T, class LMat, class RMat>
	struct cc_approx_comparer
	{
		LMAT_ENSURE_INLINE
		static bool all_approx(const LMat& a, const RMat& b, const T& tol)
		{
			return vec_approx(a.nelems(), a.ptr_data(), b.ptr_data(), tol);
		}
	};

	template<typename T, class LMat, class RMat>
	struct column_approx_comparer
	{
		LMAT_ENSURE_INLINE
		static bool all_approx(const LMat& a, const RMat& b, const T& tol)
		{
			const index_t step_s = a.row_stride();
			const index_t step_t = b.row_stride();

			if (step_s == 1)
			{
				if (step_t == 1)
				{
					return vec_approx(a.nrows(), a.ptr_data(), b.ptr_data(), tol);
				}
				else
				{
					return vec_approx(a.nrows(), a.ptr_data(), b.ptr_data(), step_t, tol);
				}
			}
			else
			{
				if (step_t == 1)
				{
					return vec_approx(a.nrows(), a.ptr_data(), step_s, b.ptr_data(), tol);
				}
				else
				{
					return vec_approx(a.nrows(), a.ptr_data(), step_s, b.ptr_data(), step_t, tol);
				}
			}
		}
	};


	template<typename T, class LMat, class RMat>
	struct row_approx_comparer
	{
		LMAT_ENSURE_INLINE
		static bool all_approx(const LMat& a, const RMat& b, const T& tol)
		{
			const index_t step_s = a.col_stride();
			const index_t step_t = b.col_stride();

			if (step_s == 1)
			{
				if (step_t == 1)
				{
					return vec_approx(a.ncolumns(), a.ptr_data(), b.ptr_data(), tol);
				}
				else
				{
					return vec_approx(a.ncolumns(), a.ptr_data(), b.ptr_data(), step_t, tol);
				}
			}
			else
			{
				if (step_t == 1)
				{
					return vec_approx(a.ncolumns(), a.ptr_data(), step_s, b.ptr_data(), tol);
				}
				else
				{
					return vec_approx(a.ncolumns(), a.ptr_data(), step_s, b.ptr_data(), step_t, tol);
				}
			}
		}
	};

	template<typename T, class LMat, class RMat>
	struct genmat_approx_comparer
	{
		static bool all_approx(const LMat& a, const RMat& b, const T& tol)
		{
			const index_t m = a.nrows();
			const index_t n = a.ncolumns();

			if (n == 1)
			{
				return column_approx_comparer<T, LMat, RMat>::all_approx(a, b, tol);
			}
			else if (m == 1)
			{
				return row_approx_comparer<T, LMat, RMat>::all_approx(a, b, tol);
			}
			else
			{
				const index_t step_s = a.row_stride();
				const index_t step_t = b.row_stride();

				if (step_s == 1)
				{
					if (step_t == 1)
					{
						if (a.col_stride() == m && b.col_stride() == m)
						{
							return vec_approx(a.nelems(), a.ptr_data(), b.ptr_data(), tol);
						}
						else
						{
							for (index_t j = 0; j < n; ++j)
							{
								if (!vec_approx(m, a.ptr_col(j), b.ptr_col(j), tol))
									return false;
							}
							return true;
						}
					}
					else
					{
						for (index_t j = 0; j < n; ++j)
						{
							if (!vec_approx(m, a.ptr_col(j), b.ptr_col(j), step_t, tol))
								return false;
						}
						return true;
					}
				}
				else
				{
					if (step_t == 1)
					{
						for (index_t j = 0; j < n; ++j)
						{
							if (!vec_approx(m, a.ptr_col(j), step_s, b.ptr_col(j), tol))
								return false;
						}
						return true;
					}
					else
					{
						for (index_t j = 0; j < n; ++j)
						{
							if (!vec_approx(m, a.ptr_col(j), step_s, b.ptr_col(j), step_t, tol))
								return false;
						}
						return true;
					}
				}
			}
		}
	};


	template<class LMat, class RMat>
	struct mat_approx_comparer_map
	{
		typedef typename binary_value_type<LMat, RMat>::type T;

		static const int M = binary_ct_rows<LMat, RMat>::value;
		static const int N = binary_ct_cols<LMat, RMat>::value;

		typedef typename
				if_c<M == 1 && N == 1,
					scalar_approx_comparer<T, LMat, RMat>,
					typename
					if_c<ct_is_continuous<LMat>::value && ct_is_continuous<RMat>::value,
						cc_approx_comparer<T, LMat, RMat>,
						typename
						if_c<N == 1,
							column_approx_comparer<T, LMat, RMat>,
							typename
							if_c<M == 1,
								row_approx_comparer<T, LMat, RMat>,
								genmat_approx_comparer<T, LMat, RMat>
							>::type
						>::type
					>::type
				>::type type;
	};

} }

#endif /* MATRIX_COMPARE_INTERNAL_H_ */
