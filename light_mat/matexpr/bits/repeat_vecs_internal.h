/**
 * @file repeat_vecs_internal.h
 *
 * Internal implementation for repeat_vecs
 *
 * @author Dahua Lin
 */

#ifdef _MSC_VER
#pragma once
#endif

#ifndef REPEAT_VECS_INTERNAL_H_
#define REPEAT_VECS_INTERNAL_H_

#include <light_mat/matrix/matrix_classes.h>

namespace lmat { namespace detail {

	template<typename T, int M, int N, class Arg, class DMat>
	inline void repcol_evaluate(const matrix_shape<M, N>& shape,
			const IDenseMatrix<Arg, T>& arg_,
			IDenseMatrix<DMat, T>& dmat_)
	{
		typedef typename colview_map<DMat, whole>::type col_t;

		const Arg& arg = arg_.derived();
		DMat& dmat = dmat_.derived();

		const index_t n = shape.ncolumns();

		if (n == 1)
		{
			col_t col(dmat.column(0));
			copy(arg, col);
		}
		else
		{
			const index_t m = shape.nrows();
			if (m == 1)
			{
				const T v = *arg.ptr_data();
				fill(dmat, v);
			}
			else
			{
				const index_t rs = arg.row_stride();
				if (rs == 1)
				{
					for (index_t j = 0; j < n; ++j)
					{
						col_t col(dmat.column(j));
						copy(arg, col);
					}
				}
				else
				{
					dense_matrix<T, M, 1> cache(arg);
					for (index_t j = 0; j < n; ++j)
					{
						col_t col(dmat.column(j));
						copy(arg, col);
					}
				}
			}
		}
	}


	template<typename T, int M, int N, class Arg, class DMat>
	inline void repcol_evaluate(const matrix_shape<M, N>& shape,
			const IMatrixXpr<Arg, T>& arg_,
			IDenseMatrix<DMat, T>& dmat_)
	{
		typedef typename colview_map<DMat, whole>::type col_t;

		const Arg& arg = arg_.derived();
		DMat& dmat = dmat_.derived();

		const index_t n = shape.ncolumns();

		if (n == 1)
		{
			col_t col(dmat.column(0));
			default_evaluate(arg, col);
		}
		else
		{
			dense_matrix<T, M, 1> cache(arg);

			const index_t m = shape.nrows();
			if (m == 1)
			{
				T v = *cache.ptr_data();
				fill(dmat, v);
			}
			else
			{
				for (index_t j = 0; j < n; ++j)
				{
					col_t col(dmat.column(j));
					copy(cache, col);
				}
			}
		}
	}


	template<typename T, int M, int N, class Arg, class DMat>
	inline void reprow_evaluate(const matrix_shape<M, N>& shape,
			const IDenseMatrix<Arg, T>& arg_,
			IDenseMatrix<DMat, T>& dmat_)
	{
		typedef typename rowview_map<DMat, whole>::type row_t;
		typedef typename colview_map<DMat, whole>::type col_t;

		const Arg& arg = arg_.derived();
		DMat& dmat = dmat_.derived();

		const index_t m = shape.nrows();

		if (m == 1)
		{
			row_t row(dmat.row(0));
			copy(arg, row);
		}
		else
		{
			const index_t n = shape.ncolumns();
			if (n == 1)
			{
				T v = *arg.ptr_data();
				fill(dmat, v);
			}
			else
			{
				for (index_t j = 0; j < n; ++j)
				{
					T v = arg(0, j);
					col_t col(dmat.column(j));
					fill(col, v);
				}
			}
		}
	}


	template<typename T, int M, int N, class Arg, class DMat>
	inline void reprow_evaluate(const matrix_shape<M, N>& shape,
			const IMatrixXpr<Arg, T>& arg_,
			IDenseMatrix<DMat, T>& dmat_)
	{
		typedef typename rowview_map<DMat, whole>::type row_t;
		typedef typename colview_map<DMat, whole>::type col_t;

		const Arg& arg = arg_.derived();
		DMat& dmat = dmat_.derived();

		const index_t m = shape.nrows();

		if (m == 1)
		{
			row_t row(dmat.row(0));
			default_evaluate(arg, row);
		}
		else
		{
			dense_matrix<T, 1, N> cache(arg);
			const index_t n = shape.ncolumns();

			if (n == 1)
			{
				T v = *cache.ptr_data();
				fill(dmat, v);
			}
			else
			{
				for (index_t j = 0; j < n; ++j)
				{
					T v = cache[j];
					col_t col(dmat.column(j));
					fill(col, v);
				}
			}
		}
	}


	template<typename T, int M, bool IsDense> struct repcol_cap;

	template<typename T, int M>
	struct repcol_cap<T, M, true>
	{
		const bool to_cache;
		dense_col<T> cache;
		const T *pdata;

		template<class Arg>
		LMAT_ENSURE_INLINE
		repcol_cap(const Arg& arg)
		: to_cache(arg.row_stride() != 1)
		, cache(to_cache ? arg.nrows() : 0)
		, pdata(to_cache ? cache.ptr_data() : arg.ptr_data())
		{
			if (to_cache)
			{
				copy(arg, cache);
			}
		}
	};

	template<typename T, int M>
	struct repcol_cap<T, M, false>
	{
		dense_col<T> cache;
		const T *pdata;

		template<class Arg>
		LMAT_ENSURE_INLINE
		repcol_cap(const Arg& arg)
		: cache(arg)
		, pdata(cache.ptr_data())
		{ }
	};

	template<class Arg, bool SupportsLinearIndex> struct reprow_cap;

	template<class Arg>
	struct reprow_cap<Arg, true>
	{
		typedef typename matrix_traits<Arg>::value_type T;
		const Arg& arg;

		LMAT_ENSURE_INLINE
		reprow_cap(const Arg& arg_)
		: arg(arg_) { }

		LMAT_ENSURE_INLINE
		T operator[] (const index_t j) const
		{
			return arg[j];
		}
	};

	template<class Arg>
	struct reprow_cap<Arg, false>
	{
		typedef typename matrix_traits<Arg>::value_type T;
		dense_row<T, ct_cols<Arg>::value> cache;

		LMAT_ENSURE_INLINE
		reprow_cap(const Arg& arg_)
		: cache(arg_) { }

		LMAT_ENSURE_INLINE
		T operator[] (const index_t j) const
		{
			return cache[j];
		}
	};


} }


#endif /* REPEAT_VECS_INTERNAL_H_ */
