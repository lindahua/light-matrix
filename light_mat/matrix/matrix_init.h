/**
 * @file matrix_init.h
 *
 * @brief Matrix initialization support
 *
 * @author Dahua Lin
 */

#ifdef _MSC_VER
#pragma once
#endif

#ifndef LIGHTMAT_MATRIX_INIT_H_
#define LIGHTMAT_MATRIX_INIT_H_

#include <light_mat/matrix/matrix_properties.h>
#include <initializer_list>

namespace lmat
{

	// initialization routines

	template<typename T, class Mat>
	LMAT_ENSURE_INLINE
	typename meta::enable_if<meta::supports_linear_index<Mat>, void>::type
	vec_initialize(IRegularMatrix<Mat, T>& mat, const std::initializer_list<T>& lst)
	{
		LMAT_CHECK_DIMS( static_cast<index_t>(lst.size()) == mat.nelems() );

		const index_t n = mat.nelems();

		auto ps = begin(lst);
		Mat& dst = mat.derived();

		for (index_t j = 0; j < n; ++j)
		{
			dst[j] = *ps++;
		}
	}

	template<typename T, class Mat>
	LMAT_ENSURE_INLINE
	typename meta::enable_if<meta::is_percol_continuous<Mat>, void>::type
	row_major_initialize(IRegularMatrix<Mat, T>& mat, const std::initializer_list<T>& lst)
	{
		LMAT_CHECK_DIMS( static_cast<index_t>(lst.size()) == mat.nelems() );

		const index_t m = mat.nrows();
		const index_t n = mat.ncolumns();

		auto ps = begin(lst);
		index_t cs = mat.col_stride();

		for (index_t i = 0; i < m; ++i)
		{
			T *pd = mat.ptr_data() + i;
			for (index_t j = 0; j < n; ++j, pd += cs)
				*pd = *ps++;
		}
	}

	template<typename T, class Mat>
	LMAT_ENSURE_INLINE
	typename meta::enable_if<meta::is_percol_continuous<Mat>, void>::type
	col_major_initialize(IRegularMatrix<Mat, T>& mat, const std::initializer_list<T>& lst)
	{
		LMAT_CHECK_DIMS( static_cast<index_t>(lst.size()) == mat.nelems() );

		const index_t m = mat.nrows();
		const index_t n = mat.ncolumns();

		auto ps = begin(lst);

		for (index_t j = 0; j < n; ++j)
		{
			T *pd = mat.ptr_col(j);
			for (index_t i = 0; i < m; ++i)
				*pd++ = *ps++;
		}
	}



	// initializer wrapper

	template<typename T>
	struct row_major_initializer
	{
		const std::initializer_list<T>& _list;

		LMAT_ENSURE_INLINE
		explicit row_major_initializer(const std::initializer_list<T>& lst)
		: _list(lst)
		{ }
	};

	template<typename T>
	struct col_major_initializer
	{
		const std::initializer_list<T>& _list;

		LMAT_ENSURE_INLINE
		explicit col_major_initializer(const std::initializer_list<T>& lst)
		: _list(lst)
		{ }
	};


	// construct initializers

	template<typename T>
	LMAT_ENSURE_INLINE
	inline row_major_initializer<T> rm_(const std::initializer_list<T>& lst)
	{
		return row_major_initializer<T>(lst);
	}

	template<typename T>
	LMAT_ENSURE_INLINE
	inline col_major_initializer<T> cm_(const std::initializer_list<T>& lst)
	{
		return col_major_initializer<T>(lst);
	}


}

#endif
