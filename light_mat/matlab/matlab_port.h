/*
 * @file matlab_port.h
 *
 * @brief
 *
 * @author Dahua Lin
 */

#ifdef _MSC_VER
#pragma once
#endif

#ifndef LIGHTMAT_MATLAB_PORT_H_
#define LIGHTMAT_MATLAB_PORT_H_

#include <light_mat/matlab/marray.h>
#include <light_mat/matlab/mdispatch.h>

#include <light_mat/matrix/matrix_eval.h>

namespace lmat { namespace matlab {

	template<typename T>
	LMAT_ENSURE_INLINE
	inline cref_matrix<T> view2d(const_marray m)
	{
		return cref_matrix<T>(m.data<T>(), m.nrows(), m.ncolumns());
	}

	template<typename T>
	LMAT_ENSURE_INLINE
	inline cref_col<T> view_as_col(const_marray m)
	{
		return cref_col<T>(m.data<T>(), m.nelems());
	}

	template<typename T>
	LMAT_ENSURE_INLINE
	inline cref_row<T> view_as_row(const_marray m)
	{
		return cref_row<T>(m.data<T>(), m.nelems());
	}

	template<class Mat, typename T>
	LMAT_ENSURE_INLINE
	inline marray to_marray(const IMatrixXpr<Mat, T>& in)
	{
		marray r = marray::numeric_matrix<T>(in.nrows(), in.ncolumns());

		typedef ref_matrix<T, ct_rows<Mat>::value, ct_cols<Mat>::value> view_t;
		view_t view(r.data<T>(), in.nrows(), in.ncolumns());
		view = in;

		return r;
	}

	template<typename T>
	LMAT_ENSURE_INLINE
	inline marray to_marray(const T *src, index_t m, index_t n)
	{
		marray r = marray::numeric_matrix<T>(m, n);
		copy_mem(m * n, src, r.data<T>());
		return r;
	}

	template<typename T>
	LMAT_ENSURE_INLINE
	inline marray to_marray(const double *src, index_t m, index_t n)
	{
		marray r = marray::double_matrix(m, n);
		copy_mem(m * n, src, r.data<T>());
		return r;
	}

} }


#endif /* MATLAB_PORT_H_ */
