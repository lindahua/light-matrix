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

namespace lmat { namespace internal {

	template<int M, int N, typename ContLevel>
	struct matrix_fill_scheme
	{
		const matrix_shape<M, N> shape;

		LMAT_ENSURE_INLINE
		matrix_fill_scheme(index_t m, index_t n)
		: shape(m, n) { }

		LMAT_ENSURE_INLINE
		index_t nrows() const { return shape.nrows(); }

		LMAT_ENSURE_INLINE
		index_t ncolumns() const { return shape.ncolumns(); }

		LMAT_ENSURE_INLINE
		index_t nelems() const { return shape.nelems(); }
	};

	template<typename T, typename DMat>
	matrix_fill_scheme<
		meta::nrows< DMat >::value,
		meta::ncols< DMat >::value,
		typename meta::contiguousness< DMat>::type >
	LMAT_ENSURE_INLINE
	inline get_fill_scheme(const IRegularMatrix<DMat, T>& dmat)
	{
		typedef matrix_fill_scheme<
			meta::nrows< DMat >::value,
			meta::ncols< DMat >::value,
			typename meta::contiguousness<DMat>::type > scheme_t;

		return scheme_t(dmat.nrows(), dmat.ncolumns());
	}


	/********************************************
	 *
	 *  Core routines
	 *
	 ********************************************/

	template<typename T>
	LMAT_ENSURE_INLINE
	inline void _fill_singlevec(index_t len, const T& v, T *pd, index_t d_step)
	{
		if (d_step == 1) fill_vec(len, pd, v);
		else fill_vec(len, step_ptr(pd, d_step), v);
	}


	template<typename T>
	LMAT_ENSURE_INLINE
	inline void _fill_multicol(index_t m, index_t n, const T& v,
			T *pd, index_t dst_cs)
	{
		for (index_t j = 0; j < n; ++j)
		{
			fill_vec(m, pd + j * dst_cs, v);
		}
	}


	template<typename T>
	LMAT_ENSURE_INLINE
	inline void _fill_multicol(index_t m, index_t n, const T& v,
			T *pd, index_t dst_rs, index_t dst_cs)
	{
		if (dst_rs == 1)
		{
			_fill_multicol(m, n, v, pd, dst_cs);
		}
		else
		{
			for (index_t j = 0; j < n; ++j)
			{
				fill_vec(m, step_ptr(pd + j * dst_cs, dst_rs), v);
			}
		}
	}


	template<typename T>
	LMAT_ENSURE_INLINE
	inline void _zero_singlevec(index_t len, T *pd, index_t d_step)
	{
		if (d_step == 1) zero_vec(len, pd);
		else zero_vec(len, step_ptr(pd, d_step));
	}


	template<typename T>
	LMAT_ENSURE_INLINE
	inline void _zero_multicol(index_t m, index_t n, T *pd, index_t dst_cs)
	{
		for (index_t j = 0; j < n; ++j)
		{
			zero_vec(m, pd + j * dst_cs);
		}
	}


	template<typename T>
	LMAT_ENSURE_INLINE
	inline void _zero_multicol(index_t m, index_t n, T *pd, index_t dst_rs, index_t dst_cs)
	{
		if (dst_rs == 1)
		{
			_zero_multicol(m, n, pd, dst_cs);
		}
		else
		{
			for (index_t j = 0; j < n; ++j)
			{
				zero_vec(m, step_ptr(pd + j * dst_cs, dst_rs));
			}
		}
	}


	/******************************************************
	 *
	 *  fill routines
	 *
	 ******************************************************/

	template<typename T, class DMat, int M, int N>
	LMAT_ENSURE_INLINE
	inline void fill(const T& v, IRegularMatrix<DMat, T>& dmat, const matrix_fill_scheme<M, N, cont_level::whole>& sch)
	{
		fill_vec(sch.nelems(), dmat.ptr_data(), v);
	}

	template<typename T, class DMat, int M, int N>
	inline void fill(const T& v, IRegularMatrix<DMat, T>& dmat,
			const matrix_fill_scheme<M, N, cont_level::percol>& sch)
	{
		const index_t m = sch.nrows();
		const index_t n = sch.ncolumns();

		T *pd = dmat.ptr_data();

		if (n == 1)
		{
			if (m == 1)
				*pd = v;
			else
				fill_vec(m, pd, v);
		}
		else
		{
			const index_t cs = dmat.col_stride();
			if (m == 1)
				_fill_singlevec(n, v, pd, cs);
			else
				_fill_multicol(m, n, v, pd, cs);
		}
	}

	template<typename T, class DMat, int M, int N>
	inline void fill(const T& v, IRegularMatrix<DMat, T>& dmat,
			const matrix_fill_scheme<1, N, cont_level::percol>& sch)
	{
		const index_t n = sch.ncolumns();

		T *pd = dmat.ptr_data();

		if (n == 1)
		{
			*pd = v;
		}
		else
		{
			_fill_singlevec(n, v, pd, dmat.col_stride());
		}
	}

	template<typename T, class DMat, int M, int N>
	inline void fill(const T& v, IRegularMatrix<DMat, T>& dmat,
			const matrix_fill_scheme<M, N, cont_level::none>& sch)
	{
		const index_t m = sch.nrows();
		const index_t n = sch.ncolumns();

		T *pd = dmat.ptr_data();

		if (n == 1)
		{
			if (m == 1)
				*pd = v;
			else
			{
				_fill_singlevec(m, v, pd, dmat.row_stride());
			}
		}
		else
		{
			const index_t cs = dmat.col_stride();
			if (m == 1)
			{
				_fill_singlevec(n, v, pd, cs);
			}
			else
			{
				const index_t rs = dmat.row_stride();
				_fill_multicol(m, n, v, pd, rs, cs);
			}
		}
	}



	/******************************************************
	 *
	 *  zero routines
	 *
	 ******************************************************/

	template<typename T, class DMat, int M, int N>
	LMAT_ENSURE_INLINE
	inline void zero(IRegularMatrix<DMat, T>& dmat, const matrix_fill_scheme<M, N, cont_level::whole>& sch)
	{
		zero_vec(sch.nelems(), dmat.ptr_data());
	}

	template<typename T, class DMat, int M, int N>
	inline void zero(IRegularMatrix<DMat, T>& dmat,
			const matrix_fill_scheme<M, N, cont_level::percol>& sch)
	{
		const index_t m = sch.nrows();
		const index_t n = sch.ncolumns();

		T *pd = dmat.ptr_data();

		if (n == 1)
		{
			if (m == 1)
				set_zero_value(*pd);
			else
				zero_vec(m, pd);
		}
		else
		{
			const index_t cs = dmat.col_stride();
			if (m == 1)
				_zero_singlevec(n, pd, cs);
			else
				_zero_multicol(m, n, pd, cs);
		}
	}

	template<typename T, class DMat, int M, int N>
	inline void zero(IRegularMatrix<DMat, T>& dmat,
			const matrix_fill_scheme<1, N, cont_level::percol>& sch)
	{
		const index_t n = sch.ncolumns();

		T *pd = dmat.ptr_data();

		if (n == 1)
		{
			set_zero_value(*pd);
		}
		else
		{
			_zero_singlevec(n, pd, dmat.col_stride());
		}
	}

	template<typename T, class DMat, int M, int N>
	inline void zero(IRegularMatrix<DMat, T>& dmat,
			const matrix_fill_scheme<M, N, cont_level::none>& sch)
	{
		const index_t m = sch.nrows();
		const index_t n = sch.ncolumns();

		T *pd = dmat.ptr_data();

		if (n == 1)
		{
			if (m == 1)
				set_zero_value(*pd);
			else
			{
				_zero_singlevec(m, pd, dmat.row_stride());
			}
		}
		else
		{
			const index_t cs = dmat.col_stride();
			if (m == 1)
			{
				_zero_singlevec(n, pd, cs);
			}
			else
			{
				const index_t rs = dmat.row_stride();
				_zero_multicol(m, n, pd, rs, cs);
			}
		}
	}


} }

#endif /* MATRIX_FILL_INTERNAL_H_ */


