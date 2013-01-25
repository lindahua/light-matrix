/*
 * @file matrix_copy_internal.h
 *
 * Internal implementation of matrix copy
 *
 * @author Dahua Lin
 */

#ifdef _MSC_VER
#pragma once
#endif

#ifndef LIGHTMAT_MATRIX_COPY_INTERNAL_H_
#define LIGHTMAT_MATRIX_COPY_INTERNAL_H_

#include <light_mat/matrix/matrix_properties.h>

namespace lmat { namespace internal {

	template<index_t M, index_t N, typename ContLevel>
	struct matrix_copy_scheme
	{
		const matrix_shape<M, N> shape;

		LMAT_ENSURE_INLINE
		matrix_copy_scheme(index_t m, index_t n)
		: shape(m, n) { }

		LMAT_ENSURE_INLINE
		index_t nrows() const { return shape.nrows(); }

		LMAT_ENSURE_INLINE
		index_t ncolumns() const { return shape.ncolumns(); }

		LMAT_ENSURE_INLINE
		index_t nelems() const { return shape.nelems(); }
	};


	template<typename T, typename DMat>
	matrix_copy_scheme<
		meta::nrows< DMat >::value,
		meta::ncols< DMat >::value,
		typename meta::contiguousness< DMat>::type >
	LMAT_ENSURE_INLINE
	inline get_copy_scheme(const IRegularMatrix<DMat, T>& dmat)
	{
		typedef matrix_copy_scheme<
			meta::nrows< DMat >::value,
			meta::ncols< DMat >::value,
			typename meta::contiguousness<DMat>::type > scheme_t;

		return scheme_t(dmat.nrows(), dmat.ncolumns());
	}

	template<typename T, typename SMat, typename DMat>
	matrix_copy_scheme<
		meta::common_nrows< SMat, DMat >::value,
		meta::common_ncols< SMat, DMat >::value,
		typename meta::contiguousness< SMat, DMat>::type >
	LMAT_ENSURE_INLINE
	inline get_copy_scheme(const IRegularMatrix<SMat, T>& smat, const IRegularMatrix<DMat, T>& dmat)
	{
		typedef matrix_copy_scheme<
			meta::common_nrows< SMat, DMat >::value,
			meta::common_ncols< SMat, DMat >::value,
			typename meta::contiguousness< SMat, DMat>::type > scheme_t;

		return scheme_t(common_nrows(smat, dmat), common_ncols(smat, dmat));
	}


	/********************************************
	 *
	 *  Core routines
	 *
	 ********************************************/

	template<typename T>
	LMAT_ENSURE_INLINE
	inline void _copy_singlevec(index_t len, const T* ps, index_t s_step, T *pd)
	{
		if (s_step == 1) copy_vec(len, ps, pd);
		else copy_vec(len, step_ptr(ps, s_step), pd);
	}

	template<typename T>
	LMAT_ENSURE_INLINE
	inline void _copy_singlevec(index_t len, const T* ps, T *pd, index_t d_step)
	{
		if (d_step == 1) copy_vec(len, ps, pd);
		else copy_vec(len, ps, step_ptr(pd, d_step));
	}

	template<typename T>
	LMAT_ENSURE_INLINE
	inline void _copy_singlevec(index_t len, const T* ps, index_t s_step, T *pd, index_t d_step)
	{
		if (s_step == 1)
		{
			if (d_step == 1) copy_vec(len, ps, pd);
			else copy_vec(len, ps, step_ptr(pd, d_step));
		}
		else
		{
			if (d_step == 1) copy_vec(len, step_ptr(ps, s_step), pd);
			else copy_vec(len, step_ptr(ps, s_step), step_ptr(pd, d_step));
		}
	}

	template<typename T>
	LMAT_ENSURE_INLINE
	inline void _copy_multicol(index_t m, index_t n,
			const T *ps, index_t src_cs, T *pd, index_t dst_cs)
	{
		for (index_t j = 0; j < n; ++j)
		{
			const T* scol = ps + j * src_cs;
			T *dcol = pd + j * dst_cs;
			copy_vec(m, scol, dcol);
		}
	}

	template<typename T>
	LMAT_ENSURE_INLINE
	inline void _copy_multicol(index_t m, index_t n,
			const T *ps, index_t src_rs, index_t src_cs, T *pd, index_t dst_cs)
	{
		for (index_t j = 0; j < n; ++j)
		{
			const T* scol = ps + j * src_cs;
			T *dcol = pd + j * dst_cs;
			copy_vec(m, step_ptr(scol, src_rs), dcol);
		}
	}

	template<typename T>
	LMAT_ENSURE_INLINE
	inline void _copy_multicol(index_t m, index_t n,
			const T *ps, index_t src_cs, T *pd, index_t dst_rs, index_t dst_cs)
	{
		for (index_t j = 0; j < n; ++j)
		{
			const T* scol = ps + j * src_cs;
			T *dcol = pd + j * dst_cs;
			copy_vec(m, scol, step_ptr(dcol, dst_rs));
		}
	}


	template<typename T>
	LMAT_ENSURE_INLINE
	inline void _copy_multicol(index_t m, index_t n,
			const T *ps, index_t src_rs, index_t src_cs, T *pd, index_t dst_rs, index_t dst_cs)
	{
		for (index_t j = 0; j < n; ++j)
		{
			const T* scol = ps + j * src_cs;
			T *dcol = pd + j * dst_cs;
			copy_vec(m, step_ptr(scol, src_rs), step_ptr(dcol, dst_rs));
		}
	}


	/********************************************
	 *
	 *  Copy from pointer to matrix
	 *
	 ********************************************/

	template<typename T, class DMat, index_t M, index_t N>
	LMAT_ENSURE_INLINE
	inline void copy(const T *ps, IRegularMatrix<DMat, T>& dmat, const matrix_copy_scheme<M, N, cont_level::whole>& sch)
	{
		copy_vec(sch.nelems(), ps, dmat.ptr_data());
	}

	template<typename T, class DMat, index_t M, index_t N>
	inline void copy(const T *ps, IRegularMatrix<DMat, T>& dmat,
			const matrix_copy_scheme<M, N, cont_level::percol>& sch)
	{
		const index_t m = sch.nrows();
		const index_t n = sch.ncolumns();

		T *pd = dmat.ptr_data();

		if (n == 1)
		{
			if (m == 1)
				*pd = *ps;
			else
				copy_vec(m, ps, pd);
		}
		else
		{
			const index_t cs = dmat.col_stride();
			if (m == 1)
				_copy_singlevec(n, ps, pd, cs);
			else
				_copy_multicol(m, n, ps, m, pd, cs);
		}
	}

	template<typename T, class DMat, index_t M, index_t N>
	inline void copy(const T *ps, IRegularMatrix<DMat, T>& dmat,
			const matrix_copy_scheme<1, N, cont_level::percol>& sch)
	{
		const index_t n = sch.ncolumns();

		T *pd = dmat.ptr_data();

		if (n == 1)
		{
			*pd = *ps;
		}
		else
		{
			_copy_singlevec(n, ps, pd, dmat.col_stride());
		}
	}

	template<typename T, class DMat, index_t M, index_t N>
	inline void copy(const T *ps, IRegularMatrix<DMat, T>& dmat,
			const matrix_copy_scheme<M, N, cont_level::none>& sch)
	{
		const index_t m = sch.nrows();
		const index_t n = sch.ncolumns();

		T *pd = dmat.ptr_data();

		if (n == 1)
		{
			if (m == 1)
				*pd = *ps;
			else
			{
				_copy_singlevec(m, ps, pd, dmat.row_stride());
			}
		}
		else
		{
			const index_t cs = dmat.col_stride();
			if (m == 1)
			{
				_copy_singlevec(n, ps, pd, cs);
			}
			else
			{
				const index_t rs = dmat.row_stride();

				if (rs == 1)
					_copy_multicol(m, n, ps, m, pd, cs);
				else
					_copy_multicol(m, n, ps, m, pd, rs, cs);
			}
		}
	}


	/********************************************
	 *
	 *  Copy from matrix to pointer
	 *
	 ********************************************/

	template<typename T, class SMat, index_t M, index_t N>
	LMAT_ENSURE_INLINE
	inline void copy(const IRegularMatrix<SMat, T>& smat, T *pd, const matrix_copy_scheme<M, N, cont_level::whole>& sch)
	{
		copy_vec(sch.nelems(), smat.ptr_data(), pd);
	}

	template<typename T, class SMat, index_t M, index_t N>
	inline void copy(const IRegularMatrix<SMat, T>& smat, T *pd,
			const matrix_copy_scheme<M, N, cont_level::percol>& sch)
	{
		const index_t m = sch.nrows();
		const index_t n = sch.ncolumns();

		const T *ps = smat.ptr_data();

		if (n == 1)
		{
			if (m == 1)
				*pd = *ps;
			else
				copy_vec(m, ps, pd);
		}
		else
		{
			const index_t cs = smat.col_stride();

			if (m == 1)
				_copy_singlevec(n, ps, cs, pd);
			else
				_copy_multicol(m, n, ps, cs, pd, m);
		}
	}

	template<typename T, class SMat, index_t M, index_t N>
	inline void copy(const IRegularMatrix<SMat, T>& smat, T *pd,
			const matrix_copy_scheme<1, N, cont_level::percol>& sch)
	{
		const index_t n = sch.ncolumns();

		const T *ps = smat.ptr_data();

		if (n == 1)
		{
			*pd = *ps;
		}
		else
		{
			const index_t cs = smat.col_stride();
			if (cs == 1)
				copy_vec(n, ps, pd);
			else
				copy_vec(n, step_ptr(ps, cs), pd);
		}
	}

	template<typename T, class SMat, index_t M, index_t N>
	inline void copy(const IRegularMatrix<SMat, T>& smat, T *pd,
			const matrix_copy_scheme<M, N, cont_level::none>& sch)
	{
		const index_t m = sch.nrows();
		const index_t n = sch.ncolumns();

		const T *ps = smat.ptr_data();

		if (n == 1)
		{
			if (m == 1)
			{
				*pd = *ps;
			}
			else
			{
				const index_t rs = smat.row_stride();
				_copy_singlevec(m, ps, rs, pd);
			}
		}
		else
		{
			const index_t cs = smat.col_stride();
			if (m == 1)
			{
				_copy_singlevec(n, ps, cs, pd);
			}
			else
			{
				const index_t rs = smat.row_stride();
				if (rs == 1)
					_copy_multicol(m, n, ps, cs, pd, m);
				else
					_copy_multicol(m, n, ps, rs, cs, pd, m);
			}
		}
	}


	/********************************************
	 *
	 *  Copy from matrix to matrix
	 *
	 ********************************************/

	template<typename T, class SMat, class DMat, index_t M, index_t N>
	LMAT_ENSURE_INLINE
	inline void copy(const IRegularMatrix<SMat, T>& smat, IRegularMatrix<DMat, T>& dmat,
			const matrix_copy_scheme<M, N, cont_level::whole>& sch)
	{
		copy_vec(sch.nelems(), smat.ptr_data(), dmat.ptr_data());
	}

	template<typename T, class SMat, class DMat, index_t M, index_t N>
	inline void copy(const IRegularMatrix<SMat, T>& smat, IRegularMatrix<DMat, T>& dmat,
			const matrix_copy_scheme<M, N, cont_level::percol>& sch)
	{
		const index_t m = sch.nrows();
		const index_t n = sch.ncolumns();

		const T *ps = smat.ptr_data();
		T *pd = dmat.ptr_data();

		if (n == 1)
		{
			if (m == 1)
				*pd = *ps;
			else
				copy_vec(m, ps, pd);
		}
		else
		{
			if (m == 1)
				_copy_singlevec(n, ps, smat.col_stride(), pd, dmat.col_stride());
			else
				_copy_multicol(m, n, ps, smat.col_stride(), pd, dmat.col_stride());
		}
	}


	template<typename T, class SMat, class DMat, index_t M, index_t N>
	inline void copy(const IRegularMatrix<SMat, T>& smat, IRegularMatrix<DMat, T>& dmat,
			const matrix_copy_scheme<1, N, cont_level::percol>& sch)
	{
		const index_t n = sch.ncolumns();

		const T *ps = smat.ptr_data();
		T *pd = dmat.ptr_data();

		if (n == 1)
			*pd = *ps;
		else
			_copy_singlevec(n, ps, smat.col_stride(), pd, dmat.col_stride());
	}


	template<typename T, class SMat, class DMat, index_t M, index_t N>
	inline void copy(const IRegularMatrix<SMat, T>& smat, IRegularMatrix<DMat, T>& dmat,
			const matrix_copy_scheme<M, N, cont_level::none>& sch)
	{
		const index_t m = sch.nrows();
		const index_t n = sch.ncolumns();

		const T *ps = smat.ptr_data();
		T *pd = dmat.ptr_data();

		if (n == 1)
		{
			if (m == 1)
				*pd = *ps;
			else
				_copy_singlevec(m, ps, smat.row_stride(), pd, dmat.row_stride());
		}
		else
		{
			if (m == 1)
				_copy_singlevec(n, ps, smat.col_stride(), pd, dmat.col_stride());
			else
				_copy_multicol(m, n,
						ps, smat.row_stride(), smat.col_stride(),
						pd, dmat.row_stride(), dmat.col_stride());
		}
	}


} }

#endif 
