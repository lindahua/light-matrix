/**
 * @file lapack_syev.h
 *
 * @brief Eigenvalues of Symmetric matrices
 *
 * @author Dahua Lin
 */

#ifdef _MSC_VER
#pragma once
#endif

#ifndef LIGHTMAT_LAPACK_SYEV_H_
#define LIGHTMAT_LAPACK_SYEV_H_

#include <light_mat/linalg/lapack_fwd.h>

extern "C"
{
	void LMAT_LAPACK_NAME(ssyev)( const char* jobz, const char* uplo, const lapack_int* n, float* a,
	            const lapack_int* lda, float* w, float* work, const lapack_int* lwork,
	            lapack_int* info );

	void LMAT_LAPACK_NAME(ssyevd)( const char* jobz, const char* uplo, const lapack_int* n, float* a,
	             const lapack_int* lda, float* w, float* work, const lapack_int* lwork,
	             lapack_int* iwork, const lapack_int* liwork, lapack_int* info );

	void LMAT_LAPACK_NAME(ssyevr)( const char* jobz, const char* range, const char* uplo,
	             const lapack_int* n, float* a, const lapack_int* lda, const float* vl,
	             const float* vu, const lapack_int* il, const lapack_int* iu,
	             const float* abstol, lapack_int* m, float* w, float* z,
	             const lapack_int* ldz, lapack_int* isuppz, float* work,
	             const lapack_int* lwork, lapack_int* iwork, const lapack_int* liwork,
	             lapack_int* info );

	void LMAT_LAPACK_NAME(ssyevx)( const char* jobz, const char* range, const char* uplo,
	             const lapack_int* n, float* a, const lapack_int* lda, const float* vl,
	             const float* vu, const lapack_int* il, const lapack_int* iu,
	             const float* abstol, lapack_int* m, float* w, float* z,
	             const lapack_int* ldz, float* work, const lapack_int* lwork,
	             lapack_int* iwork, lapack_int* ifail, lapack_int* info );

	void LMAT_LAPACK_NAME(dsyev)( const char* jobz, const char* uplo, const lapack_int* n, double* a,
	            const lapack_int* lda, double* w, double* work, const lapack_int* lwork,
	            lapack_int* info );

	void LMAT_LAPACK_NAME(dsyevd)( const char* jobz, const char* uplo, const lapack_int* n, double* a,
	             const lapack_int* lda, double* w, double* work, const lapack_int* lwork,
	             lapack_int* iwork, const lapack_int* liwork, lapack_int* info );

	void LMAT_LAPACK_NAME(dsyevr)( const char* jobz, const char* range, const char* uplo,
	             const lapack_int* n, double* a, const lapack_int* lda, const double* vl,
	             const double* vu, const lapack_int* il, const lapack_int* iu,
	             const double* abstol, lapack_int* m, double* w, double* z,
	             const lapack_int* ldz, lapack_int* isuppz, double* work,
	             const lapack_int* lwork, lapack_int* iwork, const lapack_int* liwork,
	             lapack_int* info );

	void LMAT_LAPACK_NAME(dsyevx)( const char* jobz, const char* range, const char* uplo,
	             const lapack_int* n, double* a, const lapack_int* lda, const double* vl,
	             const double* vu, const lapack_int* il, const lapack_int* iu,
	             const double* abstol, lapack_int* m, double* w, double* z,
	             const lapack_int* ldz, double* work, const lapack_int* lwork,
	             lapack_int* iwork, lapack_int* ifail, lapack_int* info );
}


namespace lmat { namespace lapack {


	/********************************************
	 *
	 *  range value range
	 *
	 ********************************************/

	struct eigval_irange
	{
		const index_t ibegin;
		const index_t iend;

		LMAT_ENSURE_INLINE
		eigval_irange(index_t b, index_t e)
		: ibegin(b), iend(e)
		{
			check_arg(b >= 0 && b <= e, "Invalid values for eigval_irange");
		}

		LMAT_ENSURE_INLINE
		index_t num() const
		{
			return iend - ibegin;
		}
	};


	LMAT_ENSURE_INLINE
	inline eigval_irange evr_I(index_t b, index_t e)
	{
		return eigval_irange(b, e);
	}


	template<typename T>
	struct eigval_vrange
	{
		const T lbound;
		const T ubound;  // lbound < e <= ubound

		LMAT_ENSURE_INLINE
		eigval_vrange(T l, T u)
		: lbound(l), ubound(u)
		{
			check_arg(l <= u, "Invalid values for eigval_vrange");
		}
	};


	template<typename T>
	LMAT_ENSURE_INLINE
	inline eigval_vrange<T> evr_V(T l, T u)
	{
		return eigval_vrange<T>(l, u);
	}


	/********************************************
	 *
	 *  SYEV
	 *
	 ********************************************/

	namespace internal
	{
		template<class A, class W>
		inline void _syev(IRegularMatrix<A, float>& a, IRegularMatrix<W, float>& w, char jobz, char uplo)
		{
			LMAT_CHECK_WHOLE_CONT(W)
			LMAT_CHECK_PERCOL_CONT(A)

			lapack_int n = (lapack_int)a.nrows();
			lapack_int lda = (lapack_int)a.col_stride();
			float lwork_opt = 0;
			lapack_int lwork = -1;
			lapack_int info = 0;

			LMAT_CALL_LAPACK(ssyev, (&jobz, &uplo, &n, a.ptr_data(), &lda,
					w.ptr_data(), &lwork_opt, &lwork, &info));

			lwork = (lapack_int)lwork_opt;

			dense_col<float> ws((index_t)lwork);

			LMAT_CALL_LAPACK(ssyev, (&jobz, &uplo, &n, a.ptr_data(), &lda,
					w.ptr_data(), ws.ptr_data(), &lwork, &info));
		}

		template<class A, class W>
		inline void _syev(IRegularMatrix<A, double>& a, IRegularMatrix<W, double>& w, char jobz, char uplo)
		{
			LMAT_CHECK_WHOLE_CONT(W)
			LMAT_CHECK_PERCOL_CONT(A)

			lapack_int n = (lapack_int)a.nrows();
			lapack_int lda = (lapack_int)a.col_stride();
			double lwork_opt = 0;
			lapack_int lwork = -1;
			lapack_int info = 0;

			LMAT_CALL_LAPACK(dsyev, (&jobz, &uplo, &n, a.ptr_data(), &lda,
					w.ptr_data(), &lwork_opt, &lwork, &info));

			lwork = (lapack_int)lwork_opt;

			dense_col<double> ws((index_t)lwork);

			LMAT_CALL_LAPACK(dsyev, (&jobz, &uplo, &n, a.ptr_data(), &lda,
					w.ptr_data(), ws.ptr_data(), &lwork, &info));
		}
	}


	template<class A, class W>
	inline void syev(const IMatrixXpr<A, float>& a, IRegularMatrix<W, float>& w, char uplo='L')
	{
		LMAT_CHECK_WHOLE_CONT(W)

		index_t n = a.nrows();
		LMAT_CHECK_DIMS(a.ncolumns() == n)
		w.require_size(n, 1);

		dense_matrix<float> a_(a);
		internal::_syev(a_, w, 'N', uplo);
	}

	template<class A, class W, class V>
	inline void syev(const IMatrixXpr<A, float>& a, IRegularMatrix<W, float>& w, IRegularMatrix<V, float>& v, char uplo='L')
	{
		LMAT_CHECK_WHOLE_CONT(W)
		LMAT_CHECK_PERCOL_CONT(V)

		index_t n = a.nrows();
		LMAT_CHECK_DIMS(a.ncolumns() == n)
		w.require_size(n, 1);
		v.require_size(n, n);

		dense_matrix<float> a_(a);
		internal::_syev(a_, w, 'V', uplo);

		copy(a_.derived(), v.derived());
	}

	template<class A, class W>
	inline void syev(const IMatrixXpr<A, double>& a, IRegularMatrix<W, double>& w, char uplo='L')
	{
		LMAT_CHECK_WHOLE_CONT(W)

		index_t n = a.nrows();
		LMAT_CHECK_DIMS(a.ncolumns() == n)
		w.require_size(n, 1);

		dense_matrix<double> a_(a);
		internal::_syev(a_, w, 'N', uplo);
	}

	template<class A, class W, class V>
	inline void syev(const IMatrixXpr<A, double>& a, IRegularMatrix<W, double>& w, IRegularMatrix<V, double>& v, char uplo='L')
	{
		LMAT_CHECK_WHOLE_CONT(W)
		LMAT_CHECK_PERCOL_CONT(V)

		index_t n = a.nrows();
		LMAT_CHECK_DIMS(a.ncolumns() == n)
		w.require_size(n, 1);
		v.require_size(n, n);

		dense_matrix<double> a_(a);
		internal::_syev(a_, w, 'V', uplo);

		copy(a_.derived(), v.derived());
	}


	/********************************************
	 *
	 *  SYEVD
	 *
	 ********************************************/

	namespace internal
	{
		template<class A, class W>
		inline void _syevd(IRegularMatrix<A, float>& a, IRegularMatrix<W, float>& w, char jobz, char uplo)
		{
			LMAT_CHECK_WHOLE_CONT(W)
			LMAT_CHECK_PERCOL_CONT(A)

			lapack_int n = (lapack_int)a.nrows();
			lapack_int lda = (lapack_int)a.col_stride();
			float lwork_opt = 0;
			lapack_int lwork = -1;
			lapack_int liwork_opt = 0;
			lapack_int liwork = -1;
			lapack_int info = 0;

			LMAT_CALL_LAPACK(ssyevd, (&jobz, &uplo, &n, a.ptr_data(), &lda, w.ptr_data(),
					&lwork_opt, &lwork, &liwork_opt, &liwork, &info));

			lwork = (lapack_int)lwork_opt;
			liwork = liwork_opt;

			dense_col<float> ws((index_t)lwork);
			dense_col<lapack_int> iws((index_t)liwork);

			LMAT_CALL_LAPACK(ssyevd, (&jobz, &uplo, &n, a.ptr_data(), &lda, w.ptr_data(),
					ws.ptr_data(), &lwork, iws.ptr_data(), &liwork, &info));
		}

		template<class A, class W>
		inline void _syevd(IRegularMatrix<A, double>& a, IRegularMatrix<W, double>& w, char jobz, char uplo)
		{
			LMAT_CHECK_WHOLE_CONT(W)
			LMAT_CHECK_PERCOL_CONT(A)

			lapack_int n = (lapack_int)a.nrows();
			lapack_int lda = (lapack_int)a.col_stride();
			double lwork_opt = 0;
			lapack_int lwork = -1;
			lapack_int liwork_opt = 0;
			lapack_int liwork = -1;
			lapack_int info = 0;

			LMAT_CALL_LAPACK(dsyevd, (&jobz, &uplo, &n, a.ptr_data(), &lda, w.ptr_data(),
					&lwork_opt, &lwork, &liwork_opt, &liwork, &info));

			lwork = (lapack_int)lwork_opt;
			liwork = liwork_opt;

			dense_col<double> ws((index_t)lwork);
			dense_col<lapack_int> iws((index_t)liwork);

			LMAT_CALL_LAPACK(dsyevd, (&jobz, &uplo, &n, a.ptr_data(), &lda, w.ptr_data(),
					ws.ptr_data(), &lwork, iws.ptr_data(), &liwork, &info));
		}
	}


	template<class A, class W>
	inline void syevd(const IMatrixXpr<A, float>& a, IRegularMatrix<W, float>& w, char uplo='L')
	{
		LMAT_CHECK_WHOLE_CONT(W)

		index_t n = a.nrows();
		LMAT_CHECK_DIMS(a.ncolumns() == n)
		w.require_size(n, 1);

		dense_matrix<float> a_(a);
		internal::_syevd(a_, w, 'N', uplo);
	}

	template<class A, class W, class V>
	inline void syevd(const IMatrixXpr<A, float>& a, IRegularMatrix<W, float>& w, IRegularMatrix<V, float>& v, char uplo='L')
	{
		LMAT_CHECK_WHOLE_CONT(W)
		LMAT_CHECK_PERCOL_CONT(V)

		index_t n = a.nrows();
		LMAT_CHECK_DIMS(a.ncolumns() == n)
		w.require_size(n, 1);
		v.require_size(n, n);

		dense_matrix<float> a_(a);
		internal::_syevd(a_, w, 'V', uplo);

		copy(a_.derived(), v.derived());
	}

	template<class A, class W>
	inline void syevd(const IMatrixXpr<A, double>& a, IRegularMatrix<W, double>& w, char uplo='L')
	{
		LMAT_CHECK_WHOLE_CONT(W)

		index_t n = a.nrows();
		LMAT_CHECK_DIMS(a.ncolumns() == n)
		w.require_size(n, 1);

		dense_matrix<double> a_(a);
		internal::_syevd(a_, w, 'N', uplo);
	}

	template<class A, class W, class V>
	inline void syevd(const IMatrixXpr<A, double>& a, IRegularMatrix<W, double>& w, IRegularMatrix<V, double>& v, char uplo='L')
	{
		LMAT_CHECK_WHOLE_CONT(W)
		LMAT_CHECK_PERCOL_CONT(V)

		index_t n = a.nrows();
		LMAT_CHECK_DIMS(a.ncolumns() == n)
		w.require_size(n, 1);
		v.require_size(n, n);

		dense_matrix<double> a_(a);
		internal::_syevd(a_, w, 'V', uplo);

		copy(a_.derived(), v.derived());
	}


	/********************************************
	 *
	 *  SYEVR
	 *
	 ********************************************/

	namespace internal
	{
		template<typename T>
		LMAT_ENSURE_INLINE
		inline void _set_eigval_range(index_t n, whole,
				char &range, T &vl, T &vu, lapack_int &il, lapack_int &iu)
		{
			range = 'A';
			vl = T(0);
			vu = T(0);
			il = 1;
			iu = (lapack_int)n;
		}

		template<typename T>
		LMAT_ENSURE_INLINE
		inline void _set_eigval_range(index_t n, const eigval_irange& irgn,
				char &range, T &vl, T &vu, lapack_int &il, lapack_int &iu)
		{
			check_arg(irgn.iend <= n, "Invalid eigenvalue index range");

			range = 'I';
			vl = T(0);
			vu = T(0);
			il = (lapack_int)(irgn.ibegin + 1);
			iu = (lapack_int)(irgn.iend);
		}

		template<typename T>
		LMAT_ENSURE_INLINE
		inline void _set_eigval_range(index_t n, const eigval_vrange<T>& vrgn,
				char &range, T &vl, T &vu, lapack_int &il, lapack_int &iu)
		{
			range = 'V';
			vl = vrgn.lbound;
			vu = vrgn.ubound;
			il = 0;
			iu = 0;
		}


		template<class A, class W, class V, typename Range>
		inline index_t _syevr(IRegularMatrix<A, float>& a, IRegularMatrix<W, float>& w, IRegularMatrix<V, float>& z,
				char jobz, char uplo, float abstol, const Range& ergn)
		{
			LMAT_CHECK_WHOLE_CONT(W)
			LMAT_CHECK_PERCOL_CONT(A)
			LMAT_CHECK_PERCOL_CONT(V)

			lapack_int m = 0;
			lapack_int n = (lapack_int)a.nrows();
			lapack_int lda = (lapack_int)a.col_stride();
			lapack_int ldz = (lapack_int)z.col_stride();
			if (ldz == 0) ldz = 1;

			float lwork_opt = 0;
			lapack_int lwork = -1;
			lapack_int liwork_opt = 0;
			lapack_int liwork = -1;
			lapack_int info = 0;

			dense_col<lapack_int> isuppz;
			if (jobz == 'V' || jobz == 'v')
				isuppz.require_size(2 * n);

			char range;
			float vl, vu;
			lapack_int il, iu;
			_set_eigval_range(n, ergn, range, vl, vu, il, iu);

			LMAT_CALL_LAPACK(ssyevr, (&jobz, &range, &uplo, &n, a.ptr_data(), &lda, &vl, &vu, &il, &iu, &abstol,
					&m, w.ptr_data(), z.ptr_data(), &ldz, isuppz.ptr_data(),
					&lwork_opt, &lwork, &liwork_opt, &liwork, &info));

			lwork = (lapack_int)lwork_opt;
			liwork = liwork_opt;

			dense_col<float> ws((index_t)lwork);
			dense_col<lapack_int> iws((index_t)liwork);

			LMAT_CALL_LAPACK(ssyevr, (&jobz, &range, &uplo, &n, a.ptr_data(), &lda, &vl, &vu, &il, &iu, &abstol,
					&m, w.ptr_data(), z.ptr_data(), &ldz, isuppz.ptr_data(),
					ws.ptr_data(), &lwork, iws.ptr_data(), &liwork, &info));

			return (index_t)m;
		}


		template<class A, class W, class V, typename Range>
		inline index_t _syevr(IRegularMatrix<A, double>& a, IRegularMatrix<W, double>& w, IRegularMatrix<V, double>& z,
				char jobz, char uplo, double abstol, const Range& ergn)
		{
			LMAT_CHECK_WHOLE_CONT(W)
			LMAT_CHECK_PERCOL_CONT(A)
			LMAT_CHECK_PERCOL_CONT(V)

			lapack_int m = 0;
			lapack_int n = (lapack_int)a.nrows();
			lapack_int lda = (lapack_int)a.col_stride();
			lapack_int ldz = (lapack_int)z.col_stride();
			if (ldz == 0) ldz = 1;

			double lwork_opt = 0;
			lapack_int lwork = -1;
			lapack_int liwork_opt = 0;
			lapack_int liwork = -1;
			lapack_int info = 0;

			dense_col<lapack_int> isuppz;
			if (jobz == 'V' || jobz == 'v')
				isuppz.require_size(2 * n);

			char range;
			double vl, vu;
			lapack_int il, iu;
			_set_eigval_range(n, ergn, range, vl, vu, il, iu);

			LMAT_CALL_LAPACK(dsyevr, (&jobz, &range, &uplo, &n, a.ptr_data(), &lda, &vl, &vu, &il, &iu, &abstol,
					&m, w.ptr_data(), z.ptr_data(), &ldz, isuppz.ptr_data(),
					&lwork_opt, &lwork, &liwork_opt, &liwork, &info));

			lwork = (lapack_int)lwork_opt;
			liwork = liwork_opt;

			dense_col<double> ws((index_t)lwork);
			dense_col<lapack_int> iws((index_t)liwork);
			m = 0;

			LMAT_CALL_LAPACK(dsyevr, (&jobz, &range, &uplo, &n, a.ptr_data(), &lda, &vl, &vu, &il, &iu, &abstol,
					&m, w.ptr_data(), z.ptr_data(), &ldz, isuppz.ptr_data(),
					ws.ptr_data(), &lwork, iws.ptr_data(), &liwork, &info));

			return (index_t)m;
		}

		template<typename T, class A, class W, typename Range>
		inline index_t _syevr_n(const IMatrixXpr<A, T>& a, IRegularMatrix<W, T>& w, index_t ns, const Range& ergn, T abstol, char uplo)
		{
			LMAT_CHECK_WHOLE_CONT(W)

			index_t n = a.nrows();
			LMAT_CHECK_DIMS(a.ncolumns() == n)
			w.require_size(ns, 1);

			dense_matrix<T> a_(a);
			dense_matrix<T> v_;

			index_t ret;

			if (ns == n)
			{
				ret = internal::_syevr(a_, w, v_, 'N', uplo, abstol, ergn);
			}
			else
			{
				dense_col<T> w_(n);
				ret = internal::_syevr(a_, w_, v_, 'N', uplo, abstol, ergn);
				copy_vec(ns, w_.ptr_data(), w.ptr_data());
			}

			return ret;
		}

		template<typename T, class A, class W, class V, typename Range>
		inline index_t _syevr_v(const IMatrixXpr<A, T>& a, IRegularMatrix<W, T>& w, IRegularMatrix<V, T>& v,
				index_t ns, const Range& ergn, T abstol, char uplo)
		{
			LMAT_CHECK_WHOLE_CONT(W)
			LMAT_CHECK_PERCOL_CONT(V)

			index_t n = a.nrows();
			LMAT_CHECK_DIMS(a.ncolumns() == n)
			w.require_size(ns, 1);
			v.require_size(n, ns);

			dense_matrix<T> a_(a);

			index_t ret;

			if (ns == n)
			{
				ret = internal::_syevr(a_, w, v, 'V', uplo, abstol, ergn);
			}
			else
			{
				dense_col<T> w_(n);
				ret = internal::_syevr(a_, w_, v, 'V', uplo, abstol, ergn);
				copy_vec(ns, w_.ptr_data(), w.ptr_data());
			}

			return ret;
		}

	}


	template<class A, class W>
	inline index_t syevr(const IMatrixXpr<A, float>& a, IRegularMatrix<W, float>& w,
			float abstol=0.0f, char uplo='L')
	{
		return internal::_syevr_n(a, w, a.nrows(), whole(), abstol, uplo);
	}

	template<class A, class W>
	inline index_t syevr(const IMatrixXpr<A, float>& a, IRegularMatrix<W, float>& w,
			const eigval_irange& ergn, float abstol=0.0f, char uplo='L')
	{
		return internal::_syevr_n(a, w, ergn.num(), ergn, abstol, uplo);
	}

	template<class A, class W>
	inline index_t syevr(const IMatrixXpr<A, float>& a, IRegularMatrix<W, float>& w,
			const eigval_vrange<float>& ergn, float abstol=0.0f, char uplo='L')
	{
		return internal::_syevr_n(a, w, a.nrows(), ergn, abstol, uplo);
	}

	template<class A, class W, class V>
	inline index_t syevr(const IMatrixXpr<A, float>& a, IRegularMatrix<W, float>& w, IRegularMatrix<V, float>& v,
			float abstol=0.0f, char uplo='L')
	{
		return internal::_syevr_v(a, w, v, a.nrows(), whole(), abstol, uplo);
	}

	template<class A, class W, class V>
	inline index_t syevr(const IMatrixXpr<A, float>& a, IRegularMatrix<W, float>& w, IRegularMatrix<V, float>& v,
			const eigval_irange& ergn, float abstol=0.0f, char uplo='L')
	{
		return internal::_syevr_v(a, w, v, ergn.num(), ergn, abstol, uplo);
	}

	template<class A, class W, class V>
	inline index_t syevr(const IMatrixXpr<A, float>& a, IRegularMatrix<W, float>& w, IRegularMatrix<V, float>& v,
			const eigval_vrange<float>& ergn, float abstol=0.0f, char uplo='L')
	{
		return internal::_syevr_v(a, w, v, a.nrows(), ergn, abstol, uplo);
	}


	template<class A, class W>
	inline index_t syevr(const IMatrixXpr<A, double>& a, IRegularMatrix<W, double>& w,
			double abstol=0.0, char uplo='L')
	{
		return internal::_syevr_n(a, w, a.nrows(), whole(), abstol, uplo);
	}

	template<class A, class W>
	inline index_t syevr(const IMatrixXpr<A, double>& a, IRegularMatrix<W, double>& w,
			const eigval_irange& ergn, double abstol=0.0, char uplo='L')
	{
		return internal::_syevr_n(a, w, ergn.num(), ergn, abstol, uplo);
	}

	template<class A, class W>
	inline index_t syevr(const IMatrixXpr<A, double>& a, IRegularMatrix<W, double>& w,
			const eigval_vrange<double>& ergn, double abstol=0.0, char uplo='L')
	{
		return internal::_syevr_n(a, w, a.nrows(), ergn, abstol, uplo);
	}

	template<class A, class W, class V>
	inline index_t syevr(const IMatrixXpr<A, double>& a, IRegularMatrix<W, double>& w, IRegularMatrix<V, double>& v,
			double abstol=0.0, char uplo='L')
	{
		return internal::_syevr_v(a, w, v, a.nrows(), whole(), abstol, uplo);
	}

	template<class A, class W, class V>
	inline index_t syevr(const IMatrixXpr<A, double>& a, IRegularMatrix<W, double>& w, IRegularMatrix<V, double>& v,
			const eigval_irange& ergn, double abstol=0.0, char uplo='L')
	{
		return internal::_syevr_v(a, w, v, ergn.num(), ergn, abstol, uplo);
	}

	template<class A, class W, class V>
	inline index_t syevr(const IMatrixXpr<A, double>& a, IRegularMatrix<W, double>& w, IRegularMatrix<V, double>& v,
			const eigval_vrange<double>& ergn, double abstol=0.0, char uplo='L')
	{
		return internal::_syevr_v(a, w, v, a.nrows(), ergn, abstol, uplo);
	}

} }


#endif /* LAPACK_SYEV_H_ */
