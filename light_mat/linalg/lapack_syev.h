/**
 * @file lapack_syev.h
 *
 * @brief Eigenvalues of Symmetric matrices
 *
 * @author Dahua Lin
 */

#ifndef LIGHTMAT_LAPACK_SYEV_H_
#define LIGHTMAT_LAPACK_SYEV_H_

#include <light_mat/linalg/lapack_fwd.h>

extern "C"
{
	void ssyev( const char* jobz, const char* uplo, const lapack_int* n, float* a,
	            const lapack_int* lda, float* w, float* work, const lapack_int* lwork,
	            lapack_int* info );

	void ssyevd( const char* jobz, const char* uplo, const lapack_int* n, float* a,
	             const lapack_int* lda, float* w, float* work, const lapack_int* lwork,
	             lapack_int* iwork, const lapack_int* liwork, lapack_int* info );

	void ssyevr( const char* jobz, const char* range, const char* uplo,
	             const lapack_int* n, float* a, const lapack_int* lda, const float* vl,
	             const float* vu, const lapack_int* il, const lapack_int* iu,
	             const float* abstol, lapack_int* m, float* w, float* z,
	             const lapack_int* ldz, lapack_int* isuppz, float* work,
	             const lapack_int* lwork, lapack_int* iwork, const lapack_int* liwork,
	             lapack_int* info );

	void ssyevx( const char* jobz, const char* range, const char* uplo,
	             const lapack_int* n, float* a, const lapack_int* lda, const float* vl,
	             const float* vu, const lapack_int* il, const lapack_int* iu,
	             const float* abstol, lapack_int* m, float* w, float* z,
	             const lapack_int* ldz, float* work, const lapack_int* lwork,
	             lapack_int* iwork, lapack_int* ifail, lapack_int* info );

	void dsyev( const char* jobz, const char* uplo, const lapack_int* n, double* a,
	            const lapack_int* lda, double* w, double* work, const lapack_int* lwork,
	            lapack_int* info );

	void dsyevd( const char* jobz, const char* uplo, const lapack_int* n, double* a,
	             const lapack_int* lda, double* w, double* work, const lapack_int* lwork,
	             lapack_int* iwork, const lapack_int* liwork, lapack_int* info );

	void dsyevr( const char* jobz, const char* range, const char* uplo,
	             const lapack_int* n, double* a, const lapack_int* lda, const double* vl,
	             const double* vu, const lapack_int* il, const lapack_int* iu,
	             const double* abstol, lapack_int* m, double* w, double* z,
	             const lapack_int* ldz, lapack_int* isuppz, double* work,
	             const lapack_int* lwork, lapack_int* iwork, const lapack_int* liwork,
	             lapack_int* info );

	void dsyevx( const char* jobz, const char* range, const char* uplo,
	             const lapack_int* n, double* a, const lapack_int* lda, const double* vl,
	             const double* vu, const lapack_int* il, const lapack_int* iu,
	             const double* abstol, lapack_int* m, double* w, double* z,
	             const lapack_int* ldz, double* work, const lapack_int* lwork,
	             lapack_int* iwork, lapack_int* ifail, lapack_int* info );
}


namespace lmat { namespace lapack {

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
		LMAT_CHECK_DIMS(a.ncolumns() == n && w.nelems() == n)

		dense_matrix<float> a_(a);
		internal::_syev(a_, w, 'N', uplo);
	}

	template<class A, class W, class V>
	inline void syev(const IMatrixXpr<A, float>& a, IRegularMatrix<W, float>& w, IRegularMatrix<V, float>& v, char uplo='L')
	{
		LMAT_CHECK_WHOLE_CONT(W)
		LMAT_CHECK_PERCOL_CONT(V)

		index_t n = a.nrows();
		LMAT_CHECK_DIMS(a.ncolumns() == n && w.nelems() == n && v.nrows() == n && v.ncolumns() == n)

		dense_matrix<float> a_(a);
		internal::_syev(a_, w, 'V', uplo);

		copy(a_.derived(), v.derived());
	}

	template<class A, class W>
	inline void syev(const IMatrixXpr<A, double>& a, IRegularMatrix<W, double>& w, char uplo='L')
	{
		LMAT_CHECK_WHOLE_CONT(W)

		index_t n = a.nrows();
		LMAT_CHECK_DIMS(a.ncolumns() == n && w.nelems() == n)

		dense_matrix<double> a_(a);
		internal::_syev(a_, w, 'N', uplo);
	}

	template<class A, class W, class V>
	inline void syev(const IMatrixXpr<A, double>& a, IRegularMatrix<W, double>& w, IRegularMatrix<V, double>& v, char uplo='L')
	{
		LMAT_CHECK_WHOLE_CONT(W)
		LMAT_CHECK_PERCOL_CONT(V)

		index_t n = a.nrows();
		LMAT_CHECK_DIMS(a.ncolumns() == n && w.nelems() == n && v.nrows() == n && v.ncolumns() == n)

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
		LMAT_CHECK_DIMS(a.ncolumns() == n && w.nelems() == n)

		dense_matrix<float> a_(a);
		internal::_syevd(a_, w, 'N', uplo);
	}

	template<class A, class W, class V>
	inline void syevd(const IMatrixXpr<A, float>& a, IRegularMatrix<W, float>& w, IRegularMatrix<V, float>& v, char uplo='L')
	{
		LMAT_CHECK_WHOLE_CONT(W)
		LMAT_CHECK_PERCOL_CONT(V)

		index_t n = a.nrows();
		LMAT_CHECK_DIMS(a.ncolumns() == n && w.nelems() == n && v.nrows() == n && v.ncolumns() == n)

		dense_matrix<float> a_(a);
		internal::_syevd(a_, w, 'V', uplo);

		copy(a_.derived(), v.derived());
	}

	template<class A, class W>
	inline void syevd(const IMatrixXpr<A, double>& a, IRegularMatrix<W, double>& w, char uplo='L')
	{
		LMAT_CHECK_WHOLE_CONT(W)

		index_t n = a.nrows();
		LMAT_CHECK_DIMS(a.ncolumns() == n && w.nelems() == n)

		dense_matrix<double> a_(a);
		internal::_syevd(a_, w, 'N', uplo);
	}

	template<class A, class W, class V>
	inline void syevd(const IMatrixXpr<A, double>& a, IRegularMatrix<W, double>& w, IRegularMatrix<V, double>& v, char uplo='L')
	{
		LMAT_CHECK_WHOLE_CONT(W)
		LMAT_CHECK_PERCOL_CONT(V)

		index_t n = a.nrows();
		LMAT_CHECK_DIMS(a.ncolumns() == n && w.nelems() == n && v.nrows() == n && v.ncolumns() == n)

		dense_matrix<double> a_(a);
		internal::_syevd(a_, w, 'V', uplo);

		copy(a_.derived(), v.derived());
	}



} }


#endif /* LAPACK_SYEV_H_ */
