/**
 * @file lapack_svd.h
 *
 * @brief Singular Value Decomposition
 *
 * @author Dahua Lin
 */

#ifndef LIGHTMAT_LAPACK_SVD_H_
#define LIGHTMAT_LAPACK_SVD_H_

#include <light_mat/linalg/lapack_fwd.h>

extern "C"
{
	void LMAT_LAPACK_NAME(sgesvd)( const char* jobu, const char* jobvt, const lapack_int* m,
	             const lapack_int* n, float* a, const lapack_int* lda, float* s,
	             float* u, const lapack_int* ldu, float* vt, const lapack_int* ldvt,
	             float* work, const lapack_int* lwork, lapack_int* info );

	void LMAT_LAPACK_NAME(sgesdd)( const char* jobz, const lapack_int* m, const lapack_int* n, float* a,
	             const lapack_int* lda, float* s, float* u, const lapack_int* ldu,
	             float* vt, const lapack_int* ldvt, float* work, const lapack_int* lwork,
	             lapack_int* iwork, lapack_int* info );

	void LMAT_LAPACK_NAME(dgesvd)( const char* jobu, const char* jobvt, const lapack_int* m,
	             const lapack_int* n, double* a, const lapack_int* lda, double* s,
	             double* u, const lapack_int* ldu, double* vt, const lapack_int* ldvt,
	             double* work, const lapack_int* lwork, lapack_int* info );

	void LMAT_LAPACK_NAME(dgesdd)( const char* jobz, const lapack_int* m, const lapack_int* n, double* a,
	             const lapack_int* lda, double* s, double* u, const lapack_int* ldu,
	             double* vt, const lapack_int* ldvt, double* work,
	             const lapack_int* lwork, lapack_int* iwork, lapack_int* info );
}


namespace lmat { namespace lapack {

	namespace internal
	{
		LMAT_ENSURE_INLINE
		inline char check_svd_jobc(char c)
		{
			if (c == 'A' || c == 'a') return 'A';
			else if (c == 'S' || c == 's') return 'S';
			else if (c == 'N' || c == 'n') return 'n';
			else
				throw invalid_argument("Invalid character for SVD job.");
		}
	}


	/********************************************
	 *
	 *  GESVD
	 *
	 ********************************************/

	namespace internal
	{
		template<class A, class S, class U, class VT>
		inline void _gesvd(IRegularMatrix<A, float>& a, IRegularMatrix<S, float>& s,
				IRegularMatrix<U, float>& u, IRegularMatrix<VT, float>& vt,
				char jobu, char jobvt)
		{
			lapack_int m = (lapack_int)a.nrows();
			lapack_int n = (lapack_int)a.ncolumns();

			lapack_int lda = (lapack_int)a.col_stride();
			lapack_int ldu = (lapack_int)u.col_stride();
			lapack_int ldvt = (lapack_int)vt.col_stride();

			if (ldu <= 0) ldu = 1;
			if (ldvt <= 0) ldvt = 1;

			float lwork_opt = 0;
			lapack_int lwork = -1;
			lapack_int info = 0;

			LMAT_CALL_LAPACK(sgesvd, (&jobu, &jobvt, &m, &n, a.ptr_data(), &lda,
					s.ptr_data(), u.ptr_data(), &ldu, vt.ptr_data(), &ldvt,
					&lwork_opt, &lwork, &info));

			lwork = (lapack_int)lwork_opt;
			dense_col<float> ws((index_t)lwork);

			LMAT_CALL_LAPACK(sgesvd, (&jobu, &jobvt, &m, &n, a.ptr_data(), &lda,
					s.ptr_data(), u.ptr_data(), &ldu, vt.ptr_data(), &ldvt,
					ws.ptr_data(), &lwork, &info));
		}

		template<class A, class S, class U, class VT>
		inline void _gesvd(IRegularMatrix<A, double>& a, IRegularMatrix<S, double>& s,
				IRegularMatrix<U, double>& u, IRegularMatrix<VT, double>& vt,
				char jobu, char jobvt)
		{
			lapack_int m = (lapack_int)a.nrows();
			lapack_int n = (lapack_int)a.ncolumns();

			lapack_int lda = (lapack_int)a.col_stride();
			lapack_int ldu = (lapack_int)u.col_stride();
			lapack_int ldvt = (lapack_int)vt.col_stride();

			if (ldu <= 0) ldu = 1;
			if (ldvt <= 0) ldvt = 1;

			double lwork_opt = 0;
			lapack_int lwork = -1;
			lapack_int info = 0;

			LMAT_CALL_LAPACK(dgesvd, (&jobu, &jobvt, &m, &n, a.ptr_data(), &lda,
					s.ptr_data(), u.ptr_data(), &ldu, vt.ptr_data(), &ldvt,
					&lwork_opt, &lwork, &info));

			lwork = (lapack_int)lwork_opt;
			dense_col<double> ws((index_t)lwork);

			LMAT_CALL_LAPACK(dgesvd, (&jobu, &jobvt, &m, &n, a.ptr_data(), &lda,
					s.ptr_data(), u.ptr_data(), &ldu, vt.ptr_data(), &ldvt,
					ws.ptr_data(), &lwork, &info));
		}


		template<typename T, class A, class S>
		inline void _gesvd_(const IMatrixXpr<A, T>& a, IRegularMatrix<S, T>& s)
		{
			LMAT_CHECK_WHOLE_CONT(S)

			dense_matrix<T> a_(a);

			index_t m = a_.nrows();
			index_t n = a_.ncolumns();
			index_t rk = math::min(m, n);

			s.require_size(rk, 1);

			dense_matrix<T> u;
			dense_matrix<T> vt;
			_gesvd(a_, s, u, vt, 'N', 'N');
		}

		template<typename T, class A, class S, class U, class VT>
		inline void _gesvd_(const IMatrixXpr<A, T>& a, IRegularMatrix<S, T>& s,
				IRegularMatrix<U, T>& u, IRegularMatrix<VT, T>& vt,
				char jobu, char jobvt)
		{
			LMAT_CHECK_PERCOL_CONT(U)
			LMAT_CHECK_PERCOL_CONT(VT)
			LMAT_CHECK_WHOLE_CONT(S)

			jobu = check_svd_jobc(jobu);
			jobvt = check_svd_jobc(jobvt);

			dense_matrix<T> a_(a);

			index_t m = a_.nrows();
			index_t n = a_.ncolumns();
			index_t rk = math::min(m, n);

			s.require_size(rk, 1);

			if (jobu == 'A')
				u.require_size(m, m);
			else if (jobu == 'S')
				u.require_size(m, rk);

			if (jobvt == 'A')
				vt.require_size(n, n);
			else if (jobvt == 'S')
				vt.require_size(rk, n);

			_gesvd(a_, s, u, vt, jobu, jobvt);
		}

	}


	template<class A, class S>
	inline void gesvd(const IMatrixXpr<A, float>& a, IRegularMatrix<S, float>& s)
	{
		internal::_gesvd_(a, s);
	}

	template<class A, class S>
	inline void gesvd(const IMatrixXpr<A, double>& a, IRegularMatrix<S, double>& s)
	{
		internal::_gesvd_(a, s);
	}

	template<class A, class S, class U, class VT>
	inline void gesvd(const IMatrixXpr<A, float>& a, IRegularMatrix<S, float>& s,
			IRegularMatrix<U, float>& u, IRegularMatrix<VT, float>& vt, char jobu ='A', char jobvt='A')
	{
		internal::_gesvd_(a, s, u, vt, jobu, jobvt);
	}

	template<class A, class S, class U, class VT>
	inline void gesvd(const IMatrixXpr<A, double>& a, IRegularMatrix<S, double>& s,
			IRegularMatrix<U, double>& u, IRegularMatrix<VT, double>& vt, char jobu ='A', char jobvt='A')
	{
		internal::_gesvd_(a, s, u, vt, jobu, jobvt);
	}


	/********************************************
	 *
	 *  GESDD
	 *
	 ********************************************/

	namespace internal
	{
		template<class A, class S, class U, class VT>
		inline void _gesdd(IRegularMatrix<A, float>& a, IRegularMatrix<S, float>& s,
				IRegularMatrix<U, float>& u, IRegularMatrix<VT, float>& vt, char jobz)
		{
			lapack_int m = (lapack_int)a.nrows();
			lapack_int n = (lapack_int)a.ncolumns();

			lapack_int lda = (lapack_int)a.col_stride();
			lapack_int ldu = (lapack_int)u.col_stride();
			lapack_int ldvt = (lapack_int)vt.col_stride();

			if (ldu <= 0) ldu = 1;
			if (ldvt <= 0) ldvt = 1;

			float lwork_opt = 0;
			lapack_int lwork = -1;
			lapack_int info = 0;

			index_t iws_len = 8 * math::max(1, math::min(m, n));
			dense_col<lapack_int> iws(iws_len);

			LMAT_CALL_LAPACK(sgesdd, (&jobz, &m, &n, a.ptr_data(), &lda, s.ptr_data(), u.ptr_data(), &ldu,
					vt.ptr_data(), &ldvt, &lwork_opt, &lwork, iws.ptr_data(), &info));

			lwork = (lapack_int)lwork_opt;
			dense_col<float> ws((index_t)lwork);

			LMAT_CALL_LAPACK(sgesdd, (&jobz, &m, &n, a.ptr_data(), &lda, s.ptr_data(), u.ptr_data(), &ldu,
					vt.ptr_data(), &ldvt, ws.ptr_data(), &lwork, iws.ptr_data(), &info));
		}


		template<class A, class S, class U, class VT>
		inline void _gesdd(IRegularMatrix<A, double>& a, IRegularMatrix<S, double>& s,
				IRegularMatrix<U, double>& u, IRegularMatrix<VT, double>& vt, char jobz)
		{
			lapack_int m = (lapack_int)a.nrows();
			lapack_int n = (lapack_int)a.ncolumns();

			lapack_int lda = (lapack_int)a.col_stride();
			lapack_int ldu = (lapack_int)u.col_stride();
			lapack_int ldvt = (lapack_int)vt.col_stride();

			if (ldu <= 0) ldu = 1;
			if (ldvt <= 0) ldvt = 1;

			double lwork_opt = 0;
			lapack_int lwork = -1;
			lapack_int info = 0;

			index_t iws_len = 8 * math::max(1, math::min(m, n));
			dense_col<lapack_int> iws(iws_len);

			LMAT_CALL_LAPACK(dgesdd, (&jobz, &m, &n, a.ptr_data(), &lda, s.ptr_data(), u.ptr_data(), &ldu,
					vt.ptr_data(), &ldvt, &lwork_opt, &lwork, iws.ptr_data(), &info));

			lwork = (lapack_int)lwork_opt;
			dense_col<double> ws((index_t)lwork);

			LMAT_CALL_LAPACK(dgesdd, (&jobz, &m, &n, a.ptr_data(), &lda, s.ptr_data(), u.ptr_data(), &ldu,
					vt.ptr_data(), &ldvt, ws.ptr_data(), &lwork, iws.ptr_data(), &info));
		}


		template<typename T, class A, class S>
		inline void _gesdd_(const IMatrixXpr<A, T>& a, IRegularMatrix<S, T>& s)
		{
			LMAT_CHECK_WHOLE_CONT(S)

			dense_matrix<T> a_(a);

			index_t m = a_.nrows();
			index_t n = a_.ncolumns();
			index_t rk = math::min(m, n);

			s.require_size(rk, 1);

			dense_matrix<T> u;
			dense_matrix<T> vt;
			_gesdd(a_, s, u, vt, 'N');
		}

		template<typename T, class A, class S, class U, class VT>
		inline void _gesdd_(const IMatrixXpr<A, T>& a, IRegularMatrix<S, T>& s,
				IRegularMatrix<U, T>& u, IRegularMatrix<VT, T>& vt, char jobz)
		{
			LMAT_CHECK_PERCOL_CONT(U)
			LMAT_CHECK_PERCOL_CONT(VT)
			LMAT_CHECK_WHOLE_CONT(S)

			jobz = check_svd_jobc(jobz);

			dense_matrix<T> a_(a);

			index_t m = a_.nrows();
			index_t n = a_.ncolumns();
			index_t rk = math::min(m, n);

			s.require_size(rk, 1);

			if (jobz == 'A')
			{
				u.require_size(m, m);
				vt.require_size(n, n);
			}
			else if (jobz == 'S')
			{
				u.require_size(m, rk);
				vt.require_size(rk, n);
			}

			_gesdd(a_, s, u, vt, jobz);
		}

	}


	template<class A, class S>
	inline void gesdd(const IMatrixXpr<A, float>& a, IRegularMatrix<S, float>& s)
	{
		internal::_gesdd_(a, s);
	}

	template<class A, class S>
	inline void gesdd(const IMatrixXpr<A, double>& a, IRegularMatrix<S, double>& s)
	{
		internal::_gesdd_(a, s);
	}

	template<class A, class S, class U, class VT>
	inline void gesdd(const IMatrixXpr<A, float>& a, IRegularMatrix<S, float>& s,
			IRegularMatrix<U, float>& u, IRegularMatrix<VT, float>& vt, char jobz ='A')
	{
		internal::_gesdd_(a, s, u, vt, jobz);
	}

	template<class A, class S, class U, class VT>
	inline void gesdd(const IMatrixXpr<A, double>& a, IRegularMatrix<S, double>& s,
			IRegularMatrix<U, double>& u, IRegularMatrix<VT, double>& vt, char jobz='A')
	{
		internal::_gesdd_(a, s, u, vt, jobz);
	}


} }


#endif /* LAPACK_SVD_H_ */
