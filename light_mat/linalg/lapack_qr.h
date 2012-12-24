/**
 * @file lapack_qr.h
 *
 * @brief QR factorization
 *
 * @author Dahua Lin
 */

#ifndef LIGHTMAT_LAPACK_QR_H_
#define LIGHTMAT_LAPACK_QR_H_

#include <light_mat/linalg/lapack_fwd.h>

extern "C"
{
	void LMAT_LAPACK_NAME(sgeqrf)( const lapack_int* m, const lapack_int* n, float* a, const lapack_int* lda,
	             float* tau, float* work, const lapack_int* lwork, lapack_int* info );

	void LMAT_LAPACK_NAME(dgeqrf)( const lapack_int* m, const lapack_int* n, double* a, const lapack_int* lda,
	             double* tau, double* work, const lapack_int* lwork, lapack_int* info );

	void LMAT_LAPACK_NAME(sgeqp3)( const lapack_int* m, const lapack_int* n, float* a, const lapack_int* lda,
	             lapack_int* jpvt, float* tau, float* work, const lapack_int* lwork,
	             lapack_int* info );

	void LMAT_LAPACK_NAME(dgeqp3)( const lapack_int* m, const lapack_int* n, double* a, const lapack_int* lda,
	             lapack_int* jpvt, double* tau, double* work, const lapack_int* lwork,
	             lapack_int* info );

	void LMAT_LAPACK_NAME(sorgqr)( const lapack_int* m, const lapack_int* n, const lapack_int* k, float* a,
	             const lapack_int* lda, const float* tau, float* work,
	             const lapack_int* lwork, lapack_int* info );

	void LMAT_LAPACK_NAME(dorgqr)( const lapack_int* m, const lapack_int* n, const lapack_int* k, double* a,
	             const lapack_int* lda, const double* tau, double* work,
	             const lapack_int* lwork, lapack_int* info );

	void LMAT_LAPACK_NAME(sormqr)( const char* side, const char* trans, const lapack_int* m,
	             const lapack_int* n, const lapack_int* k, const float* a,
	             const lapack_int* lda, const float* tau, float* c,
	             const lapack_int* ldc, float* work, const lapack_int* lwork,
	             lapack_int* info );

	void LMAT_LAPACK_NAME(dormqr)( const char* side, const char* trans, const lapack_int* m,
	             const lapack_int* n, const lapack_int* k, const double* a,
	             const lapack_int* lda, const double* tau, double* c,
	             const lapack_int* ldc, double* work, const lapack_int* lwork,
	             lapack_int* info );

	void LMAT_BLAS_NAME(strsm)(const char *side, const char *uplo, const char *transa, const char *diag,
	           const blas_int *m, const blas_int *n, const float *alpha, const float *a, const blas_int *lda,
	           float *b, const blas_int *ldb);

	void LMAT_BLAS_NAME(dtrsm)(const char *side, const char *uplo, const char *transa, const char *diag,
	           const blas_int *m, const blas_int *n, const double *alpha, const double *a, const blas_int *lda,
	           double *b, const blas_int *ldb);
}


namespace lmat { namespace lapack {

	// forward declaration

	template<typename T> class qr_fac;

	/********************************************
	 *
	 *  QR factorization classes
	 *
	 ********************************************/

	template<typename T>
	class qr_base
	{
	public:
		qr_base()
		: m_nrows(0), m_ncols(0), m_a(), m_tau()
		{ }

		bool emtpy() const
		{
			return m_nrows == 0 || m_ncols == 0;
		}

		index_t nrows() const
		{
			return m_nrows;
		}

		index_t ncolumns() const
		{
			return m_ncols;
		}

		const dense_matrix<T>& intern() const
		{
			return m_a;
		}

		const T* tau() const
		{
			return m_tau.ptr_data();
		}

		template<class R>
		void getr(IRegularMatrix<R, T>& r) const
		{
			LMAT_CHECK_PERCOL_CONT(R)
			LMAT_CHECK_DIMS( r.nrows() <= m_nrows && r.ncolumns() == m_ncols )

			index_t mr = r.nrows();

			zero(r);

			if (mr == m_nrows)
			{
				copy_triu(m_a, r);
			}
			else
			{
				copy_triu(m_a(range(0, mr), whole()), r);
			}
		}

	protected:
		bool check_multq_dims(char side, index_t mx, index_t nx) const
		{
			index_t r = (side == 'L' || side == 'l') ? mx : nx;
			return m_nrows == r;
		}

		template<class Mat>
		void set_mat(const IMatrixXpr<Mat, T>& mat)
		{
			m_a = mat.derived();
			m_nrows = m_a.nrows();
			m_ncols = m_a.ncolumns();

			index_t ltau = math::max(1, math::min(m_nrows, m_ncols));
			m_tau.require_size(ltau);
		}

	protected:
		index_t m_nrows;
		index_t m_ncols;
		dense_matrix<T> m_a;
		dense_col<T> m_tau;
	};


	template<>
	class qr_fac<float> : public qr_base<float>
	{
	public:
		explicit qr_fac() { }

		template<class Mat>
		explicit qr_fac(const IMatrixXpr<Mat, float>& mat)
		{
			set(mat);
		}

		template<class Mat>
		void set(const IMatrixXpr<Mat, float>& mat)
		{
			this->set_mat(mat);

			lapack_int m = (lapack_int)(this->m_nrows);
			lapack_int n = (lapack_int)(this->m_ncols);
			lapack_int lda = (lapack_int)(this->m_a.col_stride());

			lapack_int info = 0;
			lapack_int lwork = -1;
			float lwork_opt = 0;

			LMAT_CALL_LAPACK(sgeqrf, (&m, &n, this->m_a.ptr_data(), &lda, this->m_tau.ptr_data(),
					&lwork_opt, &lwork, &info));

			lwork = (lapack_int)lwork_opt;
			dense_col<float> ws((index_t)lwork);

			LMAT_CALL_LAPACK(sgeqrf, (&m, &n, this->m_a.ptr_data(), &lda, this->m_tau.ptr_data(),
					ws.ptr_data(), &lwork, &info));
		}

		template<class Q>
		void getq(IRegularMatrix<Q, float>& q, index_t k=-1) const  // q: m x n (with n <= m)
		{
			LMAT_CHECK_PERCOL_CONT(Q)
			LMAT_CHECK_DIMS( q.nrows() == this->m_nrows && q.ncolumns() <= this->m_nrows )

			index_t kmax = math::min(q.ncolumns(), this->m_ncols);
			LMAT_CHECK_DIMS( k <= kmax )
			if (k < 0) k = kmax;

			const dense_matrix<float>& a = m_a;

			auto av = a(whole(), range(0, k));
			auto qv = q(whole(), range(0, k));
			copy(av, qv);

			lapack_int m = (lapack_int)q.nrows();
			lapack_int n = (lapack_int)q.ncolumns();
			lapack_int k_ = (lapack_int)k;
			lapack_int ldq = (lapack_int)q.col_stride();

			lapack_int lwork = -1;
			float lwork_opt = 0;
			lapack_int info = 0;

			LMAT_CALL_LAPACK(sorgqr, (&m, &n, &k_, q.ptr_data(), &ldq,
					this->m_tau.ptr_data(), &lwork_opt, &lwork, &info));

			lwork = (lapack_int)lwork_opt;
			dense_col<float> ws((index_t)lwork);

			LMAT_CALL_LAPACK(sorgqr, (&m, &n, &k_, q.ptr_data(), &ldq,
					this->m_tau.ptr_data(), ws.ptr_data(), &lwork, &info));
		}


		template<class X>
		void multq_inplace(IRegularMatrix<X, float>& x, char trans='N', char side='L') const
		{
			LMAT_CHECK_PERCOL_CONT(X)
			LMAT_CHECK_DIMS( this->check_multq_dims(side, x.nrows(), x.ncolumns()) )

			lapack_int m = (lapack_int)x.nrows();
			lapack_int n = (lapack_int)x.ncolumns();
			lapack_int k = math::min(this->m_nrows, this->m_ncols);
			lapack_int lda = (lapack_int)this->m_a.col_stride();
			lapack_int ldx = (lapack_int)x.col_stride();

			lapack_int lwork = -1;
			float lwork_opt = 0;
			lapack_int info = 0;

			LMAT_CALL_LAPACK(sormqr, (&side, &trans, &m, &n, &k, this->m_a.ptr_data(), &lda,
					this->m_tau.ptr_data(), x.ptr_data(), &ldx, &lwork_opt, &lwork, &info));

			lwork = (lapack_int)lwork_opt;
			dense_col<float> ws((index_t)lwork);

			LMAT_CALL_LAPACK(sormqr, (&side, &trans, &m, &n, &k, this->m_a.ptr_data(), &lda,
					this->m_tau.ptr_data(), x.ptr_data(), &ldx, ws.ptr_data(), &lwork, &info));
		}

		template<class X, class Y>
		void multq(const IMatrixXpr<X, float>& x, IRegularMatrix<Y, float>& y, char trans='N', char side='L') const
		{
			LMAT_CHECK_PERCOL_CONT(Y)
			y.derived() = x.derived();
			multq_inplace(y, trans, side);
		}

		template<class X>
		void solve_inplace(IRegularMatrix<X, float>& x) const // require: m >= n
		{
			LMAT_CHECK_PERCOL_CONT(X)

			check_arg(this->m_nrows >= this->m_ncols, "QR-solve only applies when m >= n");
			LMAT_CHECK_DIMS( x.nrows() == this->m_nrows );

			multq_inplace(x, 'T', 'L');

			char side = 'L';
			char uplo = 'U';
			char transa = 'N';
			char diag = 'N';

			lapack_int m = (lapack_int)this->m_ncols;
			lapack_int n = (lapack_int)x.ncolumns();
			lapack_int lda = (lapack_int)this->m_a.col_stride();
			lapack_int ldb = (lapack_int)x.col_stride();
			float alpha = 1;

			LMAT_BLAS_NAME(strsm)(&side, &uplo, &transa, &diag,
					&m, &n, &alpha, this->m_a.ptr_data(), &lda, x.ptr_data(), &ldb);
		}

		template<class X, class B>
		void solve(const IMatrixXpr<X, float>& x, IRegularMatrix<B, float>& b) const
		{
			LMAT_CHECK_PERCOL_CONT(B)
			LMAT_CHECK_DIMS( x.nrows() == this->m_nrows && b.nrows() == this->m_ncols
					&& x.ncolumns() == b.ncolumns() );

			dense_matrix<float> x_(x);
			solve_inplace(x_);

			copy(x_(range(0, this->m_ncols), whole()), b.derived());
		}
	};


	template<>
	class qr_fac<double> : public qr_base<double>
	{
	public:
		explicit qr_fac() { }

		template<class Mat>
		explicit qr_fac(const IMatrixXpr<Mat, double>& mat)
		{
			set(mat);
		}

		template<class Mat>
		void set(const IMatrixXpr<Mat, double>& mat)
		{
			this->set_mat(mat);

			lapack_int m = (lapack_int)(this->m_nrows);
			lapack_int n = (lapack_int)(this->m_ncols);
			lapack_int lda = (lapack_int)(this->m_a.col_stride());

			lapack_int info = 0;
			lapack_int lwork = -1;
			double lwork_opt = 0;

			LMAT_CALL_LAPACK(dgeqrf, (&m, &n, this->m_a.ptr_data(), &lda, this->m_tau.ptr_data(),
					&lwork_opt, &lwork, &info));

			lwork = (lapack_int)lwork_opt;
			dense_col<double> ws((index_t)lwork);

			LMAT_CALL_LAPACK(dgeqrf, (&m, &n, this->m_a.ptr_data(), &lda, this->m_tau.ptr_data(),
					ws.ptr_data(), &lwork, &info));
		}

		template<class Q>
		void getq(IRegularMatrix<Q, double>& q, index_t k=-1) const  // q: m x n (with n <= m)
		{
			LMAT_CHECK_PERCOL_CONT(Q)
			LMAT_CHECK_DIMS( q.nrows() == this->m_nrows && q.ncolumns() <= this->m_nrows )

			index_t kmax = math::min(q.ncolumns(), this->m_ncols);
			LMAT_CHECK_DIMS( k <= kmax )
			if (k < 0) k = kmax;

			const dense_matrix<double>& a = m_a;

			auto av = a(whole(), range(0, k));
			auto qv = q(whole(), range(0, k));
			copy(av, qv);

			lapack_int m = (lapack_int)q.nrows();
			lapack_int n = (lapack_int)q.ncolumns();
			lapack_int k_ = (lapack_int)k;
			lapack_int ldq = (lapack_int)q.col_stride();

			lapack_int lwork = -1;
			double lwork_opt = 0;
			lapack_int info = 0;

			LMAT_CALL_LAPACK(dorgqr, (&m, &n, &k_, q.ptr_data(), &ldq,
					this->m_tau.ptr_data(), &lwork_opt, &lwork, &info));

			lwork = (lapack_int)lwork_opt;
			dense_col<double> ws((index_t)lwork);

			LMAT_CALL_LAPACK(dorgqr, (&m, &n, &k_, q.ptr_data(), &ldq,
					this->m_tau.ptr_data(), ws.ptr_data(), &lwork, &info));
		}

		template<class X>
		void multq_inplace(IRegularMatrix<X, double>& x, char trans='N', char side='L') const
		{
			LMAT_CHECK_PERCOL_CONT(X)
			LMAT_CHECK_DIMS( this->check_multq_dims(side, x.nrows(), x.ncolumns()) )

			lapack_int m = (lapack_int)x.nrows();
			lapack_int n = (lapack_int)x.ncolumns();
			lapack_int k = math::min(this->m_nrows, this->m_ncols);
			lapack_int lda = (lapack_int)this->m_a.col_stride();
			lapack_int ldx = (lapack_int)x.col_stride();

			lapack_int lwork = -1;
			double lwork_opt = 0;
			lapack_int info = 0;

			LMAT_CALL_LAPACK(dormqr, (&side, &trans, &m, &n, &k, this->m_a.ptr_data(), &lda,
					this->m_tau.ptr_data(), x.ptr_data(), &ldx, &lwork_opt, &lwork, &info));

			lwork = (lapack_int)lwork_opt;
			dense_col<double> ws((index_t)lwork);

			LMAT_CALL_LAPACK(dormqr, (&side, &trans, &m, &n, &k, this->m_a.ptr_data(), &lda,
					this->m_tau.ptr_data(), x.ptr_data(), &ldx, ws.ptr_data(), &lwork, &info));
		}

		template<class X, class Y>
		void multq(const IMatrixXpr<X, double>& x, IRegularMatrix<Y, double>& y, char trans='N', char side='L') const
		{
			LMAT_CHECK_PERCOL_CONT(Y)
			y.derived() = x.derived();
			multq_inplace(y, trans, side);
		}

		template<class X>
		void solve_inplace(IRegularMatrix<X, double>& x) const // require: m >= n
		{
			LMAT_CHECK_PERCOL_CONT(X)

			check_arg(this->m_nrows >= this->m_ncols, "QR-solve only applies when m >= n");
			LMAT_CHECK_DIMS( x.nrows() == this->m_nrows );

			multq_inplace(x, 'T', 'L');

			char side = 'L';
			char uplo = 'U';
			char transa = 'N';
			char diag = 'N';

			lapack_int m = (lapack_int)this->m_ncols;
			lapack_int n = (lapack_int)x.ncolumns();
			lapack_int lda = (lapack_int)this->m_a.col_stride();
			lapack_int ldb = (lapack_int)x.col_stride();
			double alpha = 1;

			LMAT_BLAS_NAME(dtrsm)(&side, &uplo, &transa, &diag,
					&m, &n, &alpha, this->m_a.ptr_data(), &lda, x.ptr_data(), &ldb);
		}

		template<class X, class B>
		void solve(const IMatrixXpr<X, double>& x, IRegularMatrix<B, double>& b) const
		{
			LMAT_CHECK_PERCOL_CONT(B)
			LMAT_CHECK_DIMS( x.nrows() == this->m_nrows && b.nrows() == this->m_ncols
					&& x.ncolumns() == b.ncolumns() );

			dense_matrix<double> x_(x);
			solve_inplace(x_);

			copy(x_(range(0, this->m_ncols), whole()), b.derived());
		}

	};


} }

#endif /* LAPACK_QR_H_ */
