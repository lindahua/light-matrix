/**
 * @file lapack_qr.h
 *
 * @brief QR factorization
 *
 * @author Dahua Lin
 */

#ifndef LIGHTMAT_LAPACK_QR_H_
#define LIGHTMAT_LAPACK_QR_H_

#include "internal/linalg_aux.h"
#include <light_mat/math/math_base.h>

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
		: m_nrows(0), m_ncols(0), m_a()
		, m_tau(), m_ws(), m_ws_size(0)
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
			LMAT_CHECK_DIMS( r.ncolumns() == m_ncols )

			copy_triu(m_a, r);
		}

	protected:
		template<class Mat>
		void set_mat(const IMatrixXpr<Mat, T>& mat)
		{
			m_a = mat.derived();
			m_nrows = m_a.nrows();
			m_ncols = m_a.ncolumns();

			index_t ltau = math::max(1, math::min(m_nrows, m_ncols));
			m_tau.require_size(ltau);

			index_t lws = 64 * math::max(m_nrows, m_ncols);
			if (lws != m_ws_size)
			{
				m_ws.require_size(lws);
				m_ws_size = lws;
			}
		}

	protected:
		index_t m_nrows;
		index_t m_ncols;
		dense_matrix<T> m_a;
		dense_col<T> m_tau;

		dense_col<T> m_ws;
		index_t m_ws_size;
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
			lapack_int lwork = (lapack_int)(this->m_ws_size);
			lapack_int info = 0;

			LMAT_LAPACK_NAME(sgeqrf)(&m, &n, this->m_a.ptr_data(), &lda, this->m_tau.ptr_data(),
					this->m_ws.ptr_data(), &lwork, &info);
			LMAT_CHECK_LAPACK_INFO(sgeqrf, info);
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
			lapack_int lwork = (lapack_int)this->m_ws_size;
			float *pws = const_cast<float*>(m_ws.ptr_data());
			lapack_int info;

			LMAT_LAPACK_NAME(sorgqr)(&m, &n, &k_, q.ptr_data(), &ldq,
					this->m_tau.ptr_data(), pws, &lwork, &info);
			LMAT_CHECK_LAPACK_INFO(sorgqr, info);
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
			lapack_int lwork = (lapack_int)(this->m_ws_size);
			lapack_int info = 0;

			LMAT_LAPACK_NAME(dgeqrf)(&m, &n, this->m_a.ptr_data(), &lda, this->m_tau.ptr_data(),
					this->m_ws.ptr_data(), &lwork, &info);
			LMAT_CHECK_LAPACK_INFO(dgeqrf, info);
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
			lapack_int lwork = (lapack_int)this->m_ws_size;
			double *pws = const_cast<double*>(m_ws.ptr_data());
			lapack_int info;

			LMAT_LAPACK_NAME(dorgqr)(&m, &n, &k_, q.ptr_data(), &ldq,
					this->m_tau.ptr_data(), pws, &lwork, &info);
			LMAT_CHECK_LAPACK_INFO(dorgqr, info);
		}

	};


} }

#endif /* LAPACK_QR_H_ */
