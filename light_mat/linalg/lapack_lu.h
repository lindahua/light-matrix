/**
 * @file lapack_lu.h
 *
 * @brief LU Factorization
 *
 * @author Dahua Lin
 */

#ifndef LIGHTMAT_LAPACK_LU_H_
#define LIGHTMAT_LAPACK_LU_H_

#include "internal/linalg_aux.h"


/************************************************
 *
 *  external LAPACK functions
 *
 ************************************************/

extern "C"
{
	void LMAT_LAPACK_NAME(sgeequ)( const lapack_int* m, const lapack_int* n, const float* a,
	             const lapack_int* lda, float* r, float* c, float* rowcnd,
	             float* colcnd, float* amax, lapack_int* info );

	void LMAT_LAPACK_NAME(dgeequ)( const lapack_int* m, const lapack_int* n, const double* a,
	             const lapack_int* lda, double* r, double* c, double* rowcnd,
	             double* colcnd, double* amax, lapack_int* info );

	void LMAT_LAPACK_NAME(sgetrf)( const lapack_int* m, const lapack_int* n, float* a, const lapack_int* lda,
	             lapack_int* ipiv, lapack_int* info );

	void LMAT_LAPACK_NAME(sgecon)( const char* norm, const lapack_int* n, const float* a,
	             const lapack_int* lda, const float* anorm, float* rcond, float* work,
	             lapack_int* iwork, lapack_int* info );

	void LMAT_LAPACK_NAME(sgetri)( const lapack_int* n, float* a, const lapack_int* lda,
	             const lapack_int* ipiv, float* work, const lapack_int* lwork,
	             lapack_int* info );

	void LMAT_LAPACK_NAME(sgetrs)( const char* trans, const lapack_int* n, const lapack_int* nrhs,
	             const float* a, const lapack_int* lda, const lapack_int* ipiv, float* b,
	             const lapack_int* ldb, lapack_int* info );

	void LMAT_LAPACK_NAME(sgerfs)( const char* trans, const lapack_int* n, const lapack_int* nrhs,
	             const float* a, const lapack_int* lda, const float* af,
	             const lapack_int* ldaf, const lapack_int* ipiv, const float* b,
	             const lapack_int* ldb, float* x, const lapack_int* ldx, float* ferr,
	             float* berr, float* work, lapack_int* iwork, lapack_int* info );

	void LMAT_LAPACK_NAME(dgetrf)( const lapack_int* m, const lapack_int* n, double* a, const lapack_int* lda,
	             lapack_int* ipiv, lapack_int* info );

	void LMAT_LAPACK_NAME(dgecon)( const char* norm, const lapack_int* n, const double* a,
	             const lapack_int* lda, const double* anorm, double* rcond,
	             double* work, lapack_int* iwork, lapack_int* info );

	void LMAT_LAPACK_NAME(dgetri)( const lapack_int* n, double* a, const lapack_int* lda,
	             const lapack_int* ipiv, double* work, const lapack_int* lwork,
	             lapack_int* info );

	void LMAT_LAPACK_NAME(dgetrs)( const char* trans, const lapack_int* n, const lapack_int* nrhs,
	             const double* a, const lapack_int* lda, const lapack_int* ipiv,
	             double* b, const lapack_int* ldb, lapack_int* info );

	void LMAT_LAPACK_NAME(dgerfs)( const char* trans, const lapack_int* n, const lapack_int* nrhs,
	             const double* a, const lapack_int* lda, const double* af,
	             const lapack_int* ldaf, const lapack_int* ipiv, const double* b,
	             const lapack_int* ldb, double* x, const lapack_int* ldx, double* ferr,
	             double* berr, double* work, lapack_int* iwork, lapack_int* info );

	void LMAT_LAPACK_NAME(sgesv)( const lapack_int* n, const lapack_int* nrhs, float* a,
	            const lapack_int* lda, lapack_int* ipiv, float* b, const lapack_int* ldb,
	            lapack_int* info );

	void LMAT_LAPACK_NAME(sgesvx)( const char* fact, const char* trans, const lapack_int* n,
	             const lapack_int* nrhs, float* a, const lapack_int* lda, float* af,
	             const lapack_int* ldaf, lapack_int* ipiv, char* equed, float* r,
	             float* c, float* b, const lapack_int* ldb, float* x,
	             const lapack_int* ldx, float* rcond, float* ferr, float* berr,
	             float* work, lapack_int* iwork, lapack_int* info );

	void LMAT_LAPACK_NAME(dgesv)( const lapack_int* n, const lapack_int* nrhs, double* a,
	            const lapack_int* lda, lapack_int* ipiv, double* b, const lapack_int* ldb,
	            lapack_int* info );

	void LMAT_LAPACK_NAME(dgesvx)( const char* fact, const char* trans, const lapack_int* n,
	             const lapack_int* nrhs, double* a, const lapack_int* lda, double* af,
	             const lapack_int* ldaf, lapack_int* ipiv, char* equed, double* r,
	             double* c, double* b, const lapack_int* ldb, double* x,
	             const lapack_int* ldx, double* rcond, double* ferr, double* berr,
	             double* work, lapack_int* iwork, lapack_int* info );
}


namespace lmat { namespace lapack {

	// forward declarations

	template<typename T> class lu_fac;

	/************************************************
	 *
	 *  LU classes
	 *
	 ************************************************/

	template<typename T>
	class lu_base
	{
	public:
		lu_base()
		: m_dim(0), m_a(), m_ipiv() { }

		bool empty() const
		{
			return m_dim == 0;
		}

		index_t dim() const
		{
			return m_dim;
		}

		const lapack_int* ipiv() const
		{
			return m_ipiv.ptr_data();
		}

		const dense_matrix<T>& intern() const
		{
			return m_a;
		}

		template<class L>
		void getl(IRegularMatrix<L, T>& mat) const
		{
			LMAT_CHECK_DIMS( mat.nrows() == m_dim && mat.ncolumns() == m_dim )

			zero(mat);
			lmat::internal::get_tril(m_dim, m_a, mat, true);
		}

		template<class L>
		void getu(IRegularMatrix<L, T>& mat) const
		{
			LMAT_CHECK_DIMS( mat.nrows() == m_dim && mat.ncolumns() == m_dim )

			zero(mat);
			lmat::internal::get_triu(m_dim, m_a, mat);
		}

	protected:
		template<class Mat>
		void set_mat(const IRegularMatrix<Mat, T>& mat)
		{
			static_assert( meta::is_percol_continuous<Mat>::value,
					"Mat must be percol continuous.");

			LMAT_CHECK_DIMS( mat.nrows() == mat.ncolumns() )

			m_dim = mat.nrows();
			m_a = mat;
			m_ipiv.require_size(m_dim);
		}

	protected:
		index_t m_dim;
		dense_matrix<T> m_a;
		dense_col<lapack_int> m_ipiv;
	};


	template<>
	class lu_fac<float> : public lu_base<float>
	{
	public:
		template<class Mat>
		void set(const IRegularMatrix<Mat, float>& mat)
		{
			this->set_mat(mat);
			trf(this->m_a, this->m_ipiv.ptr_data());
		}

		template<class X>
		void solve_inplace(IRegularMatrix<X, float>& b, char trans='N') const
		{
			lapack_int n = (lapack_int)(this->m_dim);
			lapack_int nrhs = (lapack_int)(b.ncolumns());
			lapack_int lda = (lapack_int)(this->m_a.col_stride());
			lapack_int ldb = (lapack_int)(b.col_stride());
			lapack_int info = 0;

			LMAT_LAPACK_NAME(sgetrs)(&trans, &n, &nrhs, this->m_a.ptr_data(), &lda,
					this->m_ipiv.ptr_data(), b.ptr_data(), &ldb, &info);
			LMAT_CHECK_LAPACK_INFO( sgetrs, info );
		}

		template<class B, class X>
		void solve(const IRegularMatrix<B, float>& b, IRegularMatrix<X, float>& x, char trans='N') const
		{
			x.derived() = b.derived();
			solve_inplace(x, trans);
		}

	public:
		template<class A>
		static void inv_inplace(IRegularMatrix<A, float>& a)
		{
			static_assert( meta::is_percol_continuous<A>::value,
					"a must be percol continuous.");

			LMAT_CHECK_DIMS( a.nrows() == a.ncolumns() );

			dense_col<lapack_int> ipiv(a.nrows());

			trf(a, ipiv.ptr_data());

			lapack_int n = (lapack_int)a.nrows();
			lapack_int lda = (lapack_int)a.col_stride();
			lapack_int info = 0;

			lapack_int lwork = 64 * n;
			dense_col<float> ws((index_t)lwork);

			LMAT_LAPACK_NAME(sgetri)(&n, a.ptr_data(), &lda, ipiv.ptr_data(), ws.ptr_data(), &lwork, &info);
			LMAT_CHECK_LAPACK_INFO( sgetri, info );
		}

		template<class A, class B>
		static void inv(const IRegularMatrix<A, float>& a, IRegularMatrix<B, float>& b)
		{
			b.derived() = a.derived();
			inv_inplace(b);
		}

	private:

		template<class A>
		static void trf(IRegularMatrix<A, float>& a, lapack_int* ipiv)
		{
			lapack_int n = (lapack_int)(a.nrows());
			lapack_int lda = (lapack_int)(a.col_stride());
			lapack_int info = 0;

			LMAT_LAPACK_NAME(sgetrf)(&n, &n, a.ptr_data(), &lda, ipiv, &info);
			LMAT_CHECK_LAPACK_INFO( sgetrf, info);
		}
	};



	template<>
	class lu_fac<double> : public lu_base<double>
	{
	public:
		template<class Mat>
		void set(const IRegularMatrix<Mat, double>& mat)
		{
			this->set_mat(mat);
			trf(this->m_a, this->m_ipiv.ptr_data());
		}

		template<class X>
		void solve_inplace(IRegularMatrix<X, double>& b, char trans='N') const
		{
			lapack_int n = (lapack_int)(this->m_dim);
			lapack_int nrhs = (lapack_int)(b.ncolumns());
			lapack_int lda = (lapack_int)(this->m_a.col_stride());
			lapack_int ldb = (lapack_int)(b.col_stride());
			lapack_int info = 0;

			LMAT_LAPACK_NAME(dgetrs)(&trans, &n, &nrhs, this->m_a.ptr_data(), &lda,
					this->m_ipiv.ptr_data(), b.ptr_data(), &ldb, &info);
			LMAT_CHECK_LAPACK_INFO( dgetrs, info );
		}

		template<class B, class X>
		void solve(const IRegularMatrix<B, double>& b, IRegularMatrix<X, double>& x, char trans='N') const
		{
			x.derived() = b.derived();
			solve_inplace(x, trans);
		}

	public:
		template<class A>
		static void inv_inplace(IRegularMatrix<A, double>& a)
		{
			static_assert( meta::is_percol_continuous<A>::value,
					"a must be percol continuous.");

			LMAT_CHECK_DIMS( a.nrows() == a.ncolumns() );

			dense_col<lapack_int> ipiv(a.nrows());

			trf(a, ipiv.ptr_data());

			lapack_int n = (lapack_int)a.nrows();
			lapack_int lda = (lapack_int)a.col_stride();
			lapack_int info = 0;

			lapack_int lwork = 64 * n;
			dense_col<double> ws((index_t)lwork);

			LMAT_LAPACK_NAME(dgetri)(&n, a.ptr_data(), &lda, ipiv.ptr_data(), ws.ptr_data(), &lwork, &info);
			LMAT_CHECK_LAPACK_INFO( dgetri, info );
		}

		template<class A, class B>
		static void inv(const IRegularMatrix<A, double>& a, IRegularMatrix<B, double>& b)
		{
			b.derived() = a.derived();
			inv_inplace(b);
		}

	private:

		template<class A>
		static void trf(IRegularMatrix<A, double>& a, lapack_int* ipiv)
		{
			lapack_int n = (lapack_int)(a.nrows());
			lapack_int lda = (lapack_int)(a.col_stride());
			lapack_int info = 0;

			LMAT_LAPACK_NAME(dgetrf)(&n, &n, a.ptr_data(), &lda, ipiv, &info);
			LMAT_CHECK_LAPACK_INFO( dgetrf, info);
		}
	};


	/************************************************
	 *
	 *  convenient functions
	 *
	 ************************************************/

	template<class A, class B>
	inline dense_matrix<float, meta::nrows<B>::value, meta::ncols<B>::value>
	gesv(const IRegularMatrix<A, float>& a, IRegularMatrix<B, float>& b)
	{
		LMAT_CHECK_DIMS( a.nrows() == a.ncolumns() && a.nrows() == b.nrows() );

		lapack_int n = (lapack_int)a.nrows();
		lapack_int nrhs = (lapack_int)b.ncolumns();
		lapack_int lda = (lapack_int)a.col_stride();
		lapack_int ldb = (lapack_int)b.col_stride();
		dense_col<lapack_int> ipiv(n);

		lapack_int info = 0;
		LMAT_LAPACK_NAME(sgesv)(&n, &nrhs, a.ptr_data(), &lda, ipiv.ptr_data(), b.ptr_data(), &ldb, &info);
		LMAT_CHECK_LAPACK_INFO( sgesv, info );
	}

	template<class A, class B>
	inline dense_matrix<double, meta::nrows<B>::value, meta::ncols<B>::value>
	gesv(const IRegularMatrix<A, double>& a, IRegularMatrix<B, double>& b)
	{
		LMAT_CHECK_DIMS( a.nrows() == a.ncolumns() && a.nrows() == b.nrows() );

		lapack_int n = (lapack_int)a.nrows();
		lapack_int nrhs = (lapack_int)b.ncolumns();
		lapack_int lda = (lapack_int)a.col_stride();
		lapack_int ldb = (lapack_int)b.col_stride();
		dense_col<lapack_int> ipiv(n);

		lapack_int info = 0;
		LMAT_LAPACK_NAME(dgesv)(&n, &nrhs, a.ptr_data(), &lda, ipiv.ptr_data(), b.ptr_data(), &ldb, &info);
		LMAT_CHECK_LAPACK_INFO( dgesv, info );
	}

	template<class A, class B>
	inline dense_matrix<float, meta::nrows<B>::value, meta::ncols<B>::value>
	solve(const IRegularMatrix<A, float>& a, const IRegularMatrix<B, float>& b)
	{
		typedef dense_matrix<float, meta::nrows<B>::value, meta::ncols<B>::value> rmat_t;
		rmat_t x = b;
		gesv(a, x);
		return x;
	}

	template<class A, class B>
	inline dense_matrix<double, meta::nrows<B>::value, meta::ncols<B>::value>
	solve(const IRegularMatrix<A, double>& a, const IRegularMatrix<B, double>& b)
	{
		typedef dense_matrix<double, meta::nrows<B>::value, meta::ncols<B>::value> rmat_t;
		rmat_t x = b;
		gesv(a, x);
		return x;
	}

	template<class A>
	inline dense_matrix<float, meta::nrows<A>::value, meta::ncols<A>::value>
	inv(const IRegularMatrix<A, float>& a)
	{
		typedef dense_matrix<float, meta::nrows<A>::value, meta::ncols<A>::value> rmat_t;
		rmat_t r;
		lu_fac<float>::inv(a, r);
		return r;
	}

	template<class A>
	inline dense_matrix<double, meta::nrows<A>::value, meta::ncols<A>::value>
	inv(const IRegularMatrix<A, double>& a)
	{
		typedef dense_matrix<double, meta::nrows<A>::value, meta::ncols<A>::value> rmat_t;
		rmat_t r;
		lu_fac<double>::inv(a, r);
		return r;
	}

} }

#endif /* LAPACK_LU_H_ */
