/**
 * @file lapack_lu.h
 *
 * @brief LU Factorization
 *
 * @author Dahua Lin
 */

#ifndef LIGHTMAT_LAPACK_LU_H_
#define LIGHTMAT_LAPACK_LU_H_

#include <light_mat/linalg/lapack_fwd.h>


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
			LMAT_CHECK_PERCOL_CONT(L)
			LMAT_CHECK_DIMS( mat.nrows() == m_dim && mat.ncolumns() == m_dim )

			zero(mat);
			copy_tril(m_a, mat, -1);
			for (index_t i = 0; i < m_dim; ++i) mat(i, i) = T(1);
		}

		template<class U>
		void getu(IRegularMatrix<U, T>& mat) const
		{
			LMAT_CHECK_PERCOL_CONT(U)
			LMAT_CHECK_DIMS( mat.nrows() == m_dim && mat.ncolumns() == m_dim )

			zero(mat);
			copy_triu(m_a, mat, 0);
		}

	protected:
		template<class Mat>
		void set_mat(const IMatrixXpr<Mat, T>& mat)
		{
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
		LMAT_ENSURE_INLINE
		explicit lu_fac() { }

		template<class Mat>
		LMAT_ENSURE_INLINE
		explicit lu_fac(const IMatrixXpr<Mat, float>& mat)
		{
			set(mat);
		}

		template<class Mat>
		void set(const IMatrixXpr<Mat, float>& mat)
		{
			this->set_mat(mat);
			trf(this->m_a, this->m_ipiv.ptr_data());
		}

		template<class B>
		void solve_inplace(IRegularMatrix<B, float>& b, char trans='N') const
		{
			LMAT_CHECK_PERCOL_CONT(B)

			lapack_int n = (lapack_int)(this->m_dim);
			lapack_int nrhs = (lapack_int)(b.ncolumns());
			lapack_int lda = (lapack_int)(this->m_a.col_stride());
			lapack_int ldb = (lapack_int)(b.col_stride());
			lapack_int info = 0;

			LMAT_CALL_LAPACK(sgetrs, (&trans, &n, &nrhs, this->m_a.ptr_data(), &lda,
					this->m_ipiv.ptr_data(), b.ptr_data(), &ldb, &info));
		}

		template<class B, class X>
		void solve(const IMatrixXpr<B, float>& b, IRegularMatrix<X, float>& x, char trans='N') const
		{
			LMAT_CHECK_PERCOL_CONT(X)

			x.derived() = b.derived();
			solve_inplace(x, trans);
		}

	public:
		template<class A>
		static void inv_inplace(IRegularMatrix<A, float>& a)
		{
			LMAT_CHECK_PERCOL_CONT(A)

			LMAT_CHECK_DIMS( a.nrows() == a.ncolumns() );

			dense_col<lapack_int> ipiv(a.nrows());

			trf(a, ipiv.ptr_data());

			lapack_int n = (lapack_int)a.nrows();
			lapack_int lda = (lapack_int)a.col_stride();
			lapack_int info = 0;

			lapack_int lwork = -1;
			float lwork_opt = 0;
			LMAT_CALL_LAPACK(sgetri, (&n, a.ptr_data(), &lda, ipiv.ptr_data(), &lwork_opt, &lwork, &info));

			lwork = (lapack_int)lwork_opt;
			dense_col<float> ws((index_t)lwork);

			LMAT_CALL_LAPACK(sgetri, (&n, a.ptr_data(), &lda, ipiv.ptr_data(), ws.ptr_data(), &lwork, &info));
		}

		template<class A, class B>
		static void inv(const IMatrixXpr<A, float>& a, IRegularMatrix<B, float>& b)
		{
			LMAT_CHECK_PERCOL_CONT(B)

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

			LMAT_CALL_LAPACK(sgetrf, (&n, &n, a.ptr_data(), &lda, ipiv, &info));
		}
	};


	template<>
	class lu_fac<double> : public lu_base<double>
	{
	public:
		LMAT_ENSURE_INLINE
		explicit lu_fac() { }

		template<class Mat>
		LMAT_ENSURE_INLINE
		explicit lu_fac(const IMatrixXpr<Mat, double>& mat)
		{
			set(mat);
		}

		template<class Mat>
		void set(const IMatrixXpr<Mat, double>& mat)
		{
			this->set_mat(mat);
			trf(this->m_a, this->m_ipiv.ptr_data());
		}

		template<class B>
		void solve_inplace(IRegularMatrix<B, double>& b, char trans='N') const
		{
			LMAT_CHECK_PERCOL_CONT(B)

			lapack_int n = (lapack_int)(this->m_dim);
			lapack_int nrhs = (lapack_int)(b.ncolumns());
			lapack_int lda = (lapack_int)(this->m_a.col_stride());
			lapack_int ldb = (lapack_int)(b.col_stride());
			lapack_int info = 0;

			LMAT_CALL_LAPACK(dgetrs, (&trans, &n, &nrhs, this->m_a.ptr_data(), &lda,
					this->m_ipiv.ptr_data(), b.ptr_data(), &ldb, &info));
		}

		template<class B, class X>
		void solve(const IMatrixXpr<B, double>& b, IRegularMatrix<X, double>& x, char trans='N') const
		{
			LMAT_CHECK_PERCOL_CONT(X)

			x.derived() = b.derived();
			solve_inplace(x, trans);
		}

	public:
		template<class A>
		static void inv_inplace(IRegularMatrix<A, double>& a)
		{
			LMAT_CHECK_PERCOL_CONT(A)

			LMAT_CHECK_DIMS( a.nrows() == a.ncolumns() );

			dense_col<lapack_int> ipiv(a.nrows());

			trf(a, ipiv.ptr_data());

			lapack_int n = (lapack_int)a.nrows();
			lapack_int lda = (lapack_int)a.col_stride();
			lapack_int info = 0;

			lapack_int lwork = -1;
			double lwork_opt = 0;
			LMAT_CALL_LAPACK(dgetri, (&n, a.ptr_data(), &lda, ipiv.ptr_data(), &lwork_opt, &lwork, &info));

			lwork = (lapack_int)lwork_opt;
			dense_col<double> ws((index_t)lwork);

			LMAT_CALL_LAPACK(dgetri, (&n, a.ptr_data(), &lda, ipiv.ptr_data(), ws.ptr_data(), &lwork, &info));
		}

		template<class A, class B>
		static void inv(const IMatrixXpr<A, double>& a, IRegularMatrix<B, double>& b)
		{
			LMAT_CHECK_PERCOL_CONT(B)

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

			LMAT_CALL_LAPACK(dgetrf, (&n, &n, a.ptr_data(), &lda, ipiv, &info));
		}
	};


	/************************************************
	 *
	 *  GESV
	 *
	 ************************************************/

	template<class A, class B>
	inline void gesv(IRegularMatrix<A, float>& a, IRegularMatrix<B, float>& b)
	{
		LMAT_CHECK_PERCOL_CONT(A)
		LMAT_CHECK_PERCOL_CONT(B)

		LMAT_CHECK_DIMS( a.nrows() == a.ncolumns() && a.nrows() == b.nrows() );

		lapack_int n = (lapack_int)a.nrows();
		lapack_int nrhs = (lapack_int)b.ncolumns();
		lapack_int lda = (lapack_int)a.col_stride();
		lapack_int ldb = (lapack_int)b.col_stride();
		dense_col<lapack_int> ipiv(n);

		lapack_int info = 0;
		LMAT_CALL_LAPACK(sgesv, (&n, &nrhs, a.ptr_data(), &lda, ipiv.ptr_data(), b.ptr_data(), &ldb, &info));
	}

	template<class A, class B>
	inline void gesv(IRegularMatrix<A, double>& a, IRegularMatrix<B, double>& b)
	{
		LMAT_CHECK_PERCOL_CONT(A)
		LMAT_CHECK_PERCOL_CONT(B)

		LMAT_CHECK_DIMS( a.nrows() == a.ncolumns() && a.nrows() == b.nrows() );

		lapack_int n = (lapack_int)a.nrows();
		lapack_int nrhs = (lapack_int)b.ncolumns();
		lapack_int lda = (lapack_int)a.col_stride();
		lapack_int ldb = (lapack_int)b.col_stride();
		dense_col<lapack_int> ipiv(n);

		lapack_int info = 0;
		LMAT_CALL_LAPACK(dgesv, (&n, &nrhs, a.ptr_data(), &lda, ipiv.ptr_data(), b.ptr_data(), &ldb, &info));
	}

} }


namespace lmat
{

	template<class Arg> class inv_expr;

	/************************************************
	 *
	 *  inversion expression
	 *
	 ************************************************/

	template<class Arg>
	struct matrix_traits<inv_expr<Arg> >
	{
		static const int num_dimensions = 2;
		static const int ct_num_rows = meta::sq_dim<Arg>::value;
		static const int ct_num_cols = ct_num_rows;

		static const bool is_readonly = true;

		typedef matrix_shape<ct_num_rows, ct_num_cols> shape_type;
		typedef typename matrix_traits<Arg>::value_type value_type;
		typedef typename matrix_traits<Arg>::domain domain;
	};


	namespace internal
	{
		template<class Arg>
		LMAT_ENSURE_INLINE
		inline typename matrix_traits<inv_expr<Arg> >::shape_type
		inv_shape(const Arg& a)
		{
			LMAT_CHECK_DIMS( a.nrows() == a.ncolumns() )

			typedef typename matrix_traits<inv_expr<Arg> >::shape_type shape_t;
			return shape_t(a.nrows(), a.ncolumns());
		}
	}

	template<class Arg>
	class inv_expr : public IMatrixXpr<inv_expr<Arg>, typename matrix_traits<Arg>::value_type>
	{
	public:
		static const int ct_dim = meta::sq_dim<Arg>::value;
		typedef matrix_shape<ct_dim, ct_dim> shape_type;

		LMAT_ENSURE_INLINE
		explicit inv_expr(const Arg& a)
		: m_shape(internal::inv_shape(a)), m_arg(a)
		{ }

		LMAT_ENSURE_INLINE
		index_t nrows() const { return m_shape.nrows(); }

		LMAT_ENSURE_INLINE
		index_t ncolumns() const { return m_shape.ncolumns(); }

		LMAT_ENSURE_INLINE
		index_t nelems() const { return m_shape.nelems(); }

		LMAT_ENSURE_INLINE
		shape_type shape() const { return m_shape; }

		LMAT_ENSURE_INLINE
		const Arg& arg() const
		{
			return m_arg;
		}

	private:
		shape_type m_shape;
		const Arg& m_arg;
	};


	template<class Arg, class DMat>
	inline void evaluate(const inv_expr<Arg>& expr,
			IRegularMatrix<DMat, typename matrix_traits<Arg>::value_type>& dmat)
	{
		LMAT_CHECK_PERCOL_CONT(DMat)

		typedef typename matrix_traits<Arg>::value_type T;
		lapack::lu_fac<T>::inv(expr.arg(), dmat.derived());
	}

	template<class Arg>
	LMAT_ENSURE_INLINE
	inline inv_expr<Arg> inv(const IMatrixXpr<Arg, float>& a)
	{
		return inv_expr<Arg>(a.derived());
	}

	template<class Arg>
	LMAT_ENSURE_INLINE
	inline inv_expr<Arg> inv(const IMatrixXpr<Arg, double>& a)
	{
		return inv_expr<Arg>(a.derived());
	}

}


#endif /* LAPACK_LU_H_ */
