/**
 * @file lapack_chol.h
 *
 * @brief Cholesky factorization
 *
 * @author Dahua Lin
 */


#ifdef _MSC_VER
#pragma once
#endif

#ifndef LIGHTMAT_LAPACK_CHOL_H_
#define LIGHTMAT_LAPACK_CHOL_H_

#include <light_mat/linalg/lapack_fwd.h>


/************************************************
 *
 *  external LAPACK functions
 *
 ************************************************/

extern "C"
{
	void LMAT_LAPACK_NAME(spoequ)( const lapack_int* n, const float* a, const lapack_int* lda, float* s,
	             float* scond, float* amax, lapack_int* info );

	void LMAT_LAPACK_NAME(dpoequ)( const lapack_int* n, const double* a, const lapack_int* lda, double* s,
	             double* scond, double* amax, lapack_int* info );


	void LMAT_LAPACK_NAME(spotrf)( const char* uplo, const lapack_int* n, float* a, const lapack_int* lda,
	             lapack_int* info );

	void LMAT_LAPACK_NAME(spocon)( const char* uplo, const lapack_int* n, const float* a,
	             const lapack_int* lda, const float* anorm, float* rcond, float* work,
	             lapack_int* iwork, lapack_int* info );

	void LMAT_LAPACK_NAME(spotri)( const char* uplo, const lapack_int* n, float* a, const lapack_int* lda,
	             lapack_int* info );

	void LMAT_LAPACK_NAME(spotrs)( const char* uplo, const lapack_int* n, const lapack_int* nrhs,
	             const float* a, const lapack_int* lda, float* b, const lapack_int* ldb,
	             lapack_int* info );

	void LMAT_LAPACK_NAME(sporfs)( const char* uplo, const lapack_int* n, const lapack_int* nrhs,
	             const float* a, const lapack_int* lda, const float* af,
	             const lapack_int* ldaf, const float* b, const lapack_int* ldb, float* x,
	             const lapack_int* ldx, float* ferr, float* berr, float* work,
	             lapack_int* iwork, lapack_int* info );

	void LMAT_LAPACK_NAME(dpotrf)( const char* uplo, const lapack_int* n, double* a, const lapack_int* lda,
	             lapack_int* info );

	void LMAT_LAPACK_NAME(dpocon)( const char* uplo, const lapack_int* n, const double* a,
	             const lapack_int* lda, const double* anorm, double* rcond,
	             double* work, lapack_int* iwork, lapack_int* info );

	void LMAT_LAPACK_NAME(dpotri)( const char* uplo, const lapack_int* n, double* a, const lapack_int* lda,
	             lapack_int* info );

	void LMAT_LAPACK_NAME(dpotrs)( const char* uplo, const lapack_int* n, const lapack_int* nrhs,
	             const double* a, const lapack_int* lda, double* b,
	             const lapack_int* ldb, lapack_int* info );

	void LMAT_LAPACK_NAME(dporfs)( const char* uplo, const lapack_int* n, const lapack_int* nrhs,
	             const double* a, const lapack_int* lda, const double* af,
	             const lapack_int* ldaf, const double* b, const lapack_int* ldb,
	             double* x, const lapack_int* ldx, double* ferr, double* berr,
	             double* work, lapack_int* iwork, lapack_int* info );

	void LMAT_LAPACK_NAME(sposv)( const char* uplo, const lapack_int* n, const lapack_int* nrhs, float* a,
	            const lapack_int* lda, float* b, const lapack_int* ldb, lapack_int* info );

	void LMAT_LAPACK_NAME(sposvx)( const char* fact, const char* uplo, const lapack_int* n,
	             const lapack_int* nrhs, float* a, const lapack_int* lda, float* af,
	             const lapack_int* ldaf, char* equed, float* s, float* b,
	             const lapack_int* ldb, float* x, const lapack_int* ldx, float* rcond,
	             float* ferr, float* berr, float* work, lapack_int* iwork,
	             lapack_int* info );

	void LMAT_LAPACK_NAME(dposv)( const char* uplo, const lapack_int* n, const lapack_int* nrhs, double* a,
	            const lapack_int* lda, double* b, const lapack_int* ldb, lapack_int* info );

	void LMAT_LAPACK_NAME(dposvx)( const char* fact, const char* uplo, const lapack_int* n,
	             const lapack_int* nrhs, double* a, const lapack_int* lda, double* af,
	             const lapack_int* ldaf, char* equed, double* s, double* b,
	             const lapack_int* ldb, double* x, const lapack_int* ldx, double* rcond,
	             double* ferr, double* berr, double* work, lapack_int* iwork,
	             lapack_int* info );
}


namespace lmat { namespace lapack {

	// forward declarations

	template<typename T> class chol_fac;


	/************************************************
	 *
	 *  Cholesky classes
	 *
	 ************************************************/

	namespace internal
	{
		LMAT_ENSURE_INLINE
		inline char check_chol_uplo(char c)
		{
			char r;
			if (c == 'U' || c == 'u') r = 'U';
			else if (c == 'L' || c == 'l') r = 'L';
			else
				throw invalid_argument("Invalid value for uplo");

			return r;
		}
	}


	template<typename T>
	class chol_base
	{
	public:
		explicit chol_base( char uplo )
		: m_uplo(internal::check_chol_uplo(uplo)), m_dim(0), m_a() { }

		bool empty() const
		{
			return m_dim == 0;
		}

		bool is_lower() const
		{
			return m_uplo == 'L';
		}

		bool is_upper() const
		{
			return m_uplo == 'U';
		}

		index_t dim() const
		{
			return m_dim;
		}

		const dense_matrix<T>& intern() const
		{
			return m_a;
		}

		template<class L>
		void get(IRegularMatrix<L, T>& mat) const
		{
			mat.require_size(m_dim, m_dim);
			zero(mat);

			if (is_lower())
			{
				copy_tril(m_a, mat);
			}
			else
			{
				copy_triu(m_a, mat);
			}
		}

	protected:
		template<class Mat>
		void set_mat(const IMatrixXpr<Mat, T>& mat)
		{
			LMAT_CHECK_DIMS( mat.nrows() == mat.ncolumns() )

			m_dim = mat.nrows();
			m_a = mat;
		}

	protected:
		const char m_uplo;
		index_t m_dim;
		dense_matrix<T> m_a;
	};


	template<>
	class chol_fac<float> : public chol_base<float>
	{
	public:
		explicit chol_fac(char uplo='L')
		: chol_base<float>(uplo) { }

		template<class Mat>
		explicit chol_fac(const IMatrixXpr<Mat, float>& a, char uplo='L')
		: chol_base<float>(uplo)
		{
			set(a);
		}

		template<class Mat>
		void set(const IMatrixXpr<Mat, float>& mat)
		{
			this->set_mat(mat);
			trf(this->m_a, this->m_uplo);
		}

		template<class B>
		void solve_inplace(IRegularMatrix<B, float>& b) const
		{
			LMAT_CHECK_PERCOL_CONT(B)

			lapack_int n = (lapack_int)(this->m_dim);
			lapack_int nrhs = (lapack_int)(b.ncolumns());
			lapack_int lda = (lapack_int)(this->m_a.col_stride());
			lapack_int ldb = (lapack_int)(b.col_stride());
			lapack_int info = 0;

			LMAT_CALL_LAPACK(spotrs, (&(this->m_uplo), &n, &nrhs,
					this->m_a.ptr_data(), &lda, b.ptr_data(), &ldb, &info));
		}

		template<class B, class X>
		void solve(const IMatrixXpr<B, float>& b, IRegularMatrix<X, float>& x) const
		{
			LMAT_CHECK_PERCOL_CONT(X)

			x.derived() = b.derived();
			solve_inplace(x);
		}

	public:
		template<class A>
		static void inv_inplace(IRegularMatrix<A, float>& a, char uplo='L')
		{
			LMAT_CHECK_PERCOL_CONT(A)
			LMAT_CHECK_DIMS( a.nrows() == a.ncolumns() );

			uplo = internal::check_chol_uplo(uplo);
			trf(a, uplo);

			lapack_int n = (lapack_int)a.nrows();
			lapack_int lda = (lapack_int)a.col_stride();
			lapack_int info = 0;

			LMAT_CALL_LAPACK(spotri, (&uplo, &n, a.ptr_data(), &lda, &info));

			lmat::internal::complete_sym(a.nrows(), a, uplo);
		}

		template<class A, class B>
		static void inv(const IMatrixXpr<A, float>& a, IRegularMatrix<B, float>& b, char uplo='L')
		{
			LMAT_CHECK_PERCOL_CONT(B)

			b.derived() = a.derived();
			inv_inplace(b, uplo);
		}

	private:

		template<class A>
		static void trf(IRegularMatrix<A, float>& a, char uplo)
		{
			lapack_int n = (lapack_int)(a.nrows());
			lapack_int lda = (lapack_int)(a.col_stride());
			lapack_int info = 0;

			LMAT_CALL_LAPACK(spotrf, (&uplo, &n, a.ptr_data(), &lda, &info));
		}
	};


	template<>
	class chol_fac<double> : public chol_base<double>
	{
	public:
		explicit chol_fac(char uplo='L')
		: chol_base<double>(uplo) { }

		template<class Mat>
		explicit chol_fac(const IMatrixXpr<Mat, double>& a, char uplo='L')
		: chol_base<double>(uplo)
		{
			set(a);
		}

		template<class Mat>
		void set(const IMatrixXpr<Mat, double>& mat)
		{
			this->set_mat(mat);
			trf(this->m_a, this->m_uplo);
		}

		template<class B>
		void solve_inplace(IRegularMatrix<B, double>& b) const
		{
			LMAT_CHECK_PERCOL_CONT(B)

			lapack_int n = (lapack_int)(this->m_dim);
			lapack_int nrhs = (lapack_int)(b.ncolumns());
			lapack_int lda = (lapack_int)(this->m_a.col_stride());
			lapack_int ldb = (lapack_int)(b.col_stride());
			lapack_int info = 0;

			LMAT_CALL_LAPACK(dpotrs, (&(this->m_uplo), &n, &nrhs,
					this->m_a.ptr_data(), &lda, b.ptr_data(), &ldb, &info));
		}

		template<class B, class X>
		void solve(const IMatrixXpr<B, double>& b, IRegularMatrix<X, double>& x) const
		{
			LMAT_CHECK_PERCOL_CONT(X)

			x.derived() = b.derived();
			solve_inplace(x);
		}

	public:
		template<class A>
		static void inv_inplace(IRegularMatrix<A, double>& a, char uplo='L')
		{
			LMAT_CHECK_PERCOL_CONT(A)
			LMAT_CHECK_DIMS( a.nrows() == a.ncolumns() );

			uplo = internal::check_chol_uplo(uplo);
			trf(a, uplo);

			lapack_int n = (lapack_int)a.nrows();
			lapack_int lda = (lapack_int)a.col_stride();
			lapack_int info = 0;

			LMAT_CALL_LAPACK(dpotri, (&uplo, &n, a.ptr_data(), &lda, &info));

			lmat::internal::complete_sym(a.nrows(), a, uplo);
		}

		template<class A, class B>
		static void inv(const IMatrixXpr<A, double>& a, IRegularMatrix<B, double>& b, char uplo='L')
		{
			LMAT_CHECK_PERCOL_CONT(B)

			b.derived() = a.derived();
			inv_inplace(b, uplo);
		}

	private:

		template<class A>
		static void trf(IRegularMatrix<A, double>& a, char uplo)
		{
			lapack_int n = (lapack_int)(a.nrows());
			lapack_int lda = (lapack_int)(a.col_stride());
			lapack_int info = 0;

			LMAT_CALL_LAPACK(dpotrf, (&uplo, &n, a.ptr_data(), &lda, &info));
		}
	};


	/************************************************
	 *
	 *  POSV
	 *
	 ************************************************/

	template<class A, class B>
	inline void posv(IRegularMatrix<A, float>& a, IRegularMatrix<B, float>& b, char uplo='L')
	{
		LMAT_CHECK_PERCOL_CONT(A)
		LMAT_CHECK_PERCOL_CONT(B)

		uplo = internal::check_chol_uplo(uplo);
		LMAT_CHECK_DIMS( a.nrows() == a.ncolumns() && a.nrows() == b.nrows() );

		lapack_int n = (lapack_int)a.nrows();
		lapack_int nrhs = (lapack_int)b.ncolumns();
		lapack_int lda = (lapack_int)a.col_stride();
		lapack_int ldb = (lapack_int)b.col_stride();

		lapack_int info = 0;
		LMAT_CALL_LAPACK(sposv, (&uplo, &n, &nrhs, a.ptr_data(), &lda, b.ptr_data(), &ldb, &info));
	}

	template<class A, class B>
	inline void posv(IRegularMatrix<A, double>& a, IRegularMatrix<B, double>& b, char uplo='L')
	{
		LMAT_CHECK_PERCOL_CONT(A)
		LMAT_CHECK_PERCOL_CONT(B)

		uplo = internal::check_chol_uplo(uplo);
		LMAT_CHECK_DIMS( a.nrows() == a.ncolumns() && a.nrows() == b.nrows() );

		lapack_int n = (lapack_int)a.nrows();
		lapack_int nrhs = (lapack_int)b.ncolumns();
		lapack_int lda = (lapack_int)a.col_stride();
		lapack_int ldb = (lapack_int)b.col_stride();

		lapack_int info = 0;
		LMAT_CALL_LAPACK(dposv, (&uplo, &n, &nrhs, a.ptr_data(), &lda, b.ptr_data(), &ldb, &info));
	}


} }


namespace lmat
{

	template<class Arg> class pdinv_expr;

	/************************************************
	 *
	 *  inversion expression
	 *
	 ************************************************/

	template<class Arg>
	struct matrix_traits<pdinv_expr<Arg> >
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
		inline typename matrix_traits<pdinv_expr<Arg> >::shape_type
		pdinv_shape(const Arg& a)
		{
			LMAT_CHECK_DIMS( a.nrows() == a.ncolumns() )

			typedef typename matrix_traits<pdinv_expr<Arg> >::shape_type shape_t;
			return shape_t(a.nrows(), a.ncolumns());
		}
	}

	template<class Arg>
	class pdinv_expr : public IMatrixXpr<pdinv_expr<Arg>, typename matrix_traits<Arg>::value_type>
	{
	public:
		static const int ct_dim = meta::sq_dim<Arg>::value;
		typedef matrix_shape<ct_dim, ct_dim> shape_type;

		LMAT_ENSURE_INLINE
		explicit pdinv_expr(const Arg& a, char uplo_)
		: m_uplo(lapack::internal::check_chol_uplo(uplo_))
		, m_shape(internal::pdinv_shape(a)), m_arg(a)
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
		char uplo() const
		{
			return m_uplo;
		}

		LMAT_ENSURE_INLINE
		const Arg& arg() const
		{
			return m_arg;
		}

	private:
		const char m_uplo;
		shape_type m_shape;
		const Arg& m_arg;
	};


	template<class Arg, class DMat>
	inline void evaluate(const pdinv_expr<Arg>& expr,
			IRegularMatrix<DMat, typename matrix_traits<Arg>::value_type>& dmat)
	{
		LMAT_CHECK_PERCOL_CONT(DMat)

		typedef typename matrix_traits<Arg>::value_type T;
		lapack::chol_fac<T>::inv(expr.arg(), dmat.derived(), expr.uplo());
	}

	template<class Arg>
	LMAT_ENSURE_INLINE
	inline pdinv_expr<Arg> pdinv(const IMatrixXpr<Arg, float>& a, const char uplo='L')
	{
		return pdinv_expr<Arg>(a.derived(), uplo);
	}

	template<class Arg>
	LMAT_ENSURE_INLINE
	inline pdinv_expr<Arg> pdinv(const IMatrixXpr<Arg, double>& a, const char uplo='L')
	{
		return pdinv_expr<Arg>(a.derived(), uplo);
	}

}



#endif /* LAPACK_CHOL_H_ */
