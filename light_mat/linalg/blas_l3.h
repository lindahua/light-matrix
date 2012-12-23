/**
 * @file blas_l3.h
 *
 * @brief BLAS Level 3 functions
 *
 * @author Dahua Lin
 */

#ifndef LIGHTMAT_BLAS_L3_H_
#define LIGHTMAT_BLAS_L3_H_

#include <light_mat/matrix/matrix_classes.h>
#include "internal/linalg_aux.h"

extern "C"
{
	void LMAT_BLAS_NAME(sgemm)(const char *transa, const char *transb, const blas_int *m, const blas_int *n, const blas_int *k,
	           const float *alpha, const float *a, const blas_int *lda, const float *b, const blas_int *ldb,
	           const float *beta, float *c, const blas_int *ldc);
	void LMAT_BLAS_NAME(ssymm)(const char *side, const char *uplo, const blas_int *m, const blas_int *n,
	           const float *alpha, const float *a, const blas_int *lda, const float *b, const blas_int *ldb,
	           const float *beta, float *c, const blas_int *ldc);
	void LMAT_BLAS_NAME(ssyrk)(const char *uplo, const char *trans, const blas_int *n, const blas_int *k,
	           const float *alpha, const float *a, const blas_int *lda, const float *beta,
	           float *c, const blas_int *ldc);
	void LMAT_BLAS_NAME(strmm)(const char *side, const char *uplo, const char *transa, const char *diag,
	           const blas_int *m, const blas_int *n, const float *alpha, const float *a, const blas_int *lda,
	           float *b, const blas_int *ldb);
	void LMAT_BLAS_NAME(strsm)(const char *side, const char *uplo, const char *transa, const char *diag,
	           const blas_int *m, const blas_int *n, const float *alpha, const float *a, const blas_int *lda,
	           float *b, const blas_int *ldb);

	void LMAT_BLAS_NAME(dgemm)(const char *transa, const char *transb, const blas_int *m, const blas_int *n, const blas_int *k,
	           const double *alpha, const double *a, const blas_int *lda, const double *b, const blas_int *ldb,
	           const double *beta, double *c, const blas_int *ldc);
	void LMAT_BLAS_NAME(dsymm)(const char *side, const char *uplo, const blas_int *m, const blas_int *n,
	           const double *alpha, const double *a, const blas_int *lda, const double *b, const blas_int *ldb,
	           const double *beta, double *c, const blas_int *ldc);
	void LMAT_BLAS_NAME(dsyrk)(const char *uplo, const char *trans, const blas_int *n, const blas_int *k,
	           const double *alpha, const double *a, const blas_int *lda, const double *beta,
	           double *c, const blas_int *ldc);
	void LMAT_BLAS_NAME(dtrmm)(const char *side, const char *uplo, const char *transa, const char *diag,
	           const blas_int *m, const blas_int *n, const double *alpha, const double *a, const blas_int *lda,
	           double *b, const blas_int *ldb);
	void LMAT_BLAS_NAME(dtrsm)(const char *side, const char *uplo, const char *transa, const char *diag,
	           const blas_int *m, const blas_int *n, const double *alpha, const double *a, const blas_int *lda,
	           double *b, const blas_int *ldb);
}


namespace lmat { namespace blas {

	// gemm

	namespace internal
	{
		template<class A>
		LMAT_ENSURE_INLINE
		inline void get_op_dims(const A& a, char transa, index_t& m, index_t& n)
		{
			if (transa == 'N' || transa == 'n')
			{
				m = a.nrows();
				n = a.ncolumns();
			}
			else
			{
				m = a.ncolumns();
				n = a.nrows();
			}
		}


		template<class A, class B, class C>
		LMAT_ENSURE_INLINE
		inline void gemm_get_dims(const A& a, const B& b, const C& c, char transa, char transb,
				blas_int& m, blas_int& n, blas_int& k)
		{
			index_t ma, na, mb, nb;

			get_op_dims(a, transa, ma, na);
			get_op_dims(b, transb, mb, nb);

			LMAT_CHECK_DIMS( na == mb && ma == c.nrows() && nb == c.ncolumns() );

			m = (blas_int)ma;
			n = (blas_int)nb;
			k = (blas_int)na;
		}
	}


	template<class A, class B, class C>
	LMAT_ENSURE_INLINE
	inline void gemm(float alpha, const IRegularMatrix<A, float>& a, const IRegularMatrix<B, float>& b,
	          float beta, IRegularMatrix<C, float>& c, char transa='N', char transb='N')
	{
		blas_int m, n, k;
		internal::gemm_get_dims(a, b, c, transa, transb, m, n, k);

		blas_int lda = (blas_int)a.col_stride();
		blas_int ldb = (blas_int)b.col_stride();
		blas_int ldc = (blas_int)c.col_stride();

		LMAT_BLAS_NAME(sgemm)(&transa, &transb, &m, &n, &k, &alpha,
				a.ptr_data(), &lda, b.ptr_data(), &ldb, &beta, c.ptr_data(), &ldc);
	}

	template<class A, class B, class C>
	LMAT_ENSURE_INLINE
	inline void gemm(double alpha, const IRegularMatrix<A, double>& a, const IRegularMatrix<B, double>& b,
	          double beta, IRegularMatrix<C, double>& c, char transa='N', char transb='N')
	{
		blas_int m, n, k;
		internal::gemm_get_dims(a, b, c, transa, transb, m, n, k);

		blas_int lda = (blas_int)a.col_stride();
		blas_int ldb = (blas_int)b.col_stride();
		blas_int ldc = (blas_int)c.col_stride();

		LMAT_BLAS_NAME(dgemm)(&transa, &transb, &m, &n, &k, &alpha,
				a.ptr_data(), &lda, b.ptr_data(), &ldb, &beta, c.ptr_data(), &ldc);
	}


	template<class A, class B, class C>
	LMAT_ENSURE_INLINE
	inline void gemm(const IRegularMatrix<A, float>& a, const IRegularMatrix<B, float>& b,
	          IRegularMatrix<C, float>& c, char transa='N', char transb='N')
	{
		gemm(1.0f, a, b, 0.0f, c, transa, transb);
	}

	template<class A, class B, class C>
	LMAT_ENSURE_INLINE
	inline void gemm(const IRegularMatrix<A, double>& a, const IRegularMatrix<B, double>& b,
	          IRegularMatrix<C, double>& c, char transa='N', char transb='N')
	{
		gemm(1.0, a, b, 0.0, c, transa, transb);
	}


	// symm

	namespace internal
	{
		template<class A, class B, class C>
		LMAT_ENSURE_INLINE
		inline void spmm_get_dims(const A& a, const B& b, const C& c, char side,
				blas_int& m, blas_int& n)
		{
			index_t ma = a.nrows();
			index_t na = a.ncolumns();

			index_t mb = b.nrows();
			index_t nb = b.ncolumns();

			if (side == 'L' || side == 'l')
			{
				LMAT_CHECK_DIMS( ma == na && na == mb && ma == c.nrows() && nb == c.ncolumns() );

				m = (blas_int)ma;
				n = (blas_int)nb;
			}
			else
			{
				LMAT_CHECK_DIMS( ma == na && ma == nb && mb == c.nrows() && na == c.ncolumns());

				m = mb;
				n = na;
			}
		}

		template<class A, class B>
		LMAT_ENSURE_INLINE
		inline void spmm_get_dims(const A& a, const B& b, char side,
				blas_int& m, blas_int& n)
		{
			index_t ma = a.nrows();
			index_t na = a.ncolumns();

			index_t mb = b.nrows();
			index_t nb = b.ncolumns();

			if (side == 'L' || side == 'l')
			{
				LMAT_CHECK_DIMS( ma == na && na == mb );

				m = (blas_int)ma;
				n = (blas_int)nb;
			}
			else
			{
				LMAT_CHECK_DIMS( ma == na && ma == nb );

				m = mb;
				n = na;
			}
		}
	}

	template<class A, class B, class C>
	LMAT_ENSURE_INLINE
	inline void symm(float alpha, const IRegularMatrix<A, float>& a, const IRegularMatrix<B, float>& b,
	                 float beta, IRegularMatrix<C, float>& c, char side='L', char uplo='L')
	{
		blas_int m, n;
		internal::spmm_get_dims(a, b, c, side, m, n);

		blas_int lda = (blas_int)a.col_stride();
		blas_int ldb = (blas_int)b.col_stride();
		blas_int ldc = (blas_int)c.col_stride();

		LMAT_BLAS_NAME(ssymm)(&side, &uplo, &m, &n,
				&alpha, a.ptr_data(), &lda, b.ptr_data(), &ldb, &beta, c.ptr_data(), &ldc);
	}

	template<class A, class B, class C>
	LMAT_ENSURE_INLINE
	inline void symm(double alpha, const IRegularMatrix<A, double>& a, const IRegularMatrix<B, double>& b,
	                 double beta, IRegularMatrix<C, double>& c, char side='L', char uplo='L')
	{
		blas_int m, n;
		internal::spmm_get_dims(a, b, c, side, m, n);

		blas_int lda = (blas_int)a.col_stride();
		blas_int ldb = (blas_int)b.col_stride();
		blas_int ldc = (blas_int)c.col_stride();

		LMAT_BLAS_NAME(dsymm)(&side, &uplo, &m, &n,
				&alpha, a.ptr_data(), &lda, b.ptr_data(), &ldb, &beta, c.ptr_data(), &ldc);
	}

	template<class A, class B, class C>
	LMAT_ENSURE_INLINE
	inline void symm(const IRegularMatrix<A, float>& a, const IRegularMatrix<B, float>& b,
			         IRegularMatrix<C, float>& c, char side='L', char uplo='L')
	{
		symm(1.0f, a, b, 0.0f, c, side, uplo);
	}

	template<class A, class B, class C>
	LMAT_ENSURE_INLINE
	inline void symm(const IRegularMatrix<A, double>& a, const IRegularMatrix<B, double>& b,
			         IRegularMatrix<C, double>& c, char side='L', char uplo='L')
	{
		symm(1.0f, a, b, 0.0f, c, side, uplo);
	}


	// trmm

	template<class A, class B>
	LMAT_ENSURE_INLINE
	inline void trmm(float alpha, const IRegularMatrix<A, float>& a, IRegularMatrix<B, float>& b,
	                 const trs& ts, char side='L', char uplo='L')
	{
		blas_int m, n;
		internal::spmm_get_dims(a, b, side, m, n);

		blas_int lda = (blas_int)a.col_stride();
		blas_int ldb = (blas_int)b.col_stride();

		LMAT_BLAS_NAME(strmm)(&side, &(ts.uplo), &(ts.trans), &(ts.diag),
				&m, &n, &alpha, a.ptr_data(), &lda, b.ptr_data(), &ldb);
	}

	template<class A, class B>
	LMAT_ENSURE_INLINE
	inline void trmm(double alpha, const IRegularMatrix<A, double>& a, IRegularMatrix<B, double>& b,
	                 const trs& ts, char side='L', char uplo='L')
	{
		blas_int m, n;
		internal::spmm_get_dims(a, b, side, m, n);

		blas_int lda = (blas_int)a.col_stride();
		blas_int ldb = (blas_int)b.col_stride();

		LMAT_BLAS_NAME(dtrmm)(&side, &(ts.uplo), &(ts.trans), &(ts.diag),
				&m, &n, &alpha, a.ptr_data(), &lda, b.ptr_data(), &ldb);
	}

	template<class A, class B>
	LMAT_ENSURE_INLINE
	inline void trmm(const IRegularMatrix<A, float>& a, IRegularMatrix<B, float>& b,
	                 const trs& ts, char side='L', char uplo='L')
	{
		trmm(1.0f, a, b, ts, side, uplo);
	}

	template<class A, class B>
	LMAT_ENSURE_INLINE
	inline void trmm(const IRegularMatrix<A, double>& a, IRegularMatrix<B, double>& b,
	                 const trs& ts, char side='L', char uplo='L')
	{
		trmm(1.0, a, b, ts, side, uplo);
	}


	// trsm

	template<class A, class B>
	LMAT_ENSURE_INLINE
	inline void trsm(float alpha, const IRegularMatrix<A, float>& a, IRegularMatrix<B, float>& b,
	                 const trs& ts, char side='L', char uplo='L')
	{
		blas_int m, n;
		internal::spmm_get_dims(a, b, side, m, n);

		blas_int lda = (blas_int)a.col_stride();
		blas_int ldb = (blas_int)b.col_stride();

		LMAT_BLAS_NAME(strsm)(&side, &(ts.uplo), &(ts.trans), &(ts.diag),
				&m, &n, &alpha, a.ptr_data(), &lda, b.ptr_data(), &ldb);
	}

	template<class A, class B>
	LMAT_ENSURE_INLINE
	inline void trsm(double alpha, const IRegularMatrix<A, double>& a, IRegularMatrix<B, double>& b,
	                 const trs& ts, char side='L', char uplo='L')
	{
		blas_int m, n;
		internal::spmm_get_dims(a, b, side, m, n);

		blas_int lda = (blas_int)a.col_stride();
		blas_int ldb = (blas_int)b.col_stride();

		LMAT_BLAS_NAME(dtrsm)(&side, &(ts.uplo), &(ts.trans), &(ts.diag),
				&m, &n, &alpha, a.ptr_data(), &lda, b.ptr_data(), &ldb);
	}

	template<class A, class B>
	LMAT_ENSURE_INLINE
	inline void trsm(const IRegularMatrix<A, float>& a, IRegularMatrix<B, float>& b,
	                 const trs& ts, char side='L', char uplo='L')
	{
		trsm(1.0f, a, b, ts, side, uplo);
	}

	template<class A, class B>
	LMAT_ENSURE_INLINE
	inline void trsm(const IRegularMatrix<A, double>& a, IRegularMatrix<B, double>& b,
	                 const trs& ts, char side='L', char uplo='L')
	{
		trsm(1.0, a, b, ts, side, uplo);
	}


} }

#endif /* BLAS_L3_H_ */
