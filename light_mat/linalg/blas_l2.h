/**
 * @file blas_l2.h
 *
 * @brief BLAS Level 2 functions
 *
 * @author Dahua Lin
 */

#ifndef LIGHTMAT_BLAS_L2_H_
#define LIGHTMAT_BLAS_L2_H_

#include "internal/linalg_aux.h"

extern "C"
{
	void LMAT_BLAS_NAME(sgemv)(const char *trans, const blas_int *m, const blas_int *n, const float *alpha,
	           const float *a, const blas_int *lda, const float *x, const blas_int *incx,
	           const float *beta, float *y, const blas_int *incy);
	void LMAT_BLAS_NAME(sger)(const blas_int *m, const blas_int *n, const float *alpha, const float *x, const blas_int *incx,
	          const float *y, const blas_int *incy, float *a, const blas_int *lda);
	void LMAT_BLAS_NAME(ssymv)(const char *uplo, const blas_int *n, const float *alpha, const float *a, const blas_int *lda,
	           const float *x, const blas_int *incx, const float *beta, float *y, const blas_int *incy);
	void LMAT_BLAS_NAME(ssyr)(const char *uplo, const blas_int *n, const float *alpha, const float *x, const blas_int *incx,
	          float *a, const blas_int *lda);
	void LMAT_BLAS_NAME(strmv)(const char *uplo, const char *transa, const char *diag, const blas_int *n, const float *a,
	           const blas_int *lda, float *b, const blas_int *incx);
	void LMAT_BLAS_NAME(strsv)(const char *uplo, const char *trans, const char *diag, const blas_int *n,
	           const float *a, const blas_int *lda, float *x, const blas_int *incx);

	void LMAT_BLAS_NAME(dgemv)(const char *trans, const blas_int *m, const blas_int *n, const double *alpha,
	           const double *a, const blas_int *lda, const double *x, const blas_int *incx,
	           const double *beta, double *y, const blas_int *incy);
	void LMAT_BLAS_NAME(dger)(const blas_int *m, const blas_int *n, const double *alpha, const double *x, const blas_int *incx,
	          const double *y, const blas_int *incy, double *a, const blas_int *lda);
	void LMAT_BLAS_NAME(dsymv)(const char *uplo, const blas_int *n, const double *alpha, const double *a, const blas_int *lda,
	           const double *x, const blas_int *incx, const double *beta, double *y, const blas_int *incy);
	void LMAT_BLAS_NAME(dsyr)(const char *uplo, const blas_int *n, const double *alpha, const double *x, const blas_int *incx,
	          double *a, const blas_int *lda);
	void LMAT_BLAS_NAME(dtrmv)(const char *uplo, const char *transa, const char *diag, const blas_int *n,
	           const double *a, const blas_int *lda, double *b, const blas_int *incx);
	void LMAT_BLAS_NAME(dtrsv)(const char *uplo, const char *trans, const char *diag, const blas_int *n,
	           const double *a, const blas_int *lda, double *x, const blas_int *incx);
}


namespace lmat { namespace blas {

	// gemv

	namespace internal
	{
		template<class A, class X, class Y>
		LMAT_ENSURE_INLINE
		inline bool gemv_check_dims(const A& a, const X& x, const Y& y, char trans)
		{
			if (trans == 'N' || trans == 'n')
			{
				return a.ncolumns() == x.nelems() && a.nrows() == y.nelems();
			}
			else
			{
				return a.nrows() == x.nelems() && a.ncolumns() == y.nelems();
			}
		}
	}


	template<class A, class X, class Y>
	LMAT_ENSURE_INLINE
	inline void gemv(
			float alpha, const IRegularMatrix<A, float>& a, const IRegularMatrix<X, float>& x,
			float beta, IRegularMatrix<Y, float>& y, char trans='N')
	{
		LMAT_CHECK_DIMS( internal::gemv_check_dims(a, x, y, trans) )

		blas_int lda = (blas_int)a.col_stride();
		blas_int incx = (blas_int)lmat::internal::get_vector_intv(x);
		blas_int incy = (blas_int)lmat::internal::get_vector_intv(y);

		blas_int m = (blas_int)a.nrows();
		blas_int n = (blas_int)a.ncolumns();

		LMAT_BLAS_NAME(sgemv)(&trans, &m, &n, &alpha, a.ptr_data(), &lda, x.ptr_data(), &incx, &beta, y.ptr_data(), &incy);
	}

	template<class A, class X, class Y>
	LMAT_ENSURE_INLINE
	inline void gemv(
			double alpha, const IRegularMatrix<A, double>& a, const IRegularMatrix<X, double>& x,
			double beta, IRegularMatrix<Y, double>& y, char trans='N')
	{
		LMAT_CHECK_DIMS( internal::gemv_check_dims(a, x, y, trans) )

		blas_int lda = (blas_int)a.col_stride();
		blas_int incx = (blas_int)lmat::internal::get_vector_intv(x);
		blas_int incy = (blas_int)lmat::internal::get_vector_intv(y);

		blas_int m = (blas_int)a.nrows();
		blas_int n = (blas_int)a.ncolumns();

		LMAT_BLAS_NAME(dgemv)(&trans, &m, &n, &alpha, a.ptr_data(), &lda, x.ptr_data(), &incx, &beta, y.ptr_data(), &incy);
	}


	template<class A, class X, class Y>
	LMAT_ENSURE_INLINE
	inline void gemv(
			const IRegularMatrix<A, float>& a, const IRegularMatrix<X, float>& x, IRegularMatrix<Y, float>& y, char trans='N')
	{
		gemv(1.0f, a, x, 0.0f, y, trans);
	}

	template<class A, class X, class Y>
	LMAT_ENSURE_INLINE
	inline void gemv(
			const IRegularMatrix<A, double>& a, const IRegularMatrix<X, double>& x, IRegularMatrix<Y, double>& y, char trans='N')
	{
		gemv(1.0, a, x, 0.0, y, trans);
	}


	// symv

	namespace internal
	{
		template<class A, class X>
		LMAT_ENSURE_INLINE
		inline bool sqmv_check_dims(const A& a, const X& x)
		{
			index_t n = a.nrows();
			return a.ncolumns() == n && x.nelems() == n;
		}

		template<class A, class X, class Y>
		LMAT_ENSURE_INLINE
		inline bool sqmv_check_dims(const A& a, const X& x, const Y& y)
		{
			index_t n = a.nrows();
			return a.ncolumns() == n && x.nelems() == n && y.nelems() == n;
		}
	}


	template<class A, class X, class Y>
	LMAT_ENSURE_INLINE
	inline void symv(
			float alpha, const IRegularMatrix<A, float>& a, const IRegularMatrix<X, float>& x,
			float beta, IRegularMatrix<Y, float>& y, char uplo='L')
	{
		LMAT_CHECK_DIMS( internal::sqmv_check_dims(a, x, y) )

		blas_int lda = (blas_int)a.col_stride();
		blas_int incx = (blas_int)lmat::internal::get_vector_intv(x);
		blas_int incy = (blas_int)lmat::internal::get_vector_intv(y);

		blas_int n = (blas_int)a.nrows();

		LMAT_BLAS_NAME(ssymv)(&uplo, &n, &alpha, a.ptr_data(), &lda, x.ptr_data(), &incx, &beta, y.ptr_data(), &incy);
	}

	template<class A, class X, class Y>
	LMAT_ENSURE_INLINE
	inline void symv(
			double alpha, const IRegularMatrix<A, double>& a, const IRegularMatrix<X, double>& x,
			double beta, IRegularMatrix<Y, double>& y, char uplo='L')
	{
		LMAT_CHECK_DIMS( internal::sqmv_check_dims(a, x, y) )

		blas_int lda = (blas_int)a.col_stride();
		blas_int incx = (blas_int)lmat::internal::get_vector_intv(x);
		blas_int incy = (blas_int)lmat::internal::get_vector_intv(y);

		blas_int n = (blas_int)a.nrows();

		LMAT_BLAS_NAME(dsymv)(&uplo, &n, &alpha, a.ptr_data(), &lda, x.ptr_data(), &incx, &beta, y.ptr_data(), &incy);
	}

	template<class A, class X, class Y>
	LMAT_ENSURE_INLINE
	inline void symv(
			const IRegularMatrix<A, float>& a, const IRegularMatrix<X, float>& x, IRegularMatrix<Y, float>& y, char uplo='L')
	{
		symv(1.0f, a, x, 0.0f, y, uplo);
	}

	template<class A, class X, class Y>
	LMAT_ENSURE_INLINE
	inline void symv(
			const IRegularMatrix<A, double>& a, const IRegularMatrix<X, double>& x, IRegularMatrix<Y, double>& y, char uplo='L')
	{
		symv(1.0, a, x, 0.0, y, uplo);
	}


	// trmv

	template<class A, class X>
	LMAT_ENSURE_INLINE
	inline void trmv(const IRegularMatrix<A, float>& a, IRegularMatrix<X, float>& x, const trs& ts)
	{
		LMAT_CHECK_DIMS( internal::sqmv_check_dims(a, x) )

		blas_int lda = (blas_int)a.col_stride();
		blas_int incx = (blas_int)lmat::internal::get_vector_intv(x);
		blas_int n = (blas_int)a.nrows();

		LMAT_BLAS_NAME(strmv)(&(ts.uplo), &(ts.trans), &(ts.diag), &n, a.ptr_data(), &lda, x.ptr_data(), &incx);
	}

	template<class A, class X>
	LMAT_ENSURE_INLINE
	inline void trmv(const IRegularMatrix<A, double>& a, IRegularMatrix<X, double>& x, const trs& ts)
	{
		LMAT_CHECK_DIMS( internal::sqmv_check_dims(a, x) )

		blas_int lda = (blas_int)a.col_stride();
		blas_int incx = (blas_int)lmat::internal::get_vector_intv(x);
		blas_int n = (blas_int)a.nrows();

		LMAT_BLAS_NAME(dtrmv)(&(ts.uplo), &(ts.trans), &(ts.diag), &n, a.ptr_data(), &lda, x.ptr_data(), &incx);
	}


	// trsv

	template<class A, class X>
	LMAT_ENSURE_INLINE
	inline void trsv(const IRegularMatrix<A, float>& a, IRegularMatrix<X, float>& x, const trs& ts)
	{
		LMAT_CHECK_DIMS( internal::sqmv_check_dims(a, x) )

		blas_int lda = (blas_int)a.col_stride();
		blas_int incx = (blas_int)lmat::internal::get_vector_intv(x);
		blas_int n = (blas_int)a.nrows();

		LMAT_BLAS_NAME(strsv)(&(ts.uplo), &(ts.trans), &(ts.diag), &n, a.ptr_data(), &lda, x.ptr_data(), &incx);
	}

	template<class A, class X>
	LMAT_ENSURE_INLINE
	inline void trsv(const IRegularMatrix<A, double>& a, IRegularMatrix<X, double>& x, const trs& ts)
	{
		LMAT_CHECK_DIMS( internal::sqmv_check_dims(a, x) )

		blas_int lda = (blas_int)a.col_stride();
		blas_int incx = (blas_int)lmat::internal::get_vector_intv(x);
		blas_int n = (blas_int)a.nrows();

		LMAT_BLAS_NAME(dtrsv)(&(ts.uplo), &(ts.trans), &(ts.diag), &n, a.ptr_data(), &lda, x.ptr_data(), &incx);
	}


	// ger

	namespace internal
	{
		template<class A, class X, class Y>
		LMAT_ENSURE_INLINE
		inline bool ger_check_dims(const A& a, const X& x, const Y& y)
		{
			return a.nrows() == x.nelems() && a.ncolumns() == y.nelems();
		}
	}


	template<class A, class X, class Y>
	LMAT_ENSURE_INLINE
	inline void ger(IRegularMatrix<A, float>& a, const IRegularMatrix<X, float>& x, const IRegularMatrix<Y, float>& y,
			float alpha = 1.0f)
	{
		LMAT_CHECK_DIMS( internal::ger_check_dims(a, x, y) )

		blas_int m = (blas_int)a.nrows();
		blas_int n = (blas_int)a.ncolumns();

		blas_int lda = (blas_int)a.col_stride();
		blas_int incx = (blas_int)lmat::internal::get_vector_intv(x);
		blas_int incy = (blas_int)lmat::internal::get_vector_intv(y);

		LMAT_BLAS_NAME(sger)(&m, &n, &alpha, x.ptr_data(), &incx, y.ptr_data(), &incy, a.ptr_data(), &lda);
	}


	template<class A, class X, class Y>
	LMAT_ENSURE_INLINE
	inline void ger(IRegularMatrix<A, double>& a, const IRegularMatrix<X, double>& x, const IRegularMatrix<Y, double>& y,
			double alpha = 1.0f)
	{
		LMAT_CHECK_DIMS( internal::ger_check_dims(a, x, y) )

		blas_int m = (blas_int)a.nrows();
		blas_int n = (blas_int)a.ncolumns();

		blas_int lda = (blas_int)a.col_stride();
		blas_int incx = (blas_int)lmat::internal::get_vector_intv(x);
		blas_int incy = (blas_int)lmat::internal::get_vector_intv(y);

		LMAT_BLAS_NAME(dger)(&m, &n, &alpha, x.ptr_data(), &incx, y.ptr_data(), &incy, a.ptr_data(), &lda);
	}


} }


#endif /* BLAS_L2_H_ */
