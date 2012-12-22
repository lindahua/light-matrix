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
#include <light_mat/linalg/blas_extern.h>
#include "internal/linalg_aux.h"

namespace lmat { namespace blas {

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


} }

#endif /* BLAS_L3_H_ */
