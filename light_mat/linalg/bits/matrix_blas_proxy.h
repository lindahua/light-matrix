/**
 * @file matrix_blas_proxy.h
 *
 * Matrix proxy classes for BLAS evaluation
 *
 * @author Dahua Lin
 */

#ifdef _MSC_VER
#pragma once
#endif

#ifndef LIGHTMAT_MATRIX_BLAS_PROXY_H_
#define LIGHTMAT_MATRIX_BLAS_PROXY_H_

#include <light_mat/matrix/dense_matrix.h>

namespace lmat { namespace detail {

	template<class Vec>
	class blas_densecol_proxy
	{
	public:
		typedef typename matrix_traits<Vec>::value_type T;

		LMAT_ENSURE_INLINE
		blas_densecol_proxy(const Vec& vec)
		: m_vec(vec) { }

		LMAT_ENSURE_INLINE
		lmat_blas_int length() const
		{
			return (lmat_blas_int)(m_vec.ncolumns());
		}

		LMAT_ENSURE_INLINE
		lmat_blas_int inc() const
		{
			return 1;
		}

		LMAT_ENSURE_INLINE
		const T *pdata() const
		{
			return m_vec.ptr_data();
		}

	private:
		const Vec& m_vec;
	};


	template<class Vec>
	class blas_colxpr_proxy
	{
	public:
		typedef typename matrix_traits<Vec>::value_type T;

		LMAT_ENSURE_INLINE
		blas_colxpr_proxy(const Vec& vec)
		: m_cache(vec) { }

		LMAT_ENSURE_INLINE
		lmat_blas_int length() const
		{
			return (lmat_blas_int)m_cache.ncolumns();
		}

		LMAT_ENSURE_INLINE
		lmat_blas_int inc() const
		{
			return 1;
		}

		LMAT_ENSURE_INLINE
		const T *pdata() const
		{
			return m_cache.ptr_data();
		}

	private:
		dense_matrix<T, ct_rows<Vec>::value, 1> m_cache;
	};


	template<class Vec>
	class blas_denserow_proxy
	{
	public:
		typedef typename matrix_traits<Vec>::value_type T;

		LMAT_ENSURE_INLINE
		blas_denserow_proxy(const Vec& vec)
		: m_vec(vec) { }

		LMAT_ENSURE_INLINE
		lmat_blas_int length() const
		{
			return m_vec.ncolumns();
		}

		LMAT_ENSURE_INLINE
		lmat_blas_int inc() const
		{
			return (lmat_blas_int)m_vec.lead_dim();
		}

		LMAT_ENSURE_INLINE
		const T *pdata() const
		{
			return m_vec.ptr_data();
		}

	private:
		const Vec& m_vec;
	};


	template<class Vec>
	class blas_rowxpr_proxy
	{
	public:
		typedef typename matrix_traits<Vec>::value_type T;

		LMAT_ENSURE_INLINE
		blas_rowxpr_proxy(const Vec& vec)
		: m_cache(vec)
		{ }

		LMAT_ENSURE_INLINE
		lmat_blas_int length() const
		{
			return (lmat_blas_int)m_cache.ncolumns();
		}

		LMAT_ENSURE_INLINE
		lmat_blas_int inc() const
		{
			return 1;
		}

		LMAT_ENSURE_INLINE
		const T *pdata() const
		{
			return m_cache.ptr_data();
		}

	private:
		dense_matrix<T, 1, ct_cols<Vec>::value> m_cache;
	};


	template<class Mat>
	class blas_densemat_proxy
	{
	public:
		typedef typename matrix_traits<Mat>::value_type T;

		LMAT_ENSURE_INLINE
		blas_densemat_proxy(const Mat& mat)
		: m_mat(mat) { }

		LMAT_ENSURE_INLINE
		lmat_blas_int nrows() const
		{
			return (lmat_blas_int)m_mat.nrows();
		}

		LMAT_ENSURE_INLINE
		lmat_blas_int ncolumns() const
		{
			return (lmat_blas_int)m_mat.ncolumns();
		}

		LMAT_ENSURE_INLINE
		lmat_blas_int lead_dim() const
		{
			return (lmat_blas_int)m_mat.lead_dim();
		}

		LMAT_ENSURE_INLINE
		const T *pdata() const
		{
			return m_mat.ptr_data();
		}

	private:
		const Mat& m_mat;
	};


	template<class Mat>
	class blas_matxpr_proxy
	{
	public:
		typedef typename matrix_traits<Mat>::value_type T;

		LMAT_ENSURE_INLINE
		blas_matxpr_proxy(const Mat& mat)
		: m_cache(mat) { }

		LMAT_ENSURE_INLINE
		lmat_blas_int nrows() const
		{
			return (lmat_blas_int)m_cache.nrows();
		}

		LMAT_ENSURE_INLINE
		lmat_blas_int ncolumns() const
		{
			return (lmat_blas_int)m_cache.ncolumns();
		}

		LMAT_ENSURE_INLINE
		lmat_blas_int lead_dim() const
		{
			return (lmat_blas_int)m_cache.lead_dim();
		}

		LMAT_ENSURE_INLINE
		const T *pdata() const
		{
			return m_cache.ptr_data();
		}

	private:
		dense_matrix<T, ct_rows<Mat>::value, ct_cols<Mat>::value> m_cache;
	};


	template<class Vec>
	struct blas_col_proxy
	{
		typedef typename
				if_<is_dense_mat<Vec>,
					blas_densecol_proxy<Vec>,
					blas_colxpr_proxy<Vec> >::type type;

	};

	template<class Vec>
	struct blas_row_proxy
	{
		typedef typename
				if_<is_dense_mat<Vec>,
					blas_denserow_proxy<Vec>,
					blas_rowxpr_proxy<Vec> >::type type;

	};

	template<class Mat>
	struct blas_mat_proxy
	{
		typedef typename
				if_<is_dense_mat<Mat>,
					blas_densemat_proxy<Mat>,
					blas_matxpr_proxy<Mat> >::type type;

	};

} }

#endif /* MATRIX_BLAS_WRAPPER_H_ */
