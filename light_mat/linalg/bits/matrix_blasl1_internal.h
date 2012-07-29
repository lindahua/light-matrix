/**
 * @file matrix_blasl1_internal.h
 *
 * Internal implementation of BLAS Level-1 wrappers
 *
 * @author Dahua Lin
 */

#ifdef _MSC_VER
#pragma once
#endif

#ifndef LIGHTMAT_MATRIX_BLASL1_INTERNAL_H_
#define LIGHTMAT_MATRIX_BLASL1_INTERNAL_H_

#include <light_mat/linalg/linalg_base.h>
#include <light_mat/linalg/blas_extern.h>

namespace lmat { namespace detail {

	template<class Mat>
	class contmat_blasl1_wrapper
	{
	public:
		typedef typename matrix_traits<Mat>::value_type T;

		LMAT_ENSURE_INLINE
		contmat_blasl1_wrapper(const Mat& mat)
		: m_len((lmat_blas_int)mat.nelems())
		, m_inc(1)
		, m_pdata(mat.ptr_data())
		{ }

		LMAT_ENSURE_INLINE
		const lmat_blas_int* plength() const { return &m_len; }

		LMAT_ENSURE_INLINE
		const lmat_blas_int* pinc() const { return &m_inc; }

		LMAT_ENSURE_INLINE
		const T* pdata() const { return m_pdata; }

	private:
		lmat_blas_int m_len;
		lmat_blas_int m_inc;
		const T* m_pdata;
	};

	template<class Mat>
	class regrow_blasl1_wrapper
	{
	public:
		typedef typename matrix_traits<Mat>::value_type T;

		LMAT_ENSURE_INLINE
		regrow_blasl1_wrapper(const Mat& mat)
		: m_len((lmat_blas_int)mat.ncolumns())
		, m_inc((lmat_blas_int)mat.lead_dim())
		, m_pdata(mat.ptr_data())
		{ }

		LMAT_ENSURE_INLINE
		const lmat_blas_int* plength() const { return &m_len; }

		LMAT_ENSURE_INLINE
		const lmat_blas_int* pinc() const { return &m_inc; }

		LMAT_ENSURE_INLINE
		const T* pdata() const { return m_pdata; }

	private:
		lmat_blas_int m_len;
		lmat_blas_int m_inc;
		const T* m_pdata;
	};

	template<class Mat>
	class dense_blasl1_wrapper
	{
	public:
		typedef typename matrix_traits<Mat>::value_type T;

		LMAT_ENSURE_INLINE
		dense_blasl1_wrapper(const Mat& mat)
		: m_len((lmat_blas_int)mat.nelems())
		{
			if (has_continuous_layout(mat))
			{
				m_inc = 1;
				m_own = false;
				m_pdata = const_cast<T*>(mat.ptr_data());
			}
			else if (is_row(mat))
			{
				m_inc = (lmat_blas_int)mat.lead_dim();
				m_own = false;
				m_pdata = const_cast<T*>(mat.ptr_data());
			}
			else
			{
				m_inc = 1;
				m_own = true;
				m_pdata = (T*)aligned_allocate((size_t)m_len * sizeof(T), LMAT_DEFAULT_ALIGNMENT);
				copy(mat, m_pdata);
			}
		}

		LMAT_ENSURE_INLINE
		~dense_blasl1_wrapper()
		{
			if (m_own) aligned_release(m_pdata);
		}

		LMAT_ENSURE_INLINE
		const lmat_blas_int* plength() const { return &m_len; }

		LMAT_ENSURE_INLINE
		const lmat_blas_int* pinc() const { return &m_inc; }

		LMAT_ENSURE_INLINE
		const T* pdata() const { return m_pdata; }

	private:
		lmat_blas_int m_len;
		lmat_blas_int m_inc;
		bool m_own;
		T* m_pdata;
	};

	template<class Mat>
	class generic_blasl1_wrapper
	{
	public:
		typedef typename matrix_traits<Mat>::value_type T;

		LMAT_ENSURE_INLINE
		generic_blasl1_wrapper(const Mat& mat)
		: m_cache(mat)
		, m_len((lmat_blas_int)m_cache.nelems())
		, m_inc(1)
		{ }

		LMAT_ENSURE_INLINE
		const lmat_blas_int* plength() const { return &m_len; }

		LMAT_ENSURE_INLINE
		const lmat_blas_int* pinc() const { return &m_inc; }

		LMAT_ENSURE_INLINE
		const T* pdata() const { return m_cache.ptr_data(); }

	private:
		dense_matrix<T, ct_rows<Mat>::value, ct_cols<Mat>::value> m_cache;
		lmat_blas_int m_len;
		lmat_blas_int m_inc;
	};


	template<class Mat>
	struct blasl1_wrapper_map
	{
		static const bool is_contmat = ct_is_col<Mat>::value || ct_has_continuous_layout<Mat>::value;

		typedef typename
				if_<is_dense_mat<Mat>,
					typename
					if_c<is_contmat,
						contmat_blasl1_wrapper<Mat>,
						typename
						if_<ct_is_row<Mat>,
							regrow_blasl1_wrapper<Mat>,
							dense_blasl1_wrapper<Mat>
						>::type
					>::type,
					generic_blasl1_wrapper<Mat>
				>::type type;
	};



	template<typename T, int M, int N>
	struct blasl1_internal;

	template<int M, int N>
	struct blasl1_internal<float, M, N>
	{
		template<class Mat>
		LMAT_ENSURE_INLINE
		static float asum(const IMatrixXpr<Mat, float>& x)
		{
			typename blasl1_wrapper_map<Mat>::type xw(x.derived());
			return LMAT_SASUM(xw.plength(), xw.pdata(), xw.pinc());
		}

		template<class XMat, class YMat>
		LMAT_ENSURE_INLINE
		static float axpy(const float a, const IMatrixXpr<XMat, float>& x,
				IDenseMatrix<YMat, float>& y)
		{
			check_arg(has_same_size(x, y), "blas::axpy: x and y should have the same size.");
			check_arg(has_continuous_layout(y), "blas::axpy: y should have continuous memory layout.");

			typename blasl1_wrapper_map<XMat>::type xw(x.derived());
			const lmat_blas_int incy = 1;
			LMAT_SAXPY(xw.plength(), &a, xw.pdata(), xw.pinc(), y.ptr_data(), &incy);
		}
	};

	template<int M, int N>
	struct blasl1_internal<double, M, N>
	{
		template<class Mat>
		LMAT_ENSURE_INLINE
		static double asum(const IMatrixXpr<Mat, double>& x)
		{
			typename blasl1_wrapper_map<Mat>::type xw(x.derived());
			return LMAT_DASUM(xw.plength(), xw.pdata(), xw.pinc());
		}

		template<class XMat, class YMat>
		LMAT_ENSURE_INLINE
		static double axpy(const double a, const IMatrixXpr<XMat, double>& x,
				IDenseMatrix<YMat, double>& y)
		{
			check_arg(has_same_size(x, y), "blas::axpy: x and y should have the same size.");
			check_arg(has_continuous_layout(y), "blas::axpy: y should have continuous memory layout.");

			typename blasl1_wrapper_map<XMat>::type xw(x.derived());
			const lmat_blas_int incy = 1;
			LMAT_SAXPY(xw.plength(), &a, xw.pdata(), xw.pinc(), y.ptr_data(), &incy);
		}
	};


	template<class Mat>
	struct unary_blasl1_internal_map
	{
		typedef typename matrix_traits<Mat>::value_type T;
		static const int M = ct_rows<Mat>::value;
		static const int N = ct_cols<Mat>::value;
		typedef blasl1_internal<T, M, N> type;
	};

	template<class LMat, class RMat>
	struct binary_blasl1_internal_map
	{
		typedef typename binary_value_type<LMat, RMat>::type T;
		static const int M = binary_ct_rows<LMat, RMat>::value;
		static const int N = binary_ct_cols<LMat, RMat>::value;
		typedef blasl1_internal<T, M, N> type;
	};

} }

#endif /* MATRIX_BLASL1_INTERNAL_H_ */
