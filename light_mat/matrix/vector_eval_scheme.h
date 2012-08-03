/*
 * @file vector_eval_scheme.h
 *
 * Vector evaluation schemes
 *
 * @author Dahua Lin
 */

#ifndef _MSC_VER
#pragma once
#endif

#ifndef LIGHTMAT_VECTOR_EVAL_SCHEME_H_
#define LIGHTMAT_VECTOR_EVAL_SCHEME_H_

#include <light_mat/matrix/matrix_properties.h>
#include <light_mat/matrix/matrix_classes.h>

#include <light_mat/matrix/vector_eval_kernel.h>

namespace lmat
{
	// vector evaluation scheme type

	struct as_linear_vec { };
	struct per_column { };


	/********************************************
	 *
	 *  Interfaces
	 *
	 ********************************************/

	template<class Sch> struct vector_eval_scheme_traits;

	template<class Derived, typename T>
	class IVecEvalLinearScheme
	{
	public:
		LMAT_CRTP_REF

		typedef typename vector_eval_scheme_traits<Derived>::kernel_type kernel_type;
		typedef typename vector_eval_kernel_state<kernel_type>::type kernel_state_t;

		LMAT_ENSURE_INLINE
		kernel_type kernel() const
		{
			return derived().kernel();
		}

		LMAT_ENSURE_INLINE
		kernel_state_t vec_state() const
		{
			return derived().vec_state();
		}
	};


	template<class Derived, typename T>
	class IVecEvalPerColScheme
	{
	public:
		LMAT_CRTP_REF

		typedef typename vector_eval_scheme_traits<Derived>::kernel_type kernel_type;
		typedef typename vector_eval_kernel_state<kernel_type>::type kernel_state_t;

		LMAT_ENSURE_INLINE
		kernel_type kernel() const
		{
			return derived().kernel();
		}

		LMAT_ENSURE_INLINE
		kernel_state_t col_state(const index_t j) const
		{
			return derived().col_state(j);
		}
	};


	/********************************************
	 *
	 *  Scalar schemes
	 *
	 ********************************************/

	// forward declarations

	template<typename T> class continuous_linear_veval_scheme;
	template<typename T> class dense_percol_veval_scheme;
	template<typename T> class const_linear_veval_scheme;
	template<typename T> class const_percol_veval_scheme;
	template<typename T, int M, int N> class cached_linear_veval_scheme;
	template<typename T, int M, int N> class cached_percol_veval_scheme;


	// continuous linear

	template<typename T>
	struct vector_eval_scheme_traits<continuous_linear_veval_scheme<T> >
	{
		typedef veval_memacc_kernel<T> kernel_type;
	};

	template<typename T>
	class continuous_linear_veval_scheme
	: public IVecEvalLinearScheme<continuous_linear_veval_scheme<T>, T>
	{
	public:
		typedef veval_memacc_kernel<T> kernel_type;
		typedef typename vector_eval_kernel_state<kernel_type>::type state_t;

		template<class Mat>
		LMAT_ENSURE_INLINE
		continuous_linear_veval_scheme(const IDenseMatrix<Mat, T>& X)
		{
#ifdef LMAT_USE_STATIC_ASSERT
			static_assert(ct_has_continuous_layout<Mat>::value,
					"Mat must always have continuous layout");
#endif
			m_data = X.ptr_data();
		}

		LMAT_ENSURE_INLINE kernel_type kernel() const
		{
			return kernel_type();
		}

		LMAT_ENSURE_INLINE state_t vec_state() const
		{
			return m_data;
		}

	private:
		const T *m_data;
	};


	// dense percol

	template<typename T>
	struct vector_eval_scheme_traits<dense_percol_veval_scheme<T> >
	{
		typedef veval_memacc_kernel<T> kernel_type;
	};

	template<typename T>
	class dense_percol_veval_scheme
	: public IVecEvalPerColScheme<dense_percol_veval_scheme<T>, T>
	{
	public:
		typedef veval_memacc_kernel<T> kernel_type;
		typedef typename vector_eval_kernel_state<kernel_type>::type state_t;

		template<class Mat>
		LMAT_ENSURE_INLINE
		dense_percol_veval_scheme(const IDenseMatrix<Mat, T>& X)
		: m_ldim(X.lead_dim())
		, m_data(X.ptr_data())
		{
		}

		LMAT_ENSURE_INLINE kernel_type kernel() const
		{
			return kernel_type();
		}

		LMAT_ENSURE_INLINE state_t col_state(const index_t j) const
		{
			return m_data + m_ldim * j;
		}
	private:
		const index_t m_ldim;
		const T *m_data;
	};


	// const linear

	template<typename T>
	struct vector_eval_scheme_traits<const_linear_veval_scheme<T> >
	{
		typedef veval_const_kernel<T> kernel_type;
	};

	template<typename T>
	class const_linear_veval_scheme
	: public IVecEvalLinearScheme<const_linear_veval_scheme<T>, T>
	{
	public:
		typedef veval_const_kernel<T> kernel_type;
		typedef typename vector_eval_kernel_state<kernel_type>::type state_t;

		template<int CTRows, int CTCols>
		LMAT_ENSURE_INLINE
		const_linear_veval_scheme(const const_matrix<T, CTRows, CTCols>& X)
		: m_val(X.value())
		{
		}

		LMAT_ENSURE_INLINE kernel_type kernel() const
		{
			return kernel_type();
		}

		LMAT_ENSURE_INLINE state_t vec_state() const
		{
			return m_val;
		}
	private:
		const T m_val;
	};


	// const percol

	template<typename T>
	struct vector_eval_scheme_traits<const_percol_veval_scheme<T> >
	{
		typedef veval_const_kernel<T> kernel_type;
	};

	template<typename T>
	class const_percol_veval_scheme
	: public IVecEvalPerColScheme<const_percol_veval_scheme<T>, T>
	{
	public:
		typedef veval_const_kernel<T> kernel_type;
		typedef typename vector_eval_kernel_state<kernel_type>::type state_t;

		template<int CTRows, int CTCols>
		LMAT_ENSURE_INLINE
		const_percol_veval_scheme(const const_matrix<T, CTRows, CTCols>& X)
		: m_val(X.value())
		{
		}

		LMAT_ENSURE_INLINE kernel_type kernel() const
		{
			return kernel_type();
		}

		LMAT_ENSURE_INLINE state_t col_state(const index_t j) const
		{
			return m_val;
		}

	private:
		const T m_val;
	};


	// cached linear

	template<typename T, int M, int N>
	struct vector_eval_scheme_traits<cached_linear_veval_scheme<T, M, N> >
	{
		typedef veval_memacc_kernel<T> kernel_type;
	};

	template<typename T, int M, int N>
	class cached_linear_veval_scheme
	: public IVecEvalLinearScheme<cached_linear_veval_scheme<T, M, N>, T>
	{
	public:
		typedef veval_memacc_kernel<T> kernel_type;
		typedef typename vector_eval_kernel_state<kernel_type>::type state_t;

		template<class Expr>
		LMAT_ENSURE_INLINE
		cached_linear_veval_scheme(const IMatrixXpr<Expr, T>& X)
		: m_cache(X)
		{
		}

		LMAT_ENSURE_INLINE kernel_type kernel() const
		{
			return kernel_type();
		}

		LMAT_ENSURE_INLINE state_t vec_state() const
		{
			return m_cache.ptr_data();
		}
	private:
		dense_matrix<T, M, N> m_cache;
	};


	// cached percol

	template<typename T, int M, int N>
	struct vector_eval_scheme_traits<cached_percol_veval_scheme<T, M, N> >
	{
		typedef veval_memacc_kernel<T> kernel_type;
	};

	template<typename T, int M, int N>
	class cached_percol_veval_scheme
	: public IVecEvalPerColScheme<cached_percol_veval_scheme<T, M, N>, T>
	{
	public:
		typedef veval_memacc_kernel<T> kernel_type;
		typedef typename vector_eval_kernel_state<kernel_type>::type state_t;

		template<class Expr>
		LMAT_ENSURE_INLINE
		cached_percol_veval_scheme(const IMatrixXpr<Expr, T>& X)
		: m_cache(X)
		{
		}

		LMAT_ENSURE_INLINE kernel_type kernel() const
		{
			return kernel_type();
		}

		LMAT_ENSURE_INLINE state_t col_state(const index_t j) const
		{
			return m_cache.ptr_col(j);
		}

	private:
		dense_matrix<T, M, N> m_cache;
	};



}

#endif 
