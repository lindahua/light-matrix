/**
 * @file matrix_vec_evaluators.h
 *
 * Generic vector evaluators
 *
 * @author Dahua Lin
 */

#ifdef _MSC_VER
#pragma once
#endif

#ifndef LIGHTMAT_MATRIX_VEC_EVALUATORS_H_
#define LIGHTMAT_MATRIX_VEC_EVALUATORS_H_

#include <light_mat/matrix/dense_matrix.h>
#include <light_mat/matrix/const_matrix.h>
#include <light_mat/matrix/vec_evaluator_concepts.h>

namespace lmat
{
	/********************************************
	 *
	 *  Evaluator classes
	 *
	 ********************************************/

	template<typename T>
	class continuous_linear_evaluator
	: public ILinearVectorEvaluator<continuous_linear_evaluator<T>, T>
	{
	public:

		template<class Mat>
		LMAT_ENSURE_INLINE
		continuous_linear_evaluator(const IDenseMatrix<Mat, T>& X)
		{
#ifdef LMAT_USE_STATIC_ASSERT
			static_assert(has_continuous_layout<Mat>::value,
					"Mat must always have continuous layout");
#endif
			m_data = X.ptr_data();
		}

		LMAT_ENSURE_INLINE T get_value(const index_t i) const
		{
			return m_data[i];
		}
	private:
		const T *m_data;
	};


	template<typename T>
	class dense_percol_evaluator
	: public IPerColVectorEvaluator<dense_percol_evaluator<T>, T>
	{
	public:
		template<class Mat>
		LMAT_ENSURE_INLINE
		dense_percol_evaluator(const IDenseMatrix<Mat, T>& X)
		: m_ldim(X.lead_dim())
		, m_data(X.ptr_data())
		{
		}

		LMAT_ENSURE_INLINE T get_value(const index_t i) const
		{
			return m_data[i];
		}

		LMAT_ENSURE_INLINE void next_column()
		{
			m_data += m_ldim;
		}
	private:
		const index_t m_ldim;
		const T *m_data;
	};


	template<typename T>
	class const_linear_evaluator
	: public ILinearVectorEvaluator<const_linear_evaluator<T>, T>
	{
	public:
		template<int CTRows, int CTCols>
		LMAT_ENSURE_INLINE
		const_linear_evaluator(const const_matrix<T, CTRows, CTCols>& X)
		: m_val(X.value())
		{
		}

		LMAT_ENSURE_INLINE T get_value(const index_t i) const
		{
			return m_val;
		}

	private:
		const T m_val;
	};


	template<typename T>
	class const_percol_evaluator
	: public IPerColVectorEvaluator<const_percol_evaluator<T>, T>
	{
	public:
		template<int CTRows, int CTCols>
		LMAT_ENSURE_INLINE
		const_percol_evaluator(const const_matrix<T, CTRows, CTCols>& X)
		: m_val(X.value())
		{
		}

		LMAT_ENSURE_INLINE T get_value(const index_t i) const
		{
			return m_val;
		}

		LMAT_ENSURE_INLINE void next_column() { }

	private:
		const T m_val;
	};



	template<typename T>
	class cached_linear_evaluator
	: public ILinearVectorEvaluator<cached_linear_evaluator<T>, T>
	{
	public:
		template<class Expr>
		LMAT_ENSURE_INLINE
		cached_linear_evaluator(const IMatrixXpr<Expr, T>& X)
		: m_cache(X), m_data(m_cache.ptr_data())
		{
		}

		LMAT_ENSURE_INLINE T get_value(const index_t i) const
		{
			return m_data[i];
		}
	private:
		dense_matrix<T> m_cache;
		const T *m_data;
	};


	template<typename T>
	class cached_percol_evaluator
	: public IPerColVectorEvaluator<cached_percol_evaluator<T>, T>
	{
	public:
		template<class Expr>
		LMAT_ENSURE_INLINE
		cached_percol_evaluator(const IMatrixXpr<Expr, T>& X)
		: m_cache(X), m_ldim(m_cache.lead_dim()), m_data(m_cache.ptr_data())
		{
		}

		LMAT_ENSURE_INLINE T get_value(const index_t i) const
		{
			return m_data[i];
		}

		LMAT_ENSURE_INLINE void next_column()
		{
			m_data += m_ldim;
		}

	private:
		dense_matrix<T> m_cache;
		const index_t m_ldim;
		const T *m_data;
	};


	/********************************************
	 *
	 *  Evaluator dispatch
	 *
	 ********************************************/

	template<class Expr>
	struct linear_vector_evaluator
	{
		typedef typename matrix_traits<Expr>::value_type T;

		typedef typename
				if_<and_<is_dense_mat<Expr>, has_continuous_layout<Expr> >,
					continuous_linear_evaluator<T>,
					cached_linear_evaluator<T> >::type type;
	};

	template<class Expr>
	struct percol_vector_evaluator
	{
		typedef typename matrix_traits<Expr>::value_type T;

		typedef typename
				if_<is_dense_mat<Expr>,
					dense_percol_evaluator<T>,
					cached_percol_evaluator<T> >::type type;
	};

	template<typename T, int CTRows, int CTCols>
	struct linear_vector_evaluator<const_matrix<T, CTRows, CTCols> >
	{
		typedef const_linear_evaluator<T> type;
	};

	template<typename T, int CTRows, int CTCols>
	struct percol_vector_evaluator<const_matrix<T, CTRows, CTCols> >
	{
		typedef const_percol_evaluator<T> type;
	};

}

#endif /* MATRIX_VEC_EVALUATORS_H_ */
