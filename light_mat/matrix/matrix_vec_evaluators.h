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
		continuous_linear_evaluator(const IDenseMatrix<Mat, T>& X)
		{
#ifdef LMAT_USE_STATIC_ASSERT
			static_assert(has_continuous_layout<Mat>::value,
					"Mat must always have continuous layout");
#endif
			m_data = X.ptr_data();
		}

		T get_value(const index_t i) const
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
		dense_percol_evaluator(const IDenseMatrix<Mat, T>& X)
		{
			m_ldim = X.lead_dim();
			m_data = X.ptr_data();
		}

		T get_value(const index_t i) const
		{
			return m_data[i];
		}

		void next_column()
		{
			m_data += m_ldim;
		}
	private:
		const index_t m_ldim;
		const T *m_data;
	};


	template<typename T>
	class cached_linear_evaluator
	: public ILinearVectorEvaluator<cached_linear_evaluator<T>, T>
	{
	public:
		template<class Expr>
		cached_linear_evaluator(const IMatrixExpr<Expr, T>& X)
		: m_cache(X), m_data(m_cache.ptr_data())
		{
		}

		T get_value(const index_t i) const
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
		cached_percol_evaluator(const IMatrixExpr<Expr, T>& X)
		: m_cache(X), m_ldim(m_cache.lead_dim()), m_data(m_cache.ptr_data())
		{
		}

		T get_value(const index_t i) const
		{
			return m_data[i];
		}

		void next_column()
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
		typedef typename
				if_<and_<is_dense_mat<Expr>, has_continuous_layout<Expr> >,
					continuous_linear_evaluator<Expr>,
					cached_linear_evaluator<Expr> >::type type;
	};

	template<class Expr>
	struct percol_vector_evaluator
	{
		typedef typename
				if_<is_dense_mat<Expr>,
					dense_percol_evaluator<Expr>,
					cached_percol_evaluator<Expr> >::type type;
	};


}

#endif /* MATRIX_VEC_EVALUATORS_H_ */