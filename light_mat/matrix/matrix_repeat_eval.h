/**
 * @file matrix_repeat_eval.h
 *
 * Evaluation of matrix repeating
 *
 * @author Dahua Lin
 */

#ifdef _MSC_VER
#pragma once
#endif

#ifndef LIGHTMAT_MATRIX_REPEAT_EVAL_H_
#define LIGHTMAT_MATRIX_REPEAT_EVAL_H_

#include <light_mat/matrix/matrix_repeat_expr.h>
#include <light_mat/matrix/matrix_vector_eval.h>
#include "bits/repeat_vecs_internal.h"

namespace lmat
{
	/********************************************
	 *
	 *  Expression evaluation
	 *
	 ********************************************/

	template<class Arg_Holder, int N, class Dst>
	struct default_matrix_eval_policy<horizontal_repeat_expr<Arg_Holder, N>, Dst>
	{
		typedef matrix_copy_policy type;
	};

	template<class Arg_Holder, int M, class Dst>
	struct default_matrix_eval_policy<vertical_repeat_expr<Arg_Holder, M>, Dst>
	{
		typedef matrix_copy_policy type;
	};


	template<class Arg_Holder, int N, class Dst>
	inline void evaluate(const horizontal_repeat_expr<Arg_Holder, N>& expr,
			IDenseMatrix<Dst, typename matrix_traits<typename Arg_Holder::arg_type>::value_type >& dst,
			matrix_copy_policy)
	{
		typedef typename Arg_Holder::arg_type Arg;
		Arg s = expr.arg();
		typedef typename matrix_traits<Arg>::value_type T;

		if ( is_column(expr) )
		{
			const int M = binary_ct_rows<Arg, Dst>::value;
			ref_col<T, M> dview(dst.ptr_data(), s.nrows());
			default_evaluate(s, dview);
		}
		else
		{
			typedef typename detail::repcol_ewrapper_map<Arg>::type wrapper_t;
			wrapper_t col_wrap(s);

			const index_t m = col_wrap.nrows();
			if (m == 1)
			{
				fill(dst, col_wrap[0]);
			}
			else
			{
				const index_t n = expr.ncolumns();
				for (index_t j = 0; j < n; ++j)
				{
					copy_mem(m, col_wrap.data(), dst.ptr_col(j));
				}
			}
		}
	}


	template<class Arg_Holder, class Dst>
	LMAT_ENSURE_INLINE
	inline void evaluate(const horizontal_repeat_expr<Arg_Holder, 1>& expr,
			IDenseMatrix<Dst, typename matrix_traits<typename Arg_Holder::arg_type>::value_type >& dst,
			matrix_copy_policy)
	{

		default_evaluate(expr.arg(), dst);
	}


	template<class Arg_Holder, int M, class Dst>
	inline void evaluate(const vertical_repeat_expr<Arg_Holder, M>& expr,
			IDenseMatrix<Dst, typename matrix_traits<typename Arg_Holder::arg_type>::value_type >& dst,
			matrix_copy_policy)
	{
		typedef typename Arg_Holder::arg_type Arg;
		Arg s = expr.arg();
		typedef typename matrix_traits<Arg>::value_type T;

		if ( is_row(expr) )
		{
			const int N = binary_ct_cols<Arg, Dst>::value;
			if (has_continuous_layout(dst))
			{
				ref_row<T, N> dview(dst.ptr_data(), s.ncolumns());
				default_evaluate(s, dview);
			}
			else
			{
				ref_matrix_ex<T, 1, N> dview(dst.ptr_data(), 1, s.ncolumns(), dst.lead_dim());
				default_evaluate(s, dview);
			}
		}
		else
		{
			typedef typename detail::reprow_ewrapper_map<Arg>::type wrapper_t;
			wrapper_t row_wrap(s);

			const index_t n = row_wrap.ncolumns();

			if (M == 0)
			{
				const index_t m = expr.nrows();
				if (n == 1)
				{
					fill_mem(m, dst.ptr_data(), row_wrap[0]);
				}
				else
				{
					for (index_t j = 0; j < n; ++j)
						fill_mem(m, dst.ptr_col(j), row_wrap[j]);
				}
			}
			else
			{
				if (n == 1)
				{
					fill_mem(M, dst.ptr_data(), row_wrap[0]);
				}
				else
				{
					for (index_t j = 0; j < n; ++j)
						fill_mem(M, dst.ptr_col(j), row_wrap[j]);
				}
			}
		}

	}


	template<class Arg_Holder, class Dst>
	LMAT_ENSURE_INLINE
	inline void evaluate(const vertical_repeat_expr<Arg_Holder, 1>& expr,
			IDenseMatrix<Dst, typename matrix_traits<typename Arg_Holder::arg_type>::value_type>& dst,
			matrix_copy_policy)
	{

		default_evaluate(expr.arg(), dst);
	}




	/********************************************
	 *
	 *  Vector-based Evaluators
	 *
	 ********************************************/

	template<class Arg_Holder>
	class single_vec_percol_evaluator
	: public IPerColVectorEvaluator<
	  	  single_vec_percol_evaluator<Arg_Holder>,
	  	  typename matrix_traits<typename Arg_Holder::arg_type>::value_type>
	{
	public:
		typedef typename Arg_Holder::arg_type arg_type;
		typedef typename matrix_traits<arg_type>::value_type T;

		LMAT_ENSURE_INLINE
		single_vec_percol_evaluator(const horizontal_repeat_expr<Arg_Holder, 1>& expr)
		: m_eval(expr.arg()) { }

		LMAT_ENSURE_INLINE
		single_vec_percol_evaluator(const vertical_repeat_expr<Arg_Holder, 1>& expr)
		: m_eval(expr.arg()) { }

		LMAT_ENSURE_INLINE
		T get_value(const index_t i) const { return m_eval.get_value(i); }

		LMAT_ENSURE_INLINE
		void next_column() { m_eval.next_column(); }

	private:
		typename vector_eval<arg_type, per_column, by_scalars>::evaluator_type m_eval;
	};


	template<class Arg_Holder>
	class single_vec_linear_evaluator
	: public ILinearVectorEvaluator<
	  	  single_vec_linear_evaluator<Arg_Holder>,
	  	  typename matrix_traits<typename Arg_Holder::arg_type>::value_type>
	{
	public:
		typedef typename Arg_Holder::arg_type arg_type;
		typedef typename matrix_traits<arg_type>::value_type T;

		LMAT_ENSURE_INLINE
		single_vec_linear_evaluator(const horizontal_repeat_expr<Arg_Holder, 1>& expr)
		: m_eval(expr.arg()) { }

		LMAT_ENSURE_INLINE
		single_vec_linear_evaluator(const vertical_repeat_expr<Arg_Holder, 1>& expr)
		: m_eval(expr.arg()) { }

		LMAT_ENSURE_INLINE
		T get_value(const index_t i) const { return m_eval.get_value(i); }

	private:
		typename vector_eval<arg_type, as_linear_vec, by_scalars>::evaluator_type m_eval;
	};


	template<class Arg_Holder>
	class rep_scalar_percol_evaluator
	: public IPerColVectorEvaluator<
	  	  rep_scalar_percol_evaluator<Arg_Holder>,
	  	  typename matrix_traits<typename Arg_Holder::arg_type>::value_type>
	{
	public:
		typedef typename Arg_Holder::arg_type arg_type;
		typedef typename matrix_traits<arg_type>::value_type T;

		template<int N>
		LMAT_ENSURE_INLINE
		rep_scalar_percol_evaluator(const horizontal_repeat_expr<Arg_Holder, N>& expr)
		{
			m_val = to_scalar(expr.arg());
		}

		template<int M>
		LMAT_ENSURE_INLINE
		rep_scalar_percol_evaluator(const vertical_repeat_expr<Arg_Holder, M>& expr)
		{
			m_val = to_scalar(expr.arg());
		}

		LMAT_ENSURE_INLINE
		T get_value(const index_t) const { return m_val; }

		LMAT_ENSURE_INLINE
		void next_column() { }

	private:
		T m_val;
	};


	template<class Arg_Holder>
	class rep_scalar_linear_evaluator
	: public ILinearVectorEvaluator<
	  	  rep_scalar_linear_evaluator<Arg_Holder>,
	  	  typename matrix_traits<typename Arg_Holder::arg_type>::value_type>
	{
	public:
		typedef typename Arg_Holder::arg_type arg_type;
		typedef typename matrix_traits<arg_type>::value_type T;

		template<int N>
		LMAT_ENSURE_INLINE
		rep_scalar_linear_evaluator(const horizontal_repeat_expr<Arg_Holder, N>& expr)
		{
			m_val = to_scalar(expr.arg());
		}

		template<int M>
		LMAT_ENSURE_INLINE
		rep_scalar_linear_evaluator(const vertical_repeat_expr<Arg_Holder, M>& expr)
		{
			m_val = to_scalar(expr.arg());
		}

		LMAT_ENSURE_INLINE
		T get_value(const index_t) const { return m_val; }

	private:
		T m_val;
	};


	template<class Arg_Holder>
	class repcol_percol_evaluator
	: public IPerColVectorEvaluator<
	  	  repcol_percol_evaluator<Arg_Holder>,
	  	  typename matrix_traits<typename Arg_Holder::arg_type>::value_type>
	{
	public:
		typedef typename Arg_Holder::arg_type arg_type;
		typedef typename matrix_traits<arg_type>::value_type T;

		template<int N>
		LMAT_ENSURE_INLINE
		repcol_percol_evaluator(const horizontal_repeat_expr<Arg_Holder, N>& expr)
		: m_colwrap(expr.arg())
		, m_eval(m_colwrap.ref()) { }

		LMAT_ENSURE_INLINE
		T get_value(const index_t i) const { return m_eval.get_value(i); }

		LMAT_ENSURE_INLINE
		void next_column() { }

	private:
		typename detail::repcol_ewrapper_map<arg_type>::type m_colwrap;
		typename vector_eval<arg_type, per_column, by_scalars>::evaluator_type m_eval;
	};


	template<class Arg_Holder>
	class reprow_percol_evaluator
	: public IPerColVectorEvaluator<
	  	  reprow_percol_evaluator<Arg_Holder>,
	  	  typename matrix_traits<typename Arg_Holder::arg_type>::value_type>
	{
	public:
		typedef typename Arg_Holder::arg_type arg_type;
		typedef typename matrix_traits<arg_type>::value_type T;

		template<int N>
		LMAT_ENSURE_INLINE
		reprow_percol_evaluator(const vertical_repeat_expr<Arg_Holder, N>& expr)
		: m_rowwrap(expr.arg()), m_j(0) { }

		LMAT_ENSURE_INLINE
		T get_value(const index_t ) const { return m_rowwrap[m_j]; }

		LMAT_ENSURE_INLINE
		void next_column() { ++ m_j; }

	private:
		typename detail::reprow_ewrapper_map<arg_type>::type m_rowwrap;
		index_t m_j;
	};




	/********************************************
	 *
	 *  Dispatch
	 *
	 ********************************************/

	template<class Arg_Holder, int N>
	struct vector_eval<horizontal_repeat_expr<Arg_Holder, N>, as_linear_vec, by_scalars>
	{
		typedef typename Arg_Holder::arg_type arg_type;
		typedef typename matrix_traits<arg_type>::value_type T;
		static const int M = ct_rows<arg_type>::value;

		typedef typename
				if_c<N == 1,
					single_vec_linear_evaluator<Arg_Holder>,
					typename
					if_c<M == 1,
						rep_scalar_linear_evaluator<Arg_Holder>,
						cached_linear_evaluator<T>
					>::type
				>::type evaluator_type;

		static const int cost = (N == 1 || M == 1) ? 0 : VEC_EVAL_CACHE_COST;
	};


	template<class Arg_Holder, int M>
	struct vector_eval<vertical_repeat_expr<Arg_Holder, M>, as_linear_vec, by_scalars>
	{
		typedef typename Arg_Holder::arg_type arg_type;
		typedef typename matrix_traits<arg_type>::value_type T;
		static const int N = ct_cols<arg_type>::value;

		typedef typename
				if_c<M == 1,
					single_vec_linear_evaluator<Arg_Holder>,
					typename
					if_c<N == 1,
						rep_scalar_linear_evaluator<Arg_Holder>,
						cached_linear_evaluator<typename matrix_traits<arg_type>::value_type>
					>::type
				>::type evaluator_type;

		static const int cost = (M == 1 || N == 1) ? 0 : VEC_EVAL_CACHE_COST;
	};


	template<class Arg_Holder, int N>
	struct vector_eval<horizontal_repeat_expr<Arg_Holder, N>, per_column, by_scalars>
	{
		typedef typename Arg_Holder::arg_type arg_type;
		typedef typename matrix_traits<arg_type>::value_type T;
		static const int M = ct_rows<arg_type>::value;

		typedef typename
				if_c<N == 1,
					single_vec_percol_evaluator<Arg_Holder>,
					typename
					if_c<M == 1,
						rep_scalar_percol_evaluator<Arg_Holder>,
						repcol_percol_evaluator<Arg_Holder>
					>::type
				>::type evaluator_type;

		static const int normal_cost = 0;
		static const int shortv_cost = SHORTVEC_PERCOL_COST;
		static const int cost = M < SHORTVEC_LENGTH_THRESHOLD ? shortv_cost : normal_cost;
	};


	template<class Arg_Holder, int M>
	struct vector_eval<vertical_repeat_expr<Arg_Holder, M>, per_column, by_scalars>
	{
		typedef typename Arg_Holder::arg_type arg_type;
		typedef typename matrix_traits<arg_type>::value_type T;
		static const int N = ct_cols<arg_type>::value;

		typedef typename
				if_c<M == 1,
					single_vec_percol_evaluator<Arg_Holder>,
					typename
					if_c<N == 1,
						rep_scalar_percol_evaluator<Arg_Holder>,
						reprow_percol_evaluator<Arg_Holder>
					>::type
				>::type evaluator_type;

		static const int normal_cost = 0;
		static const int shortv_cost = SHORTVEC_PERCOL_COST;
		static const int cost = M < SHORTVEC_LENGTH_THRESHOLD ? shortv_cost : normal_cost;
	};

}

#endif /* REPEAT_VECS_EVAL_H_ */
