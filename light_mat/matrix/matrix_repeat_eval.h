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

	template<typename Arg_HP, class Arg, int N, class Dst>
	struct default_matrix_eval_policy<horizontal_repeat_expr<Arg_HP, Arg, N>, Dst>
	{
		typedef matrix_copy_policy type;
	};

	template<typename Arg_HP, class Arg, int M, class Dst>
	struct default_matrix_eval_policy<vertical_repeat_expr<Arg_HP, Arg, M>, Dst>
	{
		typedef matrix_copy_policy type;
	};


	template<typename Arg_HP, class Arg, int N, class Dst>
	inline void evaluate(const horizontal_repeat_expr<Arg_HP, Arg, N>& expr,
			IDenseMatrix<Dst, typename matrix_traits<Arg>::value_type >& dst,
			matrix_copy_policy)
	{
		typedef horizontal_repeat_expr<Arg_HP, Arg, N> expr_t;
		typedef typename expr_t::arg_type arg_type;

		const arg_type& s = expr.arg();
		typedef typename matrix_traits<Arg>::value_type T;

		if ( is_column(expr) )
		{
			const int M = binary_ct_rows<arg_type, Dst>::value;
			ref_col<T, M> dview(dst.ptr_data(), s.nrows());
			default_evaluate(s, dview);
		}
		else
		{
			typedef typename detail::repcol_ewrapper_map<arg_type>::type wrapper_t;
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


	template<typename Arg_HP, class Arg, class Dst>
	LMAT_ENSURE_INLINE
	inline void evaluate(const horizontal_repeat_expr<Arg_HP, Arg, 1>& expr,
			IDenseMatrix<Dst, typename matrix_traits<Arg>::value_type >& dst,
			matrix_copy_policy)
	{

		default_evaluate(expr.arg(), dst);
	}


	template<typename Arg_HP, class Arg, int M, class Dst>
	inline void evaluate(const vertical_repeat_expr<Arg_HP, Arg, M>& expr,
			IDenseMatrix<Dst, typename matrix_traits<Arg>::value_type >& dst,
			matrix_copy_policy)
	{
		typedef vertical_repeat_expr<Arg_HP, Arg, M> expr_t;
		typedef typename expr_t::arg_type arg_type;

		const arg_type& s = expr.arg();
		typedef typename matrix_traits<Arg>::value_type T;

		if ( is_row(expr) )
		{
			const int N = binary_ct_cols<arg_type, Dst>::value;
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
			typedef typename detail::reprow_ewrapper_map<arg_type>::type wrapper_t;
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


	template<typename Arg_HP, class Arg, class Dst>
	LMAT_ENSURE_INLINE
	inline void evaluate(const vertical_repeat_expr<Arg_HP, Arg, 1>& expr,
			IDenseMatrix<Dst, typename matrix_traits<Arg>::value_type>& dst,
			matrix_copy_policy)
	{

		default_evaluate(expr.arg(), dst);
	}



	/********************************************
	 *
	 *  Vector-based Evaluators
	 *
	 ********************************************/

	// forward declarations

	template<typename Arg_HP, class Arg> class single_vec_percol_evaluator;
	template<typename Arg_HP, class Arg> class single_vec_linear_evaluator;
	template<typename Arg_HP, class Arg> class rep_scalar_percol_evaluator;
	template<typename Arg_HP, class Arg> class rep_scalar_linear_evaluator;
	template<typename Arg_HP, class Arg> class repcol_percol_evaluator;
	template<typename Arg_HP, class Arg> class reprow_percol_evaluator;

	// definitions

	template<typename Arg_HP, class Arg>
	struct percol_eval_state<single_vec_percol_evaluator<Arg_HP, Arg> >
	{
		typedef typename arg_holder<Arg_HP, Arg>::internal_arg_type arg_type;
		typedef typename vector_eval<arg_type, per_column, by_scalars>::evaluator_type arg_eval_t;
		typedef typename percol_eval_state<arg_eval_t>::type type;
	};

	template<typename Arg_HP, class Arg>
	class single_vec_percol_evaluator
	: public IPerColVectorEvaluator<
	  	  single_vec_percol_evaluator<Arg_HP, Arg>,
	  	  typename matrix_traits<Arg>::value_type>
	{
	public:
		typedef typename arg_holder<Arg_HP, Arg>::internal_arg_type arg_type;
		typedef typename matrix_traits<Arg>::value_type T;

		typedef typename vector_eval<arg_type, per_column, by_scalars>::evaluator_type arg_eval_t;
		typedef typename percol_eval_state<arg_eval_t>::type state_t;

		LMAT_ENSURE_INLINE
		single_vec_percol_evaluator(const horizontal_repeat_expr<Arg_HP, Arg, 1>& expr)
		: m_eval(expr.arg()) { }

		LMAT_ENSURE_INLINE
		single_vec_percol_evaluator(const vertical_repeat_expr<Arg_HP, Arg, 1>& expr)
		: m_eval(expr.arg()) { }

		LMAT_ENSURE_INLINE
		T get_value(const state_t& s, const index_t i) const
		{
			return m_eval.get_value(s, i);
		}

		LMAT_ENSURE_INLINE
		state_t col_state(const index_t j) const { return m_eval.col_state(j); }

	private:
		arg_eval_t m_eval;
	};


	template<typename Arg_HP, class Arg>
	class single_vec_linear_evaluator
	: public ILinearVectorEvaluator<
	  	  single_vec_linear_evaluator<Arg_HP, Arg>,
	  	  typename matrix_traits<Arg>::value_type>
	{
	public:
		typedef typename arg_holder<Arg_HP, Arg>::internal_arg_type arg_type;
		typedef typename matrix_traits<Arg>::value_type T;
		typedef typename vector_eval<arg_type, as_linear_vec, by_scalars>::evaluator_type arg_eval_t;

		LMAT_ENSURE_INLINE
		single_vec_linear_evaluator(const horizontal_repeat_expr<Arg_HP, Arg, 1>& expr)
		: m_eval(expr.arg()) { }

		LMAT_ENSURE_INLINE
		single_vec_linear_evaluator(const vertical_repeat_expr<Arg_HP, Arg, 1>& expr)
		: m_eval(expr.arg()) { }

		LMAT_ENSURE_INLINE
		T get_value(const index_t i) const { return m_eval.get_value(i); }

	private:
		arg_eval_t m_eval;
	};


	template<typename Arg_HP, class Arg>
	struct percol_eval_state<rep_scalar_percol_evaluator<Arg_HP, Arg> >
	{
		typedef nil_eval_state type;
	};

	template<typename Arg_HP, class Arg>
	class rep_scalar_percol_evaluator
	: public IPerColVectorEvaluator<
	  	  rep_scalar_percol_evaluator<Arg_HP, Arg>,
	  	  typename matrix_traits<Arg>::value_type>
	{
	public:
		typedef typename arg_holder<Arg_HP, Arg>::internal_arg_type arg_type;
		typedef typename matrix_traits<Arg>::value_type T;

		template<int N>
		LMAT_ENSURE_INLINE
		rep_scalar_percol_evaluator(const horizontal_repeat_expr<Arg_HP, Arg, N>& expr)
		{
			m_val = to_scalar(expr.arg());
		}

		template<int M>
		LMAT_ENSURE_INLINE
		rep_scalar_percol_evaluator(const vertical_repeat_expr<Arg_HP, Arg, M>& expr)
		{
			m_val = to_scalar(expr.arg());
		}

		LMAT_ENSURE_INLINE
		T get_value(const nil_eval_state&, const index_t) const { return m_val; }

		LMAT_ENSURE_INLINE
		nil_eval_state col_state(const index_t j) const { return nil_eval_state(); }

	private:
		T m_val;
	};


	template<typename Arg_HP, class Arg>
	class rep_scalar_linear_evaluator
	: public ILinearVectorEvaluator<
	  	  rep_scalar_linear_evaluator<Arg_HP, Arg>,
	  	  typename matrix_traits<Arg>::value_type>
	{
	public:
		typedef typename arg_holder<Arg_HP, Arg>::internal_arg_type arg_type;
		typedef typename matrix_traits<Arg>::value_type T;

		template<int N>
		LMAT_ENSURE_INLINE
		rep_scalar_linear_evaluator(const horizontal_repeat_expr<Arg_HP, Arg, N>& expr)
		{
			m_val = to_scalar(expr.arg());
		}

		template<int M>
		LMAT_ENSURE_INLINE
		rep_scalar_linear_evaluator(const vertical_repeat_expr<Arg_HP, Arg, M>& expr)
		{
			m_val = to_scalar(expr.arg());
		}

		LMAT_ENSURE_INLINE
		T get_value(const index_t) const { return m_val; }

	private:
		T m_val;
	};


	template<typename Arg_HP, class Arg>
	struct percol_eval_state<repcol_percol_evaluator<Arg_HP, Arg> >
	{
		typedef nil_eval_state type;
	};

	template<typename Arg_HP, class Arg>
	class repcol_percol_evaluator
	: public IPerColVectorEvaluator<
	  	  repcol_percol_evaluator<Arg_HP, Arg>,
	  	  typename matrix_traits<Arg>::value_type>
	{
	public:
		typedef typename arg_holder<Arg_HP, Arg>::internal_arg_type arg_type;
		typedef typename matrix_traits<Arg>::value_type T;

		template<int N>
		LMAT_ENSURE_INLINE
		repcol_percol_evaluator(const horizontal_repeat_expr<Arg_HP, Arg, N>& expr)
		: m_colwrap(expr.arg()) { }

		LMAT_ENSURE_INLINE
		T get_value(const nil_eval_state&, const index_t i) const { return m_colwrap[i]; }

		LMAT_ENSURE_INLINE
		nil_eval_state col_state(const index_t j) const { return nil_eval_state(); }

	private:
		typename detail::repcol_ewrapper_map<arg_type>::type m_colwrap;
	};


	template<typename Arg_HP, class Arg>
	struct percol_eval_state<reprow_percol_evaluator<Arg_HP, Arg> >
	{
		typedef typename matrix_traits<Arg>::value_type T;
		typedef T type;
	};

	template<typename Arg_HP, class Arg>
	class reprow_percol_evaluator
	: public IPerColVectorEvaluator<
	  	  reprow_percol_evaluator<Arg_HP, Arg>,
	  	  typename matrix_traits<Arg>::value_type>
	{
	public:
		typedef typename arg_holder<Arg_HP, Arg>::internal_arg_type arg_type;
		typedef typename matrix_traits<Arg>::value_type T;

		template<int N>
		LMAT_ENSURE_INLINE
		reprow_percol_evaluator(const vertical_repeat_expr<Arg_HP, Arg, N>& expr)
		: m_rowwrap(expr.arg()) { }

		LMAT_ENSURE_INLINE
		T get_value(const T& s, const index_t ) const { return s; }

		LMAT_ENSURE_INLINE
		T col_state(const index_t j) const { return m_rowwrap[j]; }

	private:
		typename detail::reprow_ewrapper_map<arg_type>::type m_rowwrap;
	};




	/********************************************
	 *
	 *  Dispatch
	 *
	 ********************************************/

	template<typename Arg_HP, class Arg, int N>
	struct vector_eval<horizontal_repeat_expr<Arg_HP, Arg, N>, as_linear_vec, by_scalars>
	{
		typedef typename matrix_traits<Arg>::value_type T;
		static const int M = ct_rows<Arg>::value;

		typedef typename
				if_c<N == 1,
					single_vec_linear_evaluator<Arg_HP, Arg>,
					typename
					if_c<M == 1,
						rep_scalar_linear_evaluator<Arg_HP, Arg>,
						cached_linear_evaluator<T>
					>::type
				>::type evaluator_type;

		static const int cost = (N == 1 || M == 1) ? 0 : VEC_EVAL_CACHE_COST;
	};


	template<typename Arg_HP, class Arg, int M>
	struct vector_eval<vertical_repeat_expr<Arg_HP, Arg, M>, as_linear_vec, by_scalars>
	{
		typedef typename matrix_traits<Arg>::value_type T;
		static const int N = ct_cols<Arg>::value;

		typedef typename
				if_c<M == 1,
					single_vec_linear_evaluator<Arg_HP, Arg>,
					typename
					if_c<N == 1,
						rep_scalar_linear_evaluator<Arg_HP, Arg>,
						cached_linear_evaluator<T>
					>::type
				>::type evaluator_type;

		static const int cost = (M == 1 || N == 1) ? 0 : VEC_EVAL_CACHE_COST;
	};


	template<typename Arg_HP, class Arg, int N>
	struct vector_eval<horizontal_repeat_expr<Arg_HP, Arg, N>, per_column, by_scalars>
	{
		typedef typename matrix_traits<Arg>::value_type T;
		static const int M = ct_rows<Arg>::value;

		typedef typename
				if_c<N == 1,
					single_vec_percol_evaluator<Arg_HP, Arg>,
					typename
					if_c<M == 1,
						rep_scalar_percol_evaluator<Arg_HP, Arg>,
						repcol_percol_evaluator<Arg_HP, Arg>
					>::type
				>::type evaluator_type;

		static const int normal_cost = 0;
		static const int shortv_cost = SHORTVEC_PERCOL_COST;
		static const int cost = M < SHORTVEC_LENGTH_THRESHOLD ? shortv_cost : normal_cost;
	};


	template<typename Arg_HP, class Arg, int M>
	struct vector_eval<vertical_repeat_expr<Arg_HP, Arg, M>, per_column, by_scalars>
	{
		typedef typename matrix_traits<Arg>::value_type T;
		static const int N = ct_cols<Arg>::value;

		typedef typename
				if_c<M == 1,
					single_vec_percol_evaluator<Arg_HP, Arg>,
					typename
					if_c<N == 1,
						rep_scalar_percol_evaluator<Arg_HP, Arg>,
						reprow_percol_evaluator<Arg_HP, Arg>
					>::type
				>::type evaluator_type;

		static const int normal_cost = 0;
		static const int shortv_cost = SHORTVEC_PERCOL_COST;
		static const int cost = M < SHORTVEC_LENGTH_THRESHOLD ? shortv_cost : normal_cost;
	};

}

#endif /* REPEAT_VECS_EVAL_H_ */
