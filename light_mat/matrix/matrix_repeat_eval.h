/**
 * @file repeat_vecs_eval.h
 *
 * Vector evaluation of repeat-vector expressions
 *
 * @author Dahua Lin
 */

#ifdef _MSC_VER
#pragma once
#endif

#ifndef LIGHTMAT_REPEAT_VECS_VEVAL_H_
#define LIGHTMAT_REPEAT_VECS_VEVAL_H_

#include <light_mat/matrix/repeat_vecs_expr.h>
#include <light_mat/matrix/matrix_veval.h>
#include "bits/repeat_vecs_internal.h"

namespace lmat
{

	/********************************************
	 *
	 *  Expression evaluation
	 *
	 ********************************************/

	template<class Col, int N, class Dst>
	struct repcols_evalctx
	{
#ifdef LMAT_USE_STATIC_ASSERT
		static_assert(N == 0 || N > 1, "N should be either 0 or greater than 1 here");
#endif

		typedef typename binary_value_type<Col, Dst>::type T;
		typedef typename unwrapped_expr<Col>::type Col_;

		inline
		static void evaluate(const repeat_col_expr<Col, N>& s,
				Dst& dst)
		{
			if ( is_column(s) )
			{
				const int M = binary_ct_rows<Col_, Dst>::value;
				ref_col<T, M> dview(dst.ptr_data(), s.nrows());
				default_evaluate(s.column(), dview);
			}
			else
			{
				typedef typename detail::repcol_ewrapper_map<Col_>::type wrapper_t;
				wrapper_t col_wrap(s.column());

				const index_t m = col_wrap.nrows();
				if (m == 1)
				{
					fill(dst, col_wrap[0]);
				}
				else
				{
					const index_t n = s.ncolumns();
					for (index_t j = 0; j < n; ++j)
					{
						copy_mem(m, col_wrap.data(), dst.ptr_col(j));
					}
				}
			}
		}
	};


	template<class Col, class Dst>
	struct repcols_evalctx<Col, 1, Dst>
	{
		LMAT_ENSURE_INLINE
		static void evaluate(const repeat_col_expr<Col, 1>& s, Dst& dst)
		{
			default_evaluate(s.column(), dst);
		}
	};


	template<class Row, int M, class Dst>
	struct reprows_evalctx
	{
#ifdef LMAT_USE_STATIC_ASSERT
		static_assert(M == 0 || M > 1, "M should be either 0 or greater than 1 here");
#endif

		typedef typename binary_value_type<Row, Dst>::type T;
		typedef typename unwrapped_expr<Row>::type Row_;

		inline
		static void evaluate(const repeat_row_expr<Row, M>& s,
				Dst& dst)
		{
			if ( is_row(s) )
			{
				const int N = binary_ct_cols<Row_, Dst>::value;
				if (has_continuous_layout(dst))
				{
					ref_row<T, N> dview(dst.ptr_data(), s.ncolumns());
					default_evaluate(s.row(), dview);
				}
				else
				{
					ref_matrix_ex<T, 1, N> dview(dst.ptr_data(), 1, s.ncolumns(), dst.lead_dim());
					default_evaluate(s.row(), dview);
				}
			}
			else
			{
				typedef typename detail::reprow_ewrapper_map<Row_>::type wrapper_t;
				wrapper_t row_wrap(s.row());

				const index_t n = row_wrap.ncolumns();

				if (M == 0)
				{
					const index_t m = s.nrows();
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
	};


	template<class Row, class Dst>
	struct reprows_evalctx<Row, 1, Dst>
	{
		typedef typename binary_value_type<Row, Dst>::type T;

		inline
		static void evaluate(const repeat_row_expr<Row, 1>& s, Dst& dst)
		{
			default_evaluate(s.row(), dst);
		}
	};



	/********************************************
	 *
	 *  Vector-based Evaluators
	 *
	 ********************************************/

	template<class Vec>
	class single_vec_percol_evaluator
	: public IPerColVectorEvaluator<
	  	  single_vec_percol_evaluator<Vec>,
	  	  typename matrix_traits<Vec>::value_type>
	{
	public:
		typedef typename unwrapped_expr<Vec>::type Vec_;
		typedef typename matrix_traits<Vec_>::value_type T;

		LMAT_ENSURE_INLINE
		single_vec_percol_evaluator(const repeat_col_expr<Vec, 1>& expr)
		: m_eval(expr.column()) { }

		LMAT_ENSURE_INLINE
		single_vec_percol_evaluator(const repeat_row_expr<Vec, 1>& expr)
		: m_eval(expr.row()) { }

		LMAT_ENSURE_INLINE
		T get_value(const index_t i) const { return m_eval.get_value(i); }

		LMAT_ENSURE_INLINE
		void next_column() { m_eval.next_column(); }

	private:
		typename percol_eval<Vec_>::evaluator_type m_eval;
	};


	template<class Vec>
	class single_vec_linear_evaluator
	: public ILinearVectorEvaluator<
	  	  single_vec_linear_evaluator<Vec>,
	  	  typename matrix_traits<Vec>::value_type>
	{
	public:
		typedef typename unwrapped_expr<Vec>::type Vec_;
		typedef typename matrix_traits<Vec_>::value_type T;

		LMAT_ENSURE_INLINE
		single_vec_linear_evaluator(const repeat_col_expr<Vec, 1>& expr)
		: m_eval(expr.column()) { }

		LMAT_ENSURE_INLINE
		single_vec_linear_evaluator(const repeat_row_expr<Vec, 1>& expr)
		: m_eval(expr.row()) { }

		LMAT_ENSURE_INLINE
		T get_value(const index_t i) const { return m_eval.get_value(i); }

	private:
		typename linear_eval<Vec_>::evaluator_type m_eval;
	};


	template<class Vec>
	class rep_scalar_percol_evaluator
	: public IPerColVectorEvaluator<
	  	  rep_scalar_percol_evaluator<Vec>,
	  	  typename matrix_traits<Vec>::value_type>
	{
	public:
		typedef typename matrix_traits<Vec>::value_type T;

		template<int N>
		LMAT_ENSURE_INLINE
		rep_scalar_percol_evaluator(const repeat_col_expr<Vec, N>& expr)
		{
			m_val = to_scalar(expr.column());
		}

		template<int M>
		LMAT_ENSURE_INLINE
		rep_scalar_percol_evaluator(const repeat_row_expr<Vec, M>& expr)
		{
			m_val = to_scalar(expr.row());
		}

		LMAT_ENSURE_INLINE
		T get_value(const index_t) const { return m_val; }

		LMAT_ENSURE_INLINE
		void next_column() { }

	private:
		typename matrix_traits<Vec>::value_type m_val;
	};


	template<class Vec>
	class rep_scalar_linear_evaluator
	: public ILinearVectorEvaluator<
	  	  rep_scalar_linear_evaluator<Vec>,
	  	  typename matrix_traits<Vec>::value_type>
	{
	public:
		typedef typename matrix_traits<Vec>::value_type T;

		template<int N>
		LMAT_ENSURE_INLINE
		rep_scalar_linear_evaluator(const repeat_col_expr<Vec, N>& expr)
		{
			m_val = to_scalar(expr.column());
		}

		template<int M>
		LMAT_ENSURE_INLINE
		rep_scalar_linear_evaluator(const repeat_row_expr<Vec, M>& expr)
		{
			m_val = to_scalar(expr.row());
		}

		LMAT_ENSURE_INLINE
		T get_value(const index_t) const { return m_val; }

	private:
		typename matrix_traits<Vec>::value_type m_val;
	};


	template<class Col>
	class repcol_percol_evaluator
	: public IPerColVectorEvaluator<
	  	  repcol_percol_evaluator<Col>,
	  	  typename matrix_traits<Col>::value_type>
	{
	public:
		typedef typename unwrapped_expr<Col>::type Col_;
		typedef typename matrix_traits<Col_>::value_type T;

		template<int N>
		LMAT_ENSURE_INLINE
		repcol_percol_evaluator(const repeat_col_expr<Col, N>& expr)
		: m_colwrap(expr.column())
		, m_eval(m_colwrap.ref()) { }

		LMAT_ENSURE_INLINE
		T get_value(const index_t i) const { return m_eval.get_value(i); }

		LMAT_ENSURE_INLINE
		void next_column() { }

	private:
		typedef typename detail::repcol_ewrapper_map<Col_>::type wrapper_t;
		wrapper_t m_colwrap;

		typename percol_eval<typename wrapper_t::col_t>::evaluator_type m_eval;
	};


	template<class Row>
	class reprow_percol_evaluator
	: public IPerColVectorEvaluator<
	  	  reprow_percol_evaluator<Row>,
	  	  typename matrix_traits<Row>::value_type>
	{
	public:
		typedef typename unwrapped_expr<Row>::type Row_;
		typedef typename matrix_traits<Row_>::value_type T;

		template<int N>
		LMAT_ENSURE_INLINE
		reprow_percol_evaluator(const repeat_row_expr<Row, N>& expr)
		: m_rowwrap(expr.row()), m_j(0) { }

		LMAT_ENSURE_INLINE
		T get_value(const index_t ) const { return m_rowwrap[m_j]; }

		LMAT_ENSURE_INLINE
		void next_column() { ++ m_j; }

	private:
		typedef typename detail::reprow_ewrapper_map<Row_>::type wrapper_t;
		wrapper_t m_rowwrap;
		index_t m_j;
	};




	/********************************************
	 *
	 *  Dispatch
	 *
	 ********************************************/

	template<class Col, int N>
	struct linear_eval<repeat_col_expr<Col, N> >
	{
		static const int M = ct_rows<Col>::value;

		typedef typename
				if_c<N == 1,
					single_vec_linear_evaluator<Col>,
					typename
					if_c<M == 1,
						rep_scalar_linear_evaluator<Col>,
						cached_linear_evaluator<typename matrix_traits<Col>::value_type>
					>::type
				>::type evaluator_type;

		static const int cost = (N == 1 || M == 1) ? 0 : VEC_EVAL_CACHE_COST;
	};


	template<class Row, int M>
	struct linear_eval<repeat_row_expr<Row, M> >
	{
		static const int N = ct_cols<Row>::value;

		typedef typename
				if_c<M == 1,
					single_vec_linear_evaluator<Row>,
					typename
					if_c<N == 1,
						rep_scalar_linear_evaluator<Row>,
						cached_linear_evaluator<typename matrix_traits<Row>::value_type>
					>::type
				>::type evaluator_type;

		static const int cost = (M == 1 || N == 1) ? 0 : VEC_EVAL_CACHE_COST;
	};


	template<class Col, int N>
	struct percol_eval<repeat_col_expr<Col, N> >
	{
		static const int M = ct_rows<Col>::value;

		typedef typename
				if_c<N == 1,
					single_vec_percol_evaluator<Col>,
					typename
					if_c<M == 1,
						rep_scalar_percol_evaluator<Col>,
						repcol_percol_evaluator<Col>
					>::type
				>::type evaluator_type;

		static const int normal_cost = 0;
		static const int shortv_cost = SHORTVEC_PERCOL_COST;
	};


	template<class Row, int M>
	struct percol_eval<repeat_row_expr<Row, M> >
	{
		static const int N = ct_cols<Row>::value;

		typedef typename
				if_c<M == 1,
					single_vec_percol_evaluator<Row>,
					typename
					if_c<N == 1,
						rep_scalar_percol_evaluator<Row>,
						reprow_percol_evaluator<Row>
					>::type
				>::type evaluator_type;

		static const int normal_cost = 0;
		static const int shortv_cost = SHORTVEC_PERCOL_COST;
	};

}

#endif /* REPEAT_VECS_EVAL_H_ */
