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

#include <light_mat/matrix/repeat_vecs.h>
#include <light_mat/matrix/generic_matrix_eval.h>

namespace lmat
{

	/********************************************
	 *
	 *  Evaluators
	 *
	 ********************************************/

	template<class Vec>
	class single_vec_percol_evaluator
	: public IPerColVectorEvaluator<
	  	  single_vec_percol_evaluator<Vec>,
	  	  typename matrix_traits<Vec>::value_type>
	{
	public:
		typedef typename matrix_traits<Vec>::value_type T;

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
		typename percol_vector_evaluator<Vec>::type m_eval;
	};


	template<class Vec>
	class single_vec_linear_evaluator
	: public ILinearVectorEvaluator<
	  	  single_vec_linear_evaluator<Vec>,
	  	  typename matrix_traits<Vec>::value_type>
	{
	public:
		typedef typename matrix_traits<Vec>::value_type T;

		LMAT_ENSURE_INLINE
		single_vec_linear_evaluator(const repeat_col_expr<Vec, 1>& expr)
		: m_eval(expr.column()) { }

		LMAT_ENSURE_INLINE
		single_vec_linear_evaluator(const repeat_row_expr<Vec, 1>& expr)
		: m_eval(expr.row()) { }

		LMAT_ENSURE_INLINE
		T get_value(const index_t i) const { return m_eval.get_value(i); }

	private:
		typename linear_vector_evaluator<Vec>::type m_eval;
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
		typedef typename matrix_traits<Col>::value_type T;

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
		typedef typename detail::repcol_ewrapper_map<Col>::type wrapper_t;
		wrapper_t m_colwrap;

		typename percol_vector_evaluator<typename wrapper_t::col_t>::type m_eval;
	};


	template<class Row>
	class reprow_percol_evaluator
	: public IPerColVectorEvaluator<
	  	  reprow_percol_evaluator<Row>,
	  	  typename matrix_traits<Row>::value_type>
	{
	public:
		typedef typename matrix_traits<Row>::value_type T;

		template<int N>
		LMAT_ENSURE_INLINE
		reprow_percol_evaluator(const repeat_row_expr<Row, N>& expr)
		: m_rowwrap(expr.row()), m_j(0) { }

		LMAT_ENSURE_INLINE
		T get_value(const index_t ) const { return m_rowwrap[m_j]; }

		LMAT_ENSURE_INLINE
		void next_column() { ++ m_j; }

	private:
		typedef typename detail::reprow_ewrapper_map<Row>::type wrapper_t;
		wrapper_t m_rowwrap;
		index_t m_j;
	};




	/********************************************
	 *
	 *  Dispatch
	 *
	 ********************************************/

	template<class Col, int N>
	struct linear_vector_evaluator<repeat_col_expr<Col, N> >
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
				>::type type;
	};


	template<class Row, int M>
	struct linear_vector_evaluator<repeat_row_expr<Row, M> >
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
				>::type type;
	};


	template<class Col, int N>
	struct percol_vector_evaluator<repeat_col_expr<Col, N> >
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
				>::type type;
	};


	template<class Row, int M>
	struct percol_vector_evaluator<repeat_row_expr<Row, M> >
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
				>::type type;
	};


	/********************************************
	 *
	 *  Cost model
	 *
	 ********************************************/

	template<class Col, int N>
	struct linear_eval_cost<repeat_col_expr<Col, N> >
	{
		static const int M = ct_rows<Col>::value;
		static const int _cost = (M == 1 || N == 1) ? 0 : VEC_EVAL_CACHE_COST;

		LMAT_ENSURE_INLINE
		static int of(const repeat_col_expr<Col, N>& )
		{
			return _cost;
		}
	};

	template<class Row, int M>
	struct linear_eval_cost<repeat_row_expr<Row, M> >
	{
		static const int N = ct_cols<Row>::value;
		static const int _cost = (M == 1 || N == 1) ? 0 : VEC_EVAL_CACHE_COST;

		LMAT_ENSURE_INLINE
		static int of(const repeat_row_expr<Row, M>& )
		{
			return _cost;
		}
	};


	template<class Col, int N>
	struct percol_eval_cost<repeat_col_expr<Col, N> >
	{
		LMAT_ENSURE_INLINE
		static int of(const repeat_col_expr<Col, N>& )
		{
			return 0;
		}
	};

	template<class Row, int M>
	struct percol_eval_cost<repeat_row_expr<Row, M> >
	{
		LMAT_ENSURE_INLINE
		static int of(const repeat_row_expr<Row, M>& )
		{
			return 0;
		}
	};


}

#endif /* REPEAT_VECS_EVAL_H_ */
