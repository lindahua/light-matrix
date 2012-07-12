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

	template<class Vec, bool IsEmbed>
	class single_vec_percol_evaluator
	: public IPerColVectorEvaluator<
	  	  single_vec_percol_evaluator<Vec, IsEmbed>,
	  	  typename matrix_traits<Vec>::value_type>
	{
	public:
		typedef typename matrix_traits<Vec>::value_type T;

		LMAT_ENSURE_INLINE
		single_vec_percol_evaluator(const repeat_col_expr<Vec, 1, IsEmbed>& expr)
		: m_eval(expr.column()) { }

		LMAT_ENSURE_INLINE
		single_vec_percol_evaluator(const repeat_row_expr<Vec, 1, IsEmbed>& expr)
		: m_eval(expr.row()) { }

		LMAT_ENSURE_INLINE
		T get_value(const index_t i) const { return m_eval.get_value(i); }

		LMAT_ENSURE_INLINE
		void next_column() { m_eval.next_column(); }

	private:
		typename percol_vector_evaluator<Vec>::type m_eval;
	};


	template<class Vec, bool IsEmbed>
	class single_vec_linear_evaluator
	: public ILinearVectorEvaluator<
	  	  single_vec_linear_evaluator<Vec, IsEmbed>,
	  	  typename matrix_traits<Vec>::value_type>
	{
	public:
		typedef typename matrix_traits<Vec>::value_type T;

		LMAT_ENSURE_INLINE
		single_vec_linear_evaluator(const repeat_col_expr<Vec, 1, IsEmbed>& expr)
		: m_eval(expr.column()) { }

		LMAT_ENSURE_INLINE
		single_vec_linear_evaluator(const repeat_row_expr<Vec, 1, IsEmbed>& expr)
		: m_eval(expr.row()) { }

		LMAT_ENSURE_INLINE
		T get_value(const index_t i) const { return m_eval.get_value(i); }

	private:
		typename linear_vector_evaluator<Vec>::type m_eval;
	};


	template<class Vec, bool IsEmbed>
	class rep_scalar_percol_evaluator
	: public IPerColVectorEvaluator<
	  	  rep_scalar_percol_evaluator<Vec, IsEmbed>,
	  	  typename matrix_traits<Vec>::value_type>
	{
	public:
		typedef typename matrix_traits<Vec>::value_type T;

		template<int N>
		LMAT_ENSURE_INLINE
		rep_scalar_percol_evaluator(const repeat_col_expr<Vec, N, IsEmbed>& expr)
		{
			m_val = to_scalar(expr.column());
		}

		template<int M>
		LMAT_ENSURE_INLINE
		rep_scalar_percol_evaluator(const repeat_row_expr<Vec, M, IsEmbed>& expr)
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


	template<class Vec, bool IsEmbed>
	class rep_scalar_linear_evaluator
	: public ILinearVectorEvaluator<
	  	  rep_scalar_linear_evaluator<Vec, IsEmbed>,
	  	  typename matrix_traits<Vec>::value_type>
	{
	public:
		typedef typename matrix_traits<Vec>::value_type T;

		template<int N>
		LMAT_ENSURE_INLINE
		rep_scalar_linear_evaluator(const repeat_col_expr<Vec, N, IsEmbed>& expr)
		{
			m_val = to_scalar(expr.column());
		}

		template<int M>
		LMAT_ENSURE_INLINE
		rep_scalar_linear_evaluator(const repeat_row_expr<Vec, M, IsEmbed>& expr)
		{
			m_val = to_scalar(expr.row());
		}

		LMAT_ENSURE_INLINE
		T get_value(const index_t) const { return m_val; }

	private:
		typename matrix_traits<Vec>::value_type m_val;
	};


	template<class Col, bool IsEmbed>
	class repcol_percol_evaluator
	: public IPerColVectorEvaluator<
	  	  repcol_percol_evaluator<Col, IsEmbed>,
	  	  typename matrix_traits<Col>::value_type>
	{
	public:
		typedef typename matrix_traits<Col>::value_type T;

		template<int N>
		LMAT_ENSURE_INLINE
		repcol_percol_evaluator(const repeat_col_expr<Col, N, IsEmbed>& expr)
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


	template<class Row, bool IsEmbed>
	class reprow_percol_evaluator
	: public IPerColVectorEvaluator<
	  	  reprow_percol_evaluator<Row, IsEmbed>,
	  	  typename matrix_traits<Row>::value_type>
	{
	public:
		typedef typename matrix_traits<Row>::value_type T;

		template<int N>
		LMAT_ENSURE_INLINE
		reprow_percol_evaluator(const repeat_row_expr<Row, N, IsEmbed>& expr)
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

	template<class Col, int N, bool IsEmbed>
	struct linear_vector_evaluator<repeat_col_expr<Col, N, IsEmbed> >
	{
		static const int M = ct_rows<Col>::value;

		typedef typename
				if_c<N == 1,
					single_vec_linear_evaluator<Col, IsEmbed>,
					typename
					if_c<M == 1,
						rep_scalar_linear_evaluator<Col, IsEmbed>,
						cached_linear_evaluator<typename matrix_traits<Col>::value_type>
					>::type
				>::type type;
	};


	template<class Row, int M, bool IsEmbed>
	struct linear_vector_evaluator<repeat_row_expr<Row, M, IsEmbed> >
	{
		static const int N = ct_cols<Row>::value;

		typedef typename
				if_c<M == 1,
					single_vec_linear_evaluator<Row, IsEmbed>,
					typename
					if_c<N == 1,
						rep_scalar_linear_evaluator<Row, IsEmbed>,
						cached_linear_evaluator<typename matrix_traits<Row>::value_type>
					>::type
				>::type type;
	};


	template<class Col, int N, bool IsEmbed>
	struct percol_vector_evaluator<repeat_col_expr<Col, N, IsEmbed> >
	{
		static const int M = ct_rows<Col>::value;

		typedef typename
				if_c<N == 1,
					single_vec_percol_evaluator<Col, IsEmbed>,
					typename
					if_c<M == 1,
						rep_scalar_percol_evaluator<Col, IsEmbed>,
						repcol_percol_evaluator<Col, IsEmbed>
					>::type
				>::type type;
	};


	template<class Row, int M, bool IsEmbed>
	struct percol_vector_evaluator<repeat_row_expr<Row, M, IsEmbed> >
	{
		static const int N = ct_cols<Row>::value;

		typedef typename
				if_c<M == 1,
					single_vec_percol_evaluator<Row, IsEmbed>,
					typename
					if_c<N == 1,
						rep_scalar_percol_evaluator<Row, IsEmbed>,
						reprow_percol_evaluator<Row, IsEmbed>
					>::type
				>::type type;
	};


	/********************************************
	 *
	 *  Cost model
	 *
	 ********************************************/

	template<class Col, int N, bool IsEmbed>
	struct linear_eval_cost<repeat_col_expr<Col, N, IsEmbed> >
	{
		static const int M = ct_rows<Col>::value;
		static const int _cost = (M == 1 || N == 1) ? 0 : VEC_EVAL_CACHE_COST;

		LMAT_ENSURE_INLINE
		static int of(const repeat_col_expr<Col, N>& )
		{
			return _cost;
		}
	};

	template<class Row, int M, bool IsEmbed>
	struct linear_eval_cost<repeat_row_expr<Row, M, IsEmbed> >
	{
		static const int N = ct_cols<Row>::value;
		static const int _cost = (M == 1 || N == 1) ? 0 : VEC_EVAL_CACHE_COST;

		LMAT_ENSURE_INLINE
		static int of(const repeat_row_expr<Row, M>& )
		{
			return _cost;
		}
	};


	template<class Col, int N, bool IsEmbed>
	struct percol_eval_cost<repeat_col_expr<Col, N, IsEmbed> >
	{
		LMAT_ENSURE_INLINE
		static int of(const repeat_col_expr<Col, N, IsEmbed>& )
		{
			return 0;
		}
	};

	template<class Row, int M, bool IsEmbed>
	struct percol_eval_cost<repeat_row_expr<Row, M, IsEmbed> >
	{
		LMAT_ENSURE_INLINE
		static int of(const repeat_row_expr<Row, M, IsEmbed>& )
		{
			return 0;
		}
	};


}

#endif /* REPEAT_VECS_EVAL_H_ */
