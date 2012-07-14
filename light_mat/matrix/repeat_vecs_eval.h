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
#include <light_mat/matrix/generic_matrix_eval.h>
#include "bits/repeat_vecs_internal.h"

namespace lmat
{

	/********************************************
	 *
	 *  Expression evaluation
	 *
	 ********************************************/

	template<class Col, class DMat, bool IsEmbed>
	inline void evaluate_to(const repeat_col_expr<Col, 1, IsEmbed>& s,
			IDenseMatrix<DMat, typename matrix_traits<Col>::value_type>& dst)
	{
		evaluate_to(s.column(), dst.derived());
	}


	template<class Col, int N, class DMat, bool IsEmbed>
	inline void evaluate_to(const repeat_col_expr<Col, N, IsEmbed>& s,
			IDenseMatrix<DMat, typename matrix_traits<Col>::value_type>& dst)
	{
		typedef typename matrix_traits<Col>::value_type T;

		if ( is_column(s) )
		{
			const int M = binary_ct_rows<Col, DMat>::value;
			ref_col<T, M> dview(dst.ptr_data(), s.nrows());
			evaluate_to(s.column(), dview);
		}
		else
		{
			typedef typename detail::repcol_ewrapper_map<Col>::type wrapper_t;
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

	template<class Row, class DMat, bool IsEmbed>
	inline void evaluate_to(const repeat_row_expr<Row, 1, IsEmbed>& s,
			IDenseMatrix<DMat, typename matrix_traits<Row>::value_type>& dst)
	{
		evaluate_to( s.row(), dst.derived() );
	}

	template<class Row, int M, class DMat, bool IsEmbed>
	inline void evaluate_to(const repeat_row_expr<Row, M, IsEmbed>& s,
			IDenseMatrix<DMat, typename matrix_traits<Row>::value_type>& dst)
	{
		typedef typename matrix_traits<Row>::value_type T;

		if ( is_row(s) )
		{
			const int N = binary_ct_cols<Row, DMat>::value;
			if (has_continuous_layout(dst))
			{
				ref_row<T, N> dview(dst.ptr_data(), s.ncolumns());
				evaluate_to(s.row(), dview);
			}
			else
			{
				ref_matrix_ex<T, 1, N> dview(dst.ptr_data(), 1, s.ncolumns(), dst.lead_dim());
				evaluate_to(s.row(), dview);
			}
		}
		else
		{
			typedef typename detail::reprow_ewrapper_map<Row>::type wrapper_t;
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



	/********************************************
	 *
	 *  Vector-based Evaluators
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
		typename percol_eval<Vec>::evaluator_type m_eval;
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
		typename linear_eval<Vec>::evaluator_type m_eval;
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

		typename percol_eval<typename wrapper_t::col_t>::evaluator_type m_eval;
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
	struct linear_eval<repeat_col_expr<Col, N, IsEmbed> >
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
				>::type evaluator_type;

		static const int cost = (N == 1 || M == 1) ? 0 : VEC_EVAL_CACHE_COST;

		LMAT_ENSURE_INLINE
		static int cost_of(const repeat_col_expr<Col, N, IsEmbed>& )
		{
			return cost;
		}
	};


	template<class Row, int M, bool IsEmbed>
	struct linear_eval<repeat_row_expr<Row, M, IsEmbed> >
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
				>::type evaluator_type;

		static const int cost = (M == 1 || N == 1) ? 0 : VEC_EVAL_CACHE_COST;

		LMAT_ENSURE_INLINE
		static int cost_of(const repeat_row_expr<Row, M, IsEmbed>& )
		{
			return cost;
		}
	};


	template<class Col, int N, bool IsEmbed>
	struct percol_eval<repeat_col_expr<Col, N, IsEmbed> >
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
				>::type evaluator_type;

		static const int cost = M < SHORTVEC_LENGTH_THRESHOLD ? SHORTVEC_PERCOL_COST : 0;

		LMAT_ENSURE_INLINE
		static int cost_of(const repeat_col_expr<Col, N, IsEmbed>& )
		{
			return cost;
		}
	};


	template<class Row, int M, bool IsEmbed>
	struct percol_eval<repeat_row_expr<Row, M, IsEmbed> >
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
				>::type evaluator_type;

		static const int cost = M < SHORTVEC_LENGTH_THRESHOLD ? SHORTVEC_PERCOL_COST : 0;

		LMAT_ENSURE_INLINE
		static int cost_of(const repeat_row_expr<Row, M, IsEmbed>& )
		{
			return cost;
		}
	};

}

#endif /* REPEAT_VECS_EVAL_H_ */
