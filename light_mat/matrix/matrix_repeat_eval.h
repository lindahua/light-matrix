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
#include <light_mat/matrix/matrix_visit_eval.h>
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
		const Arg& s = expr.arg();
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
		const Arg& s = expr.arg();
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
	 *  Vector evaluation schemes
	 *
	 ********************************************/

	// forward declarations

	template<class Arg, int ME, int NE> class single_vec_percol_mvisitor;
	template<class Arg, int ME, int NE> class single_vec_linear_mvisitor;
	template<class Arg> class rep_scalar_percol_mvisitor;
	template<class Arg> class rep_scalar_linear_mvisitor;
	template<class Arg> class repcol_percol_mvisitor;
	template<class Arg> class reprow_percol_mvisitor;

	// schemes

	// single_vec_percol_mvisitor

	template<class Arg, int ME, int NE>
	struct matrix_visitor_state<single_vec_percol_mvisitor<Arg, ME, NE> >
	{
		typedef matrix_visit_setting<percol_vis, scalar_kernel_t, ME, NE> setting_t;

		typedef typename matrix_visitor_state<
					typename matrix_vismap<Arg, setting_t>::type>::type type;
	};

	template<class Arg, int ME, int NE>
	class single_vec_percol_mvisitor
	: public IPerColMatrixScalarVisitor<
	  	  single_vec_percol_mvisitor<Arg, ME, NE>,
	  	  typename matrix_traits<Arg>::value_type>
	{
	public:
		typedef typename matrix_traits<Arg>::value_type T;

		typedef matrix_visit_setting<percol_vis, scalar_kernel_t, ME, NE> setting_t;
		typedef typename matrix_vismap<Arg, setting_t>::type arg_visitor_t;
		typedef typename matrix_visitor_state<arg_visitor_t>::type state_t;

		template<typename Arg_HP>
		LMAT_ENSURE_INLINE
		single_vec_percol_mvisitor(const horizontal_repeat_expr<Arg_HP, Arg, 1>& expr)
		: m_arg_visitor(expr.arg()) { }

		template<typename Arg_HP>
		LMAT_ENSURE_INLINE
		single_vec_percol_mvisitor(const vertical_repeat_expr<Arg_HP, Arg, 1>& expr)
		: m_arg_visitor(expr.arg()) { }

		LMAT_ENSURE_INLINE
		T get_scalar(const index_t i, const state_t& s) const
		{
			return m_arg_visitor.get_scalar(i, s);
		}

		LMAT_ENSURE_INLINE
		state_t col_state(const index_t j) const
		{
			return m_arg_visitor.col_state(j);
		}

	private:
		arg_visitor_t m_arg_visitor;
	};


	// single_vec_linear_mvisitor

	template<class Arg, int ME, int NE>
	class single_vec_linear_mvisitor
	: public ILinearMatrixScalarVisitor<
	  	  single_vec_linear_mvisitor<Arg, ME, NE>,
	  	  typename matrix_traits<Arg>::value_type>
	{
	public:
		typedef typename matrix_traits<Arg>::value_type T;

		typedef matrix_visit_setting<linear_vis, scalar_kernel_t, ME, NE> setting_t;
		typedef typename matrix_vismap<Arg, setting_t>::type arg_visitor_t;

		template<typename Arg_HP>
		LMAT_ENSURE_INLINE
		single_vec_linear_mvisitor(const horizontal_repeat_expr<Arg_HP, Arg, 1>& expr)
		: m_arg_visitor(expr.arg()) { }

		template<typename Arg_HP>
		LMAT_ENSURE_INLINE
		single_vec_linear_mvisitor(const vertical_repeat_expr<Arg_HP, Arg, 1>& expr)
		: m_arg_visitor(expr.arg()) { }

		LMAT_ENSURE_INLINE
		T get_scalar(const index_t i) const
		{
			return m_arg_visitor.get_scalar(i);
		}

	private:
		arg_visitor_t m_arg_visitor;
	};


	// rep_scalar_percol_mvisitor

	template<class Arg>
	struct matrix_visitor_state<rep_scalar_percol_mvisitor<Arg> >
	{
		typedef nil_type type;
	};

	template<class Arg>
	class rep_scalar_percol_mvisitor
	: public IPerColMatrixScalarVisitor<
	  	  rep_scalar_percol_mvisitor<Arg>,
	  	  typename matrix_traits<Arg>::value_type>
	{
	public:
		typedef typename matrix_traits<Arg>::value_type T;

		template<typename Arg_HP, int N>
		LMAT_ENSURE_INLINE
		rep_scalar_percol_mvisitor(const horizontal_repeat_expr<Arg_HP, Arg, N>& expr)
		{
			m_val = to_scalar(expr.arg());
		}

		template<typename Arg_HP, int M>
		LMAT_ENSURE_INLINE
		rep_scalar_percol_mvisitor(const vertical_repeat_expr<Arg_HP, Arg, M>& expr)
		{
			m_val = to_scalar(expr.arg());
		}

		LMAT_ENSURE_INLINE
		T get_scalar(const index_t, nil_type) const
		{
			return m_val;
		}

		LMAT_ENSURE_INLINE
		nil_type col_state(const index_t j) const
		{
			return nil_type();
		}

	private:
		T m_val;
	};


	// rep_scalar_linear_mvisitor

	template<class Arg>
	class rep_scalar_linear_mvisitor
	: public ILinearMatrixScalarVisitor<
	  	  rep_scalar_linear_mvisitor<Arg>,
	  	  typename matrix_traits<Arg>::value_type>
	{
	public:
		typedef typename matrix_traits<Arg>::value_type T;

		template<typename Arg_HP, int N>
		LMAT_ENSURE_INLINE
		rep_scalar_linear_mvisitor(const horizontal_repeat_expr<Arg_HP, Arg, N>& expr)
		{
			m_val = to_scalar(expr.arg());
		}

		template<typename Arg_HP, int M>
		LMAT_ENSURE_INLINE
		rep_scalar_linear_mvisitor(const vertical_repeat_expr<Arg_HP, Arg, M>& expr)
		{
			m_val = to_scalar(expr.arg());
		}

		LMAT_ENSURE_INLINE
		T get_scalar(const index_t ) const
		{
			return m_val;
		}

	private:
		T m_val;
	};


	// repcol_percol_mvisitor

	template<class Arg>
	struct matrix_visitor_state<repcol_percol_mvisitor<Arg> >
	{
		typedef nil_type type;
	};

	template<class Arg>
	class repcol_percol_mvisitor
	: public IPerColMatrixScalarVisitor<
	  	  repcol_percol_mvisitor<Arg>,
	  	  typename matrix_traits<Arg>::value_type>
	{
	public:
		typedef typename matrix_traits<Arg>::value_type T;

		template<typename Arg_HP, int N>
		LMAT_ENSURE_INLINE
		repcol_percol_mvisitor(const horizontal_repeat_expr<Arg_HP, Arg, N>& expr)
		: m_colwrap(expr.arg()) { }

		LMAT_ENSURE_INLINE
		T get_scalar(const index_t i, nil_type) const
		{
			return m_colwrap[i];
		}

		LMAT_ENSURE_INLINE
		nil_type col_state(const index_t ) const
		{
			return nil_type();
		}

	private:
		typename detail::repcol_ewrapper_map<Arg>::type m_colwrap;
	};


	// reprow_percol_mvisitor

	template<class Arg>
	struct matrix_visitor_state<reprow_percol_mvisitor<Arg> >
	{
		typedef typename matrix_traits<Arg>::value_type type;
	};

	template<class Arg>
	class reprow_percol_mvisitor
	: public IPerColMatrixScalarVisitor<
	  	  reprow_percol_mvisitor<Arg>,
	  	  typename matrix_traits<Arg>::value_type>
	{
	public:
		typedef typename matrix_traits<Arg>::value_type T;

		template<typename Arg_HP, int N>
		LMAT_ENSURE_INLINE
		reprow_percol_mvisitor(const vertical_repeat_expr<Arg_HP, Arg, N>& expr)
		: m_rowwrap(expr.arg()) { }

		LMAT_ENSURE_INLINE
		T get_scalar(const index_t i, const T& s) const
		{
			return s;
		}

		LMAT_ENSURE_INLINE
		T col_state(const index_t j) const
		{
			return m_rowwrap[j];
		}

	private:
		typename detail::reprow_ewrapper_map<Arg>::type m_rowwrap;
	};




	/********************************************
	 *
	 *  Dispatch
	 *
	 ********************************************/

	template<typename Arg_HP, class Arg, int N, int ME, int NE>
	struct matrix_vismap<
		horizontal_repeat_expr<Arg_HP, Arg, N>,
		matrix_visit_setting<linear_vis, scalar_kernel_t, ME, NE> >
	{
		typedef typename matrix_traits<Arg>::value_type T;
		static const int M = ct_rows<Arg>::value;

		typedef typename
				if_c<N == 1,
					single_vec_linear_mvisitor<Arg, ME, NE>,
					typename
					if_c<M == 1,
						rep_scalar_linear_mvisitor<Arg>,
						cached_linear_mvisitor<T, M, N>
					>::type
				>::type type;

		static const int cost = (N == 1 || M == 1) ? 0 : MVISIT_CACHE_COST;
	};

	template<typename Arg_HP, class Arg, int M, int ME, int NE>
	struct matrix_vismap<
		vertical_repeat_expr<Arg_HP, Arg, M>,
		matrix_visit_setting<linear_vis, scalar_kernel_t, ME, NE> >
	{
		typedef typename matrix_traits<Arg>::value_type T;
		static const int N = ct_cols<Arg>::value;

		typedef typename
				if_c<M == 1,
					single_vec_linear_mvisitor<Arg, ME, NE>,
					typename
					if_c<N == 1,
						rep_scalar_linear_mvisitor<Arg>,
						cached_linear_mvisitor<T, M, N>
					>::type
				>::type type;

		static const int cost = (N == 1 || M == 1) ? 0 : MVISIT_CACHE_COST;
	};


	template<typename Arg_HP, class Arg, int N, int ME, int NE>
	struct matrix_vismap<
		horizontal_repeat_expr<Arg_HP, Arg, N>,
		matrix_visit_setting<percol_vis, scalar_kernel_t, ME, NE> >
	{
		typedef typename matrix_traits<Arg>::value_type T;
		static const int M = ct_rows<Arg>::value;

		typedef typename
				if_c<N == 1,
					single_vec_percol_mvisitor<Arg, ME, NE>,
					typename
					if_c<M == 1,
						rep_scalar_percol_mvisitor<Arg>,
						repcol_percol_mvisitor<Arg>
					>::type
				>::type type;

		static const int normal_cost = 0;
		static const int shortv_cost = SHORTVEC_PERCOL_COST;
		static const int cost = ME < SHORTVEC_LENGTH_THRESHOLD ? shortv_cost : normal_cost;
	};


	template<typename Arg_HP, class Arg, int M, int ME, int NE>
	struct matrix_vismap<
		vertical_repeat_expr<Arg_HP, Arg, M>,
		matrix_visit_setting<percol_vis, scalar_kernel_t, ME, NE> >
	{
		typedef typename matrix_traits<Arg>::value_type T;
		static const int N = ct_cols<Arg>::value;

		typedef typename
				if_c<M == 1,
					single_vec_percol_mvisitor<Arg, ME, NE>,
					typename
					if_c<N == 1,
						rep_scalar_percol_mvisitor<Arg>,
						reprow_percol_mvisitor<Arg>
					>::type
				>::type type;

		static const int normal_cost = 0;
		static const int shortv_cost = SHORTVEC_PERCOL_COST;
		static const int cost = ME < SHORTVEC_LENGTH_THRESHOLD ? shortv_cost : normal_cost;
	};



}

#endif /* REPEAT_VECS_EVAL_H_ */
