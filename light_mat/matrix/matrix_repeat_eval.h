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
	 *  Vector evaluation schemes
	 *
	 ********************************************/

	// forward declarations

	template<typename Arg_HP, class Arg> class single_vec_percol_veval_scheme;
	template<typename Arg_HP, class Arg> class single_vec_linear_veval_scheme;
	template<typename Arg_HP, class Arg> class rep_scalar_percol_veval_scheme;
	template<typename Arg_HP, class Arg> class rep_scalar_linear_veval_scheme;
	template<typename Arg_HP, class Arg> class repcol_percol_veval_scheme;
	template<typename Arg_HP, class Arg> class reprow_percol_veval_scheme;

	// schemes

	// single_vec_percol_veval_scheme

	template<typename Arg_HP, class Arg>
	struct vector_eval_scheme_traits<single_vec_percol_veval_scheme<Arg_HP, Arg> >
	{
		typedef typename arg_holder<Arg_HP, Arg>::internal_arg_type arg_type;
		typedef typename vector_eval<arg_type, per_column, scalar_kernel_t>::scheme_type arg_scheme_type;
		typedef typename vector_eval_scheme_traits<arg_scheme_type>::kernel_type kernel_type;
	};

	template<typename Arg_HP, class Arg>
	class single_vec_percol_veval_scheme
	: public IVecEvalPerColScheme<
	  	  single_vec_percol_veval_scheme<Arg_HP, Arg>,
	  	  typename matrix_traits<Arg>::value_type>
	{
	public:
		typedef typename arg_holder<Arg_HP, Arg>::internal_arg_type arg_type;
		typedef typename vector_eval<arg_type, per_column, scalar_kernel_t>::scheme_type arg_scheme_type;
		typedef typename vector_eval_scheme_traits<arg_scheme_type>::kernel_type kernel_type;
		typedef typename vector_eval_kernel_state<kernel_type>::type state_t;

		LMAT_ENSURE_INLINE
		single_vec_percol_veval_scheme(const horizontal_repeat_expr<Arg_HP, Arg, 1>& expr)
		: m_arg_sch(expr.arg()) { }

		LMAT_ENSURE_INLINE
		single_vec_percol_veval_scheme(const vertical_repeat_expr<Arg_HP, Arg, 1>& expr)
		: m_arg_sch(expr.arg()) { }

		LMAT_ENSURE_INLINE
		kernel_type kernel() const
		{
			return m_arg_sch.kernel();
		}

		LMAT_ENSURE_INLINE
		state_t col_state(const index_t j) const { return m_arg_sch.col_state(j); }

	private:
		arg_scheme_type m_arg_sch;
	};


	// single_vec_linear_veval_scheme

	template<typename Arg_HP, class Arg>
	struct vector_eval_scheme_traits<single_vec_linear_veval_scheme<Arg_HP, Arg> >
	{
		typedef typename arg_holder<Arg_HP, Arg>::internal_arg_type arg_type;
		typedef typename vector_eval<arg_type, as_linear_vec, scalar_kernel_t>::scheme_type arg_scheme_type;
		typedef typename vector_eval_scheme_traits<arg_scheme_type>::kernel_type kernel_type;
	};

	template<typename Arg_HP, class Arg>
	class single_vec_linear_veval_scheme
	: public IVecEvalLinearScheme<
	  	  single_vec_linear_veval_scheme<Arg_HP, Arg>,
	  	  typename matrix_traits<Arg>::value_type>
	{
	public:
		typedef typename arg_holder<Arg_HP, Arg>::internal_arg_type arg_type;
		typedef typename vector_eval<arg_type, as_linear_vec, scalar_kernel_t>::scheme_type arg_scheme_type;
		typedef typename vector_eval_scheme_traits<arg_scheme_type>::kernel_type kernel_type;
		typedef typename vector_eval_kernel_state<kernel_type>::type state_t;

		LMAT_ENSURE_INLINE
		single_vec_linear_veval_scheme(const horizontal_repeat_expr<Arg_HP, Arg, 1>& expr)
		: m_arg_sch(expr.arg()) { }

		LMAT_ENSURE_INLINE
		single_vec_linear_veval_scheme(const vertical_repeat_expr<Arg_HP, Arg, 1>& expr)
		: m_arg_sch(expr.arg()) { }

		LMAT_ENSURE_INLINE
		kernel_type kernel() const
		{
			return m_arg_sch.kernel();
		}

		LMAT_ENSURE_INLINE
		state_t vec_state() const
		{
			return m_arg_sch.vec_state();
		}

	private:
		arg_scheme_type m_arg_sch;
	};


	// rep_scalar_percol_veval_scheme

	template<typename Arg_HP, class Arg>
	struct vector_eval_scheme_traits<rep_scalar_percol_veval_scheme<Arg_HP, Arg> >
	{
		typedef typename matrix_traits<Arg>::value_type T;
		typedef veval_const_kernel<T> kernel_type;
	};


	template<typename Arg_HP, class Arg>
	class rep_scalar_percol_veval_scheme
	: public IVecEvalPerColScheme<
	  	  rep_scalar_percol_veval_scheme<Arg_HP, Arg>,
	  	  typename matrix_traits<Arg>::value_type>
	{
	public:
		typedef typename arg_holder<Arg_HP, Arg>::internal_arg_type arg_type;
		typedef typename matrix_traits<Arg>::value_type T;
		typedef veval_const_kernel<T> kernel_type;

		template<int N>
		LMAT_ENSURE_INLINE
		rep_scalar_percol_veval_scheme(const horizontal_repeat_expr<Arg_HP, Arg, N>& expr)
		{
			m_val = to_scalar(expr.arg());
		}

		template<int M>
		LMAT_ENSURE_INLINE
		rep_scalar_percol_veval_scheme(const vertical_repeat_expr<Arg_HP, Arg, M>& expr)
		{
			m_val = to_scalar(expr.arg());
		}

		LMAT_ENSURE_INLINE
		kernel_type kernel() const
		{
			return kernel_type();
		}

		LMAT_ENSURE_INLINE
		T col_state(const index_t j) const { return m_val; }

	private:
		T m_val;
	};


	// rep_scalar_linear_veval_scheme

	template<typename Arg_HP, class Arg>
	struct vector_eval_scheme_traits<rep_scalar_linear_veval_scheme<Arg_HP, Arg> >
	{
		typedef typename matrix_traits<Arg>::value_type T;
		typedef veval_const_kernel<T> kernel_type;
	};

	template<typename Arg_HP, class Arg>
	class rep_scalar_linear_veval_scheme
	: public IVecEvalLinearScheme<
	  	  rep_scalar_linear_veval_scheme<Arg_HP, Arg>,
	  	  typename matrix_traits<Arg>::value_type>
	{
	public:
		typedef typename arg_holder<Arg_HP, Arg>::internal_arg_type arg_type;
		typedef typename matrix_traits<Arg>::value_type T;
		typedef veval_const_kernel<T> kernel_type;

		template<int N>
		LMAT_ENSURE_INLINE
		rep_scalar_linear_veval_scheme(const horizontal_repeat_expr<Arg_HP, Arg, N>& expr)
		{
			m_val = to_scalar(expr.arg());
		}

		template<int M>
		LMAT_ENSURE_INLINE
		rep_scalar_linear_veval_scheme(const vertical_repeat_expr<Arg_HP, Arg, M>& expr)
		{
			m_val = to_scalar(expr.arg());
		}

		LMAT_ENSURE_INLINE
		kernel_type kernel() const
		{
			return kernel_type();
		}

		LMAT_ENSURE_INLINE
		T vec_state() const { return m_val; }

	private:
		T m_val;
	};

	// repcol_percol_veval_scheme

	template<typename Arg_HP, class Arg>
	struct vector_eval_scheme_traits<repcol_percol_veval_scheme<Arg_HP, Arg> >
	{
		typedef typename matrix_traits<Arg>::value_type T;
		typedef veval_memacc_kernel<T> kernel_type;
	};


	template<typename Arg_HP, class Arg>
	class repcol_percol_veval_scheme
	: public IVecEvalPerColScheme<
	  	  repcol_percol_veval_scheme<Arg_HP, Arg>,
	  	  typename matrix_traits<Arg>::value_type>
	{
	public:
		typedef typename arg_holder<Arg_HP, Arg>::internal_arg_type arg_type;
		typedef typename matrix_traits<Arg>::value_type T;
		typedef veval_memacc_kernel<T> kernel_type;

		template<int N>
		LMAT_ENSURE_INLINE
		repcol_percol_veval_scheme(const horizontal_repeat_expr<Arg_HP, Arg, N>& expr)
		: m_colwrap(expr.arg()) { }

		LMAT_ENSURE_INLINE
		kernel_type kernel() const
		{
			return kernel_type();
		}

		LMAT_ENSURE_INLINE
		const T* col_state(const index_t ) const { return m_colwrap.data(); }

	private:
		typename detail::repcol_ewrapper_map<arg_type>::type m_colwrap;
	};


	// reprow_percol_veval_scheme

	template<typename Arg_HP, class Arg>
	struct vector_eval_scheme_traits<reprow_percol_veval_scheme<Arg_HP, Arg> >
	{
		typedef typename matrix_traits<Arg>::value_type T;
		typedef veval_const_kernel<T> kernel_type;
	};

	template<typename Arg_HP, class Arg>
	class reprow_percol_veval_scheme
	: public IVecEvalPerColScheme<
	  	  reprow_percol_veval_scheme<Arg_HP, Arg>,
	  	  typename matrix_traits<Arg>::value_type>
	{
	public:
		typedef typename arg_holder<Arg_HP, Arg>::internal_arg_type arg_type;
		typedef typename matrix_traits<Arg>::value_type T;
		typedef veval_const_kernel<T> kernel_type;

		template<int N>
		LMAT_ENSURE_INLINE
		reprow_percol_veval_scheme(const vertical_repeat_expr<Arg_HP, Arg, N>& expr)
		: m_rowwrap(expr.arg()) { }

		LMAT_ENSURE_INLINE
		kernel_type kernel() const
		{
			return kernel_type();
		}

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
	struct vector_eval<horizontal_repeat_expr<Arg_HP, Arg, N>, as_linear_vec, scalar_kernel_t>
	{
		typedef typename matrix_traits<Arg>::value_type T;
		static const int M = ct_rows<Arg>::value;

		typedef typename
				if_c<N == 1,
					single_vec_linear_veval_scheme<Arg_HP, Arg>,
					typename
					if_c<M == 1,
						rep_scalar_linear_veval_scheme<Arg_HP, Arg>,
						cached_linear_veval_scheme<T, M, N>
					>::type
				>::type scheme_type;

		static const int cost = (N == 1 || M == 1) ? 0 : VEC_EVAL_CACHE_COST;
	};


	template<typename Arg_HP, class Arg, int M>
	struct vector_eval<vertical_repeat_expr<Arg_HP, Arg, M>, as_linear_vec, scalar_kernel_t>
	{
		typedef typename matrix_traits<Arg>::value_type T;
		static const int N = ct_cols<Arg>::value;

		typedef typename
				if_c<M == 1,
					single_vec_linear_veval_scheme<Arg_HP, Arg>,
					typename
					if_c<N == 1,
						rep_scalar_linear_veval_scheme<Arg_HP, Arg>,
						cached_linear_veval_scheme<T, M, N>
					>::type
				>::type scheme_type;

		static const int cost = (M == 1 || N == 1) ? 0 : VEC_EVAL_CACHE_COST;
	};


	template<typename Arg_HP, class Arg, int N>
	struct vector_eval<horizontal_repeat_expr<Arg_HP, Arg, N>, per_column, scalar_kernel_t>
	{
		typedef typename matrix_traits<Arg>::value_type T;
		static const int M = ct_rows<Arg>::value;

		typedef typename
				if_c<N == 1,
					single_vec_percol_veval_scheme<Arg_HP, Arg>,
					typename
					if_c<M == 1,
						rep_scalar_percol_veval_scheme<Arg_HP, Arg>,
						repcol_percol_veval_scheme<Arg_HP, Arg>
					>::type
				>::type scheme_type;

		static const int normal_cost = 0;
		static const int shortv_cost = SHORTVEC_PERCOL_COST;
		static const int cost = M < SHORTVEC_LENGTH_THRESHOLD ? shortv_cost : normal_cost;
	};


	template<typename Arg_HP, class Arg, int M>
	struct vector_eval<vertical_repeat_expr<Arg_HP, Arg, M>, per_column, scalar_kernel_t>
	{
		typedef typename matrix_traits<Arg>::value_type T;
		static const int N = ct_cols<Arg>::value;

		typedef typename
				if_c<M == 1,
					single_vec_percol_veval_scheme<Arg_HP, Arg>,
					typename
					if_c<N == 1,
						rep_scalar_percol_veval_scheme<Arg_HP, Arg>,
						reprow_percol_veval_scheme<Arg_HP, Arg>
					>::type
				>::type scheme_type;

		static const int normal_cost = 0;
		static const int shortv_cost = SHORTVEC_PERCOL_COST;
		static const int cost = M < SHORTVEC_LENGTH_THRESHOLD ? shortv_cost : normal_cost;
	};

}

#endif /* REPEAT_VECS_EVAL_H_ */
