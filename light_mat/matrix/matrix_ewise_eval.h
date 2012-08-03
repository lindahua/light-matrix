/**
 * @file matrix_ewise_eval.h
 *
 * Evaluation of element-wise matrix expressions
 *
 * @author Dahua Lin
 */

#ifdef _MSC_VER
#pragma once
#endif

#ifndef LIGHTMAT_MATRIX_EWISE_EVAL_H_
#define LIGHTMAT_MATRIX_EWISE_EVAL_H_

#include <light_mat/matrix/matrix_ewise_expr.h>
#include <light_mat/matrix/matrix_vector_eval.h>
#include <utility>

namespace lmat
{

	// forward

	template<typename Ker, class Fun, typename Arg_Kernel>
	class unary_ewise_kernel;

	template<typename Ker, class Fun, typename Arg1_Kernel, typename Arg2_Kernel>
	class binary_ewise_kernel;

	template<typename Sch, typename Ker, class Fun, typename Arg_HP, class Arg>
	class unary_ewise_veval_scheme;

	template<typename Sch, typename Ker, class Fun, typename Arg1_HP, class Arg1, typename Arg2_HP, class Arg2>
	class binary_ewise_veval_scheme;


	/********************************************
	 *
	 *  Kernels
	 *
	 ********************************************/

	template<typename Ker, class Fun, typename ArgKernel>
	struct vector_eval_kernel_state<unary_ewise_kernel<Ker, Fun, ArgKernel> >
	{
		typedef typename vector_eval_kernel_state<ArgKernel>::type type;
	};

	template<class Fun, typename Arg_Kernel>
	class unary_ewise_kernel<scalar_kernel_t, Fun, Arg_Kernel>
	: public IVecEvalScalarKernel<
	  	  unary_ewise_kernel<scalar_kernel_t, Fun, Arg_Kernel>,
	  	  typename Fun::result_type>
	{
	public:
		typedef typename Fun::result_type T;
		typedef typename vector_eval_kernel_state<Arg_Kernel>::type state_t;

		unary_ewise_kernel(const Fun& f, const Arg_Kernel& aker)
		: m_fun(f), arg_ker(aker)
		{
		}

		LMAT_ENSURE_INLINE
		T get_value(const index_t i, const state_t& s) const
		{
			return m_fun(arg_ker.get_value(i, s));
		}

	private:
		Fun m_fun;
		Arg_Kernel arg_ker;
	};


	template<typename Ker, class Fun, typename Arg1_Kernel, typename Arg2_Kernel>
	struct vector_eval_kernel_state<binary_ewise_kernel<Ker, Fun, Arg1_Kernel, Arg2_Kernel> >
	{
		typedef typename vector_eval_kernel_state<Arg1_Kernel>::type state1_t;
		typedef typename vector_eval_kernel_state<Arg2_Kernel>::type state2_t;
		typedef std::pair<state1_t, state2_t> type;
	};

	template<class Fun, typename Arg1_Kernel, typename Arg2_Kernel>
	class binary_ewise_kernel<scalar_kernel_t, Fun, Arg1_Kernel, Arg2_Kernel>
	: public IVecEvalScalarKernel<
	  	  binary_ewise_kernel<scalar_kernel_t, Fun, Arg1_Kernel, Arg2_Kernel>,
	  	  typename Fun::result_type>
	{
	public:
		typedef typename Fun::result_type T;
		typedef typename vector_eval_kernel_state<Arg1_Kernel>::type state1_t;
		typedef typename vector_eval_kernel_state<Arg2_Kernel>::type state2_t;

		binary_ewise_kernel(const Fun& f, const Arg1_Kernel& aker1, const Arg2_Kernel& aker2)
		: m_fun(f), arg1_ker(aker1), arg2_ker(aker2)
		{
		}

		LMAT_ENSURE_INLINE
		T get_value(const index_t i, const std::pair<state1_t, state2_t>& s) const
		{
			return m_fun(
					arg1_ker.get_value(i, s.first),
					arg2_ker.get_value(i, s.second));
		}

	private:
		Fun m_fun;
		Arg1_Kernel arg1_ker;
		Arg2_Kernel arg2_ker;
	};


	/********************************************
	 *
	 *  Schemes
	 *
	 ********************************************/

	// unary

	template<typename Sch, typename Ker, class Fun, typename Arg_HP, class Arg>
	struct unary_ewise_scheme_traits
	{
		typedef typename arg_holder<Arg_HP, Arg>::internal_arg_type arg_type;

		typedef typename vector_eval<arg_type, Sch, Ker>::scheme_type arg_scheme_type;

		typedef typename vector_eval_scheme_traits<arg_scheme_type>::kernel_type arg_kernel_type;

		typedef unary_ewise_kernel<Ker, Fun, arg_kernel_type> kernel_type;
	};


	template<typename Sch, typename Ker, class Fun, typename Arg_HP, class Arg>
	struct vector_eval_scheme_traits<
		unary_ewise_veval_scheme<Sch, Ker, Fun, Arg_HP, Arg> >
	: public unary_ewise_scheme_traits<Sch, Ker, Fun, Arg_HP, Arg>
	{
		typedef unary_ewise_scheme_traits<Sch, Ker, Fun, Arg_HP, Arg> _base_t;
		typedef typename _base_t::kernel_type kernel_type;
	};

	template<class Fun, typename Arg_HP, class Arg>
	class unary_ewise_veval_scheme<as_linear_vec, scalar_kernel_t, Fun, Arg_HP, Arg>
	: public unary_ewise_scheme_traits<as_linear_vec, scalar_kernel_t, Fun, Arg_HP, Arg>
	, public IVecEvalLinearScheme<
	  	  unary_ewise_veval_scheme<as_linear_vec, scalar_kernel_t, Fun, Arg_HP, Arg>,
	  	  typename Fun::result_type>
	{
	public:
		typedef unary_ewise_scheme_traits<as_linear_vec, scalar_kernel_t, Fun, Arg_HP, Arg> base_t;

		typedef typename base_t::arg_scheme_type arg_scheme_type;
		typedef typename base_t::kernel_type kernel_type;
		typedef typename vector_eval_kernel_state<kernel_type>::type kernel_state_t;

		LMAT_ENSURE_INLINE
		unary_ewise_veval_scheme(const unary_ewise_expr<Fun, Arg_HP, Arg>& expr)
		: m_fun(expr.fun()), m_arg_sch(expr.arg())
		{
		}

		LMAT_ENSURE_INLINE
		kernel_type kernel() const
		{
			return kernel_type(m_fun, m_arg_sch.kernel());
		}

		LMAT_ENSURE_INLINE
		kernel_state_t vec_state() const
		{
			return m_arg_sch.vec_state();
		}

	private:
		Fun m_fun;
		arg_scheme_type m_arg_sch;
	};


	template<class Fun, typename Arg_HP, class Arg>
	class unary_ewise_veval_scheme<per_column, scalar_kernel_t, Fun, Arg_HP, Arg>
	: public unary_ewise_scheme_traits<per_column, scalar_kernel_t, Fun, Arg_HP, Arg>
	, public IVecEvalPerColScheme<
	  	  unary_ewise_veval_scheme<per_column, scalar_kernel_t, Fun, Arg_HP, Arg>,
	  	  typename Fun::result_type>
	{
	public:
		typedef unary_ewise_scheme_traits<per_column, scalar_kernel_t, Fun, Arg_HP, Arg> base_t;

		typedef typename base_t::arg_scheme_type arg_scheme_type;
		typedef typename base_t::kernel_type kernel_type;
		typedef typename vector_eval_kernel_state<kernel_type>::type kernel_state_t;

		LMAT_ENSURE_INLINE
		unary_ewise_veval_scheme(const unary_ewise_expr<Fun, Arg_HP, Arg>& expr)
		: m_fun(expr.fun()), m_arg_sch(expr.arg())
		{
		}

		LMAT_ENSURE_INLINE
		kernel_type kernel() const
		{
			return kernel_type(m_fun, m_arg_sch.kernel());
		}

		LMAT_ENSURE_INLINE
		kernel_state_t col_state(const index_t j) const
		{
			return m_arg_sch.col_state(j);
		}

	private:
		Fun m_fun;
		arg_scheme_type m_arg_sch;
	};


	// binary

	template<typename Sch, typename Ker, class Fun, typename Arg1_HP, class Arg1, typename Arg2_HP, class Arg2>
	struct binary_ewise_scheme_traits
	{
		typedef typename arg_holder<Arg1_HP, Arg1>::internal_arg_type arg1_type;
		typedef typename arg_holder<Arg2_HP, Arg2>::internal_arg_type arg2_type;

		typedef typename vector_eval<arg1_type, Sch, Ker>::scheme_type arg1_scheme_type;
		typedef typename vector_eval<arg2_type, Sch, Ker>::scheme_type arg2_scheme_type;

		typedef typename vector_eval_scheme_traits<arg1_scheme_type>::kernel_type arg1_kernel_type;
		typedef typename vector_eval_scheme_traits<arg2_scheme_type>::kernel_type arg2_kernel_type;

		typedef binary_ewise_kernel<Ker, Fun, arg1_kernel_type, arg2_kernel_type> kernel_type;
	};

	template<typename Sch, typename Ker, class Fun, typename Arg1_HP, class Arg1, typename Arg2_HP, class Arg2>
	struct vector_eval_scheme_traits<
		binary_ewise_veval_scheme<Sch, Ker, Fun, Arg1_HP, Arg1, Arg2_HP, Arg2> >
	: public binary_ewise_scheme_traits<Sch, Ker, Fun, Arg1_HP, Arg1, Arg2_HP, Arg2>
	{
		typedef binary_ewise_scheme_traits<Sch, Ker, Fun, Arg1_HP, Arg1, Arg2_HP, Arg2> _base_t;
		typedef typename _base_t::kernel_type kernel_type;
	};

	template<class Fun, typename Arg1_HP, class Arg1, typename Arg2_HP, class Arg2>
	class binary_ewise_veval_scheme<as_linear_vec, scalar_kernel_t, Fun, Arg1_HP, Arg1, Arg2_HP, Arg2>
	: public binary_ewise_scheme_traits<as_linear_vec, scalar_kernel_t, Fun, Arg1_HP, Arg1, Arg2_HP, Arg2>
	, public IVecEvalLinearScheme<
	  	  binary_ewise_veval_scheme<as_linear_vec, scalar_kernel_t, Fun, Arg1_HP, Arg1, Arg2_HP, Arg2>,
	  	  typename Fun::result_type>
	{
	public:
		typedef binary_ewise_scheme_traits<as_linear_vec, scalar_kernel_t, Fun,
				Arg1_HP, Arg1, Arg2_HP, Arg2> base_t;

		typedef typename base_t::arg1_scheme_type arg1_scheme_type;
		typedef typename base_t::arg2_scheme_type arg2_scheme_type;
		typedef typename base_t::kernel_type kernel_type;
		typedef typename vector_eval_kernel_state<kernel_type>::type kernel_state_t;

		LMAT_ENSURE_INLINE
		binary_ewise_veval_scheme(const binary_ewise_expr<Fun, Arg1_HP, Arg1, Arg2_HP, Arg2>& expr)
		: m_fun(expr.fun()), m_arg1_sch(expr.first_arg()), m_arg2_sch(expr.second_arg())
		{
		}

		LMAT_ENSURE_INLINE
		kernel_type kernel() const
		{
			return kernel_type(m_fun, m_arg1_sch.kernel(), m_arg2_sch.kernel());
		}

		LMAT_ENSURE_INLINE
		kernel_state_t vec_state() const
		{
			return kernel_state_t(m_arg1_sch.vec_state(), m_arg2_sch.vec_state());
		}

	private:
		Fun m_fun;
		arg1_scheme_type m_arg1_sch;
		arg2_scheme_type m_arg2_sch;
	};


	template<class Fun, typename Arg1_HP, class Arg1, typename Arg2_HP, class Arg2>
	class binary_ewise_veval_scheme<per_column, scalar_kernel_t, Fun, Arg1_HP, Arg1, Arg2_HP, Arg2>
	: public binary_ewise_scheme_traits<per_column, scalar_kernel_t, Fun, Arg1_HP, Arg1, Arg2_HP, Arg2>
	, public IVecEvalPerColScheme<
	  	  binary_ewise_veval_scheme<per_column, scalar_kernel_t, Fun, Arg1_HP, Arg1, Arg2_HP, Arg2>,
	  	  typename Fun::result_type>
	{
	public:
		typedef binary_ewise_scheme_traits<per_column, scalar_kernel_t, Fun,
				Arg1_HP, Arg1, Arg2_HP, Arg2> base_t;

		typedef typename base_t::arg1_scheme_type arg1_scheme_type;
		typedef typename base_t::arg2_scheme_type arg2_scheme_type;
		typedef typename base_t::kernel_type kernel_type;
		typedef typename vector_eval_kernel_state<kernel_type>::type kernel_state_t;

		LMAT_ENSURE_INLINE
		binary_ewise_veval_scheme(const binary_ewise_expr<Fun, Arg1_HP, Arg1, Arg2_HP, Arg2>& expr)
		: m_fun(expr.fun()), m_arg1_sch(expr.first_arg()), m_arg2_sch(expr.second_arg())
		{
		}

		LMAT_ENSURE_INLINE
		kernel_type kernel() const
		{
			return kernel_type(m_fun, m_arg1_sch.kernel(), m_arg2_sch.kernel());
		}

		LMAT_ENSURE_INLINE
		kernel_state_t col_state(const index_t j) const
		{
			return kernel_state_t(m_arg1_sch.col_state(j), m_arg2_sch.col_state(j));
		}

	private:
		Fun m_fun;
		arg1_scheme_type m_arg1_sch;
		arg2_scheme_type m_arg2_sch;
	};



	/********************************************
	 *
	 *  Evaluator dispatch
	 *
	 ********************************************/

	template<typename Fun, typename Arg_HP, class Arg, typename Ker>
	struct vector_eval<unary_ewise_expr<Fun, Arg_HP, Arg>, as_linear_vec, Ker>
	{
		typedef unary_ewise_veval_scheme<as_linear_vec, Ker, Fun, Arg_HP, Arg> scheme_type;

		typedef unary_ewise_expr<Fun, Arg_HP, Arg> expr_t;
		typedef typename expr_t::arg_type arg_type;
		static const int cost = vector_eval<arg_type, as_linear_vec, Ker>::cost;
	};

	template<typename Fun, typename Arg_HP, class Arg, typename Ker>
	struct vector_eval<unary_ewise_expr<Fun, Arg_HP, Arg>, per_column, Ker>
	{
		typedef unary_ewise_veval_scheme<per_column, Ker, Fun, Arg_HP, Arg> scheme_type;

		typedef unary_ewise_expr<Fun, Arg_HP, Arg> expr_t;
		typedef typename expr_t::arg_type arg_type;
		static const int normal_cost = vector_eval<arg_type, per_column, Ker>::normal_cost;
		static const int shortv_cost = vector_eval<arg_type, per_column, Ker>::shortv_cost;

		static const bool has_short_col = ct_rows<arg_type>::value < SHORTVEC_LENGTH_THRESHOLD;
		static const int cost = has_short_col ? shortv_cost : normal_cost;
	};

	template<typename Fun, typename Arg1_HP, class Arg1, typename Arg2_HP, class Arg2, typename Ker>
	struct vector_eval<binary_ewise_expr<Fun, Arg1_HP, Arg1, Arg2_HP, Arg2>, as_linear_vec, Ker>
	{
		typedef binary_ewise_veval_scheme<as_linear_vec, Ker, Fun, Arg1_HP, Arg1, Arg2_HP, Arg2> scheme_type;

		typedef binary_ewise_expr<Fun, Arg1_HP, Arg1, Arg2_HP, Arg2> expr_t;
		typedef typename expr_t::arg1_type arg1_type;
		typedef typename expr_t::arg2_type arg2_type;

		static const int cost =
				vector_eval<arg1_type, as_linear_vec, Ker>::cost +
				vector_eval<arg2_type, as_linear_vec, Ker>::cost;
	};

	template<typename Fun, typename Arg1_HP, class Arg1, typename Arg2_HP, class Arg2, typename Ker>
	struct vector_eval<binary_ewise_expr<Fun, Arg1_HP, Arg1, Arg2_HP, Arg2>, per_column, Ker>
	{
		typedef binary_ewise_veval_scheme<per_column, Ker, Fun, Arg1_HP, Arg1, Arg2_HP, Arg2> scheme_type;

		typedef binary_ewise_expr<Fun, Arg1_HP, Arg1, Arg2_HP, Arg2> expr_t;
		typedef typename expr_t::arg1_type arg1_type;
		typedef typename expr_t::arg2_type arg2_type;

		static const int normal_cost =
				vector_eval<arg1_type, per_column, Ker>::normal_cost +
				vector_eval<arg2_type, per_column, Ker>::normal_cost;

		static const int shortv_cost =
				vector_eval<arg1_type, per_column, Ker>::shortv_cost +
				vector_eval<arg2_type, per_column, Ker>::shortv_cost;

		static const bool has_short_col = binary_ct_rows<arg1_type, arg2_type>::value < SHORTVEC_LENGTH_THRESHOLD;
		static const int cost = has_short_col ? shortv_cost : normal_cost;
	};


	template<typename Fun, typename Arg_HP, class Arg, class Dst>
	struct default_matrix_eval_policy<unary_ewise_expr<Fun, Arg_HP, Arg>, Dst>
	{
		typedef typename vector_eval_default_policy<
				unary_ewise_expr<Fun, Arg_HP, Arg> >::type type;
	};

	template<typename Fun, typename Arg1_HP, class Arg1, typename Arg2_HP, class Arg2, class Dst>
	struct default_matrix_eval_policy<binary_ewise_expr<Fun, Arg1_HP, Arg1, Arg2_HP, Arg2>, Dst>
	{
		typedef typename vector_eval_default_policy<
				binary_ewise_expr<Fun, Arg1_HP, Arg1, Arg2_HP, Arg2> >::type type;
	};

}

#endif /* MATRIX_EWISE_EVAL_H_ */
