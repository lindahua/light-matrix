/**
 * @file partial_reduce.h
 *
 * Partial reduction expression and evaluation
 *
 * @author Dahua Lin
 */

#ifdef _MSC_VER
#pragma once
#endif

#ifndef LIGHTMAT_PARTIAL_REDUCE_H_
#define LIGHTMAT_PARTIAL_REDUCE_H_

#include <light_mat/math/reduction_functors.h>
#include <light_mat/matrix/matrix_arith.h>
#include <light_mat/matrix/matrix_ewise_eval.h>

#include "bits/partial_reduce_internal.h"

namespace lmat
{
	/********************************************
	 *
	 *  Expression classes
	 *
	 ********************************************/

	template<class Fun, class Arg>
	struct matrix_traits<colwise_reduce_expr<Fun, Arg> >
	{
		static const int num_dimensions = 2;
		static const int compile_time_num_rows = 1;
		static const int compile_time_num_cols = ct_cols<Arg>::value;

		static const bool is_readonly = true;

		typedef typename Fun::result_type value_type;
	};


	template<class Fun, class Arg>
	struct matrix_traits<rowwise_reduce_expr<Fun, Arg> >
	{
		static const int num_dimensions = 2;
		static const int compile_time_num_rows = ct_rows<Arg>::value;
		static const int compile_time_num_cols = 1;

		static const bool is_readonly = true;

		typedef typename Fun::result_type value_type;
	};


	template<class Fun, class Arg>
	class colwise_reduce_expr
	: public IMatrixXpr<colwise_reduce_expr<Fun, Arg>, typename Fun::result_type>
	{
#ifdef LMAT_USE_STATIC_ASSERT
		static_assert(is_reduction_functor<Fun>::value, "Fun must be a reduction_functor");
		static_assert(is_mat_xpr<Arg>::value, "Arg must be a matrix expression class.");
#endif

		typedef typename unwrapped_expr<Arg>::type arg_expr_t;

	public:
		typedef typename Fun::result_type value_type;

		LMAT_ENSURE_INLINE
		colwise_reduce_expr(const Fun& fun, const Arg& a)
		: m_fun(fun), m_arg(a) { }

		LMAT_ENSURE_INLINE const Fun& fun() const
		{
			return m_fun;
		}

		LMAT_ENSURE_INLINE const arg_expr_t& arg() const
		{
			return m_arg.get();
		}

		LMAT_ENSURE_INLINE index_t nelems() const
		{
			return arg().ncolumns();
		}

		LMAT_ENSURE_INLINE size_t size() const
		{
			return static_cast<size_t>(nelems());
		}

		LMAT_ENSURE_INLINE index_t nrows() const
		{
			return 1;
		}

		LMAT_ENSURE_INLINE index_t ncolumns() const
		{
			return arg().ncolumns();
		}

	private:
		Fun m_fun;
		obj_wrapper<Arg> m_arg;
	};


	template<class Fun, class Arg>
	class rowwise_reduce_expr
	: public IMatrixXpr<rowwise_reduce_expr<Fun, Arg>, typename Fun::result_type>
	{
#ifdef LMAT_USE_STATIC_ASSERT
		static_assert(is_reduction_functor<Fun>::value, "Fun must be a reduction_functor");
		static_assert(is_mat_xpr<Arg>::value, "Arg must be a matrix expression class.");
#endif

		typedef typename unwrapped_expr<Arg>::type arg_expr_t;

	public:
		typedef typename Fun::result_type value_type;

		LMAT_ENSURE_INLINE
		rowwise_reduce_expr(const Fun& fun, const Arg& a)
		: m_fun(fun), m_arg(a) { }

		LMAT_ENSURE_INLINE const Fun& fun() const
		{
			return m_fun;
		}

		LMAT_ENSURE_INLINE const arg_expr_t& arg() const
		{
			return m_arg.get();
		}

		LMAT_ENSURE_INLINE index_t nelems() const
		{
			return arg().nrows();
		}

		LMAT_ENSURE_INLINE size_t size() const
		{
			return static_cast<size_t>(nelems());
		}

		LMAT_ENSURE_INLINE index_t nrows() const
		{
			return arg().nrows();
		}

		LMAT_ENSURE_INLINE index_t ncolumns() const
		{
			return 1;
		}

	private:
		Fun m_fun;
		obj_wrapper<Arg> m_arg;
	};


	/********************************************
	 *
	 *  Generic expressions
	 *
	 ********************************************/

	template<class Fun, class Arg>
	struct colwise_reduce_expr_map
	{
		typedef colwise_reduce_expr<Fun, Arg> type;
	};

	template<class Fun, class Arg>
	struct rowwise_reduce_expr_map
	{
		typedef rowwise_reduce_expr<Fun, Arg> type;
	};

	template<class Fun, class Arg>
	LMAT_ENSURE_INLINE
	inline typename colwise_reduce_expr_map<Fun, Arg>::type
	reduce( const Fun& fun,
			const IMatrixXpr<Arg, typename Fun::arg_type>& arg,
			colwise )
	{
		return colwise_reduce_expr<Fun, Arg>(fun, arg.derived());
	}

	template<class Fun, class Arg>
	LMAT_ENSURE_INLINE
	inline typename rowwise_reduce_expr_map<Fun, Arg>::type
	reduce( const Fun& fun,
			const IMatrixXpr<Arg, typename Fun::arg_type>& arg,
			rowwise )
	{
		return rowwise_reduce_expr<Fun, Arg>(fun, arg.derived());
	}

	template<class Fun, class Arg, class DMat>
	inline void evaluate_to(const colwise_reduce_expr<Fun, Arg>& expr,
			IDenseMatrix<DMat, typename Fun::result_type>& dst)
	{
		detail::colwise_reduce_internal::eval(expr.fun(), expr.arg(), dst.derived());
	}

	template<class Fun, class Arg, class DMat>
	inline void evaluate_to(const rowwise_reduce_expr<Fun, Arg>& expr,
			IDenseMatrix<DMat, typename Fun::result_type>& dst)
	{
		detail::rowwise_reduce_internal::eval(expr.fun(), expr.arg(), dst.derived());
	}


	/********************************************
	 *
	 *  Specific expressions
	 *
	 ********************************************/

	// sum

	template<class Arg>
	struct colwise_sum_expr_map
	{
		typedef typename colwise_reduce_expr_map<
				sum_fun<typename matrix_traits<Arg>::value_type>,
				Arg>::type type;
	};

	template<class Arg>
	struct rowwise_sum_expr_map
	{
		typedef typename rowwise_reduce_expr_map<
				sum_fun<typename matrix_traits<Arg>::value_type>,
				Arg>::type type;
	};

	template<typename T, class Arg>
	LMAT_ENSURE_INLINE
	inline typename colwise_sum_expr_map<Arg>::type
	sum(const IMatrixXpr<Arg, T>& arg, colwise)
	{
		return reduce(sum_fun<T>(), arg.derived(), colwise());
	}

	template<typename T, class Arg>
	LMAT_ENSURE_INLINE
	inline typename rowwise_sum_expr_map<Arg>::type
	sum(const IMatrixXpr<Arg, T>& arg, rowwise)
	{
		return reduce(sum_fun<T>(), arg.derived(), rowwise());
	}


	// mean

	template<class Arg>
	struct colwise_mean_expr_map
	{
		typedef typename mul_fix2_expr_map<
					embed_mat<typename colwise_sum_expr_map<Arg>::type>
				>::type type;
	};

	template<class Arg>
	struct rowwise_mean_expr_map
	{
		typedef typename mul_fix2_expr_map<
					embed_mat<typename rowwise_sum_expr_map<Arg>::type>
				>::type type;
	};

	template<typename T, class Arg>
	LMAT_ENSURE_INLINE
	inline typename colwise_mean_expr_map<Arg>::type
	mean(const IMatrixXpr<Arg, T>& arg, colwise)
	{
		return embed(sum(arg.derived(), colwise())) *
				math::rcp(T(arg.nrows()));
	}

	template<typename T, class Arg>
	LMAT_ENSURE_INLINE
	inline typename rowwise_sum_expr_map<Arg>::type
	mean(const IMatrixXpr<Arg, T>& arg, rowwise)
	{
		return embed(sum(arg.derived(), rowwise())) *
				math::rcp(T(arg.ncolumns()));
	}


	// prod

	template<class Arg>
	struct colwise_prod_expr_map
	{
		typedef typename colwise_reduce_expr_map<
				prod_fun<typename matrix_traits<Arg>::value_type>,
				Arg>::type type;
	};

	template<class Arg>
	struct rowwise_prod_expr_map
	{
		typedef typename rowwise_reduce_expr_map<
				prod_fun<typename matrix_traits<Arg>::value_type>,
				Arg>::type type;
	};

	template<typename T, class Arg>
	LMAT_ENSURE_INLINE
	inline typename colwise_prod_expr_map<Arg>::type
	prod(const IMatrixXpr<Arg, T>& arg, colwise)
	{
		return reduce(prod_fun<T>(), arg.derived(), colwise());
	}

	template<typename T, class Arg>
	LMAT_ENSURE_INLINE
	inline typename rowwise_prod_expr_map<Arg>::type
	prod(const IMatrixXpr<Arg, T>& arg, rowwise)
	{
		return reduce(prod_fun<T>(), arg.derived(), rowwise());
	}


	// maximum

	template<class Arg>
	struct colwise_maximum_expr_map
	{
		typedef typename colwise_reduce_expr_map<
				maximum_fun<typename matrix_traits<Arg>::value_type>,
				Arg>::type type;
	};

	template<class Arg>
	struct rowwise_maximum_expr_map
	{
		typedef typename rowwise_reduce_expr_map<
				maximum_fun<typename matrix_traits<Arg>::value_type>,
				Arg>::type type;
	};

	template<typename T, class Arg>
	LMAT_ENSURE_INLINE
	inline typename colwise_maximum_expr_map<Arg>::type
	maximum(const IMatrixXpr<Arg, T>& arg, colwise)
	{
		return reduce(maximum_fun<T>(), arg.derived(), colwise());
	}

	template<typename T, class Arg>
	LMAT_ENSURE_INLINE
	inline typename rowwise_maximum_expr_map<Arg>::type
	maximum(const IMatrixXpr<Arg, T>& arg, rowwise)
	{
		return reduce(maximum_fun<T>(), arg.derived(), rowwise());
	}


	// minimum

	template<class Arg>
	struct colwise_minimum_expr_map
	{
		typedef typename colwise_reduce_expr_map<
				minimum_fun<typename matrix_traits<Arg>::value_type>,
				Arg>::type type;
	};

	template<class Arg>
	struct rowwise_minimum_expr_map
	{
		typedef typename rowwise_reduce_expr_map<
				minimum_fun<typename matrix_traits<Arg>::value_type>,
				Arg>::type type;
	};

	template<typename T, class Arg>
	LMAT_ENSURE_INLINE
	inline typename colwise_minimum_expr_map<Arg>::type
	minimum(const IMatrixXpr<Arg, T>& arg, colwise)
	{
		return reduce(minimum_fun<T>(), arg.derived(), colwise());
	}

	template<typename T, class Arg>
	LMAT_ENSURE_INLINE
	inline typename rowwise_minimum_expr_map<Arg>::type
	minimum(const IMatrixXpr<Arg, T>& arg, rowwise)
	{
		return reduce(minimum_fun<T>(), arg.derived(), rowwise());
	}


	// dot

	template<class LArg, class RArg>
	struct colwise_dot_expr_map
	{
		typedef typename colwise_sum_expr_map<
					embed_mat<typename mul_expr_map<LArg, RArg>::type>
				>::type type;
	};

	template<class LArg, class RArg>
	struct rowwise_dot_expr_map
	{
		typedef typename rowwise_sum_expr_map<
					embed_mat<typename mul_expr_map<LArg, RArg>::type>
				>::type type;
	};

	template<typename T, class LArg, class RArg>
	LMAT_ENSURE_INLINE
	inline typename colwise_dot_expr_map<LArg, RArg>::type
	dot(const IMatrixXpr<LArg, T>& a, const IMatrixXpr<RArg, T>& b, colwise)
	{
		return sum(embed(a.derived() * b.derived()), colwise());
	}

	template<typename T, class LArg, class RArg>
	LMAT_ENSURE_INLINE
	inline typename rowwise_dot_expr_map<LArg, RArg>::type
	dot(const IMatrixXpr<LArg, T>& a, const IMatrixXpr<RArg, T>& b, rowwise)
	{
		return sum(embed(a.derived() * b.derived()), rowwise());
	}


	// L1norm

	template<class Arg>
	struct colwise_L1norm_expr_map
	{
		typedef typename colwise_sum_expr_map<
					embed_mat<typename abs_expr_map<Arg>::type>
				>::type type;
	};

	template<class Arg>
	struct rowwise_L1norm_expr_map
	{
		typedef typename rowwise_sum_expr_map<
					embed_mat<typename abs_expr_map<Arg>::type>
				>::type type;
	};

	template<typename T, class Arg>
	LMAT_ENSURE_INLINE
	inline typename colwise_L1norm_expr_map<Arg>::type
	L1norm(const IMatrixXpr<Arg, T>& arg, colwise)
	{
		return sum(embed(abs(arg.derived())), colwise());
	}

	template<typename T, class Arg>
	LMAT_ENSURE_INLINE
	inline typename rowwise_L1norm_expr_map<Arg>::type
	L1norm(const IMatrixXpr<Arg, T>& arg, rowwise)
	{
		return sum(embed(abs(arg.derived())), rowwise());
	}


	// sqL2norm

	template<class Arg>
	struct colwise_sqL2norm_expr_map
	{
		typedef typename colwise_sum_expr_map<
					embed_mat<typename sqr_expr_map<Arg>::type>
				>::type type;
	};

	template<class Arg>
	struct rowwise_sqL2norm_expr_map
	{
		typedef typename rowwise_sum_expr_map<
					embed_mat<typename sqr_expr_map<Arg>::type>
				>::type type;
	};

	template<typename T, class Arg>
	LMAT_ENSURE_INLINE
	inline typename colwise_sqL2norm_expr_map<Arg>::type
	sqL2norm(const IMatrixXpr<Arg, T>& arg, colwise)
	{
		return sum(embed(sqr(arg.derived())), colwise());
	}

	template<typename T, class Arg>
	LMAT_ENSURE_INLINE
	inline typename rowwise_sqL2norm_expr_map<Arg>::type
	sqL2norm(const IMatrixXpr<Arg, T>& arg, rowwise)
	{
		return sum(embed(sqr(arg.derived())), rowwise());
	}


	// L2norm

	template<class Arg>
	struct colwise_L2norm_expr_map
	{
		typedef typename sqrt_expr_map<
					embed_mat<typename colwise_sqL2norm_expr_map<Arg>::type>
				>::type type;
	};

	template<class Arg>
	struct rowwise_L2norm_expr_map
	{
		typedef typename sqrt_expr_map<
					embed_mat<typename rowwise_sqL2norm_expr_map<Arg>::type>
				>::type type;
	};

	template<typename T, class Arg>
	LMAT_ENSURE_INLINE
	inline typename colwise_L2norm_expr_map<Arg>::type
	L2norm(const IMatrixXpr<Arg, T>& arg, colwise)
	{
		return sqrt(embed(sqL2norm(arg.derived(), colwise())));
	}

	template<typename T, class Arg>
	LMAT_ENSURE_INLINE
	inline typename rowwise_L2norm_expr_map<Arg>::type
	L2norm(const IMatrixXpr<Arg, T>& arg, rowwise)
	{
		return sqrt(embed(sqL2norm(arg.derived(), rowwise())));
	}


	// Linfnorm

	template<class Arg>
	struct colwise_Linfnorm_expr_map
	{
		typedef typename colwise_maximum_expr_map<
					embed_mat<typename abs_expr_map<Arg>::type>
				>::type type;
	};

	template<class Arg>
	struct rowwise_Linfnorm_expr_map
	{
		typedef typename rowwise_maximum_expr_map<
					embed_mat<typename abs_expr_map<Arg>::type>
				>::type type;
	};

	template<typename T, class Arg>
	LMAT_ENSURE_INLINE
	inline typename colwise_Linfnorm_expr_map<Arg>::type
	Linfnorm(const IMatrixXpr<Arg, T>& arg, colwise)
	{
		return maximum(embed(abs(arg.derived())), colwise());
	}

	template<typename T, class Arg>
	LMAT_ENSURE_INLINE
	inline typename rowwise_Linfnorm_expr_map<Arg>::type
	Linfnorm(const IMatrixXpr<Arg, T>& arg, rowwise)
	{
		return maximum(embed(abs(arg.derived())), rowwise());
	}


	// nrmdot

	template<class LArg, class RArg>
	struct colwise_nrmdot_expr_map
	{
		typedef typename div_expr_map<
					embed_mat<typename colwise_dot_expr_map<LArg, RArg>::type>,
					embed_mat<typename mul_expr_map<
						embed_mat<typename colwise_L2norm_expr_map<LArg>::type>,
						embed_mat<typename colwise_L2norm_expr_map<RArg>::type>
					>::type>
				>::type type;
	};

	template<class LArg, class RArg>
	struct rowwise_nrmdot_expr_map
	{
		typedef typename div_expr_map<
					embed_mat<typename rowwise_dot_expr_map<LArg, RArg>::type>,
					embed_mat<typename mul_expr_map<
						embed_mat<typename rowwise_L2norm_expr_map<LArg>::type>,
						embed_mat<typename rowwise_L2norm_expr_map<RArg>::type>
					>::type>
				>::type type;
	};

	template<typename T, class LArg, class RArg>
	LMAT_ENSURE_INLINE
	inline typename colwise_nrmdot_expr_map<LArg, RArg>::type
	nrmdot(const IMatrixXpr<LArg, T>& a, const IMatrixXpr<RArg, T>& b, colwise)
	{
		return embed(dot(a.derived(), b.derived(), colwise())) /
				embed(
						embed(L2norm(a.derived(), colwise())) *
						embed(L2norm(b.derived(), colwise())) );
	}

	template<typename T, class LArg, class RArg>
	LMAT_ENSURE_INLINE
	inline typename rowwise_nrmdot_expr_map<LArg, RArg>::type
	nrmdot(const IMatrixXpr<LArg, T>& a, const IMatrixXpr<RArg, T>& b, rowwise)
	{
		return embed(dot(a.derived(), b.derived(), rowwise())) /
				embed(
						embed(L2norm(a.derived(), rowwise())) *
						embed(L2norm(b.derived(), rowwise())) );
	}

}

#endif
