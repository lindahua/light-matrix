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

	template<class Fun, typename Arg_HP, class Arg>
	struct matrix_traits<colwise_reduce_expr<Fun, Arg_HP, Arg> >
	{
		static const int num_dimensions = 2;
		static const int compile_time_num_rows = 1;
		static const int compile_time_num_cols = ct_cols<Arg>::value;

		static const bool is_readonly = true;

		typedef typename Fun::result_type value_type;
		typedef typename matrix_traits<Arg>::domain domain;
	};

	template<class Fun, typename Arg_HP, class Arg>
	struct matrix_traits<rowwise_reduce_expr<Fun, Arg_HP, Arg> >
	{
		static const int num_dimensions = 2;
		static const int compile_time_num_rows = ct_rows<Arg>::value;
		static const int compile_time_num_cols = 1;

		static const bool is_readonly = true;

		typedef typename Fun::result_type value_type;
		typedef typename matrix_traits<Arg>::domain domain;
	};


	template<class Fun, typename Arg_HP, class Arg>
	class colwise_reduce_expr
	: public unary_expr_base<Arg_HP, Arg>
	, public IMatrixXpr<
	  	  colwise_reduce_expr<Fun, Arg_HP, Arg>,
	  	  typename Fun::result_type>
	{
#ifdef LMAT_USE_STATIC_ASSERT
		static_assert(is_reduction_functor<Fun>::value, "Fun must be a reduction_functor");
		static_assert(is_mat_xpr<Arg>::value, "Arg must be a matrix expression class.");
#endif

	public:
		typedef unary_expr_base<Arg_HP, Arg> base_t;
		typedef typename Fun::result_type value_type;

		LMAT_ENSURE_INLINE
		colwise_reduce_expr(const Fun& fun, const arg_forwarder<Arg_HP, Arg>& arg_fwd)
		: base_t(arg_fwd), m_fun(fun) { }

		LMAT_ENSURE_INLINE const Fun& fun() const
		{
			return m_fun;
		}

		LMAT_ENSURE_INLINE index_t nelems() const
		{
			return this->arg().ncolumns();
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
			return this->arg().ncolumns();
		}

	private:
		Fun m_fun;
	};


	template<class Fun, typename Arg_HP, class Arg>
	class rowwise_reduce_expr
	: public unary_expr_base<Arg_HP, Arg>
	, public IMatrixXpr<rowwise_reduce_expr<Fun, Arg_HP, Arg>, typename Fun::result_type>
	{
#ifdef LMAT_USE_STATIC_ASSERT
		static_assert(is_reduction_functor<Fun>::value, "Fun must be a reduction_functor");
		static_assert(is_mat_xpr<Arg>::value, "Arg must be a matrix expression class.");
#endif

	public:
		typedef unary_expr_base<Arg_HP, Arg> base_t;
		typedef typename Fun::result_type value_type;

		LMAT_ENSURE_INLINE
		rowwise_reduce_expr(const Fun& fun, const arg_forwarder<Arg_HP, Arg>& arg_fwd)
		: base_t(arg_fwd), m_fun(fun) { }

		LMAT_ENSURE_INLINE const Fun& fun() const
		{
			return m_fun;
		}

		LMAT_ENSURE_INLINE index_t nelems() const
		{
			return this->arg().nrows();
		}

		LMAT_ENSURE_INLINE size_t size() const
		{
			return static_cast<size_t>(nelems());
		}

		LMAT_ENSURE_INLINE index_t nrows() const
		{
			return this->arg().nrows();
		}

		LMAT_ENSURE_INLINE index_t ncolumns() const
		{
			return 1;
		}

	private:
		Fun m_fun;
	};


	/********************************************
	 *
	 *  Expression Maps
	 *
	 ********************************************/

	template<class Fun>
	struct colwise_reduce_t
	{
		const Fun& fun;

		LMAT_ENSURE_INLINE
		explicit colwise_reduce_t(const Fun& f)
		: fun(f) { }
	};

	template<class Fun>
	struct rowwise_reduce_t
	{
		const Fun& fun;

		LMAT_ENSURE_INLINE
		explicit rowwise_reduce_t(const Fun& f)
		: fun(f) { }
	};

	template<class Fun>
	inline colwise_reduce_t<Fun> colwise_reduce(const Fun& f)
	{
		return colwise_reduce_t<Fun>(f);
	}

	template<class Fun>
	inline rowwise_reduce_t<Fun> rowwise_reduce(const Fun& f)
	{
		return rowwise_reduce_t<Fun>(f);
	}


	template<class Fun, class Arg>
	struct unary_expr_verifier<colwise_reduce_t<Fun>, Arg>
	{
		static const bool value = is_mat_xpr<Arg>::value;
	};

	template<class Fun, class Arg>
	struct unary_expr_verifier<rowwise_reduce_t<Fun>, Arg>
	{
		static const bool value = is_mat_xpr<Arg>::value;
	};

	template<class Fun, typename Arg_HP, class Arg>
	struct unary_expr_map<colwise_reduce_t<Fun>, Arg_HP, Arg>
	{
		typedef colwise_reduce_expr<Fun, Arg_HP, Arg> type;

		LMAT_ENSURE_INLINE
		static type get(const colwise_reduce_t<Fun>& spec,
				const arg_forwarder<Arg_HP, Arg>& arg_fwd)
		{
			return type(spec.fun, arg_fwd);
		}
	};

	template<class Fun, typename Arg_HP, class Arg>
	struct unary_expr_map<rowwise_reduce_t<Fun>, Arg_HP, Arg>
	{
		typedef rowwise_reduce_expr<Fun, Arg_HP, Arg> type;

		LMAT_ENSURE_INLINE
		static type get(const rowwise_reduce_t<Fun>& spec,
				const arg_forwarder<Arg_HP, Arg>& arg_fwd)
		{
			return type(spec.fun, arg_fwd);
		}
	};


	/********************************************
	 *
	 *  Evaluation
	 *
	 ********************************************/

	struct partial_reduce_policy { };

	template<class Fun, typename Arg_HP, class Arg, class Dst>
	struct default_matrix_eval_policy< colwise_reduce_expr<Fun, Arg_HP, Arg>, Dst >
	{
		typedef partial_reduce_policy type;
	};

	template<class Fun, typename Arg_HP, class Arg, class Dst>
	struct default_matrix_eval_policy< rowwise_reduce_expr<Fun, Arg_HP, Arg>, Dst >
	{
		typedef partial_reduce_policy type;
	};

	template<class Fun, typename Arg_HP, class Arg, class Dst>
	LMAT_ENSURE_INLINE
	inline void evaluate(const colwise_reduce_expr<Fun, Arg_HP, Arg>& expr,
			IDenseMatrix<Dst, typename Fun::result_type>& dst,
			partial_reduce_policy)
	{
		detail::colwise_reduce_internal<scalar_kernel_t>::eval(
				expr.fun(), expr.arg(), dst.derived());
	}

	template<class Fun, typename Arg_HP, class Arg, class Dst>
	LMAT_ENSURE_INLINE
	inline void evaluate(const rowwise_reduce_expr<Fun, Arg_HP, Arg>& expr,
			IDenseMatrix<Dst, typename Fun::result_type>& dst,
			partial_reduce_policy)
	{
		detail::rowwise_reduce_internal<scalar_kernel_t>::eval(
				expr.fun(), expr.arg(), dst.derived());
	}



	/********************************************
	 *
	 *  Generic partial reduction
	 *
	 ********************************************/

	template<class Fun, class Arg>
	LMAT_ENSURE_INLINE
	inline typename unary_expr_map<colwise_reduce_t<Fun>, ref_arg_t, Arg>::type
	reduce( const Fun& fun,
			const IMatrixXpr<Arg, typename Fun::arg_type>& arg,
			colwise )
	{
		return unary_expr_map<colwise_reduce_t<Fun>, ref_arg_t, Arg>::get(
				colwise_reduce(fun),
				ref_arg(arg.derived()) );
	}

	template<class Fun, class Arg>
	LMAT_ENSURE_INLINE
	inline typename unary_expr_map<rowwise_reduce_t<Fun>, ref_arg_t, Arg>::type
	reduce( const Fun& fun,
			const IMatrixXpr<Arg, typename Fun::arg_type>& arg,
			rowwise )
	{
		return unary_expr_map<rowwise_reduce_t<Fun>, ref_arg_t, Arg>::get(
				rowwise_reduce(fun),
				ref_arg(arg.derived()) );
	}


	/********************************************
	 *
	 *  Specific expressions
	 *
	 ********************************************/

	// sum

	template<typename Arg_HP, class Arg>
	struct colwise_sum_expr_map
	{
		typedef typename matrix_traits<Arg>::value_type T;

		typedef typename unary_expr_map<
					colwise_reduce_t<sum_fun<T> >, Arg_HP, Arg
				>::type type;
	};

	template<typename Arg_HP, class Arg>
	struct rowwise_sum_expr_map
	{
		typedef typename matrix_traits<Arg>::value_type T;

		typedef typename unary_expr_map<
					rowwise_reduce_t<sum_fun<T> >, Arg_HP, Arg
				>::type type;
	};

	template<typename T, class Arg>
	LMAT_ENSURE_INLINE
	inline typename colwise_sum_expr_map<ref_arg_t, Arg>::type
	sum(const IMatrixXpr<Arg, T>& arg, colwise)
	{
		return reduce(sum_fun<T>(), arg.derived(), colwise());
	}

	template<typename T, class Arg>
	LMAT_ENSURE_INLINE
	inline typename rowwise_sum_expr_map<ref_arg_t, Arg>::type
	sum(const IMatrixXpr<Arg, T>& arg, rowwise)
	{
		return reduce(sum_fun<T>(), arg.derived(), rowwise());
	}


	// mean

	template<typename Arg_HP, class Arg>
	struct colwise_mean_expr_map
	{
		typedef typename matrix_traits<Arg>::value_type T;

		typedef typename binary_fix2_ewise_expr_map<
					mul_op<T>,
					copy_arg_t,
					typename colwise_sum_expr_map<Arg_HP, Arg>::type
				>::type type;
	};

	template<typename Arg_HP, class Arg>
	struct rowwise_mean_expr_map
	{
		typedef typename matrix_traits<Arg>::value_type T;

		typedef typename binary_fix2_ewise_expr_map<
					mul_op<T>,
					copy_arg_t,
					typename rowwise_sum_expr_map<Arg_HP, Arg>::type
				>::type type;
	};

	template<typename T, class Arg>
	LMAT_ENSURE_INLINE
	inline typename colwise_mean_expr_map<ref_arg_t, Arg>::type
	mean(const IMatrixXpr<Arg, T>& arg, colwise)
	{
		const_matrix<T, 1, ct_cols<Arg>::value>
		s( 1, arg.ncolumns(), math::rcp(T(arg.nrows())) );

		return make_expr(
				ewise(mul_op<T>()),
				copy_arg(
					make_expr(
						colwise_reduce(sum_fun<T>()),
						ref_arg(arg.derived())
					)
				),
				copy_arg(s)
		);
	}

	template<typename T, class Arg>
	LMAT_ENSURE_INLINE
	inline typename rowwise_mean_expr_map<ref_arg_t, Arg>::type
	mean(const IMatrixXpr<Arg, T>& arg, rowwise)
	{
		const_matrix<T, ct_rows<Arg>::value, 1>
		s( arg.nrows(), 1, math::rcp(T(arg.ncolumns())) );

		return make_expr(
				ewise(mul_op<T>()),
				copy_arg(
					make_expr(
						rowwise_reduce(sum_fun<T>()),
						ref_arg(arg.derived())
					)
				),
				copy_arg(s)
		);
	}


	// prod

	template<typename Arg_HP, class Arg>
	struct colwise_prod_expr_map
	{
		typedef typename unary_expr_map<
				colwise_reduce_t<
					prod_fun<typename matrix_traits<Arg>::value_type>
				>, Arg_HP, Arg>::type type;
	};

	template<typename Arg_HP, class Arg>
	struct rowwise_prod_expr_map
	{
		typedef typename unary_expr_map<
				rowwise_reduce_t<
					prod_fun<typename matrix_traits<Arg>::value_type>
				>, Arg_HP, Arg>::type type;
	};

	template<typename T, class Arg>
	LMAT_ENSURE_INLINE
	inline typename colwise_prod_expr_map<ref_arg_t, Arg>::type
	prod(const IMatrixXpr<Arg, T>& arg, colwise)
	{
		return reduce(prod_fun<T>(), arg.derived(), colwise());
	}

	template<typename T, class Arg>
	LMAT_ENSURE_INLINE
	inline typename rowwise_prod_expr_map<ref_arg_t, Arg>::type
	prod(const IMatrixXpr<Arg, T>& arg, rowwise)
	{
		return reduce(prod_fun<T>(), arg.derived(), rowwise());
	}


	// maximum

	template<typename Arg_HP, class Arg>
	struct colwise_maximum_expr_map
	{
		typedef typename matrix_traits<Arg>::value_type T;

		typedef typename unary_expr_map<
				colwise_reduce_t<
					maximum_fun<T>
				>, Arg_HP, Arg>::type type;
	};

	template<typename Arg_HP, class Arg>
	struct rowwise_maximum_expr_map
	{
		typedef typename matrix_traits<Arg>::value_type T;

		typedef typename unary_expr_map<
				rowwise_reduce_t<
					maximum_fun<T>
				>, Arg_HP, Arg>::type type;
	};

	template<typename T, class Arg>
	LMAT_ENSURE_INLINE
	inline typename colwise_maximum_expr_map<ref_arg_t, Arg>::type
	maximum(const IMatrixXpr<Arg, T>& arg, colwise)
	{
		return reduce(maximum_fun<T>(), arg.derived(), colwise());
	}

	template<typename T, class Arg>
	LMAT_ENSURE_INLINE
	inline typename rowwise_maximum_expr_map<ref_arg_t, Arg>::type
	maximum(const IMatrixXpr<Arg, T>& arg, rowwise)
	{
		return reduce(maximum_fun<T>(), arg.derived(), rowwise());
	}


	// minimum

	template<typename Arg_HP, class Arg>
	struct colwise_minimum_expr_map
	{
		typedef typename matrix_traits<Arg>::value_type T;

		typedef typename unary_expr_map<
				colwise_reduce_t<
					minimum_fun<T>
				>, Arg_HP, Arg>::type type;
	};

	template<typename Arg_HP, class Arg>
	struct rowwise_minimum_expr_map
	{
		typedef typename matrix_traits<Arg>::value_type T;

		typedef typename unary_expr_map<
				rowwise_reduce_t<
					minimum_fun<T>
				>, Arg_HP, Arg>::type type;
	};

	template<typename T, class Arg>
	LMAT_ENSURE_INLINE
	inline typename colwise_minimum_expr_map<ref_arg_t, Arg>::type
	minimum(const IMatrixXpr<Arg, T>& arg, colwise)
	{
		return reduce(minimum_fun<T>(), arg.derived(), colwise());
	}

	template<typename T, class Arg>
	LMAT_ENSURE_INLINE
	inline typename rowwise_minimum_expr_map<ref_arg_t, Arg>::type
	minimum(const IMatrixXpr<Arg, T>& arg, rowwise)
	{
		return reduce(minimum_fun<T>(), arg.derived(), rowwise());
	}


	// dot

	template<typename LArg_HP, class LArg, typename RArg_HP, class RArg>
	struct colwise_dot_expr_map
	{
		typedef typename binary_value_type<LArg, RArg>::type T;

		typedef typename colwise_sum_expr_map<
					copy_arg_t,
					typename binary_expr_map<
						ewise_t< mul_op<T> >,
						LArg_HP, LArg,
						RArg_HP, RArg>::type
				>::type type;
	};

	template<typename LArg_HP, class LArg, typename RArg_HP, class RArg>
	struct rowwise_dot_expr_map
	{
		typedef typename binary_value_type<LArg, RArg>::type T;

		typedef typename rowwise_sum_expr_map<
					copy_arg_t,
					typename binary_expr_map<
						ewise_t< mul_op<T> >,
						LArg_HP, LArg,
						RArg_HP, RArg>::type
				>::type type;
	};

	template<typename T, class LArg, class RArg>
	LMAT_ENSURE_INLINE
	inline typename colwise_dot_expr_map<ref_arg_t, LArg, ref_arg_t, RArg>::type
	dot(const IMatrixXpr<LArg, T>& a, const IMatrixXpr<RArg, T>& b, colwise)
	{
		return make_expr(colwise_reduce(sum_fun<T>()),
				copy_arg(a.derived() * b.derived()));
	}

	template<typename T, class LArg, class RArg>
	LMAT_ENSURE_INLINE
	inline typename rowwise_dot_expr_map<ref_arg_t, LArg, ref_arg_t, RArg>::type
	dot(const IMatrixXpr<LArg, T>& a, const IMatrixXpr<RArg, T>& b, rowwise)
	{
		return make_expr(rowwise_reduce(sum_fun<T>()),
				copy_arg(a.derived() * b.derived()));
	}


	// L1norm

	template<typename Arg_HP, class Arg>
	struct colwise_L1norm_expr_map
	{
		typedef typename matrix_traits<Arg>::value_type T;

		typedef typename colwise_sum_expr_map<
					copy_arg_t,
					typename unary_expr_map<
						ewise_t<abs_fun<T> >,
						Arg_HP, Arg
					>::type
				>::type type;
	};

	template<typename Arg_HP, class Arg>
	struct rowwise_L1norm_expr_map
	{
		typedef typename matrix_traits<Arg>::value_type T;

		typedef typename rowwise_sum_expr_map<
					copy_arg_t,
					typename unary_expr_map<
						ewise_t<abs_fun<T> >,
						Arg_HP, Arg
					>::type
				>::type type;
	};

	template<typename T, class Arg>
	LMAT_ENSURE_INLINE
	inline typename colwise_L1norm_expr_map<ref_arg_t, Arg>::type
	L1norm(const IMatrixXpr<Arg, T>& arg, colwise)
	{
		return make_expr(
					colwise_reduce(sum_fun<T>()),
					copy_arg(abs(arg.derived()))
				);
	}

	template<typename T, class Arg>
	LMAT_ENSURE_INLINE
	inline typename rowwise_L1norm_expr_map<ref_arg_t, Arg>::type
	L1norm(const IMatrixXpr<Arg, T>& arg, rowwise)
	{
		return make_expr(
					rowwise_reduce(sum_fun<T>()),
					copy_arg(abs(arg.derived()))
				);
	}


	// sqL2norm

	template<typename Arg_HP, class Arg>
	struct colwise_sqL2norm_expr_map
	{
		typedef typename matrix_traits<Arg>::value_type T;

		typedef typename colwise_sum_expr_map<
					copy_arg_t,
					typename unary_expr_map<
						ewise_t<sqr_fun<T> >,
						Arg_HP, Arg
					>::type
				>::type type;
	};

	template<typename Arg_HP, class Arg>
	struct rowwise_sqL2norm_expr_map
	{
		typedef typename matrix_traits<Arg>::value_type T;

		typedef typename rowwise_sum_expr_map<
					copy_arg_t,
					typename unary_expr_map<
						ewise_t<sqr_fun<T> >,
						Arg_HP, Arg
					>::type
				>::type type;
	};

	template<typename T, class Arg>
	LMAT_ENSURE_INLINE
	inline typename colwise_sqL2norm_expr_map<ref_arg_t, Arg>::type
	sqL2norm(const IMatrixXpr<Arg, T>& arg, colwise)
	{
		return make_expr(
					colwise_reduce(sum_fun<T>()),
					copy_arg(sqr(arg.derived()))
				);
	}

	template<typename T, class Arg>
	LMAT_ENSURE_INLINE
	inline typename rowwise_sqL2norm_expr_map<ref_arg_t, Arg>::type
	sqL2norm(const IMatrixXpr<Arg, T>& arg, rowwise)
	{
		return make_expr(
					rowwise_reduce(sum_fun<T>()),
					copy_arg(sqr(arg.derived()))
				);
	}


	// L2norm

	template<typename Arg_HP, class Arg>
	struct colwise_L2norm_expr_map
	{
		typedef typename matrix_traits<Arg>::value_type T;

		typedef typename unary_expr_map<
					ewise_t< sqrt_fun<T> >,
					copy_arg_t,
					typename colwise_sqL2norm_expr_map<Arg_HP, Arg>::type
				>::type type;
	};

	template<typename Arg_HP, class Arg>
	struct rowwise_L2norm_expr_map
	{
		typedef typename matrix_traits<Arg>::value_type T;

		typedef typename unary_expr_map<
					ewise_t< sqrt_fun<T> >,
					copy_arg_t,
					typename rowwise_sqL2norm_expr_map<Arg_HP, Arg>::type
				>::type type;
	};

	template<typename T, class Arg>
	LMAT_ENSURE_INLINE
	inline typename colwise_L2norm_expr_map<ref_arg_t, Arg>::type
	L2norm(const IMatrixXpr<Arg, T>& arg, colwise)
	{
		return make_expr(
					ewise(sqrt_fun<T>()),
					copy_arg(
						make_expr(
							colwise_reduce(sum_fun<T>()),
							copy_arg( sqr(arg.derived()))
						)
					)
				);
	}

	template<typename T, class Arg>
	LMAT_ENSURE_INLINE
	inline typename rowwise_L2norm_expr_map<ref_arg_t, Arg>::type
	L2norm(const IMatrixXpr<Arg, T>& arg, rowwise)
	{
		return make_expr(
					ewise(sqrt_fun<T>()),
					copy_arg(
						make_expr(
							rowwise_reduce(sum_fun<T>()),
							copy_arg( sqr(arg.derived()))
						)
					)
				);
	}


	// Linfnorm

	template<typename Arg_HP, class Arg>
	struct colwise_Linfnorm_expr_map
	{
		typedef typename matrix_traits<Arg>::value_type T;

		typedef typename colwise_maximum_expr_map<
					copy_arg_t,
					typename unary_expr_map<ewise_t<abs_fun<T> >, Arg_HP, Arg>::type
				>::type type;
	};

	template<typename Arg_HP, class Arg>
	struct rowwise_Linfnorm_expr_map
	{
		typedef typename matrix_traits<Arg>::value_type T;

		typedef typename rowwise_maximum_expr_map<
					copy_arg_t,
					typename unary_expr_map<ewise_t<abs_fun<T> >, Arg_HP, Arg>::type
				>::type type;
	};

	template<typename T, class Arg>
	LMAT_ENSURE_INLINE
	inline typename colwise_Linfnorm_expr_map<ref_arg_t, Arg>::type
	Linfnorm(const IMatrixXpr<Arg, T>& arg, colwise)
	{
		return make_expr(
					colwise_reduce(maximum_fun<T>()),
					copy_arg(abs(arg.derived()))
				);
	}

	template<typename T, class Arg>
	LMAT_ENSURE_INLINE
	inline typename rowwise_Linfnorm_expr_map<ref_arg_t, Arg>::type
	Linfnorm(const IMatrixXpr<Arg, T>& arg, rowwise)
	{
		return make_expr(
					rowwise_reduce(maximum_fun<T>()),
					copy_arg(abs(arg.derived()))
				);
	}


	// nrmdot

	template<typename LArg_HP, class LArg, typename RArg_HP, class RArg>
	struct colwise_nrmdot_expr_map
	{
		typedef typename binary_value_type<LArg, RArg>::type T;

		typedef typename binary_expr_map<
					ewise_t<div_op<T> >,
					copy_arg_t,
					typename colwise_dot_expr_map<LArg_HP, LArg, RArg_HP, RArg>::type,
					copy_arg_t,
					typename binary_expr_map<
						ewise_t<mul_op<T> >,
						copy_arg_t, typename colwise_L2norm_expr_map<LArg_HP, LArg>::type,
						copy_arg_t, typename colwise_L2norm_expr_map<RArg_HP, RArg>::type
					>::type
				>::type type;
	};

	template<typename LArg_HP, class LArg, typename RArg_HP, class RArg>
	struct rowwise_nrmdot_expr_map
	{
		typedef typename binary_value_type<LArg, RArg>::type T;

		typedef typename binary_expr_map<
					ewise_t<div_op<T> >,
					copy_arg_t,
					typename rowwise_dot_expr_map<LArg_HP, LArg, RArg_HP, RArg>::type,
					copy_arg_t,
					typename binary_expr_map<
						ewise_t<mul_op<T> >,
						copy_arg_t, typename rowwise_L2norm_expr_map<LArg_HP, LArg>::type,
						copy_arg_t, typename rowwise_L2norm_expr_map<RArg_HP, RArg>::type
					>::type
				>::type type;
	};

	template<typename T, class LArg, class RArg>
	LMAT_ENSURE_INLINE
	inline typename colwise_nrmdot_expr_map<ref_arg_t, LArg, ref_arg_t, RArg>::type
	nrmdot(const IMatrixXpr<LArg, T>& a, const IMatrixXpr<RArg, T>& b, colwise)
	{
		return make_expr(
					ewise(div_op<T>()),
					copy_arg( dot(a.derived(), b.derived(), colwise()) ),
					copy_arg(
						make_expr(
							ewise(mul_op<T>()),
							copy_arg( L2norm(a.derived(), colwise()) ),
							copy_arg( L2norm(b.derived(), colwise()) )
						)
					)
				);
	}

	template<typename T, class LArg, class RArg>
	LMAT_ENSURE_INLINE
	inline typename rowwise_nrmdot_expr_map<ref_arg_t, LArg, ref_arg_t, RArg>::type
	nrmdot(const IMatrixXpr<LArg, T>& a, const IMatrixXpr<RArg, T>& b, rowwise)
	{
		return make_expr(
					ewise(div_op<T>()),
					copy_arg( dot(a.derived(), b.derived(), rowwise()) ),
					copy_arg(
						make_expr(
							ewise(mul_op<T>()),
							copy_arg( L2norm(a.derived(), rowwise()) ),
							copy_arg( L2norm(b.derived(), rowwise()) )
						)
					)
				);
	}


}

#endif
