/**
 * @file matrix_eval_cost.h
 *
 * The cost model for matrix evaluation
 *
 * @author Dahua Lin
 */

#ifdef _MSC_VER
#pragma once
#endif

#ifndef LIGHTMAT_MATRIX_EVAL_COST_H_
#define LIGHTMAT_MATRIX_EVAL_COST_H_

#include <light_mat/matrix/matrix_properties.h>
#include <light_mat/matrix/const_matrix.h>

namespace lmat
{
	const int VEC_EVAL_CACHE_COST = 1000;
	const int SHORTVEC_LENGTH_THRESHOLD = 4;
	const int SHORTVEC_PERCOL_COST = 200;

	namespace detail
	{
		template<class Expr>
		struct continuous_linear_eval_cost
		{
			LMAT_ENSURE_INLINE
			static int of(const Expr& )
			{
				return 0;
			}
		};

		template<class Expr>
		struct generic_linear_eval_cost
		{
			LMAT_ENSURE_INLINE
			static int of(const Expr& )
			{
				return VEC_EVAL_CACHE_COST;
			}
		};

		template<class Expr, int CTRows>
		struct dense_percol_eval_cost
		{
			static const int _cost =
					(CTRows < SHORTVEC_LENGTH_THRESHOLD ?
							SHORTVEC_PERCOL_COST : 0);

			LMAT_ENSURE_INLINE
			static int of(const Expr& )
			{
				return _cost;
			}
		};

		template<class Expr>
		struct dense_percol_eval_cost<Expr, DynamicDim>
		{
			LMAT_ENSURE_INLINE
			static int of(const Expr& expr)
			{
				return expr.nrows() < SHORTVEC_LENGTH_THRESHOLD ?
						SHORTVEC_PERCOL_COST : 0;
			}
		};

		template<class Expr, int CTRows>
		struct generic_percol_eval_cost
		{
			LMAT_ENSURE_INLINE
			static int of(const Expr& expr)
			{
				return VEC_EVAL_CACHE_COST +
						dense_percol_eval_cost<Expr, CTRows>::of(expr);
			}
		};
	}

	template<class Expr>
	struct linear_eval_cost
	{
		typedef typename
				if_<and_<is_dense_mat<Expr>, has_continuous_layout<Expr> >,
					detail::continuous_linear_eval_cost<Expr>,
					detail::generic_linear_eval_cost<Expr>
				>::type internal_t;

		LMAT_ENSURE_INLINE
		static int of(const Expr& expr)
		{
			return internal_t::of(expr);
		}
	};


	template<class Expr>
	struct percol_eval_cost
	{
		typedef typename
				if_<and_<is_dense_mat<Expr>, has_continuous_layout<Expr> >,
					detail::dense_percol_eval_cost<Expr, ct_rows<Expr>::value>,
					detail::generic_percol_eval_cost<Expr, ct_rows<Expr>::value>
				>::type internal_t;

		LMAT_ENSURE_INLINE
		static int of(const Expr& expr)
		{
			return internal_t::of(expr);
		}
	};


	template<typename T, int M, int N>
	struct linear_eval_cost<const_matrix<T, M, N> >
	{
		LMAT_ENSURE_INLINE
		static int of(const const_matrix<T, M, N>& ) { return 0; }
	};

	template<typename T, int M, int N>
	struct percol_eval_cost<const_matrix<T, M, N> >
	{
		LMAT_ENSURE_INLINE
		static int of(const const_matrix<T, M, N>& ) { return 0; }
	};

}

#endif /* MATRIX_EVAL_COST_H_ */
