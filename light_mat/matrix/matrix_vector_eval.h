/**
 * @file matrix_vector_eval.h
 *
 * Generic vector-based matrix evaluation
 *
 * @author Dahua Lin
 */

#ifdef _MSC_VER
#pragma once
#endif

#ifndef LIGHTMAT_MATRIX_VECTOR_EVAL_H_
#define LIGHTMAT_MATRIX_VECTOR_EVAL_H_

#include <light_mat/matrix/vector_eval_scheme.h>

namespace lmat
{

	/********************************************
	 *
	 *  implementation of vector evaluation
	 *
	 ********************************************/

	template<typename T, int CTSize, typename Ker> struct linear_eval_impl;
	template<typename T, int CTRows, int CTCols, typename Ker> struct percol_eval_impl;

	template<class Expr, class Dst, typename Sch, typename Ker> struct vector_eval_impl_map;

	template<class Expr, class Dst, typename Ker>
	struct vector_eval_impl_map<Expr, Dst, as_linear_vec, Ker>
	{
		typedef linear_eval_impl<
			typename matrix_traits<Expr>::value_type,
			binary_ct_size<Expr, Dst>::value,
			Ker> type;
	};

	template<class Expr, class Dst, typename Ker>
	struct vector_eval_impl_map<Expr, Dst, per_column, Ker>
	{
		typedef percol_eval_impl<
			typename matrix_traits<Expr>::value_type,
			binary_ct_rows<Expr, Dst>::value,
			binary_ct_cols<Expr, Dst>::value,
			Ker> type;
	};



	template<typename T, int CTSize>
	struct linear_eval_impl<T, CTSize, scalar_kernel_t>
	{
#ifdef LMAT_USE_STATIC_ASSERT
		static_assert(CTSize > 0, "CTSize must be positive.");
#endif

		template<class Sch, class Mat>
		LMAT_ENSURE_INLINE
		static void evaluate(
				const IVecEvalLinearScheme<Sch, T>& sch,
				IDenseMatrix<Mat, T>& dst)
		{
			typedef typename vector_eval_scheme_traits<Sch>::kernel_type kernel_t;
			typedef typename vector_eval_kernel_state<kernel_t>::type state_t;

			kernel_t ker = sch.kernel();
			state_t s = sch.vec_state();
			T* pd = dst.ptr_data();

			for (index_t i = 0; i < CTSize; ++i)
			{
				pd[i] = ker.get_value(i, s);
			}
		}
	};

	template<typename T>
	struct linear_eval_impl<T, DynamicDim, scalar_kernel_t>
	{
		template<class Sch, class Mat>
		LMAT_ENSURE_INLINE
		static void evaluate(
				const IVecEvalLinearScheme<Sch, T>& sch,
				IDenseMatrix<Mat, T>& dst)
		{
			typedef typename vector_eval_scheme_traits<Sch>::kernel_type kernel_t;
			typedef typename vector_eval_kernel_state<kernel_t>::type state_t;

			kernel_t ker = sch.kernel();
			state_t s = sch.vec_state();
			T* pd = dst.ptr_data();
			const index_t len = dst.nelems();

			for (index_t i = 0; i < len; ++i)
			{
				pd[i] = ker.get_value(i, s);
			}
		}
	};


	template<typename T, int CTRows, int CTCols>
	struct percol_eval_impl<T, CTRows, CTCols, scalar_kernel_t>
	{
#ifdef LMAT_USE_STATIC_ASSERT
		static_assert(CTRows > 0, "CTSize must be positive.");
#endif

		template<class Sch, class Mat>
		LMAT_ENSURE_INLINE
		static void evaluate(
				IVecEvalPerColScheme<Sch, T>& sch,
				IDenseMatrix<Mat, T>& dst)
		{
			typedef typename vector_eval_scheme_traits<Sch>::kernel_type kernel_t;
			typedef typename vector_eval_kernel_state<kernel_t>::type state_t;

			kernel_t ker = sch.kernel();
			const index_t ncols = dst.ncolumns();
			const index_t ldim = dst.lead_dim();
			T *pd = dst.ptr_data();

			for (index_t j = 0; j < ncols; ++j, pd += ldim)
			{
				state_t s = sch.col_state(j);

				for (index_t i = 0; i < CTRows; ++i)
				{
					pd[i] = ker.get_value(i, s);
				}
			}
		}
	};


	template<typename T, int CTCols>
	struct percol_eval_impl<T, DynamicDim, CTCols, scalar_kernel_t>
	{
		template<class Sch, class Mat>
		LMAT_ENSURE_INLINE
		static void evaluate(
				IVecEvalPerColScheme<Sch, T>& sch,
				IDenseMatrix<Mat, T>& dst)
		{
			typedef typename vector_eval_scheme_traits<Sch>::kernel_type kernel_t;
			typedef typename vector_eval_kernel_state<kernel_t>::type state_t;

			kernel_t ker = sch.kernel();
			const index_t nrows = dst.nrows();
			const index_t ncols = dst.ncolumns();
			const index_t ldim = dst.lead_dim();
			T *pd = dst.ptr_data();

			for (index_t j = 0; j < ncols; ++j, pd += ldim)
			{
				state_t s = sch.col_state(j);

				for (index_t i = 0; i < nrows; ++i)
				{
					pd[i] = ker.get_value(i, s);
				}
			}
		}
	};


	/********************************************
	 *
	 *  Evaluator Map and Cost Model
	 *
	 ********************************************/

	// Ker can be either scalar_kernel_t or simd_kernel_t
	// Sch can be either as_linear_vec or per_column

	template<typename Sch, typename Ker> struct vector_eval_policy { };

	// generic vector evaluation map

	const int VEC_EVAL_CACHE_COST = 1000;
	const int SHORTVEC_LENGTH_THRESHOLD = 4;
	const int SHORTVEC_PERCOL_COST = 200;

	template<class Expr, typename Sch, typename Ker> struct generic_vector_eval;

	template<class Expr>
	struct generic_vector_eval<Expr, as_linear_vec, scalar_kernel_t>
	{
		static const bool can_direct = is_dense_mat<Expr>::value && ct_has_continuous_layout<Expr>::value;

		typedef typename matrix_traits<Expr>::value_type T;
		typedef cached_linear_veval_scheme<T, ct_rows<Expr>::value, ct_cols<Expr>::value> cached_mat_type;

		typedef typename
				if_c<can_direct,
					continuous_linear_veval_scheme<T>,
					cached_mat_type>::type scheme_type;

		static const int cost = can_direct ? 0 : VEC_EVAL_CACHE_COST;
	};

	template<class Expr>
	struct generic_vector_eval<Expr, per_column, scalar_kernel_t>
	{
		static const bool can_direct = is_dense_mat<Expr>::value;
		static const bool has_short_col = ct_rows<Expr>::value < SHORTVEC_LENGTH_THRESHOLD;

		typedef typename matrix_traits<Expr>::value_type T;
		typedef cached_percol_veval_scheme<T, ct_rows<Expr>::value, ct_cols<Expr>::value> cached_mat_type;

		typedef typename
				if_c<can_direct,
					dense_percol_veval_scheme<T>,
					cached_mat_type>::type scheme_type;

		static const int normal_cost = can_direct ? 0 : VEC_EVAL_CACHE_COST;
		static const int shortv_cost = SHORTVEC_PERCOL_COST + normal_cost;

		static const int cost = has_short_col ? shortv_cost : normal_cost;
	};


	// vector evaluation map

	template<class Expr, typename Sch, typename Ker>
	struct vector_eval
	{
		typedef typename generic_vector_eval<Expr, Sch, Ker>::scheme_type scheme_type;
		static const int cost = generic_vector_eval<Expr, Sch, Ker>::cost;
	};

	template<class Expr, typename Ker>
	struct vector_eval<Expr, per_column, Ker>
	{
		typedef typename generic_vector_eval<Expr, per_column, Ker>::scheme_type scheme_type;
		static const int normal_cost = generic_vector_eval<Expr, per_column, Ker>::normal_cost;
		static const int shortv_cost = generic_vector_eval<Expr, per_column, Ker>::shortv_cost;
		static const int cost = generic_vector_eval<Expr, per_column, Ker>::cost;
	};


	template<typename T, int CTRows, int CTCols>
	struct vector_eval<const_matrix<T, CTRows, CTCols>, as_linear_vec, scalar_kernel_t>
	{
		typedef const_linear_veval_scheme<T> scheme_type;

		static const int cost = 0;
	};

	template<typename T, int CTRows, int CTCols>
	struct vector_eval<const_matrix<T, CTRows, CTCols>, per_column, scalar_kernel_t>
	{
		typedef const_percol_veval_scheme<T> scheme_type;

		static const int normal_cost = 0;
		static const int shortv_cost = 0;
		static const int cost = 0;
	};


	// default policy

	template<class Expr>
	struct vector_eval_default_policy
	{
		// Note: SIMD has not been implemented

		typedef scalar_kernel_t kernel_t;

		static const int linear_cost = vector_eval<Expr, as_linear_vec, kernel_t>::cost;
		static const int percol_cost = vector_eval<Expr, per_column, kernel_t>::cost;

		static const bool choose_linear = (linear_cost <= percol_cost);

		typedef typename if_c<choose_linear, as_linear_vec, per_column>::type sch_t;

		typedef vector_eval_policy<sch_t, kernel_t> type;
	};


	/********************************************
	 *
	 *  Evaluation functions
	 *
	 ********************************************/

	template<typename T, class Expr, class Dst, typename Sch, typename Ker>
	LMAT_ENSURE_INLINE
	void evaluate(const IMatrixXpr<Expr, T>& src, IDenseMatrix<Dst, T>& dst, vector_eval_policy<Sch, Ker>)
	{
		typedef typename vector_eval_impl_map<Expr, Dst, Sch, Ker>::type impl_t;
		typename vector_eval<Expr, Sch, Ker>::scheme_type scheme(src.derived());
		impl_t::evaluate(scheme, dst.derived());
	}

	template<typename T, class Expr, class Dst>
	LMAT_ENSURE_INLINE
	void linear_by_scalars_evaluate(const IMatrixXpr<Expr, T>& expr, IDenseMatrix<Dst, T>& dst)
	{
		evaluate(expr.derived(), dst.derived(), vector_eval_policy<as_linear_vec, scalar_kernel_t>());
	}

	template<typename T, class Expr, class Dst>
	LMAT_ENSURE_INLINE
	void percol_by_scalars_evaluate(const IMatrixXpr<Expr, T>& expr, IDenseMatrix<Dst, T>& dst)
	{
		evaluate(expr.derived(), dst.derived(), vector_eval_policy<per_column, scalar_kernel_t>());
	}

}

#endif
