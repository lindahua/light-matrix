/**
 * @file vec_eval_concepts.h
 *
 * Generic concepts for vector-based evaluation
 *
 * @author Dahua Lin
 */

#ifdef _MSC_VER
#pragma once
#endif

#ifndef LIGHTMAT_VEC_EVAL_CONCEPTS_H_
#define LIGHTMAT_VEC_EVAL_CONCEPTS_H_

#include <light_mat/matrix/matrix_properties.h>

namespace lmat
{
	/********************************************
	 *
	 *  vector evaluator interfaces
	 *
	 ********************************************/

	template<class Derived, typename T>
	class ILinearVectorEvaluator
	{
	public:
		LMAT_CRTP_REF

		LMAT_ENSURE_INLINE
		T get_value(const index_t i) const
		{
			return derived().get_value(i);
		}
	};

	template<class Derived, typename T>
	class IPerColVectorEvaluator
	{
	public:
		LMAT_CRTP_REF

		LMAT_ENSURE_INLINE
		T get_value(const index_t i) const
		{
			return derived().get_value(i);
		}

		LMAT_ENSURE_INLINE
		void next_column()
		{
			derived().next_column();
		}
	};


	template<class Eva, typename T>
	struct is_linear_vector_evaluator
	{
		static const bool value = is_base_of<
				ILinearVectorEvaluator<Eva, T>,
				Eva>::value;
	};

	template<class Eva, typename T>
	struct is_percol_vector_evaluator
	{
		static const bool value = is_base_of<
				IPerColVectorEvaluator<Eva, T>,
				Eva>::value;
	};



	/********************************************
	 *
	 *  evaluation context
	 *
	 ********************************************/

	template<typename T, int CTSize>
	struct linear_byscalars_eval_impl
	{
#ifdef LMAT_USE_STATIC_ASSERT
		static_assert(CTSize > 0, "CTSize must be positive.");
#endif

		template<class Eva, class Mat>
		LMAT_ENSURE_INLINE
		static void evaluate(
				const ILinearVectorEvaluator<Eva, T>& evaluator,
				IDenseMatrix<Mat, T>& dst)
		{
			T* pd = dst.ptr_data();
			for (index_t i = 0; i < CTSize; ++i)
			{
				pd[i] = evaluator.get_value(i);
			}
		}
	};

	template<typename T>
	struct linear_byscalars_eval_impl<T, DynamicDim>
	{
		template<class Eva, class Mat>
		LMAT_ENSURE_INLINE
		static void evaluate(
				const ILinearVectorEvaluator<Eva, T>& evaluator,
				IDenseMatrix<Mat, T>& dst)
		{
			T* pd = dst.ptr_data();
			const index_t len = dst.nelems();

			for (index_t i = 0; i < len; ++i)
			{
				pd[i] = evaluator.get_value(i);
			}
		}
	};


	template<class Expr, class Dst>
	struct linear_scalar_evalctx
	{
		typedef typename binary_value_type<Expr, Dst>::type T;
		typedef linear_byscalars_eval_impl<T, binary_ct_size<Expr, Dst>::value> impl_t;

		LMAT_ENSURE_INLINE
		static void evaluate(const Expr& src, Dst& dst)
		{
			typename linear_eval<Expr>::evaluator_type evaluator(src);
			impl_t::evaluate(evaluator, dst);
		}
	};



	template<typename T, int CTRows, int CTCols>
	struct percol_byscalars_eval_impl
	{
#ifdef LMAT_USE_STATIC_ASSERT
		static_assert(CTRows > 0, "CTSize must be positive.");
#endif

		template<class Eva, class Mat>
		LMAT_ENSURE_INLINE
		static void evaluate(
				IPerColVectorEvaluator<Eva, T>& evaluator,
				IDenseMatrix<Mat, T>& dst)
		{
			const index_t ncols = dst.ncolumns();
			const index_t ldim = dst.lead_dim();
			T *pd = dst.ptr_data();

			for (index_t j = 0; j < ncols; ++j,
				pd += ldim, evaluator.next_column())
			{
				for (index_t i = 0; i < CTRows; ++i)
				{
					pd[i] = evaluator.get_value(i);
				}
			}
		}
	};


	template<typename T, int CTCols>
	struct percol_byscalars_eval_impl<T, DynamicDim, CTCols>
	{
		template<class Eva, class Mat>
		LMAT_ENSURE_INLINE
		static void evaluate(
				IPerColVectorEvaluator<Eva, T>& evaluator,
				IDenseMatrix<Mat, T>& dst)
		{
			const index_t nrows = dst.nrows();
			const index_t ncols = dst.ncolumns();
			const index_t ldim = dst.lead_dim();
			T *pd = dst.ptr_data();

			for (index_t j = 0; j < ncols; ++j,
				pd += ldim, evaluator.next_column())
			{
				for (index_t i = 0; i < nrows; ++i)
				{
					pd[i] = evaluator.get_value(i);
				}
			}
		}
	};

	template<class Expr, class Dst>
	struct percol_scalar_evalctx
	{
		typedef typename binary_value_type<Expr, Dst>::type T;
		typedef percol_byscalars_eval_impl<T,
				binary_ct_rows<Expr, Dst>::value,
				binary_ct_cols<Expr, Dst>::value> impl_t;

		LMAT_ENSURE_INLINE
		static void evaluate(const Expr& src, Dst& dst)
		{
			typename percol_eval<Expr>::evaluator_type evaluator(src);
			impl_t::evaluate(evaluator, dst);
		}
	};

}

#endif


