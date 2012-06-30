/**
 * @file vec_evaluator_concepts.h
 *
 * Generic concepts for vector-based evaluation
 *
 * @author Dahua Lin
 */

#ifdef _MSC_VER
#pragma once
#endif

#ifndef LIGHTMAT_VEC_EVALUATOR_CONCEPTS_H_
#define LIGHTMAT_VEC_EVALUATOR_CONCEPTS_H_

#include <light_mat/matrix/matrix_concepts.h>

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

		T get_value(const index_t i) const
		{
			return derived().get_value(i);
		}

		void next_column()
		{
			derived().next_column();
		}
	};


	/********************************************
	 *
	 *  evaluation context
	 *
	 ********************************************/

	template<typename T, int CTSize>
	class linear_eval_context
	{
#ifdef LMAT_USE_STATIC_ASSERT
		static_assert(CTSize > 0, "CTSize must be positive.");
#endif

	private:
		T *dst;

	public:
		template<class Mat>
		LMAT_ENSURE_INLINE
		linear_eval_context(IDenseMatrix<Mat, T>& mat)
		{
#ifdef LMAT_USE_STATIC_ASSERT
			static_assert(has_continuous_layout<Mat>::value,
					"Mat must have continuous layout");
#endif
			dst = mat.ptr_data();
		}

		template<class Eva>
		LMAT_ENSURE_INLINE
		void eval_by_scalars(const ILinearVectorEvaluator<Eva, T>& evaluator)
		{
			for (index_t i = 0; i < CTSize; ++i)
			{
				dst[i] = evaluator.get_value(i);
			}
		}
	};

	template<typename T>
	class linear_eval_context<T, DynamicDim>
	{
	private:
		const index_t nelems;
		T *dst;

	public:
		template<class Mat>
		LMAT_ENSURE_INLINE
		linear_eval_context(IDenseMatrix<Mat, T>& mat)
		{
#ifdef LMAT_USE_STATIC_ASSERT
			static_assert(has_continuous_layout<Mat>::value,
					"Mat must have continuous layout");
#endif
			nelems = mat.nelems();
			dst = mat.ptr_data();
		}

		template<class Eva>
		LMAT_ENSURE_INLINE
		void eval_by_scalars(const ILinearVectorEvaluator<Eva, T>& evaluator) const
		{
			for (index_t i = 0; i < nelems; ++i)
			{
				dst[i] = evaluator.get_value(i);
			}
		}
	};


	template<typename T, int CTRows, int CTCols>
	class percol_eval_context
	{
#ifdef LMAT_USE_STATIC_ASSERT
		static_assert(CTRows > 0, "CTRows must be positive.");
#endif

	private:
		const index_t ncols;
		const index_t ldim;
		T *dst;

	public:
		template<class Mat>
		LMAT_ENSURE_INLINE
		percol_eval_context(IDenseMatrix<Mat, T>& mat)
		{
			ncols = mat.ncolumns();
			ldim = mat.lead_dim();
			dst = mat.ptr_data();
		}

		template<class Eva>
		LMAT_ENSURE_INLINE
		void eval_by_scalars(IPerColVectorEvaluator<Eva, T>& evaluator) const
		{
			T *pd = dst;
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
	class percol_eval_context<T, DynamicDim, CTCols>
	{
	private:
		const index_t nrows;
		const index_t ncols;
		const index_t ldim;
		T *dst;

	public:
		template<class Mat>
		LMAT_ENSURE_INLINE
		percol_eval_context(IDenseMatrix<Mat, T>& mat)
		{
			nrows = mat.nrows();
			ncols = mat.ncolumns();
			ldim = mat.lead_dim();
			dst = mat.ptr_data();
		}

		template<class Eva>
		LMAT_ENSURE_INLINE
		void eval_by_scalars(IPerColVectorEvaluator<Eva, T>& evaluator) const
		{
			T *pd = dst;
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



}

#endif /* VEC_EVALUATOR_CONCEPTS_H_ */
