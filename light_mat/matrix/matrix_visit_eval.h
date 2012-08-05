/**
 * @file matrix_visit_eval.h
 *
 * Generic visitor-based matrix evaluation
 *
 * @author Dahua Lin
 */

#ifdef _MSC_VER
#pragma once
#endif

#ifndef LIGHTMAT_MATRIX_VISIT_EVAL_H_
#define LIGHTMAT_MATRIX_VISIT_EVAL_H_

#include <light_mat/matrix/matrix_visitors.h>

namespace lmat
{

	/********************************************
	 *
	 *  implementation of vector evaluation
	 *
	 ********************************************/

	template<int CTSize, typename KerCate> struct linear_viseval_impl;
	template<int CTRows, int CTCols, typename KerCate> struct percol_viseval_impl;

	template<typename VisSetting> struct matrix_viseval_impl_map;

	template<typename KerCate, int ME, int NE>
	struct matrix_viseval_impl_map<matrix_visit_setting<linear_vis, KerCate, ME, NE> >
	{
		typedef linear_viseval_impl<ME * NE, KerCate> type;
	};

	template<typename KerCate, int ME, int NE>
	struct matrix_viseval_impl_map<matrix_visit_setting<percol_vis, KerCate, ME, NE> >
	{
		typedef percol_viseval_impl<ME, NE, KerCate> type;
	};



	template<int CTSize>
	struct linear_viseval_impl<CTSize, scalar_kernel_t>
	{
#ifdef LMAT_USE_STATIC_ASSERT
		static_assert(CTSize > 0, "CTSize must be positive.");
#endif

		template<typename T, class Visitor, class Mat>
		LMAT_ENSURE_INLINE
		static void evaluate(
				const ILinearMatrixScalarVisitor<Visitor, T>& visitor,
				IDenseMatrix<Mat, T>& dst)
		{
			T* pd = dst.ptr_data();

			for (index_t i = 0; i < CTSize; ++i)
			{
				pd[i] = visitor.get_scalar(i);
			}
		}
	};

	template<>
	struct linear_viseval_impl<DynamicDim, scalar_kernel_t>
	{
		template<typename T, class Visitor, class Mat>
		LMAT_ENSURE_INLINE
		static void evaluate(
				const ILinearMatrixScalarVisitor<Visitor, T>& visitor,
				IDenseMatrix<Mat, T>& dst)
		{

			T* pd = dst.ptr_data();
			const index_t len = dst.nelems();

			for (index_t i = 0; i < len; ++i)
			{
				pd[i] = visitor.get_scalar(i);
			}
		}
	};


	template<int CTRows, int CTCols>
	struct percol_viseval_impl<CTRows, CTCols, scalar_kernel_t>
	{
#ifdef LMAT_USE_STATIC_ASSERT
		static_assert(CTRows > 0, "CTSize must be positive.");
#endif

		template<typename T, class Visitor, class Mat>
		LMAT_ENSURE_INLINE
		static void evaluate(
				IPerColMatrixScalarVisitor<Visitor, T>& visitor,
				IDenseMatrix<Mat, T>& dst)
		{
			const index_t ncols = dst.ncolumns();
			const index_t ldim = dst.lead_dim();
			T *pd = dst.ptr_data();

			for (index_t j = 0; j < ncols; ++j, pd += ldim)
			{
				typename matrix_visitor_state<Visitor>::type s = visitor.col_state(j);

				for (index_t i = 0; i < CTRows; ++i)
				{
					pd[i] = visitor.get_scalar(i, s);
				}
			}
		}
	};


	template<int CTCols>
	struct percol_viseval_impl<DynamicDim, CTCols, scalar_kernel_t>
	{
		template<typename T, class Visitor, class Mat>
		LMAT_ENSURE_INLINE
		static void evaluate(
				IPerColMatrixScalarVisitor<Visitor, T>& visitor,
				IDenseMatrix<Mat, T>& dst)
		{
			const index_t nrows = dst.nrows();
			const index_t ncols = dst.ncolumns();
			const index_t ldim = dst.lead_dim();
			T *pd = dst.ptr_data();

			for (index_t j = 0; j < ncols; ++j, pd += ldim)
			{
				typename matrix_visitor_state<Visitor>::type s = visitor.col_state(j);

				for (index_t i = 0; i < nrows; ++i)
				{
					pd[i] = visitor.get_scalar(i, s);
				}
			}
		}
	};


	/********************************************
	 *
	 *  Evaluation functions
	 *
	 ********************************************/

	template<typename T, class Expr, class Dst, typename VisCate, typename KerCate>
	LMAT_ENSURE_INLINE
	void evaluate(const IMatrixXpr<Expr, T>& src, IDenseMatrix<Dst, T>& dst,
			matrix_visit_policy<VisCate, KerCate>)
	{
		typedef matrix_visit_setting<VisCate, KerCate,
				binary_ct_rows<Expr, Dst>::value,
				binary_ct_cols<Expr, Dst>::value> setting_t;

		typedef typename matrix_viseval_impl_map<setting_t>::type impl_t;
		typename matrix_vismap<Expr, setting_t>::type visitor(src.derived());
		impl_t::evaluate(visitor, dst.derived());
	}

	template<typename T, class Expr, class Dst>
	LMAT_ENSURE_INLINE
	void linear_by_scalars_evaluate(const IMatrixXpr<Expr, T>& expr, IDenseMatrix<Dst, T>& dst)
	{
		evaluate(expr.derived(), dst.derived(), matrix_visit_policy<linear_vis, scalar_kernel_t>());
	}

	template<typename T, class Expr, class Dst>
	LMAT_ENSURE_INLINE
	void percol_by_scalars_evaluate(const IMatrixXpr<Expr, T>& expr, IDenseMatrix<Dst, T>& dst)
	{
		evaluate(expr.derived(), dst.derived(), matrix_visit_policy<percol_vis, scalar_kernel_t>());
	}

}

#endif
