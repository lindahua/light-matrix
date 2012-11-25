/**
 * @file matrix_acc_eval.h
 *
 * Generic accessor-based matrix evaluation
 *
 * @author Dahua Lin
 */

#ifdef _MSC_VER
#pragma once
#endif

#ifndef LIGHTMAT_MATRIX_ACC_EVAL_H_
#define LIGHTMAT_MATRIX_ACC_EVAL_H_

#include <light_mat/matexpr/dense_accessors.h>

namespace lmat
{

	template<class SExpr, class DMat, class Setting> struct macc_eval_scheme;

	template<class SExpr, class DMat>
	struct default_macc_setting
	{
		static const int M = binary_ct_rows<SExpr, DMat>::value;
		static const int N = binary_ct_cols<SExpr, DMat>::value;

		static const int linacc_cost = macc_cost<SExpr, M, N, linear_macc, scalar_kernel_t>::value;
		static const int percol_cost = macc_cost<SExpr, M, N, percol_macc, scalar_kernel_t>::value;

		static const int use_linear =
				ct_supports_linear_index<DMat>::value &&
				(linacc_cost <= percol_cost);

		typedef scalar_kernel_t ker_category;
		typedef typename if_c<use_linear, linear_macc, percol_macc>::type acc_category;

		typedef macc_setting<acc_category, ker_category, M, N> type;
	};

	template<class SExpr, class DMat>
	struct default_macc_eval_scheme
	{
		typedef typename default_macc_setting<SExpr, DMat>::type setting_t;
		typedef macc_eval_scheme<SExpr, DMat, setting_t> type;
	};


	/********************************************
	 *
	 *  scheme implementation
	 *
	 ********************************************/

	template<class SExpr, class DMat, int M, int N>
	struct macc_eval_scheme<SExpr, DMat, macc_setting<linear_macc, scalar_kernel_t, M, N> >
	{
#ifdef LMAT_USE_STATIC_ASSERT
		static_assert(ct_supports_linear_index<DMat>::value, "DMat must support linear indexing.");
#endif

		typedef typename
				macc_accessor_map<SExpr, linear_macc, scalar_kernel_t>::type
				accessor_t;

		const SExpr& sexpr;
		DMat& dmat;

		LMAT_ENSURE_INLINE
		macc_eval_scheme(const SExpr& s, DMat& t)
		: sexpr(s), dmat(t)
		{ }

		void evaluate() const
		{
			accessor_t a(sexpr);
			const index_t ne = dmat.nelems();

			for (index_t i = 0; i < ne; ++i)
			{
				dmat[i] = a.get_scalar(i);
			}
		}
	};


	template<class SExpr, class DMat, int M, int N>
	struct macc_eval_scheme<SExpr, DMat, macc_setting<percol_macc, scalar_kernel_t, M, N> >
	{
		typedef typename binary_value_type<SExpr, DMat>::type T;

		typedef typename
				macc_accessor_map<SExpr, percol_macc, scalar_kernel_t>::type
				accessor_t;

		typedef typename percol_macc_state_map<accessor_t>::type col_state_t;

		const SExpr& sexpr;
		DMat& dmat;

		LMAT_ENSURE_INLINE
		macc_eval_scheme(const SExpr& s, DMat& t)
		: sexpr(s), dmat(t)
		{ }

		void evaluate() const
		{
			accessor_t a(sexpr);
			const index_t m = dmat.nrows();
			const index_t n = dmat.ncolumns();

			if (m == 1)
			{
				T *pd = dmat.ptr_data();
				const index_t cs = dmat.col_stride();

				if (cs == 1)
				{
					for (index_t j = 0; j < n; ++j)
						pd[j] = a.get_scalar(0, a.col_state(j));
				}
				else
				{
					for (index_t j = 0; j < n; ++j)
						pd[j * cs] = a.get_scalar(0, a.col_state(j));
				}
			}
			else
			{
				const index_t rs = dmat.row_stride();

				if (rs == 1)
				{
					for (index_t j = 0; j < n; ++j)
					{
						T *pd = dmat.ptr_col(j);
						col_state_t s = a.col_state(j);

						for (index_t i = 0; i < m; ++i)
							pd[i] = a.get_scalar(i, s);
					}
				}
				else
				{
					for (index_t j = 0; j < n; ++j)
					{
						T *pd = dmat.ptr_col(j);
						col_state_t s = a.col_state(j);

						for (index_t i = 0; i < m; ++i)
							pd[i * rs] = a.get_scalar(i, s);
					}
				}
			}
		}
	};

}

#endif


