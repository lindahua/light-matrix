/**
 * @file matrix_access_eval.h
 *
 * @brief Evaluation based on matrix access
 *
 * @author Dahua Lin
 */

#ifndef LIGHTMAT_MATRIX_ACCESS_EVAL_H_
#define LIGHTMAT_MATRIX_ACCESS_EVAL_H_

#include "bits/macc_eval_impl.h"

namespace lmat
{
	/********************************************
	 *
	 *  macc schemes
	 *
	 ********************************************/

	template<class Xpr, typename AccCate, typename KerCate>
	struct macc_cost;

	const int MACC_CACHE_COST = 1200;
	const int MACC_SHORT_PERCOL_COST = 100;
	const int MACC_SHORTCOL_UBOUND = 4;

	// matrix access categories

	struct any_macc { };
	struct linear_macc { };
	struct percol_macc { };

	// matrix access setting

	template<typename AccCate, typename KerCate, int M, int N>
	struct macc_scheme;

	template<int M, int N>
	struct macc_scheme<any_macc, scalar_ker, M, N>
	{
		const matrix_shape<M, N> shape;
		const bool _use_linear;

		LMAT_ENSURE_INLINE
		macc_scheme(index_t m, index_t n, bool lin)
		: shape(m, n), _use_linear(lin) { }

		LMAT_ENSURE_INLINE bool use_linear() const { return _use_linear; }
		LMAT_ENSURE_INLINE bool use_simd() const { return false; }

		template<class SExpr, class DMat>
		LMAT_ENSURE_INLINE
		void evaluate(const SExpr& sexpr, DMat& dmat)
		{
			if (_use_linear)
			{
				internal::macc_eval_impl<linear_macc, scalar_ker>::evaluate(
						shape.nelems(), sexpr, dmat);
			}
			else
			{
				internal::macc_eval_impl<percol_macc, scalar_ker>::evaluate(
						shape.nrows(), shape.ncolumns(), sexpr, dmat);
			}
		}

		template<class SExpr, class DMat>
		LMAT_ENSURE_INLINE
		static macc_scheme get_default(const SExpr& expr, const DMat& dmat)
		{
			const index_t m = dmat.nrows();

			int macc_linear_cost = macc_cost<SExpr, linear_macc, scalar_ker>::value;
			int macc_percol_cost = macc_cost<SExpr, percol_macc, scalar_ker>::value;

			if (m <= MACC_SHORT_PERCOL_COST)
				macc_percol_cost += MACC_SHORT_PERCOL_COST;

			return macc_scheme(m, dmat.ncolumns(),
					macc_linear_cost <= macc_percol_cost);
		}
	};

	template<int M, int N>
	struct macc_scheme<linear_macc, scalar_ker, M, N>
	{
		const matrix_shape<M, N> shape;

		LMAT_ENSURE_INLINE
		macc_scheme(index_t m, index_t n)
		: shape(m, n) { }

		LMAT_ENSURE_INLINE bool use_linear() const { return true; }
		LMAT_ENSURE_INLINE bool use_simd() const { return false; }

		template<class SExpr, class DMat>
		LMAT_ENSURE_INLINE
		void evaluate(const SExpr& sexpr, DMat& dmat)
		{
			internal::macc_eval_impl<linear_macc, scalar_ker>::evaluate(
					shape.nelems(), sexpr, dmat);
		}

		template<class SExpr, class DMat>
		LMAT_ENSURE_INLINE
		static macc_scheme get_default(const SExpr& expr, const DMat& dmat)
		{
			return macc_scheme(dmat.nrows(), dmat.ncolumns());
		}
	};


	template<int M, int N>
	struct macc_scheme<percol_macc, scalar_ker, M, N>
	{
		const matrix_shape<M, N> shape;

		LMAT_ENSURE_INLINE
		macc_scheme(index_t m, index_t n)
		: shape(m, n) { }

		LMAT_ENSURE_INLINE bool use_linear() const { return false; }
		LMAT_ENSURE_INLINE bool use_simd() const { return false; }

		template<class SExpr, class DMat>
		LMAT_ENSURE_INLINE
		void evaluate(const SExpr& sexpr, DMat& dmat)
		{
			internal::macc_eval_impl<percol_macc, scalar_ker>::evaluate(
					shape.nrows(), shape.ncolumns(), sexpr, dmat);
		}

		template<class SExpr, class DMat>
		LMAT_ENSURE_INLINE
		static macc_scheme get_default(const SExpr& expr, const DMat& dmat)
		{
			return macc_scheme(dmat.nrows(), dmat.ncolumns());
		}
	};


	template<class SExpr, class DMat>
	struct default_macc_scheme
	{
		typedef typename
				meta::if_<meta::supports_linear_index<DMat>,
				any_macc,
				percol_macc>::type access_category;

		typedef scalar_ker kernel_category;

		static const int M = meta::common_nrows< LMAT_TYPELIST_2(SExpr, DMat) >::value;
		static const int N = meta::common_ncols< LMAT_TYPELIST_2(SExpr, DMat) >::value;

		typedef macc_scheme<access_category, kernel_category, M, N> type;

		static type get(const SExpr& sexpr, const DMat& dmat)
		{
			return type::get_default(sexpr, dmat);
		}
	};
}

#endif /* MATRIX_ACCESS_EVAL_H_ */
