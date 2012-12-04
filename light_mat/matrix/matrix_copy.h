/*
 * @file matrix_copy.h
 *
 * Functions to copy or cast matrices
 *
 * @author Dahua Lin
 */

#ifndef _MSC_VER
#pragma once
#endif

#ifndef LIGHTMAT_MATRIX_COPY_H_
#define LIGHTMAT_MATRIX_COPY_H_

#include "bits/matrix_copy_internal.h"

namespace lmat
{
	template<typename T, class RMat>
	LMAT_ENSURE_INLINE
	inline void copy(const T *ps, IDenseMatrix<RMat, T>& dst)
	{
		typedef typename internal::mat_copier_p_map<RMat>::type copier_t;
		copier_t::copy(ps, dst.derived());
	}

	template<typename T, class LMat>
	LMAT_ENSURE_INLINE
	inline void copy(const IDenseMatrix<LMat, T>& src, T* pd)
	{
		typedef typename internal::mat_copier_p_map<LMat>::type copier_t;
		copier_t::copy(src.derived(), pd);
	}


	template<typename T, class LMat, class RMat>
	LMAT_ENSURE_INLINE
	inline void copy(const IDenseMatrix<LMat, T>& src, IDenseMatrix<RMat, T>& dst)
	{
		typedef typename internal::mat_copier_map<LMat, RMat>::type copier_t;
		copier_t::copy(src.derived(), dst.derived());
	}


	template<typename T, class DMat>
	LMAT_ENSURE_INLINE
	inline DMat& operator << (IDenseMatrix<DMat, T>& dmat, const T *src)
	{
		copy(src, dmat.derived());
		return dmat.derived();
	}


	template<int M, int N>
	struct matrix_copy_scheme
	{
		matrix_shape<M, N> m_shape;

		matrix_copy_scheme(const index_t m, const index_t n)
		: m_shape(m, n) { }

		template<class SExpr, class DMat>
		LMAT_ENSURE_INLINE
		void evaluate(const SExpr& sexpr, DMat& dmat)
		{
			copy(sexpr, dmat);
		}
	};


	template<typename T, class SExpr, class DMat>
	LMAT_ENSURE_INLINE
	inline matrix_copy_scheme<
		meta::common_nrows< LMAT_TYPELIST_2(SExpr, DMat) >::value,
		meta::common_ncols< LMAT_TYPELIST_2(SExpr, DMat) >::value>
	get_default_eval_scheme(
			const IDenseMatrix<SExpr, T>& sexpr,
			IDenseMatrix<DMat, T>& dmat)
	{
		const int M = meta::common_nrows< LMAT_TYPELIST_2(SExpr, DMat) >::value;
		const int N = meta::common_ncols< LMAT_TYPELIST_2(SExpr, DMat) >::value;

		return matrix_copy_scheme<M, N>(
				common_nrows(sexpr, dmat),
				common_ncols(sexpr, dmat));
	}

	template<typename T, class SExpr, class DMat>
	LMAT_ENSURE_INLINE
	inline typename meta::enable_if< meta::is_mat_assignable<SExpr, DMat>, void>::type
	default_evaluate(const IMatrixXpr<SExpr, T>& sexpr, IDenseMatrix<DMat, T>& dmat)
	{
		const SExpr& s = sexpr.derived();
		DMat& t = dmat.derived();

		get_default_eval_scheme(s, t).evaluate(s, t);
	}

}

#endif 
