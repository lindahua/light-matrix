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
	using internal::matrix_copy_scheme;

	template<typename T, class RMat>
	LMAT_ENSURE_INLINE
	inline void copy(const T *ps, IDenseMatrix<RMat, T>& dst)
	{
		matrix_copy_scheme<
			meta::nrows<RMat>::value,
			meta::ncols<RMat>::value,
			typename meta::continuous_level<RMat>::type>
		copy_sch(dst.nrows(), dst.ncolumns());

		internal::copy(ps, dst.derived(), copy_sch);
	}

	template<typename T, class LMat>
	LMAT_ENSURE_INLINE
	inline void copy(const IDenseMatrix<LMat, T>& src, T* pd)
	{
		matrix_copy_scheme<
			meta::nrows<LMat>::value,
			meta::ncols<LMat>::value,
			typename meta::continuous_level<LMat>::type>
		copy_sch(src.nrows(), src.ncolumns());

		internal::copy(src.derived(), pd, copy_sch);
	}


	template<class SMat, class DMat>
	struct default_copy_scheme
	{
		typedef typename meta::binary_continuous_level<SMat, DMat>::type cont_level;

		typedef matrix_copy_scheme<
			meta::common_nrows< LMAT_TYPELIST_2(SMat, DMat) >::value,
			meta::common_ncols< LMAT_TYPELIST_2(SMat, DMat) >::value,
			cont_level> type;

		LMAT_ENSURE_INLINE
		static type get(const SMat& smat, const DMat& dmat)
		{
			return type(common_nrows(smat, dmat), common_ncols(smat, dmat));
		}
	};

	template<typename T, class LMat, class RMat>
	LMAT_ENSURE_INLINE
	inline void copy(const IDenseMatrix<LMat, T>& src, IDenseMatrix<RMat, T>& dst)
	{
		const LMat& s = src.derived();
		RMat& d = dst.derived();

		internal::copy(s, d, default_copy_scheme<LMat, RMat>::get(s, d));
	}


	template<typename T, class DMat>
	LMAT_ENSURE_INLINE
	inline DMat& operator << (IDenseMatrix<DMat, T>& dmat, const T *src)
	{
		copy(src, dmat.derived());
		return dmat.derived();
	}

	template<typename T, class SExpr, class DMat>
	LMAT_ENSURE_INLINE
	inline typename meta::enable_if< meta::is_mat_assignable<SExpr, DMat>, void>::type
	evaluate(const IDenseMatrix<SExpr, T>& sexpr, IDenseMatrix<DMat, T>& dmat)
	{
		copy(sexpr.derived(), dmat.derived());
	}

}

#endif 
