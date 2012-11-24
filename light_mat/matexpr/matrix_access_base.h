/**
 * @file matrix_access_base.h
 *
 * @brief
 *
 * @author Dahua Lin
 */

#ifdef _MSC_VER
#pragma once
#endif

#ifndef LIGHTMAT_MATRIX_ACCESS_BASE_H_
#define LIGHTMAT_MATRIX_ACCESS_BASE_H_

#include <light_mat/matrix/matrix_meta.h>

namespace lmat
{
	/********************************************
	 *
	 *  Types for Matrix Access
	 *
	 ********************************************/

	// kernel categories

	struct scalar_kernel_t { };
	struct simd_kernel_t { };

	// matrix access categories

	struct linear_macc { };
	struct percol_macc { };

	// matrix access setting

	template<typename AccCate, typename KerCate, int ME, int NE>
	struct macc_setting
	{
		typedef AccCate access_category;
		typedef KerCate kernel_category;

		 static const int ct_rows = ME;
		 static const int ct_cols = NE;
	};

	template<class SExpr, class AccCate, class KerCate>
	struct macc_accessor_map;

	template<class SExpr, class DMat>
	struct default_macc_setting;


	/********************************************
	 *
	 *  Accessor Interfaces
	 *
	 ********************************************/

	template<class Derived, typename T>
	class ILinearMatrixScalarAccessor
	{
	public:
		LMAT_CRTP_REF

		LMAT_ENSURE_INLINE
		T get_scalar(const index_t i) const
		{
			return derived().get_scalar(i);
		}
	};


	template<class Visitor>
	struct percol_macc_state_map { };

	template<class Derived, typename T>
	class IPerColMatrixScalarAccessor
	{
	public:
		typedef typename percol_macc_state_map<Derived>::type col_state_t;

		LMAT_CRTP_REF

		LMAT_ENSURE_INLINE
		T get_scalar(const index_t i, const col_state_t& s) const
		{
			return derived().get_scalar(i, s);
		}

		LMAT_ENSURE_INLINE
		col_state_t col_state(const index_t j) const
		{
			return derived().col_state(j);
		}
	};


}

#endif /* MATRIX_ACCESS_BASE_H_ */





