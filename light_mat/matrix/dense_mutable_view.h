/**
 * @file dense_mutable_view.h
 *
 * The dense_mutable_view classes
 *
 * @author Dahua Lin
 */

#ifdef _MSC_VER
#pragma once
#endif

#ifndef LIGHTMAT_DENSE_MUTABLE_VIEW_H_
#define LIGHTMAT_DENSE_MUTABLE_VIEW_H_

#include <light_mat/matrix/matrix_concepts.h>

namespace lmat
{

	template<class Mat>
	class dense_mutable_view : public Mat
	{
		static_assert(meta::is_regular_mat<Mat>::value, "Mat should be a dense matrix class");

	public:
		typedef typename meta::value_type_of<Mat>::type value_type;

		LMAT_ENSURE_INLINE
		dense_mutable_view(const Mat& base_mat) : Mat(base_mat)
		{
		}

	public:
		template<class Expr>
		LMAT_ENSURE_INLINE const dense_mutable_view& operator = (const IMatrixXpr<Expr, value_type>& r) const
		{
			(const_cast<dense_mutable_view*>(this))->assign(r);
			return *this;
		}

	private:
		template<class Expr>
		LMAT_ENSURE_INLINE void assign(const IMatrixXpr<Expr, value_type>& r)
		{
			Mat::operator = (r);
		}


	}; // end dense_wref_mat
}

#endif /* DENSE_WREF_MAT_H_ */
