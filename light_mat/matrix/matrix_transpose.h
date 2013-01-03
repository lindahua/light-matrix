/**
 * @file matrix_transpose.h
 *
 * @brief Matrix transposition
 *
 * @author Dahua Lin
 */

#ifdef _MSC_VER
#pragma once
#endif

#ifndef LIGHTMAT_MATRIX_TRANSPOSE_H_
#define LIGHTMAT_MATRIX_TRANSPOSE_H_

#include "internal/matrix_transpose_internal.h"

namespace lmat
{
	// forward declaration

	template<class Arg> class transpose_expr;

	// direct transpose

	template<typename T, class SMat, class DMat>
	LMAT_ENSURE_INLINE
	inline void transpose(const IRegularMatrix<SMat, T>& smat, IRegularMatrix<DMat, T>& dmat)
	{
		index_t m = smat.nrows();
		index_t n = smat.ncolumns();

		LMAT_CHECK_DIMS( dmat.nrows() == n && dmat.ncolumns() == m );

		internal::direct_transpose(m, n, smat, dmat);
	}


	/******************************************************
	 *
	 *  transpose expression class
	 *
	 ******************************************************/

	template<class Arg>
	struct matrix_traits<transpose_expr<Arg> >
	{
		static const int num_dimensions = 2;
		static const int ct_num_rows = meta::ncols<Arg>::value;
		static const int ct_num_cols = meta::nrows<Arg>::value;

		static const bool is_readonly = true;

		typedef typename matrix_traits<Arg>::value_type value_type;
		typedef typename matrix_traits<Arg>::domain domain;
	};


	template<class Arg>
	class transpose_expr
	: public IMatrixXpr<transpose_expr<Arg>, typename matrix_traits<Arg>::value_type>
	{
		static_assert(meta::is_regular_mat<Arg>::value, "Arg must be a regular matrix.");
		typedef matrix_shape<meta::ncols<Arg>::value, meta::nrows<Arg>::value> shape_type;

	public:
		LMAT_ENSURE_INLINE
		explicit transpose_expr(const Arg& arg) : m_arg(arg) { }

		LMAT_ENSURE_INLINE
		const Arg& arg() const
		{
			return m_arg;
		}

		LMAT_ENSURE_INLINE index_t nelems() const
		{
			return m_arg.nelems();
		}

		LMAT_ENSURE_INLINE index_t nrows() const
		{
			return m_arg.ncolumns();
		}

		LMAT_ENSURE_INLINE index_t ncolumns() const
		{
			return m_arg.nrows();
		}

		LMAT_ENSURE_INLINE shape_type shape() const
		{
			return shape_type(nrows(), ncolumns());
		}

	private:
		const Arg& m_arg;
	};


	template<typename T, class Arg>
	LMAT_ENSURE_INLINE
	inline transpose_expr<Arg> transpose(const IRegularMatrix<Arg, T>& arg)
	{
		return transpose_expr<Arg>(arg.derived());
	}


	/******************************************************
	 *
	 *  evaluation
	 *
	 ******************************************************/

	template<class Arg, class DMat>
	LMAT_ENSURE_INLINE
	inline void evaluate(const transpose_expr<Arg>& expr,
			IRegularMatrix<DMat, typename matrix_traits<Arg>::value_type>& dmat)
	{
		transpose(expr.arg(), dmat);
	}


}

#endif
