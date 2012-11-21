/**
 * @file const_matrix.h
 *
 * Constant matrix class
 *
 * @author Dahua Lin
 */

#ifndef _MSC_VER
#pragma once
#endif

#ifndef LIGHTMAT_CONST_MATRIX_H_
#define LIGHTMAT_CONST_MATRIX_H_

#include <light_mat/matrix/matrix_fill.h>
#include <light_mat/matrix/matrix_shape.h>

namespace lmat
{

	template<typename T, int CTRows, int CTCols>
	struct matrix_traits<const_matrix<T, CTRows, CTCols> >
	{
		static const int num_dimensions = 2;
		static const int ct_num_rows = CTRows;
		static const int ct_num_cols = CTCols;

		static const bool is_readonly = true;

		typedef T value_type;
		typedef cpu_domain domain;
	};


	template<typename T, int CTRows, int CTCols>
	class const_matrix : public IMatrixXpr<const_matrix<T, CTRows, CTCols>, T>
	{
#ifdef LMAT_USE_STATIC_ASSERT
		static_assert(is_supported_matrix_value_type<T>::value,
				"T must be a supported matrix value type");
#endif

	public:
		LMAT_MAT_TRAITS_DEFS(T)

		LMAT_ENSURE_INLINE const_matrix(const index_type m, const index_type n, const T& val)
		: m_shape(m, n), m_val(val) { }

		LMAT_ENSURE_INLINE index_type nelems() const
		{
			return m_shape.nelems();
		}

		LMAT_ENSURE_INLINE index_type nrows() const
		{
			return m_shape.nrows();
		}

		LMAT_ENSURE_INLINE index_type ncolumns() const
		{
			return m_shape.ncolumns();
		}

		LMAT_ENSURE_INLINE T value() const
		{
			return m_val;
		}

	private:
		matrix_shape<CTRows, CTCols> m_shape;
		const T m_val;
	};


	// Evaluation

	template<int M0, int N0, class DMat>
	struct const_matrix_eval_scheme
	{
		typedef typename matrix_traits<DMat>::value_type T;

		T val;
		DMat& dmat;

		LMAT_ENSURE_INLINE
		const_matrix_eval_scheme(const const_matrix<T, M0, N0>& sm, DMat& dm)
		: val(sm.value()), dmat(dm) { }

		LMAT_ENSURE_INLINE
		void evaluate()
		{
			fill(dmat, val);
		}
	};


	template<typename T, int CTRows, int CTCols, class DMat>
	struct default_matrix_eval_scheme<const_matrix<T, CTRows, CTCols>, DMat>
	{
		typedef const_matrix_eval_scheme<CTRows, CTCols, DMat> type;
	};


	// Transpose

	template<typename T, int M, int N, typename Arg_HP>
	struct unary_expr_map<transpose_t, Arg_HP, const_matrix<T, M, N> >
	{
		typedef const_matrix<T, N, M> type;

		LMAT_ENSURE_INLINE
		static type get(transpose_t, const arg_forwarder<Arg_HP, const_matrix<T, M, N> >& arg_fwd)
		{
			const const_matrix<T, M, N>& arg = arg_fwd.arg;
			return type(arg.ncolumns(), arg.nrows(), arg.value());
		}
	};

}


#endif 
