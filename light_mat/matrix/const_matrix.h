/*
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

#include <light_mat/matrix/matrix_base.h>

namespace lmat
{

	template<typename T, int CTRows, int CTCols>
	struct matrix_traits<const_matrix<T, CTRows, CTCols> >
	{
		static const int num_dimensions = 2;
		static const int compile_time_num_rows = CTRows;
		static const int compile_time_num_cols = CTCols;

		static const bool is_readonly = true;

		typedef T value_type;
	};

	template<typename T, int CTRows, int CTCols>
	struct is_linear_accessible<const_matrix<T, CTRows, CTCols> >
	{
		static const bool value = true;
	};

	template<typename T, int CTRows, int CTCols>
	class const_matrix : public IMatrixView<const_matrix<T, CTRows, CTCols>, T>
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

		LMAT_ENSURE_INLINE size_type size() const
		{
			return static_cast<size_type>(nelems());
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

		LMAT_ENSURE_INLINE T elem(const index_type, const index_type) const
		{
			return m_val;
		}

		LMAT_ENSURE_INLINE T operator[] (const index_type) const
		{
			return m_val;
		}

	private:
		matrix_shape<CTRows, CTCols> m_shape;
		const T m_val;
	};


	template<typename T, int CTRows, int CTCols, class DMat>
	LMAT_ENSURE_INLINE
	inline void evaluate_to(const const_matrix<T, CTRows, CTCols>& s, IDenseMatrix<DMat, T>& dst)
	{
		fill(dst, s.value());
	}


	template<class Mat>
	struct const_mat_same_form
	{
		typedef typename matrix_traits<Mat>::value_type value_type;

		typedef const_matrix<value_type,
				ct_rows<Mat>::value,
				ct_cols<Mat>::value> type;

		LMAT_ENSURE_INLINE
		static type get(const value_type& v)
		{
			return type(v);
		}
	};

}


#endif 
