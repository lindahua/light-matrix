/**
 * @file ref_matrix.h
 *
 * ref_matrix classes
 *
 * @author Dahua Lin
 */

#ifdef _MSC_VER
#pragma once
#endif

#ifndef LIGHTMAT_REF_MATRIX_H_
#define LIGHTMAT_REF_MATRIX_H_

#include <light_mat/matrix/matrix_base.h>

namespace lmat
{
	/********************************************
	 *
	 *  cref_matrix
	 *
	 ********************************************/

	template<typename T, int CTRows, int CTCols, typename Align>
	struct matrix_traits<cref_matrix<T, CTRows, CTCols, Align> >
	{
		static const int num_dimensions = 2;
		static const int compile_time_num_rows = CTRows;
		static const int compile_time_num_cols = CTCols;

		static const bool is_readonly = true;

		typedef T value_type;
		typedef cpu_domain domain;
	};

	template<typename T, int CTRows, int CTCols, typename Align>
	struct ct_has_continuous_layout<cref_matrix<T, CTRows, CTCols, Align> >
	{
		static const bool value = true;
	};

	template<typename T, int CTRows, int CTCols, typename Align>
	struct is_base_aligned<cref_matrix<T, CTRows, CTCols, Align> >
	{
		static const bool value = is_base_aligned_from_tag<Align>::value;
	};

	template<typename T, int CTRows, int CTCols, typename Align>
	struct is_percol_aligned<cref_matrix<T, CTRows, CTCols, Align> >
	{
		static const bool value = is_percol_aligned_from_tag<Align>::value;
	};


	template<typename T, int CTRows, int CTCols, typename Align>
	struct is_linear_accessible<cref_matrix<T, CTRows, CTCols, Align> >
	{
		static const bool value = true;
	};

	template<typename T, int CTRows, int CTCols, typename Align, class DMat>
	struct default_matrix_eval_policy<cref_matrix<T, CTRows, CTCols, Align>, DMat>
	{
		typedef matrix_copy_policy type;
	};


	template<typename T, int CTRows, int CTCols, typename Align>
	class cref_matrix : public IDenseMatrix<cref_matrix<T, CTRows, CTCols, Align>, T>
	{
#ifdef LMAT_USE_STATIC_ASSERT
		static_assert(is_supported_matrix_value_type<T>::value,
				"T must be a supported matrix value type");
#endif

	public:
		LMAT_MAT_TRAITS_CDEFS(T)

	private:
		static const int CTSize = CTRows * CTCols;
		static const bool IsDynamic = (CTSize == 0);

	public:

		LMAT_ENSURE_INLINE
		cref_matrix(const T* pdata, index_type m, index_type n)
		: m_data(pdata), m_shape(m, n)
		{
		}

	private:
		cref_matrix& operator = (const cref_matrix& );  // no assignment

	public:
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

		LMAT_ENSURE_INLINE index_type lead_dim() const
		{
			return m_shape.nrows();
		}

		LMAT_ENSURE_INLINE const_pointer ptr_data() const
		{
			return m_data;
		}

		LMAT_ENSURE_INLINE const_pointer ptr_col(const index_type j) const
		{
			return m_data + j * lead_dim();
		}

		LMAT_ENSURE_INLINE index_type offset(const index_type i, const index_type j) const
		{
			return m_shape.offset_c(i, j);
		}

		LMAT_ENSURE_INLINE const_reference elem(const index_type i, const index_type j) const
		{
			return m_data[offset(i, j)];
		}

		LMAT_ENSURE_INLINE const_reference operator[] (const index_type i) const
		{
			return m_data[i];
		}

	private:
		const T *m_data;
		matrix_shape<CTRows, CTCols> m_shape;

	}; // end class cref_matrix


	/********************************************
	 *
	 *  ref_matrix
	 *
	 ********************************************/


	template<typename T, int CTRows, int CTCols, typename Align>
	struct matrix_traits<ref_matrix<T, CTRows, CTCols, Align> >
	{
		static const int num_dimensions = 2;
		static const int compile_time_num_rows = CTRows;
		static const int compile_time_num_cols = CTCols;

		static const bool is_readonly = false;

		typedef T value_type;
		typedef cpu_domain domain;
	};

	template<typename T, int CTRows, int CTCols, typename Align>
	struct ct_has_continuous_layout<ref_matrix<T, CTRows, CTCols, Align> >
	{
		static const bool value = true;
	};

	template<typename T, int CTRows, int CTCols, typename Align>
	struct is_base_aligned<ref_matrix<T, CTRows, CTCols, Align> >
	{
		static const bool value = is_base_aligned_from_tag<Align>::value;
	};

	template<typename T, int CTRows, int CTCols, typename Align>
	struct is_percol_aligned<ref_matrix<T, CTRows, CTCols, Align> >
	{
		static const bool value = is_percol_aligned_from_tag<Align>::value;
	};

	template<typename T, int CTRows, int CTCols, typename Align>
	struct is_linear_accessible<ref_matrix<T, CTRows, CTCols, Align> >
	{
		static const bool value = true;
	};

	template<typename T, int CTRows, int CTCols, typename Align, class DMat>
	struct default_matrix_eval_policy<ref_matrix<T, CTRows, CTCols, Align>, DMat>
	{
		typedef matrix_copy_policy type;
	};


	template<typename T, int CTRows, int CTCols, typename Align>
	class ref_matrix : public IDenseMatrix<ref_matrix<T, CTRows, CTCols, Align>, T>
	{
#ifdef LMAT_USE_STATIC_ASSERT
		static_assert(is_supported_matrix_value_type<T>::value,
				"T must be a supported matrix value type");
#endif

	public:
		LMAT_MAT_TRAITS_DEFS(T)

	public:
		LMAT_ENSURE_INLINE
		ref_matrix(T* pdata, index_type m, index_type n)
		: m_data(pdata), m_shape(m, n)
		{
		}

	public:
		LMAT_ENSURE_INLINE ref_matrix& operator = (const ref_matrix& r)
		{
			if (this != &r)
			{
				copy(r, *this);
			}
			return *this;
		}

		template<class Expr>
		LMAT_ENSURE_INLINE ref_matrix& operator = (const IMatrixXpr<Expr, T>& r)
		{
			default_assign(*this, r);
			return *this;
		}

		template<class Gen>
		LMAT_ENSURE_INLINE ref_matrix& operator = (const IMatrixGenerator<Gen, T>& gen)
		{
			gen.generate_to(*this);
			return *this;
		}

	public:
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

		LMAT_ENSURE_INLINE index_type lead_dim() const
		{
			return m_shape.nrows();
		}

		LMAT_ENSURE_INLINE const_pointer ptr_data() const
		{
			return m_data;
		}

		LMAT_ENSURE_INLINE pointer ptr_data()
		{
			return m_data;
		}

		LMAT_ENSURE_INLINE const_pointer ptr_col(const index_type j) const
		{
			return m_data + j * lead_dim();
		}

		LMAT_ENSURE_INLINE pointer ptr_col(const index_type j)
		{
			return m_data + j * lead_dim();
		}

		LMAT_ENSURE_INLINE index_type offset(const index_type i, const index_type j) const
		{
			return m_shape.offset_c(i, j);
		}

		LMAT_ENSURE_INLINE const_reference elem(const index_type i, const index_type j) const
		{
			return m_data[offset(i, j)];
		}

		LMAT_ENSURE_INLINE reference elem(const index_type i, const index_type j)
		{
			return m_data[offset(i, j)];
		}

		LMAT_ENSURE_INLINE const_reference operator[] (const index_type i) const
		{
			return m_data[i];
		}

		LMAT_ENSURE_INLINE reference operator[] (const index_type i)
		{
			return m_data[i];
		}

		LMAT_ENSURE_INLINE void require_size(const index_type m, const index_type n)
		{
			check_arg(nrows() == m && ncolumns() == n,
					"ref_matrix::require_size: The required size is invalid.");
		}

	private:
		T *m_data;
		matrix_shape<CTRows, CTCols> m_shape;

	}; // end ref_matrix


	/********************************************
	 *
	 *  vectors derived from (c)ref_matrix
	 *
	 ********************************************/

	template<typename T, int CTRows, typename Align>
	class cref_col: public cref_matrix<T, CTRows, 1, Align>
	{
		typedef cref_matrix<T, CTRows, 1> base_mat_t;

	public:
		typedef index_t index_type;

		LMAT_ENSURE_INLINE
		cref_col(const T* pdata, index_type m)
		: base_mat_t(pdata, m, 1) { }
	};

	template<typename T, int CTRows, typename Align>
	class ref_col: public ref_matrix<T, CTRows, 1, Align>
	{
		typedef ref_matrix<T, CTRows, 1> base_mat_t;

	public:
		typedef index_t index_type;

		LMAT_ENSURE_INLINE
		ref_col(T* pdata, index_type m)
		: base_mat_t(pdata, m, 1) { }

		LMAT_ENSURE_INLINE ref_col& operator = (const base_mat_t& r)
		{
			base_mat_t::operator = (r);
			return *this;
		}

		template<class Expr>
		LMAT_ENSURE_INLINE ref_col& operator = (const IMatrixXpr<Expr, T>& r)
		{
			base_mat_t::operator = (r);
			return *this;
		}

		template<class Gen>
		LMAT_ENSURE_INLINE ref_col& operator = (const IMatrixGenerator<Gen, T>& gen)
		{
			base_mat_t::operator = (gen);
			return *this;
		}
	};


	template<typename T, int CTCols, typename Align>
	class cref_row: public cref_matrix<T, 1, CTCols, Align>
	{
		typedef cref_matrix<T, 1, CTCols> base_mat_t;

	public:
		typedef index_t index_type;

		LMAT_ENSURE_INLINE
		cref_row(const T* pdata, index_type n)
		: base_mat_t(pdata, 1, n) { }
	};

	template<typename T, int CTCols, typename Align>
	class ref_row: public ref_matrix<T, 1, CTCols, Align>
	{
		typedef ref_matrix<T, 1, CTCols> base_mat_t;

	public:
		typedef index_t index_type;

		LMAT_ENSURE_INLINE
		ref_row(T* pdata, index_type n)
		: base_mat_t(pdata, 1, n) { }

		LMAT_ENSURE_INLINE ref_row& operator = (const base_mat_t& r)
		{
			base_mat_t::operator = (r);
			return *this;
		}

		template<class Expr>
		LMAT_ENSURE_INLINE ref_row& operator = (const IMatrixXpr<Expr, T>& r)
		{
			base_mat_t::operator = (r);
			return *this;
		}

		template<class Gen>
		LMAT_ENSURE_INLINE ref_row& operator = (const IMatrixGenerator<Gen, T>& gen)
		{
			base_mat_t::operator = (gen);
			return *this;
		}
	};

}

#endif /* REF_MATRIX_H_ */
