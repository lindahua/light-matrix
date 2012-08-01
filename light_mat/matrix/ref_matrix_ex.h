/**
 * @file ref_matrix_ex.h
 *
 * ref_matrix_ex classes
 *
 * @author Dahua Lin
 */

#ifdef _MSC_VER
#pragma once
#endif

#ifndef LIGHTMAT_REF_MATRIX_EX_H_
#define LIGHTMAT_REF_MATRIX_EX_H_

#include <light_mat/matrix/matrix_base.h>

namespace lmat
{
	namespace detail
	{
		/********************************************
		 *
		 *  linear offset helper
		 *
		 ********************************************/

		template<class Mat, bool IsRow, bool IsCol>
		struct ref_matrix_ex_linear_offset_helper;

		template<class Mat>
		struct ref_matrix_ex_linear_offset_helper<Mat, false, false>
		{
			LMAT_ENSURE_INLINE static index_t get(const Mat& a, const index_t)
			{
				throw invalid_operation(
						"Accessing a (c)ref_matrix_ex object that is not a compile-time vector with linear index is not allowed.");
			}
		};

		template<class Mat>
		struct ref_matrix_ex_linear_offset_helper<Mat, false, true>
		{
			LMAT_ENSURE_INLINE static index_t get(const Mat& a, const index_t i)
			{
				return i;
			}
		};

		template<class Mat>
		struct ref_matrix_ex_linear_offset_helper<Mat, true, false>
		{
			LMAT_ENSURE_INLINE static index_t get(const Mat& a, const index_t i)
			{
				return i * a.lead_dim();
			}
		};

		template<class Mat>
		struct ref_matrix_ex_linear_offset_helper<Mat, true, true>
		{
			LMAT_ENSURE_INLINE static index_t get(const Mat& a, const index_t)
			{
				return 0;
			}
		};


		template<class Mat>
		index_t ref_ex_linear_offset(const Mat& a, const index_t i)
		{
			return ref_matrix_ex_linear_offset_helper<Mat,
					ct_is_row<Mat>::value,
					ct_is_col<Mat>::value>::get(a, i);
		}
	}

	/********************************************
	 *
	 *  cref_matrix_ex
	 *
	 ********************************************/

	template<typename T, int CTRows, int CTCols, typename Align>
	struct matrix_traits<cref_matrix_ex<T, CTRows, CTCols, Align> >
	{
		static const int num_dimensions = 2;
		static const int compile_time_num_rows = CTRows;
		static const int compile_time_num_cols = CTCols;

		static const bool is_readonly = true;

		typedef T value_type;
		typedef cpu_domain domain;
	};

	template<typename T, int CTRows, int CTCols, typename Align>
	struct ct_has_continuous_layout<cref_matrix_ex<T, CTRows, CTCols, Align> >
	{
		static const bool value = (CTCols == 1);
	};

	template<typename T, int CTRows, int CTCols, typename Align>
	struct is_base_aligned<cref_matrix_ex<T, CTRows, CTCols, Align> >
	{
		static const bool value = is_base_aligned_from_tag<Align>::value;
	};

	template<typename T, int CTRows, int CTCols, typename Align>
	struct is_percol_aligned<cref_matrix_ex<T, CTRows, CTCols, Align> >
	{
		static const bool value = is_percol_aligned_from_tag<Align>::value;
	};

	template<typename T, int CTRows, int CTCols, typename Align>
	struct is_linear_accessible<cref_matrix_ex<T, CTRows, CTCols, Align> >
	{
		static const bool value = (CTRows == 1 || CTCols == 1);
	};

	template<typename T, int CTRows, int CTCols, typename Align, class DMat>
	struct default_matrix_eval_policy<cref_matrix_ex<T, CTRows, CTCols, Align>, DMat>
	{
		typedef matrix_copy_policy type;
	};


	template<typename T, int CTRows, int CTCols, typename Align>
	class cref_matrix_ex : public IDenseMatrix<cref_matrix_ex<T, CTRows, CTCols, Align>, T>
	{
#ifdef LMAT_USE_STATIC_ASSERT
		static_assert(is_supported_matrix_value_type<T>::value,
				"T must be a supported matrix value type");
#endif

	public:
		LMAT_MAT_TRAITS_CDEFS(T)

	public:

		LMAT_ENSURE_INLINE
		cref_matrix_ex(const T* pdata, index_type m, index_type n, index_type ldim)
		: m_data(pdata), m_shape(m, n), m_ldim(ldim)
		{
		}

	private:
		cref_matrix_ex& operator = (const cref_matrix_ex& );  // no assignment

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
			return m_ldim;
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
			return m_shape.offset(m_ldim, i, j);
		}

		LMAT_ENSURE_INLINE const_reference elem(const index_type i, const index_type j) const
		{
			return m_data[offset(i, j)];
		}

		LMAT_ENSURE_INLINE const_reference operator[] (const index_type i) const
		{
			return m_data[detail::ref_ex_linear_offset(*this, i)];
		}

	private:
		const T *m_data;
		matrix_shape<CTRows, CTCols> m_shape;
		const index_type m_ldim;

	}; // end class cref_matrix_ex



	/********************************************
	 *
	 *  ref_matrix_ex
	 *
	 ********************************************/


	template<typename T, int CTRows, int CTCols, typename Align>
	struct matrix_traits<ref_matrix_ex<T, CTRows, CTCols, Align> >
	{
		static const int num_dimensions = 2;
		static const int compile_time_num_rows = CTRows;
		static const int compile_time_num_cols = CTCols;

		static const bool is_readonly = false;

		typedef T value_type;
		typedef cpu_domain domain;
	};

	template<typename T, int CTRows, int CTCols, typename Align>
	struct ct_has_continuous_layout<ref_matrix_ex<T, CTRows, CTCols, Align> >
	{
		static const bool value = (CTCols == 1);
	};

	template<typename T, int CTRows, int CTCols, typename Align>
	struct is_base_aligned<ref_matrix_ex<T, CTRows, CTCols, Align> >
	{
		static const bool value = is_base_aligned_from_tag<Align>::value;
	};

	template<typename T, int CTRows, int CTCols, typename Align>
	struct is_percol_aligned<ref_matrix_ex<T, CTRows, CTCols, Align> >
	{
		static const bool value = is_percol_aligned_from_tag<Align>::value;
	};

	template<typename T, int CTRows, int CTCols, typename Align>
	struct is_linear_accessible<ref_matrix_ex<T, CTRows, CTCols, Align> >
	{
		static const bool value = (CTRows == 1 || CTCols == 1);
	};

	template<typename T, int CTRows, int CTCols, typename Align, class DMat>
	struct default_matrix_eval_policy<ref_matrix_ex<T, CTRows, CTCols, Align>, DMat>
	{
		typedef matrix_copy_policy type;
	};


	template<typename T, int CTRows, int CTCols, typename Align>
	class ref_matrix_ex : public IDenseMatrix<ref_matrix_ex<T, CTRows, CTCols, Align>, T>
	{
#ifdef LMAT_USE_STATIC_ASSERT
		static_assert(is_supported_matrix_value_type<T>::value,
				"T must be a supported matrix value type");
#endif

	public:
		LMAT_MAT_TRAITS_DEFS(T)

	public:
		LMAT_ENSURE_INLINE
		ref_matrix_ex(T* pdata, index_type m, index_type n, index_type ldim)
		: m_data(pdata), m_shape(m, n), m_ldim(ldim)
		{
		}

	public:
		LMAT_ENSURE_INLINE ref_matrix_ex& operator = (const ref_matrix_ex& r)
		{
			if (this != &r)
			{
				copy(r, *this);
			}
			return *this;
		}

		template<class Expr>
		LMAT_ENSURE_INLINE ref_matrix_ex& operator = (const IMatrixXpr<Expr, T>& r)
		{
			default_assign(*this, r);
			return *this;
		}

		template<class Gen>
		LMAT_ENSURE_INLINE ref_matrix_ex& operator = (const IMatrixGenerator<Gen, T>& gen)
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
			return m_ldim;
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
			return m_shape.offset(m_ldim, i, j);
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
			return m_data[detail::ref_ex_linear_offset(*this, i)];
		}

		LMAT_ENSURE_INLINE reference operator[] (const index_type i)
		{
			return m_data[detail::ref_ex_linear_offset(*this, i)];
		}

		LMAT_ENSURE_INLINE void require_size(const index_type m, const index_type n)
		{
			check_arg(nrows() == m && ncolumns() == n,
					"ref_matrix::require_size: The required size is invalid.");
		}

	private:
		T *m_data;
		matrix_shape<CTRows, CTCols> m_shape;
		const index_type m_ldim;

	}; // end ref_matrix_ex

}

#endif /* REF_MATRIX_EX_H_ */
