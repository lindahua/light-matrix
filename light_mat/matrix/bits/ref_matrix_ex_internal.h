/**
 * @file ref_matrix_ex_internal.h
 *
 * Internal implementation of ref_matrix_ex
 *
 * @author Dahua Lin
 */

#ifdef _MSC_VER
#pragma once
#endif

#ifndef LIGHTMAT_REF_MATRIX_EX_INTERNAL_H_
#define LIGHTMAT_REF_MATRIX_EX_INTERNAL_H_

#include <light_mat/matrix/matrix_base.h>

namespace lmat { namespace detail {

	/********************************************
	 *
	 *  ref_matrix_ex_internal
	 *
	 ********************************************/

	template<class Mat, bool IsRow, bool IsCol>
	struct ref_matrix_ex_internal_linear_offset_helper;

	template<class Mat>
	struct ref_matrix_ex_internal_linear_offset_helper<Mat, false, false>
	{
		LMAT_ENSURE_INLINE static index_t get(const Mat& a, const index_t)
		{
			throw invalid_operation(
					"Accessing a (c)ref_matrix_ex object that is not a compile-time vector with linear index is not allowed.");
		}
	};

	template<class Mat>
	struct ref_matrix_ex_internal_linear_offset_helper<Mat, false, true>
	{
		LMAT_ENSURE_INLINE static index_t get(const Mat& a, const index_t i)
		{
			return i;
		}
	};

	template<class Mat>
	struct ref_matrix_ex_internal_linear_offset_helper<Mat, true, false>
	{
		LMAT_ENSURE_INLINE static index_t get(const Mat& a, const index_t i)
		{
			return i * a.lead_dim();
		}
	};

	template<class Mat>
	struct ref_matrix_ex_internal_linear_offset_helper<Mat, true, true>
	{
		LMAT_ENSURE_INLINE static index_t get(const Mat& a, const index_t)
		{
			return 0;
		}
	};


	template<class Mat>
	index_t ref_ex_linear_offset(const Mat& a, const index_t i)
	{
		return ref_matrix_ex_internal_linear_offset_helper<Mat,
				ct_is_row<Mat>::value,
				ct_is_col<Mat>::value>::get(a, i);
	}



	template<typename T, int CTRows, int CTCols>
	class ref_matrix_ex_internal
	{
	public:
		LMAT_ENSURE_INLINE
		ref_matrix_ex_internal(T *pdata, index_t m, index_t n, index_t ldim)
		: m_pdata(pdata), m_leaddim(ldim)
		{
			check_arg(m == CTRows && n == CTCols,
					"Attempted to construct a ref_matrix_ex with incorrect dimensions.");
		}

		LMAT_ENSURE_INLINE index_t nelems() const { return CTRows * CTCols; }

		LMAT_ENSURE_INLINE index_t nrows() const { return CTRows; }

		LMAT_ENSURE_INLINE index_t ncolumns() const { return CTCols; }

		LMAT_ENSURE_INLINE index_t lead_dim() const { return m_leaddim; }

		LMAT_ENSURE_INLINE T* ptr_data() const { return m_pdata; }

	private:
		T* m_pdata;
		index_t m_leaddim;
	};

	template<typename T, int CTCols>
	class ref_matrix_ex_internal<T, DynamicDim, CTCols>
	{
	public:
		LMAT_ENSURE_INLINE
		ref_matrix_ex_internal(T *pdata, index_t m, index_t n, index_t ldim)
		: m_pdata(pdata), m_nrows(m), m_leaddim(ldim)
		{
			check_arg(n == CTCols,
					"Attempted to construct a ref_matrix_ex with incorrect dimension (n != CTCols)");
		}

		LMAT_ENSURE_INLINE index_t nelems() const { return m_nrows * CTCols; }

		LMAT_ENSURE_INLINE index_t nrows() const { return m_nrows; }

		LMAT_ENSURE_INLINE index_t ncolumns() const { return CTCols; }

		LMAT_ENSURE_INLINE index_t lead_dim() const { return m_leaddim; }

		LMAT_ENSURE_INLINE T* ptr_data() const { return m_pdata; }

	private:
		T* m_pdata;
		index_t m_nrows;
		index_t m_leaddim;
	};


	template<typename T, int CTRows>
	class ref_matrix_ex_internal<T, CTRows, DynamicDim>
	{
	public:
		LMAT_ENSURE_INLINE
		ref_matrix_ex_internal(T *pdata, index_t m, index_t n, index_t ldim)
		: m_pdata(pdata), m_ncols(n), m_leaddim(ldim)
		{
			check_arg(m == CTRows,
					"Attempted to construct a ref_matrix_ex with incorrect dimension (m != CTRows)");
		}

		LMAT_ENSURE_INLINE index_t nelems() const { return CTRows * m_ncols; }

		LMAT_ENSURE_INLINE index_t nrows() const { return CTRows; }

		LMAT_ENSURE_INLINE index_t ncolumns() const { return m_ncols; }

		LMAT_ENSURE_INLINE index_t lead_dim() const { return m_leaddim; }

		LMAT_ENSURE_INLINE T* ptr_data() const { return m_pdata; }

	private:
		T* m_pdata;
		index_t m_ncols;
		index_t m_leaddim;
	};


	template<typename T>
	class ref_matrix_ex_internal<T, DynamicDim, DynamicDim>
	{
	public:
		LMAT_ENSURE_INLINE
		ref_matrix_ex_internal(T *pdata, index_t m, index_t n, index_t ldim)
		: m_pdata(pdata), m_nrows(m), m_ncols(n), m_leaddim(ldim)
		{
		}

		LMAT_ENSURE_INLINE index_t nelems() const { return m_nrows * m_ncols; }

		LMAT_ENSURE_INLINE index_t nrows() const { return m_nrows; }

		LMAT_ENSURE_INLINE index_t ncolumns() const { return m_ncols; }

		LMAT_ENSURE_INLINE index_t lead_dim() const { return m_leaddim; }

		LMAT_ENSURE_INLINE T* ptr_data() const { return m_pdata; }

	private:
		T* m_pdata;
		index_t m_nrows;
		index_t m_ncols;
		index_t m_leaddim;
	};

} }


#endif /* REF_MATRIX_EX_INTERNAL_H_ */
