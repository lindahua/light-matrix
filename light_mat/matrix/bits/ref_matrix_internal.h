/**
 * @file ref_matrix_internal.h
 *
 * Internal implementation of ref_matrix classes
 *
 * @author Dahua Lin
 */

#ifdef _MSC_VER
#pragma once
#endif

#ifndef LIGHTMAT_REF_MATRIX_INTERNAL_H_
#define LIGHTMAT_REF_MATRIX_INTERNAL_H_

#include <light_mat/matrix/matrix_base.h>

namespace lmat { namespace detail {

	/********************************************
	 *
	 *  ref_matrix_internal
	 *
	 ********************************************/

	template<typename T, int CTRows, int CTCols>
	class ref_matrix_internal
	{
	public:
		LMAT_ENSURE_INLINE
		ref_matrix_internal(T *pdata, index_t m, index_t n)
		: m_pdata(pdata)
		{
			check_arg(m == CTRows && n == CTCols,
					"Attempted to construct a static ref_matrix with incorrect dimensions.");
		}

		LMAT_ENSURE_INLINE index_t nelems() const { return CTRows * CTCols; }

		LMAT_ENSURE_INLINE index_t nrows() const { return CTRows; }

		LMAT_ENSURE_INLINE index_t ncolumns() const { return CTCols; }

		LMAT_ENSURE_INLINE index_t lead_dim() const { return CTRows; }

		LMAT_ENSURE_INLINE T* ptr_data() const { return m_pdata; }

	private:
		T* m_pdata;
	};

	template<typename T, int CTCols>
	class ref_matrix_internal<T, DynamicDim, CTCols>
	{
	public:
		LMAT_ENSURE_INLINE
		ref_matrix_internal(T *pdata, index_t m, index_t n)
		: m_pdata(pdata), m_nrows(m)
		{
			check_arg(n == CTCols,
					"Attempted to construct a ref_matrix with incorrect dimension (n != CTCols).");
		}

		LMAT_ENSURE_INLINE index_t nelems() const { return m_nrows * CTCols; }

		LMAT_ENSURE_INLINE index_t nrows() const { return m_nrows; }

		LMAT_ENSURE_INLINE index_t ncolumns() const { return CTCols; }

		LMAT_ENSURE_INLINE index_t lead_dim() const { return m_nrows; }

		LMAT_ENSURE_INLINE T* ptr_data() const { return m_pdata; }

	private:
		T* m_pdata;
		index_t m_nrows;
	};


	template<typename T, int CTRows>
	class ref_matrix_internal<T, CTRows, DynamicDim>
	{
	public:
		LMAT_ENSURE_INLINE
		ref_matrix_internal(T *pdata, index_t m, index_t n)
		: m_pdata(pdata), m_ncols(n)
		{
			check_arg(m == CTRows,
					"Attempted to construct a ref_matrix with incorrect dimension (m != CTRows).");
		}

		LMAT_ENSURE_INLINE index_t nelems() const { return CTRows * m_ncols; }

		LMAT_ENSURE_INLINE index_t nrows() const { return CTRows; }

		LMAT_ENSURE_INLINE index_t ncolumns() const { return m_ncols; }

		LMAT_ENSURE_INLINE index_t lead_dim() const { return CTRows; }

		LMAT_ENSURE_INLINE T* ptr_data() const { return m_pdata; }

	private:
		T* m_pdata;
		index_t m_ncols;
	};


	template<typename T>
	class ref_matrix_internal<T, DynamicDim, DynamicDim>
	{
	public:
		LMAT_ENSURE_INLINE
		ref_matrix_internal(T *pdata, index_t m, index_t n)
		: m_pdata(pdata), m_nrows(m), m_ncols(n)
		{
		}

		LMAT_ENSURE_INLINE index_t nelems() const { return m_nrows * m_ncols; }

		LMAT_ENSURE_INLINE index_t nrows() const { return m_nrows; }

		LMAT_ENSURE_INLINE index_t ncolumns() const { return m_ncols; }

		LMAT_ENSURE_INLINE index_t lead_dim() const { return m_nrows; }

		LMAT_ENSURE_INLINE T* ptr_data() const { return m_pdata; }

	private:
		T* m_pdata;
		index_t m_nrows;
		index_t m_ncols;
	};




} }

#endif /* REF_MATRIX_INTERNAL_H_ */
