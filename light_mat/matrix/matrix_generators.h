/*
 * @file matrix_generators.h
 *
 * Basic matrix generators
 *
 * @author Dahua Lin
 */

#ifndef _MSC_VER
#pragma once
#endif

#ifndef LIGHTMAT_MATRIX_GENERATORS_H_
#define LIGHTMAT_MATRIX_GENERATORS_H_

#include <light_mat/matrix/matrix_concepts.h>
#include <light_mat/core/mem_op.h>

namespace lmat
{

	template<typename T>
	class zero_gen : public IMatrixGenerator<zero_gen<T>, T>
	{
	public:
		LMAT_ENSURE_INLINE
		void generate_to(const index_t m, const index_t n,
				const index_t ldim, T *dst) const
		{
			if (n == 1 || ldim == m)
			{
				zero_mem(m * n, dst);
			}
			else
			{
				for (index_t j = 0; j < n; ++j, dst += ldim)
				{
					zero_mem(m, dst);
				}
			}
		}

	}; // end class zero_gen

	template<typename T>
	LMAT_ENSURE_INLINE
	inline zero_gen<T> zeros()
	{
		return zero_gen<T>();
	}


	template<typename T>
	class fill_gen : public IMatrixGenerator<fill_gen<T>, T>
	{
	public:
		LMAT_ENSURE_INLINE
		explicit fill_gen(const T& v) : m_val(v) { }

		LMAT_ENSURE_INLINE
		void generate_to(const index_t m, const index_t n,
				const index_t ldim, T *dst) const
		{
			if (n == 1 || ldim == m)
			{
				fill_mem(m * n, dst, m_val);
			}
			else
			{
				for (index_t j = 0; j < n; ++j, dst += ldim)
				{
					fill_mem(m, dst, m_val);
				}
			}
		}

	private:
		const T m_val;

	}; // end class fill_gen

	template<typename T>
	LMAT_ENSURE_INLINE
	inline fill_gen<T> fill_value(const T& v)
	{
		return fill_gen<T>(v);
	}


	template<typename T>
	class copy_gen : public IMatrixGenerator<copy_gen<T>, T>
	{
	public:
		LMAT_ENSURE_INLINE
		explicit copy_gen(const T *src) : m_src(src) { }

		LMAT_ENSURE_INLINE
		void generate_to(const index_t m, const index_t n,
				const index_t ldim, T *dst) const
		{
			if (n == 1 || m == ldim)
			{
				copy_mem(m * n, m_src, dst);
			}
			else if (m == 1)
			{
				unpack_vec(n, m_src, dst, ldim);
			}
			else
			{
				const T *s = m_src;
				for (index_t j = 0; j < n; ++j, s += m, dst += ldim)
				{
					copy_mem(m, s, dst);
				}
			}
		}

	private:
		const T *m_src;

	}; // end class copy_gen


	template<typename T>
	LMAT_ENSURE_INLINE
	inline copy_gen<T> copy_from(const T* src)
	{
		return copy_gen<T>(src);
	}


	// Convenient functions

	template<typename T, class Mat>
	LMAT_ENSURE_INLINE
	inline void zero(IDenseMatrix<Mat, T>& X)
	{
		zero_gen<T>().generate_to(X.nrows(), X.ncolumns(), X.lead_dim(), X.ptr_data());
	}

	template<typename T, class Mat>
	LMAT_ENSURE_INLINE
	inline void fill(IDenseMatrix<Mat, T>& X, const T& val)
	{
		fill_gen<T>(val).generate_to(X.nrows(), X.ncolumns(), X.lead_dim(), X.ptr_data());
	}

	template<typename T, class Mat>
	LMAT_ENSURE_INLINE
	inline void copy_to(const T *src, IDenseMatrix<Mat, T>& X)
	{
		copy_gen<T>(src).generate_to(X.nrows(), X.ncolumns(), X.lead_dim(), X.ptr_data());
	}


}


#endif 
