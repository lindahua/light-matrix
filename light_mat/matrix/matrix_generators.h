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

#include <light_mat/matrix/matrix_copy.h>
#include <light_mat/matrix/matrix_fill.h>

namespace lmat
{
	template<class Derived, typename T>
	class IMatrixGenerator
	{
	public:
		LMAT_CRTP_REF

		template<class Mat>
		LMAT_ENSURE_INLINE void generate_to(IDenseMatrix<Mat, T>& mat) const
		{
			derived().generate_to(mat);
		}
	};


	template<typename T>
	class zero_gen : public IMatrixGenerator<zero_gen<T>, T>
	{
	public:
		template<class Mat>
		LMAT_ENSURE_INLINE
		void generate_to(IDenseMatrix<Mat, T>& mat) const
		{
			zero(mat);
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

		template<class Mat>
		LMAT_ENSURE_INLINE
		void generate_to(IDenseMatrix<Mat, T>& mat) const
		{
			fill(mat, m_val);
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

		template<class Mat>
		LMAT_ENSURE_INLINE
		void generate_to(IDenseMatrix<Mat, T>& mat) const
		{
			copy(m_src, mat.derived());
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

}


#endif 
