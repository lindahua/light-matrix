/*
 * @file vector_eval_kernel.h
 *
 * Vector evaluation kernels
 *
 * @author Dahua Lin
 */

#ifndef _MSC_VER
#pragma once
#endif

#ifndef LIGHTMAT_VECTOR_EVAL_KERNEL_H_
#define LIGHTMAT_VECTOR_EVAL_KERNEL_H_

#include <light_mat/core/basic_defs.h>

namespace lmat
{
	// Kernel types

	struct scalar_kernel_t { };
	struct simd_kernel_t { };

	/********************************************
	 *
	 *  Interfaces
	 *
	 ********************************************/

	template<class Ker> struct vector_eval_kernel_state;  // kernel --> state
	struct veval_nil_state { };

	template<class Derived, typename T>
	class IVecEvalScalarKernel
	{
	public:
		LMAT_CRTP_REF

		LMAT_ENSURE_INLINE
		T get_value(const index_t i,
				const typename vector_eval_kernel_state<Derived>::type& s) const
		{
			return derived().get_value(i, s);
		}
	};


	/********************************************
	 *
	 *  Scalar kernels
	 *
	 ********************************************/

	template<typename T> struct veval_memacc_kernel;
	template<typename T> struct veval_const_kernel;

	template<typename T>
	struct vector_eval_kernel_state<veval_memacc_kernel<T> >
	{
		typedef const T* type;
	};

	template<typename T>
	struct vector_eval_kernel_state<veval_const_kernel<T> >
	{
		typedef T type;
	};

	template<typename T>
	struct veval_memacc_kernel
	{
	public:
		LMAT_ENSURE_INLINE
		T get_value(const index_t i, const T* s) const
		{
			return s[i];
		}
	};

	template<typename T>
	struct veval_const_kernel
	{
	public:
		LMAT_ENSURE_INLINE
		T get_value(const index_t, const T& s) const
		{
			return s;
		}
	};
}

#endif 
