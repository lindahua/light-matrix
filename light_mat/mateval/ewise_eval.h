/**
 * @file ewise_eval.h
 *
 * @brief Element-wise evaluation
 *
 * @author Dahua Lin
 */

#ifdef _MSC_VER
#pragma once
#endif

#ifndef LIGHTMAT_EWISE_EVAL_H_
#define LIGHTMAT_EWISE_EVAL_H_

#include <light_mat/mateval/common_kernels.h>
#include <light_mat/mateval/macc_policy.h>
#include "internal/ewise_eval_internal.h"

namespace lmat
{

	/********************************************
	 *
	 *  vectorized kernel
	 *
	 ********************************************/

	template<class Kernel>
	struct ewise_kernel
	{
	public:
		typedef Kernel kernel_type;

		LMAT_ENSURE_INLINE
		ewise_kernel(const Kernel& kernel)
		: m_kernel(kernel) { }

		LMAT_ENSURE_INLINE
		const Kernel& kernel() const
		{
			return m_kernel;
		}

		template<typename U, typename... Wraps>
		LMAT_ENSURE_INLINE
		void eval(macc_<linear_, U>, index_t m, index_t n, const Wraps&... wraps) const
		{
			dimension<0> dim(m * n);
			internal::_linear_ewise_eval(dim, U(), m_kernel, make_vec_accessor(U(), wraps)...);
		}

		template<typename U, index_t CM, index_t CN, typename... Wraps>
		LMAT_ENSURE_INLINE
		void eval(macc_<linear_, U>, const matrix_shape<CM, CN>& shape, const Wraps&... wraps) const
		{
			dimension<CM * CN> dim(shape.nelems());
			internal::_linear_ewise_eval(dim, U(), m_kernel, make_vec_accessor(U(), wraps)...);
		}

		template<typename U, typename... Wraps>
		LMAT_ENSURE_INLINE
		void eval(macc_<percol_, U>, index_t m, index_t n, const Wraps&... wraps) const
		{
			matrix_shape<0, 0> shape(m, n);
			internal::_percol_ewise_eval(shape, U(), m_kernel, make_multicol_accessor(U(), wraps)...);
		}

		template<typename U, index_t CM, index_t CN, typename... Wraps>
		LMAT_ENSURE_INLINE
		void eval(macc_<percol_, U>, const matrix_shape<CM, CN>& shape, const Wraps&... wraps) const
		{
			internal::_percol_ewise_eval(shape, U(), m_kernel, make_multicol_accessor(U(), wraps)...);
		}

		template<typename... Wraps>
		LMAT_ENSURE_INLINE
		void operator() (index_t m, index_t n, const Wraps&... wraps) const
		{
			eval(get_preferred_macc_policy(m, n, m_kernel, wraps...), m, n, wraps...);
		}

		template<index_t CM, index_t CN, typename... Wraps>
		LMAT_ENSURE_INLINE
		void operator() (const matrix_shape<CM, CN>& shape, const Wraps&... wraps) const
		{
			eval(get_preferred_macc_policy(shape, m_kernel, wraps...), shape, wraps...);
		}

	private:
		const Kernel& m_kernel;
	};


	template<class Kernel>
	LMAT_ENSURE_INLINE
	inline ewise_kernel<Kernel> ewise(const Kernel& kernel)
	{
		return ewise_kernel<Kernel>(kernel);
	}


	/********************************************
	 *
	 *  more friendly syntax construction
	 *
	 ********************************************/


	// map

	template<class Fun>
	LMAT_ENSURE_INLINE
	inline ewise_kernel<map_kernel<Fun> >
	map(const Fun& fun)
	{
		return ewise(map_kernel<Fun>(fun));
	}

	template<typename T, class DMat, class Fun, typename... Wraps>
	LMAT_ENSURE_INLINE
	inline void map_to(IRegularMatrix<DMat, T>& dmat, const Fun& fun, const Wraps&... wraps)
	{
		map(fun)(dmat.shape(), out_(dmat), wraps...);
	}


	// accum

	template<typename T, class DMat, class SMat>
	LMAT_ENSURE_INLINE
	inline void accum_to(IRegularMatrix<DMat, T>& dmat, const IEWiseMatrix<SMat, T>& smat)
	{
		return ewise(accum_kernel<T>())(dmat.shape(), in_out_(dmat), in_(smat));
	}

	template<typename T, class DMat, class SMat>
	LMAT_ENSURE_INLINE
	inline void accum_to(IRegularMatrix<DMat, T>& dmat, const T& c, const IEWiseMatrix<SMat, T>& smat)
	{
		return ewise(accumx_kernel<T>())(dmat.shape(), in_out_(dmat), const_(c), in_(smat));
	}

	template<typename T, class DMat, class CMat, class SMat>
	LMAT_ENSURE_INLINE
	inline void accum_to(IRegularMatrix<DMat, T>& dmat, const IEWiseMatrix<CMat, T>& cmat, const IEWiseMatrix<SMat, T>& smat)
	{
		return ewise(accumx_kernel<T>())(dmat.shape(), in_out_(dmat), in_(cmat), in_(smat));
	}


	/********************************************
	 *
	 *  Generic by-map evaluation
	 *
	 ********************************************/

	template<typename T, typename Acc, typename U, class Expr, class DMat>
	LMAT_ENSURE_INLINE
	inline void macc_evaluate(const IEWiseMatrix<Expr, T>& s, IRegularMatrix<DMat, T>& d, macc_<Acc, U> policy)
	{
		ewise(copy_kernel<T>()).eval(policy, common_shape(s.derived(), d.derived()), in_(s), out_(d));
	}

	template<typename T, class Expr, class DMat>
	LMAT_ENSURE_INLINE
	inline void macc_evaluate(const IEWiseMatrix<Expr, T>& s, IRegularMatrix<DMat, T>& d)
	{
		ewise(copy_kernel<T>())(common_shape(s.derived(), d.derived()), in_(s), out_(d));
	}

}

#endif /* EWISE_EVAL_H_ */
