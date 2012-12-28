/**
 * @file ewise_eval.h
 *
 * @brief Element-wise evaluation
 *
 * @author Dahua Lin
 */

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

	template<class Kernel, typename U>
	struct vectorized_kernel
	{
	public:
		typedef Kernel kernel_type;
		typedef typename meta::kernel_value_type<Kernel>::type value_type;
		typedef U unit_type;

		LMAT_ENSURE_INLINE
		vectorized_kernel(const Kernel& kernel)
		: m_kernel(kernel) { }

		LMAT_ENSURE_INLINE
		const Kernel& kernel() const
		{
			return m_kernel;
		}

		template<int L, typename... Accessors>
		LMAT_ENSURE_INLINE
		void apply(const dimension<L>& dim, const Accessors&... accs) const
		{
			internal::linear_ewise_eval(dim, U(), m_kernel, accs...);
		}

		template<typename... Accessors>
		LMAT_ENSURE_INLINE
		void apply(index_t len, const Accessors&... accs) const
		{
			apply(dimension<0>(len), accs...);
		}

		template<int M, int N, typename... Accessors>
		LMAT_ENSURE_INLINE
		void apply (const matrix_shape<M, N>& shape, const Accessors&... accs) const
		{
			apply(dimension<M * N>(shape.nelems()), accs...);
		}

		template<typename... Wraps>
		LMAT_ENSURE_INLINE
		void operator() (index_t len, const Wraps&... wraps) const
		{
			apply(len, make_vec_accessor(U(), wraps)...);
		}

		template<int L, typename... Wraps>
		LMAT_ENSURE_INLINE
		void operator() (const dimension<L>& dim, const Wraps&... wraps) const
		{
			apply(dim, make_vec_accessor(U(), wraps)...);
		}

		template<int M, int N, typename... Wraps>
		LMAT_ENSURE_INLINE
		void operator() (const matrix_shape<M, N>& shape, const Wraps&... wraps) const
		{
			apply(shape, make_vec_accessor(U(), wraps)...);
		}

	private:
		const Kernel& m_kernel;
	};


	/********************************************
	 *
	 *  more friendly syntax construction
	 *
	 ********************************************/

	// ewise

	template<class Kernel, typename U>
	LMAT_ENSURE_INLINE
	inline vectorized_kernel<Kernel, U>
	ewise(const Kernel& kernel, U)
	{
		return vectorized_kernel<Kernel, U>(kernel);
	}


	template<class Kernel>
	LMAT_ENSURE_INLINE
	inline vectorized_kernel<Kernel, default_access_unit_t>
	ewise(const Kernel& kernel)
	{
		return vectorized_kernel<Kernel, default_access_unit_t>(kernel);
	}


	// map

	template<class Fun, typename U>
	LMAT_ENSURE_INLINE
	inline vectorized_kernel<map_kernel<Fun>, U>
	map(const Fun& fun, U)
	{
		return ewise(map_kernel<Fun>(fun), U());
	}

	template<class Fun>
	LMAT_ENSURE_INLINE
	inline vectorized_kernel<map_kernel<Fun>, default_access_unit_t>
	map(const Fun& fun)
	{
		return ewise(map_kernel<Fun>(fun), default_access_unit_t());
	}

	template<typename T, class DMat, class Fun, typename... Wraps>
	LMAT_ENSURE_INLINE
	inline void map_to(IRegularMatrix<DMat, T>& dmat, const Fun& fun, const Wraps&... wraps)
	{
		map(fun)(dmat.shape(), out_(dmat.derived()), wraps...);
	}


	// accum

	template<typename T, class DMat, class SMat>
	LMAT_ENSURE_INLINE
	inline void accum_to(IRegularMatrix<DMat, T>& dmat, const IEWiseMatrix<SMat, T>& smat)
	{
		return ewise(accum_kernel<T>())(dmat.shape(), in_out_(dmat.derived()), in_(smat.derived()));
	}

	template<typename T, class DMat, class SMat>
	LMAT_ENSURE_INLINE
	inline void accum_to(IRegularMatrix<DMat, T>& dmat, const T& c, const IEWiseMatrix<SMat, T>& smat)
	{
		return ewise(accumx_kernel<T>())(dmat.shape(), in_out_(dmat.derived()), in_(c, atags::single()), in_(smat.derived()));
	}

	template<typename T, class DMat, class CMat, class SMat>
	LMAT_ENSURE_INLINE
	inline void accum_to(IRegularMatrix<DMat, T>& dmat, const IEWiseMatrix<CMat, T>& cmat, const IEWiseMatrix<SMat, T>& smat)
	{
		return ewise(accumx_kernel<T>())(dmat.shape(), in_out_(dmat.derived()), in_(cmat.derived()), in_(smat.derived()));
	}


	template<typename T, class DMat, typename Fun, typename... Wraps>
	LMAT_ENSURE_INLINE
	inline void accumf_to(IRegularMatrix<DMat, T>& dmat, const Fun& fun, const Wraps&... wraps)
	{
		return ewise(accumf_kernel<Fun>(fun))(dmat.shape(), in_out_(dmat.derived()), wraps...);
	}

	template<typename T, class DMat, typename Fun, typename... Wraps>
	LMAT_ENSURE_INLINE
	inline void accumf_to(IRegularMatrix<DMat, T>& dmat, const T& c, const Fun& fun, const Wraps&... wraps)
	{
		return ewise(accumfx_kernel<Fun>(fun))(dmat.shape(), in_out_(dmat.derived()), in_(c, atags::single()), wraps...);
	}


	// percol

	template<class VecKernel, typename... Wraps>
	LMAT_ENSURE_INLINE
	inline void percol(const VecKernel& veckernel, index_t m, index_t n, const Wraps&... wraps)
	{
		typedef typename VecKernel::unit_type U;
		dimension<0> coldim(m);
		internal::percol_eval_a(coldim, n, veckernel, make_multicol_accessor(U(), wraps)...);
	}

	template<class VecKernel, int M, int N, typename... Wraps>
	LMAT_ENSURE_INLINE
	inline void percol(const VecKernel& veckernel, const matrix_shape<M, N>& shape, const Wraps&... wraps)
	{
		typedef typename VecKernel::unit_type U;
		dimension<M> coldim(shape.nrows());
		internal::percol_eval_a(coldim, shape.ncolumns(), veckernel, make_multicol_accessor(U(), wraps)...);
	}


}

#endif /* EWISE_EVAL_H_ */
