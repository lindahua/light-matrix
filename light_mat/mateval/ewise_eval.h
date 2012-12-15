/**
 * @file ewise_eval.h
 *
 * @brief Element-wise evaluation
 *
 * @author Dahua Lin
 */

#ifndef LIGHTMAT_EWISE_EVAL_H_
#define LIGHTMAT_EWISE_EVAL_H_

#include <light_mat/mateval/mateval_fwd.h>
#include <light_mat/matrix/matrix_concepts.h>

namespace lmat
{
	// overall ewise evaluation function

	template<class Kernel, class Scheme, typename... Wraps>
	LMAT_ENSURE_INLINE
	inline void ewise_eval(const Kernel& kernel, const Scheme& scheme, Wraps... wraps)
	{
		scheme.apply(kernel, wraps...);
	}


	/********************************************
	 *
	 *  ewise schemes
	 *
	 ********************************************/

	// linear scheme

	template<typename U, int Len=0>
	class linear_ewise_scheme
	{
	public:
		LMAT_ENSURE_INLINE
		linear_ewise_scheme(index_t len)
		: m_dim(len) { }

		LMAT_ENSURE_INLINE
		index_t length() const
		{
			return m_dim.value();
		}

		template<class Kernel, typename... Wraps>
		LMAT_ENSURE_INLINE
		void apply(const Kernel& kernel, Wraps... wraps) const
		{
			internal::linear_ewise_eval(m_dim, U(), kernel, wraps...);
		}

	private:
		dimension<Len> m_dim;
	};


	template<typename U, int Len>
	LMAT_ENSURE_INLINE
	inline linear_ewise_scheme<U, Len> linear_ewise(U, fix_int<Len>)
	{
		return linear_ewise_scheme<U, Len>(Len);
	}

	template<typename U>
	LMAT_ENSURE_INLINE
	inline linear_ewise_scheme<U, 0> linear_ewise(U, index_t len)
	{
		return linear_ewise_scheme<U, 0>(len);
	}

	template<typename U, int M, int N>
	LMAT_ENSURE_INLINE
	inline linear_ewise_scheme<U, M * N> linear_ewise(U, const matrix_shape<M, N>& s)
	{
		return linear_ewise_scheme<U, M * N>(s.nrows() * s.ncolumns());
	}


	// per column scheme

	template<typename U, int M=0, int N=0>
	struct percol_ewise_scheme
	{
	public:
		LMAT_ENSURE_INLINE
		percol_ewise_scheme(index_t m, index_t n)
		: m_shape(m, n) { }

		LMAT_ENSURE_INLINE
		percol_ewise_scheme(const matrix_shape<M, N>& s)
		: m_shape(s) { }

		LMAT_ENSURE_INLINE
		index_t nrows() const
		{
			return m_shape.nrows();
		}

		LMAT_ENSURE_INLINE
		index_t ncolumns() const
		{
			return m_shape.ncolumns();
		}

		template<class Kernel, typename... Wraps>
		LMAT_ENSURE_INLINE
		void apply(const Kernel& kernel, Wraps... wraps) const
		{
			internal::percol_ewise_eval(m_shape, U(), kernel, wraps...);
		}

	private:
		matrix_shape<M, N> m_shape;
	};


	template<typename U, int M, int N>
	LMAT_ENSURE_INLINE
	inline percol_ewise_scheme<U, M, N> percol_ewise(U, fix_int<M>, fix_int<N>)
	{
		return percol_ewise_scheme<U, M, N>(M, N);
	}

	template<typename U, int M>
	LMAT_ENSURE_INLINE
	inline percol_ewise_scheme<U, M, 0> percol_ewise(U, fix_int<M>, index_t n)
	{
		return percol_ewise_scheme<U, M, 0>(M, n);
	}

	template<typename U, int N>
	LMAT_ENSURE_INLINE
	inline percol_ewise_scheme<U, 0, N> percol_ewise(U, index_t m, fix_int<N>)
	{
		return percol_ewise_scheme<U, 0, N>(m, N);
	}

	template<typename U>
	LMAT_ENSURE_INLINE
	inline percol_ewise_scheme<U, 0, 0> percol_ewise(U, index_t m, index_t n)
	{
		return percol_ewise_scheme<U, 0, 0>(m, n);
	}

	template<typename U, int M, int N>
	LMAT_ENSURE_INLINE
	inline percol_ewise_scheme<U, M * N> percol_ewise(U, const matrix_shape<M, N>& s)
	{
		return percol_ewise_scheme<U, M, N>(s);
	}

}

#endif /* EWISE_EVAL_H_ */
