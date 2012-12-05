/**
 * @file macc_eval_core.h
 *
 * @brief The core of access-based evaluation
 *
 * @author Dahua Lin
 */

#ifndef LIGHTMAT_MACC_EVAL_CORE_H_
#define LIGHTMAT_MACC_EVAL_CORE_H_

#include <light_mat/common/basic_defs.h>

namespace lmat {

	struct linear_macc;
	struct percol_macc;

	template<class SExpr, class Policy>
	struct macc_accessor_map;

	template<class Accessor>
	struct percol_macc_state_map;

	namespace internal {

		/****************************************
		 *
		 *  vector reader/writer
		 *
		 ****************************************/

		template<typename Ker, typename T> class ContVecRW;
		template<typename Ker, typename T> class StepVecRW;

		template<typename T>
		class ContVecRW<scalar_ker, T>
		{
		public:
			LMAT_ENSURE_INLINE
			ContVecRW(T* p) : m_p(p) { }

			LMAT_ENSURE_INLINE
			T operator[] (index_t i) const { return m_p[i]; }

			LMAT_ENSURE_INLINE
			T& operator[] (index_t i) { return m_p[i]; }

		private:
			T *m_p;
		};

		template<typename T>
		class StepVecRW<scalar_ker, T>
		{
		public:
			LMAT_ENSURE_INLINE
			StepVecRW(T* p, index_t step) : m_p(p), m_step(step) { }

			LMAT_ENSURE_INLINE
			T operator[] (index_t i) const { return m_p[i * m_step]; }

			LMAT_ENSURE_INLINE
			T& operator[] (index_t i) { return m_p[i * m_step]; }

		private:
			T *m_p;
			const index_t m_step;
		};


		/****************************************
		 *
		 *  copy
		 *
		 ****************************************/

		template<typename Ker, int Len>
		struct macc_vec_copy;

		template<int Len>
		struct macc_vec_copy<scalar_ker, Len>
		{
			template<class Acc, class Vec>
			LMAT_ENSURE_INLINE
			static void eval(index_t, const Acc& in, Vec& out)
			{
				for (index_t i = 0; i < Len; ++i)
				{
					out[i] = in.get_scalar(i);
				}
			}
		};


		template<>
		struct macc_vec_copy<scalar_ker, 0>
		{
			template<class Acc, class Vec>
			LMAT_ENSURE_INLINE
			static void eval(index_t len, const Acc& in, Vec& out)
			{
				for (index_t i = 0; i < len; ++i)
				{
					out[i] = in.get_scalar(i);
				}
			}
		};



	}
}

#endif /* MACC_EVAL_CORE_H_ */
