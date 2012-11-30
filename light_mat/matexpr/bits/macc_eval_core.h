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

	template<class SExpr, class AccCate, class KerCate>
	struct macc_accessor_map;

	template<class Accessor>
	struct percol_macc_state_map;

	namespace internal {

		/****************************************
		 *
		 *  vector reader/writer
		 *
		 ****************************************/

		template<class Vec>
		class LinVecRW
		{
		public:
			typedef typename matrix_traits<Vec>::value_type T;

			LMAT_ENSURE_INLINE
			LinVecRW(Vec& vec) : m_vec(vec) { }

			LMAT_ENSURE_INLINE
			T operator[] (index_t i) const { return m_vec[i]; }

			LMAT_ENSURE_INLINE
			T& operator[] (index_t i) { return m_vec[i]; }

		private:
			Vec& m_vec;
		};

		template<typename T>
		class ContVecRW
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
		class StepVecRW
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
		struct macc_vec_copy<scalar_kernel_t, Len>
		{
			template<class Acc, class Vec>
			static void eval(index_t, const Acc& in, Vec& out)
			{
				for (index_t i = 0; i < Len; ++i)
				{
					out[i] = in.get_scalar(i);
				}
			}
		};


		template<>
		struct macc_vec_copy<scalar_kernel_t, 0>
		{
			template<class Acc, class Vec>
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
