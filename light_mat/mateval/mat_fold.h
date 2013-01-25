/**
 * @file mat_fold.h
 *
 * Matrix folding (this is the basis for matrix reduction)
 * 
 * @author Dahua Lin 
 */

#ifdef _MSC_VER
#pragma once
#endif

#ifdef _MSC_VER
#pragma once
#endif

#ifndef LIGHTMAT_MAT_FOLD_H_
#define LIGHTMAT_MAT_FOLD_H_

#include <light_mat/mateval/ewise_eval.h>
#include "internal/mat_fold_internal.h"

#include <utility>


#define LMAT_DEFINE_SIMPLE_FOLD_KERNEL( Name, InitExpr, FoldExpr, ReducExpr ) \
	template<typename T> \
	struct Name##_kernel { \
		typedef T value_type; \
		typedef T accumulated_type; \
		LMAT_ENSURE_INLINE \
		T init(const T& x) const { return InitExpr; } \
		LMAT_ENSURE_INLINE \
		void operator()(T& a, const T& x) const { FoldExpr; } \
	}; \
	template<typename T, typename Kind> \
	struct Name##_kernel<simd_pack<T, Kind> > { \
		typedef simd_pack<T, Kind> value_type; \
		typedef simd_pack<T, Kind> accumulated_type; \
		LMAT_ENSURE_INLINE \
		accumulated_type init(const value_type& x) const { return InitExpr; } \
		LMAT_ENSURE_INLINE \
		void operator()(accumulated_type& a, const value_type& x) const { FoldExpr; } \
		LMAT_ENSURE_INLINE \
		T reduce(const accumulated_type& a) const { return ReducExpr; } \
	}; \
	LMAT_DECL_SIMDIZABLE_ON_REAL( Name##_kernel ) \
	LMAT_DEF_TRIVIAL_SIMDIZE_MAP( Name##_kernel )


#define LMAT_DEFINE_AGGREG_SIMD_FOLDKERNEL(StatT, Kernel, NA) \
		template<typename T> struct StatT##_; \
		template<typename T, typename Kind> \
		struct StatT##_<lmat::simd_pack<T, Kind> > : public StatT<lmat::simd_pack<T, Kind> > \
		{ \
			static const unsigned int pack_width = lmat::simd_traits<T, Kind>::pack_width; \
			LMAT_ENSURE_INLINE \
			StatT##_() { } \
			LMAT_ENSURE_INLINE \
			StatT##_(const StatT<lmat::simd_pack<T, Kind> >& s) \
			: StatT<lmat::simd_pack<T, Kind> >(s) { } \
		}; \
		template<typename T> \
		struct Kernel { \
			typedef StatT<T> accumulated_type; \
			typedef T value_type; \
			LMAT_ENSURE_INLINE \
			accumulated_type init(LMAT_REPEAT_ARGS_P##NA(const T& x)) const \
			{ \
				return accumulated_type(LMAT_REPEAT_ARGS_P##NA(x)); \
			} \
			LMAT_ENSURE_INLINE \
			void operator() (accumulated_type& a, LMAT_REPEAT_ARGS_P##NA(const T& x)) const \
			{ \
				a.update(LMAT_REPEAT_ARGS_P##NA(x)); \
			} \
			LMAT_ENSURE_INLINE \
			void operator() (accumulated_type& a, const accumulated_type& b) const \
			{ \
				a.update(b); \
			} \
		}; \
		template<typename T, typename Kind> \
		struct Kernel<lmat::simd_pack<T, Kind> > { \
			typedef StatT##_<lmat::simd_pack<T, Kind> > accumulated_type; \
			typedef lmat::simd_pack<T, Kind> value_type; \
			LMAT_ENSURE_INLINE \
			accumulated_type init(LMAT_REPEAT_ARGS_P##NA(const value_type& x)) const \
			{ \
				return StatT<lmat::simd_pack<T, Kind> >(LMAT_REPEAT_ARGS_P##NA(x)); \
			} \
			LMAT_ENSURE_INLINE \
			void operator() (accumulated_type& a, LMAT_REPEAT_ARGS_P##NA(const value_type& x)) const \
			{ \
				a.update(LMAT_REPEAT_ARGS_P##NA(x)); \
			} \
			LMAT_ENSURE_INLINE \
			void operator() (accumulated_type& a, const accumulated_type& b) const \
			{ \
				a.update(b); \
			} \
			LMAT_ENSURE_INLINE \
			StatT<T> reduce(const accumulated_type& a) const \
			{ \
				return reduce_impl(a); \
			} \
		}; \


namespace lmat
{

	/********************************************
	 *
	 *  common folders
	 *
	 ********************************************/

	LMAT_DEFINE_SIMPLE_FOLD_KERNEL( sum, x, a += x, sum(a) )

	LMAT_DEFINE_SIMPLE_FOLD_KERNEL( maximum, x, a = math::max(a, x), maximum(a) )

	LMAT_DEFINE_SIMPLE_FOLD_KERNEL( minimum, x, a = math::min(a, x), minimum(a) )



	/********************************************
	 *
	 *  folder on matrices
	 *
	 ********************************************/

	template<class FoldKernel>
	class matrix_folder
	{
	public:
		typedef typename FoldKernel::accumulated_type result_type;

		LMAT_ENSURE_INLINE
		matrix_folder(const FoldKernel& folder)
		: m_kernel(folder) { }

		template<typename U, index_t CM, index_t CN, typename... Wrap>
		LMAT_ENSURE_INLINE
		result_type eval(macc_<linear_, U>, const matrix_shape<CM, CN>& shape, const Wrap&... wrap) const
		{
			dimension<CM * CN> dim(shape.nelems());
			return internal::linear_fold_impl(dim, U(), m_kernel, make_vec_accessor(U(), wrap)...);
		}

		template<typename U, typename... Wrap>
		LMAT_ENSURE_INLINE
		result_type eval(macc_<linear_, U>, index_t m, index_t n, const Wrap&... wrap) const
		{
			dimension<0> dim(m * n);
			return internal::linear_fold_impl(dim, U(), m_kernel, make_vec_accessor(U(), wrap)...);
		}

		template<typename U, index_t CM, index_t CN, typename... Wrap>
		LMAT_ENSURE_INLINE
		result_type eval(macc_<percol_, U>, const matrix_shape<CM, CN>& shape, const Wrap&... wrap) const
		{
			return internal::percol_fold_impl(shape, U(), m_kernel, make_multicol_accessor(U(), wrap)...);
		}

		template<typename U, typename... Wrap>
		LMAT_ENSURE_INLINE
		result_type eval(macc_<percol_, U>, index_t m, index_t n, const Wrap&... wrap) const
		{
			matrix_shape<0,0> shape(m, n);
			return internal::percol_fold_impl(shape, U(), m_kernel, make_multicol_accessor(U(), wrap)...);
		}

		template<index_t CM, index_t CN, typename... Wrap>
		LMAT_ENSURE_INLINE
		result_type operator() (const matrix_shape<CM, CN>& shape, const Wrap&... wrap) const
		{
			typedef typename internal::fold_policy<FoldKernel, matrix_shape<CM, CN>, Wrap...>::type policy_t;
			return eval(policy_t(), shape, wrap...);
		}

		template<typename... Wrap>
		LMAT_ENSURE_INLINE
		result_type operator() (index_t m, index_t n, const Wrap&... wrap) const
		{
			typedef typename internal::fold_policy<FoldKernel, matrix_shape<0, 0>, Wrap...>::type policy_t;
			return eval(policy_t(), m, n, wrap...);
		}

	private:
		FoldKernel m_kernel;
	};


	/********************************************
	 *
	 *  convenient function
	 *
	 ********************************************/

	template<class Folder>
	LMAT_ENSURE_INLINE
	inline matrix_folder<Folder> fold(const Folder& folder)
	{
		return matrix_folder<Folder>(folder);
	}

}

#endif 
