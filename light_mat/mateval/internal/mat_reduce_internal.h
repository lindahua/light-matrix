/**
 * @file mat_reduce_internal.h
 *
 * Internal implementation of matrix reduction
 * 
 * @author Dahua Lin 
 */

#ifdef _MSC_VER
#pragma once
#endif

#ifndef LIGHTMAT_MAT_REDUCE_INTERNAL_H_
#define LIGHTMAT_MAT_REDUCE_INTERNAL_H_

#include <light_mat/mateval/mat_fold.h>
#include <light_mat/mateval/common_kernels.h>
#include <light_mat/math/math_functors.h>
#include <light_mat/mateval/mat_arith.h>
#include <light_mat/mateval/mat_emath.h>

namespace lmat { namespace internal {


	template<class RFun>
	struct _fold_kernel
	{
		typedef typename RFun::value_type value_type;
		RFun rfun;

		LMAT_ENSURE_INLINE
		_fold_kernel(const RFun& rf) : rfun(rf) { }

		LMAT_ENSURE_INLINE
		void operator() (value_type& a, const value_type& x) const
		{
			rfun.fold(a, x);
		}

		template<typename Kind>
		LMAT_ENSURE_INLINE
		void operator() (math::simd_pack<value_type, Kind>& a,
				const math::simd_pack<value_type, Kind>& x) const
		{
			rfun.fold(a, x);
		}
	};


	/********************************************
	 *
	 *  helpers on atag
	 *
	 ********************************************/

	template<class Folder, class TExpr>
	struct full_reduc_policy
	{
		static_assert(prefers_linear<TExpr>::value,
				"TExpr should allow linear access.");

		static const bool use_linear = prefers_linear<TExpr>::value;

		typedef typename matrix_traits<TExpr>::value_type vtype;
		typedef default_simd_kind simd_kind;

		static const bool use_simd =
				folder_supports_simd<Folder>::value &&
				prefers_simd<TExpr, vtype, simd_kind, use_linear>::value;

		typedef typename meta::if_c<use_simd,
				atags::simd<simd_kind>,
				atags::scalar>::type atag;
	};


	template<class Folder, class TExpr, class DMat>
	struct colwise_reduc_policy
	{
		static_assert(meta::supports_linear_index<DMat>::value,
				"DMat should support linear indexing.");

		static const bool use_linear = false;

		typedef typename matrix_traits<TExpr>::value_type vtype;
		typedef default_simd_kind simd_kind;

		static const bool use_simd = folder_supports_simd<Folder>::value &&
			prefers_simd<TExpr, vtype, simd_kind, use_linear>::value;

		typedef typename meta::if_c<use_simd,
				atags::simd<simd_kind>,
				atags::scalar>::type atag;
	};


	template<class Folder, class TExpr, class DMat>
	struct rowwise_reduc_policy
	{
		static_assert(meta::is_continuous<DMat>::value,
				"DMat should be continuous.");

		static const bool use_linear = false;

		typedef typename matrix_traits<TExpr>::value_type vtype;
		typedef default_simd_kind simd_kind;

		static const bool use_simd = folder_supports_simd<Folder>::value &&
			prefers_simd<TExpr, vtype, simd_kind, use_linear>::value;

		typedef typename meta::if_c<use_simd,
				atags::simd<simd_kind>,
				atags::scalar>::type atag;
	};


	/********************************************
	 *
	 *  helpers on shape and empty value
	 *
	 ********************************************/

	template<typename T, class Mat>
	LMAT_ENSURE_INLINE
	inline index_t reduc_get_length(const IEWiseMatrix<Mat, T>& a)
	{
		return a.nelems();
	}

	template<typename T, class Mat1, class Mat2>
	LMAT_ENSURE_INLINE
	inline index_t reduc_get_length(const IEWiseMatrix<Mat1, T>& a, const IEWiseMatrix<Mat2, T>& b)
	{
		return common_shape(a.derived(), b.derived()).nelems();
	}


	template<typename T, class Mat>
	LMAT_ENSURE_INLINE
	inline typename meta::shape<Mat>::type
	reduc_get_shape(const IEWiseMatrix<Mat, T>& a)
	{
		return a.shape();
	}

	template<typename T, class Mat1, class Mat2>
	LMAT_ENSURE_INLINE
	inline typename meta::common_shape<Mat1, Mat2>::type
	reduc_get_shape(const IEWiseMatrix<Mat1, T>& a, const IEWiseMatrix<Mat2, T>& b)
	{
		return common_shape(a.derived(), b.derived());
	}

	template<typename T>
	struct empty_values
	{
		LMAT_ENSURE_INLINE
		static T sum() { return T(0); }

		LMAT_ENSURE_INLINE
		static T mean() { return std::numeric_limits<T>::quiet_NaN(); }

		LMAT_ENSURE_INLINE
		static T maximum() { return - std::numeric_limits<T>::infinity(); }

		LMAT_ENSURE_INLINE
		static T minimum() { return std::numeric_limits<T>::infinity(); }
	};



	/********************************************
	 *
	 *  colwise reduction
	 *
	 ********************************************/

	template<int M, int N, typename U, class Folder, class DMat, class MultiColReader>
	inline void colwise_fold_impl(const matrix_shape<M, N>& shape, U u,
			const Folder& folder, DMat& dmat, const MultiColReader& rd)
	{
		dimension<M> col_dim(shape.nrows());
		const index_t n = shape.ncolumns();

		LMAT_CHECK_DIMS( n == dmat.nelems() )

		vecfold_kernel<Folder, U> fker = fold(folder, u);

		for (index_t j = 0; j < n; ++j)
		{
			dmat[j] = fker.apply(col_dim, rd.col(j));
		}
	}


	/********************************************
	 *
	 *  rowwise reduction
	 *
	 ********************************************/

	template<int M, int N, typename U, class RFun, class DMat, class MultiColReader>
	inline void rowwise_fold_impl(const matrix_shape<M, N>& shape, U u,
			RFun rfun, DMat& dmat, const MultiColReader& rd)
	{
		typedef typename matrix_traits<DMat>::value_type T;

		dimension<M> col_dim(shape.nrows());
		const index_t n = shape.ncolumns();

		LMAT_CHECK_DIMS( col_dim.value() == dmat.nelems() )

		auto a = make_vec_accessor(u, in_out_(dmat));

		internal::linear_ewise_eval(col_dim, u, copy_kernel<T>(), rd.col(0), a);

		_fold_kernel<RFun> fker(rfun);
		for (index_t j = 1; j < n; ++j)
		{
			internal::linear_ewise_eval(col_dim, u, fker, a, rd.col(j));
		}
	}


} }

#endif 
