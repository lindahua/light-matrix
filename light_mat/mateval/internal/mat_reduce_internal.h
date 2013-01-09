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

namespace lmat {


	template<class Folder>
	struct _parfold_kernel
	{
		typedef typename Folder::value_type value_type;
		Folder folder;

		LMAT_ENSURE_INLINE
		explicit _parfold_kernel(const Folder& rf) : folder(rf) { }

		LMAT_ENSURE_INLINE
		void operator() (value_type& a, const value_type& x) const
		{
			folder.fold(a, x);
		}
	};

	template<class Folder>
	struct is_simdizable<_parfold_kernel<Folder> >
	{
		static const bool value = is_simdizable<Folder>::value;
	};

	template<class Folder, typename Kind>
	struct simdize_map<_parfold_kernel<Folder>, Kind>
	{
		typedef typename simdize_map<Folder, Kind>::type simd_folder_t;
		typedef _parfold_kernel<simd_folder_t> type;

		LMAT_ENSURE_INLINE
		static type get(const _parfold_kernel<Folder>& k)
		{
			return type(simdize_map<Folder, Kind>::get(k.folder));
		}
	};


namespace internal {

	/********************************************
	 *
	 *  helpers on atag
	 *
	 ********************************************/

	template<class Folder, class TExpr>
	struct full_reduc_policy
	{
		static const bool use_linear = supports_linear_macc<TExpr>::value;

		typedef typename matrix_traits<TExpr>::value_type vtype;
		typedef default_simd_kind simd_kind;

		static const bool use_simd =
				is_simdizable<Folder>::value &&
				supports_simd<TExpr, vtype, simd_kind, use_linear>::value;

		typedef typename std::conditional<use_simd,
				atags::simd<simd_kind>,
				atags::scalar>::type atag;

		typedef typename std::conditional<use_linear,
				linear_macc<atag>,
				percol_macc<atag> >::type type;
	};


	template<class Folder, class TExpr, class DMat>
	struct colwise_reduc_policy
	{
		static_assert(meta::supports_linear_index<DMat>::value,
				"DMat should support linear indexing.");

		typedef typename matrix_traits<TExpr>::value_type vtype;
		typedef default_simd_kind simd_kind;

		static const bool use_simd = is_simdizable<Folder>::value &&
			supports_simd<TExpr, vtype, simd_kind, false>::value;

		typedef typename std::conditional<use_simd,
				atags::simd<simd_kind>,
				atags::scalar>::type atag;
	};


	template<class Folder, class TExpr, class DMat>
	struct rowwise_reduc_policy
	{
		static_assert(supports_linear_macc<DMat>::value,
				"DMat should support linear access.");

		typedef typename matrix_traits<TExpr>::value_type vtype;
		typedef default_simd_kind simd_kind;

		static const bool use_simd = is_simdizable<Folder>::value &&
			supports_simd<TExpr, vtype, simd_kind, false>::value;

		typedef typename std::conditional<use_simd,
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
	 *  full reduction implementation
	 *
	 ********************************************/

	template<int M, int N, class Folder, class A, typename U>
	LMAT_ENSURE_INLINE
	inline typename Folder::value_type
	_full_reduce(const matrix_shape<M, N>& shape, const Folder& folder,
			const A& a, linear_macc<U>)
	{
		dimension<meta::nelems<A>::value> dim = a.nelems();
		return fold_impl( dim, U(), folder, make_vec_accessor(U(), in_(a)) );
	}

	template<int M, int N, class Folder, class A, typename U>
	inline typename Folder::value_type
	_full_reduce(const matrix_shape<M, N>& shape, const Folder& folder,
			const A& a, percol_macc<U>)
	{
		typedef typename Folder::value_type T;

		dimension<meta::nrows<A>::value> col_dim = a.nrows();
		auto rd = make_multicol_accessor(U(), in_(a));

		T r = fold_impl(col_dim, U(), folder, rd.col(0));

		const index_t n = shape.ncolumns();
		for (index_t j = 1; j < n; ++j)
		{
			T rj = fold_impl(col_dim, U(), folder, rd.col(j));
			folder.fold(r, rj);
		}

		return r;
	}

	template<int M, int N, class Folder, class A>
	LMAT_ENSURE_INLINE
	inline typename Folder::value_type
	_full_reduce(const matrix_shape<M, N>& shape, const Folder& folder, const A& a)
	{
		typedef typename full_reduc_policy<Folder, A>::type policy_t;
		return _full_reduce(shape, folder, a, policy_t());
	}


	/********************************************
	 *
	 *  vector-wise reduction implementation
	 *
	 ********************************************/

	// column wise reduction

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


	// row wise reduction

	template<int M, int N, typename U, class Folder, class DMat, class MultiColReader>
	inline void rowwise_fold_impl(const matrix_shape<M, N>& shape, U u,
			const Folder& folder, DMat& dmat, const MultiColReader& rd)
	{
		typedef typename matrix_traits<DMat>::value_type T;

		dimension<M> col_dim(shape.nrows());
		const index_t n = shape.ncolumns();

		LMAT_CHECK_DIMS( col_dim.value() == dmat.nelems() )

		auto a = make_vec_accessor(u, in_out_(dmat));

		internal::linear_ewise_eval(col_dim, u, copy_kernel<T>(), rd.col(0), a);

		_parfold_kernel<Folder> pfker(folder);
		for (index_t j = 1; j < n; ++j)
		{
			internal::linear_ewise_eval(col_dim, u, pfker, a, rd.col(j));
		}
	}


} }

#endif 
