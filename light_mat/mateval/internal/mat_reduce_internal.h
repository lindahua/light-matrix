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


	template<class RFun, class TFun>
	struct _foldx_kernel
	{
		typedef typename RFun::value_type value_type;
		RFun rfun;
		TFun tfun;

		LMAT_ENSURE_INLINE
		_foldx_kernel(const RFun& rf, const TFun& tf)
		: rfun(rf), tfun(tf) { }

		template<typename... A>
		LMAT_ENSURE_INLINE
		void operator() (value_type& a, const A&... x) const
		{
			rfun.fold(a, tfun(x...));
		}

		template<typename Kind, typename... A>
		LMAT_ENSURE_INLINE
		void operator() (math::simd_pack<value_type, Kind>& a,
				const math::simd_pack<A, Kind>&... x) const
		{
			rfun.fold(a, tfun(x...));
		}
	};


	/********************************************
	 *
	 *  full reduction
	 *
	 ********************************************/

	template<typename T, int N, typename U, class Wrap>
	inline T sum_(type_<T>, const dimension<N>& dim, U u, const Wrap& wrap)
	{
		return fold(sum_folder<T>(), u)(dim, wrap);
	}

	template<int N, typename T, typename U, class Wrap>
	inline T mean_(type_<T>, const dimension<N>& dim, U u, const Wrap& wrap)
	{
		T r = sum_(type_<T>(), dim, u, wrap);
		return r / T(dim.value());
	}

	template<int N, typename T, typename U, class Wrap>
	inline T maximum_(type_<T>, const dimension<N>& dim, U u, const Wrap& wrap)
	{
		return fold(maximum_folder<T>(), u)(dim, wrap);
	}

	template<int N, typename T, typename U, class Wrap>
	inline T minimum_(type_<T>, const dimension<N>& dim, U u, const Wrap& wrap)
	{
		return fold(minimum_folder<T>(), u)(dim, wrap);
	}


	template<int N, typename T, typename U, typename TFun, typename... Wrap>
	inline T sumx_(type_<T>, const dimension<N>& dim, U u, const TFun& tfun, const Wrap&... wraps)
	{
		return foldf(sum_folder<T>(), tfun, u)(dim, wraps...);
	}

	template<int N, typename T, typename U, typename TFun, typename... Wrap>
	inline T meanx_(type_<T>, const dimension<N>& dim, U u, const TFun& tfun, const Wrap&... wraps)
	{
		T r = sumx_(type_<T>(), dim, u, tfun, wraps...);
		return r / T(dim.value());
	}

	template<int N, typename T, typename U, typename TFun, typename... Wrap>
	inline T maximumx_(type_<T>, const dimension<N>& dim, U u, const TFun& tfun, const Wrap&... wraps)
	{
		return foldf(maximum_folder<T>(), tfun, u)(dim, wraps...);
	}

	template<int N, typename T, typename U, typename TFun, typename... Wrap>
	inline T minimumx_(type_<T>, const dimension<N>& dim, U u, const TFun& tfun, const Wrap&... wraps)
	{
		return foldf(minimum_folder<T>(), tfun, u)(dim, wraps...);
	}

	/********************************************
	 *
	 *  colwise reduction
	 *
	 ********************************************/

	template<typename T>
	struct div_by_dim
	{
		LMAT_ENSURE_INLINE
		div_by_dim(index_t n)
		: coeff(T(1) / T(n)) { }

		T operator() (const T& x) const
		{
			return x * coeff;
		}

	private:
		T coeff;
	};


	template<typename Kind, class DMat, class TFun>
	LMAT_ENSURE_INLINE
	inline void colwise_post(atags::simd<Kind>, index_t n, DMat& dmat, const TFun& tfun)
	{
		for (index_t j = 0; j < n; ++j) dmat[j] = tfun(dmat[j]);
	}


	template<int M, int N, typename U, class Folder, class DMat, class MultiColReader>
	inline void colwise_fold_impl(const matrix_shape<M, N>& shape, U u,
			const Folder& folder, DMat& dmat, const MultiColReader& rd)
	{
		dimension<M> col_dim(shape.nrows());
		const index_t n = shape.ncolumns();

		vecfold_kernel<Folder, U> fker = fold(folder, u);

		for (index_t j = 0; j < n; ++j)
		{
			dmat[j] = fker.apply(col_dim, rd.col(j));
		}
	}

	template<int M, int N, typename U, class Folder, class DMat, typename TFun, typename... MultiColReader>
	inline void colwise_foldx_impl(const matrix_shape<M, N>& shape, U u,
			const Folder& folder, DMat& dmat, const TFun& tfun, const MultiColReader&... rds)
	{
		dimension<M> col_dim(shape.nrows());
		const index_t n = shape.ncolumns();

		vecfoldf_kernel<Folder, TFun, U> fker = foldf(folder, tfun, u);

		for (index_t j = 0; j < n; ++j)
		{
			dmat[j] = fker.apply(col_dim, rds.col(j)...);
		}
	}


	template<int M, int N, typename Kind, class DMat, class Wrap>
	inline void colwise_sum_(const matrix_shape<M, N>& shape, atags::simd<Kind> u, DMat& dmat, const Wrap& wrap)
	{
		typedef typename matrix_traits<DMat>::value_type T;
		colwise_fold_impl(shape, u, sum_folder<T>(), dmat, make_multicol_accessor(u, wrap));
	}

	template<int M, int N, typename Kind, class DMat, class Wrap>
	inline void colwise_mean_(const matrix_shape<M, N>& shape, atags::simd<Kind> u, DMat& dmat, const Wrap& wrap)
	{
		typedef typename matrix_traits<DMat>::value_type T;
		const index_t m = shape.nrows();
		const index_t n = shape.ncolumns();

		colwise_fold_impl(shape, u, sum_folder<T>(), dmat, make_multicol_accessor(u, wrap));
		colwise_post(u, n, dmat, div_by_dim<T>(m));
	}

	template<int M, int N, typename Kind, class DMat, class Wrap>
	inline void colwise_maximum_(const matrix_shape<M, N>& shape, atags::simd<Kind> u, DMat& dmat, const Wrap& wrap)
	{
		typedef typename matrix_traits<DMat>::value_type T;
		colwise_fold_impl(shape, u, maximum_folder<T>(), dmat, make_multicol_accessor(u, wrap));
	}

	template<int M, int N, typename Kind, class DMat, class Wrap>
	inline void colwise_minimum_(const matrix_shape<M, N>& shape, atags::simd<Kind> u, DMat& dmat, const Wrap& wrap)
	{
		typedef typename matrix_traits<DMat>::value_type T;
		colwise_fold_impl(shape, u, minimum_folder<T>(), dmat, make_multicol_accessor(u, wrap));
	}


	template<int M, int N, typename Kind, class DMat, typename TFun, typename... Wrap>
	inline void colwise_sumx_(const matrix_shape<M, N>& shape, atags::simd<Kind> u, DMat& dmat,
			const TFun& tfun, const Wrap&... wraps)
	{
		typedef typename matrix_traits<DMat>::value_type T;
		colwise_foldx_impl(shape, u, sum_folder<T>(), dmat, tfun, make_multicol_accessor(u, wraps)...);
	}

	template<int M, int N, typename Kind, class DMat, typename TFun, typename... Wrap>
	inline void colwise_meanx_(const matrix_shape<M, N>& shape, atags::simd<Kind> u, DMat& dmat,
			const TFun& tfun, const Wrap&... wraps)
	{
		typedef typename matrix_traits<DMat>::value_type T;
		const index_t m = shape.nrows();
		const index_t n = shape.ncolumns();

		colwise_foldx_impl(shape, u, sum_folder<T>(), dmat, tfun, make_multicol_accessor(u, wraps)...);
		colwise_post(u, n, dmat, div_by_dim<T>(m));
	}

	template<int M, int N, typename Kind, class DMat, typename TFun, typename... Wrap>
	inline void colwise_maximumx_(const matrix_shape<M, N>& shape, atags::simd<Kind> u, DMat& dmat,
			const TFun& tfun, const Wrap&... wraps)
	{
		typedef typename matrix_traits<DMat>::value_type T;
		colwise_foldx_impl(shape, u, maximum_folder<T>(), dmat, tfun, make_multicol_accessor(u, wraps)...);
	}

	template<int M, int N, typename Kind, class DMat, typename TFun, typename... Wrap>
	inline void colwise_minimumx_(const matrix_shape<M, N>& shape, atags::simd<Kind> u, DMat& dmat,
			const TFun& tfun, const Wrap&... wraps)
	{
		typedef typename matrix_traits<DMat>::value_type T;
		colwise_foldx_impl(shape, u, minimum_folder<T>(), dmat, tfun, make_multicol_accessor(u, wraps)...);
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
		dimension<M> col_dim(shape.nrows());
		const index_t n = shape.ncolumns();
		typedef typename matrix_traits<DMat>::value_type T;

		auto a = make_vec_accessor(u, in_out_(dmat));

		internal::linear_ewise_eval(col_dim, u, copy_kernel<T>(), rd.col(0), a);

		_fold_kernel<RFun> fker(rfun);
		for (index_t j = 1; j < n; ++j)
		{
			internal::linear_ewise_eval(col_dim, u, fker, a, rd.col(j));
		}
	}

	template<int M, int N, typename U, class RFun, class DMat, typename TFun, typename... MultiColReader>
	inline void rowwise_foldx_impl(const matrix_shape<M, N>& shape, U u,
			RFun rfun, DMat& dmat, const TFun& tfun, const MultiColReader&... rds)
	{
		dimension<M> col_dim(shape.nrows());
		const index_t n = shape.ncolumns();

		auto a = make_vec_accessor(u, in_out_(dmat));

		internal::linear_ewise_eval(col_dim, u, map_kernel<TFun>(tfun), a, rds.col(0)...);

		_foldx_kernel<RFun, TFun> fker(rfun, tfun);
		for (index_t j = 1; j < n; ++j)
		{
			internal::linear_ewise_eval(col_dim, u, fker, a, rds.col(j)...);
		}
	}


	template<int M, int N, typename Kind, class DMat, class Wrap>
	inline void rowwise_sum_(const matrix_shape<M, N>& shape, atags::simd<Kind> u, DMat& dmat, const Wrap& wrap)
	{
		typedef typename matrix_traits<DMat>::value_type T;
		rowwise_fold_impl(shape, u, sum_folder<T>(), dmat, make_multicol_accessor(u, wrap));
	}

	template<int M, int N, typename Kind, class DMat, class Wrap>
	inline void rowwise_mean_(const matrix_shape<M, N>& shape, atags::simd<Kind> u, DMat& dmat, const Wrap& wrap)
	{
		typedef typename matrix_traits<DMat>::value_type T;
		rowwise_fold_impl(shape, u, sum_folder<T>(), dmat, make_multicol_accessor(u, wrap));

		T c = T(1) / T(shape.ncolumns());
		dimension<M> dim(shape.nrows());
		map(math::mul_fun<T>(), atags::simd<Kind>())(dim, out_(dmat), in_(dmat), in_(c, atags::single()));
	}

	template<int M, int N, typename Kind, class DMat, class Wrap>
	inline void rowwise_maximum_(const matrix_shape<M, N>& shape, atags::simd<Kind> u, DMat& dmat, const Wrap& wrap)
	{
		typedef typename matrix_traits<DMat>::value_type T;
		rowwise_fold_impl(shape, u, maximum_folder<T>(), dmat, make_multicol_accessor(u, wrap));
	}

	template<int M, int N, typename Kind, class DMat, class Wrap>
	inline void rowwise_minimum_(const matrix_shape<M, N>& shape, atags::simd<Kind> u, DMat& dmat, const Wrap& wrap)
	{
		typedef typename matrix_traits<DMat>::value_type T;
		rowwise_fold_impl(shape, u, minimum_folder<T>(), dmat, make_multicol_accessor(u, wrap));
	}


	template<int M, int N, typename Kind, class DMat, class TFun, typename... Wrap>
	inline void rowwise_sumx_(const matrix_shape<M, N>& shape, atags::simd<Kind> u, DMat& dmat,
			const TFun& tfun, const Wrap&... wraps)
	{
		typedef typename matrix_traits<DMat>::value_type T;
		rowwise_foldx_impl(shape, u, sum_folder<T>(), dmat, tfun, make_multicol_accessor(u, wraps)...);
	}

	template<int M, int N, typename Kind, class DMat, class TFun, typename... Wrap>
	inline void rowwise_meanx_(const matrix_shape<M, N>& shape, atags::simd<Kind> u, DMat& dmat,
			const TFun& tfun, const Wrap&... wraps)
	{
		typedef typename matrix_traits<DMat>::value_type T;
		rowwise_foldx_impl(shape, u, sum_folder<T>(), dmat, tfun, make_multicol_accessor(u, wraps)...);

		T c = T(1) / T(shape.ncolumns());
		dimension<M> dim(shape.nrows());
		map(math::mul_fun<T>(), atags::simd<Kind>())(dim, out_(dmat), in_(dmat), in_(c, atags::single()));
	}

	template<int M, int N, typename Kind, class DMat, class TFun, typename... Wrap>
	inline void rowwise_maximumx_(const matrix_shape<M, N>& shape, atags::simd<Kind> u, DMat& dmat,
			const TFun& tfun, const Wrap&... wraps)
	{
		typedef typename matrix_traits<DMat>::value_type T;
		rowwise_foldx_impl(shape, u, maximum_folder<T>(), dmat, tfun, make_multicol_accessor(u, wraps)...);
	}

	template<int M, int N, typename Kind, class DMat, class TFun, typename... Wrap>
	inline void rowwise_minimumx_(const matrix_shape<M, N>& shape, atags::simd<Kind> u, DMat& dmat,
			const TFun& tfun, const Wrap&... wraps)
	{
		typedef typename matrix_traits<DMat>::value_type T;
		rowwise_foldx_impl(shape, u, minimum_folder<T>(), dmat, tfun, make_multicol_accessor(u, wraps)...);
	}

} }

#endif 
