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
#include <light_mat/matexpr/mat_arith.h>

namespace lmat { namespace internal {


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
	 *  vector-wise reduction implementation
	 *
	 ********************************************/

	// column wise reduction

	template<index_t CM, index_t CN, class FoldKernel, typename... Args>
	class colwise_fold_getter;

	template<index_t CM, index_t CN, class FoldKernel, typename Arg1>
	class colwise_fold_getter<CM, CN, FoldKernel, Arg1>
	{
	public:
		typedef fold_policy<FoldKernel, matrix_shape<CM, 1>, Arg1> pmap;
		typedef typename pmap::unit U;

		LMAT_ENSURE_INLINE
		colwise_fold_getter(const FoldKernel& kernel, const matrix_shape<CM, CN>& shape,
				const Arg1& arg1)
		: m_kernel(kernel)
		, m_coldim(shape.nrows())
		, m_rd1(make_multicol_accessor(U(), in_(arg1)))
		{ }

		LMAT_ENSURE_INLINE
		typename FoldKernel::accumulated_type operator[] (index_t j) const
		{
			return linear_fold_impl(m_coldim, U(), m_kernel, m_rd1.col(j));
		}

	private:
		const FoldKernel& m_kernel;
		dimension<CM> m_coldim;
		typename multicol_reader_map<Arg1, U>::type m_rd1;
	};


	template<index_t CM, index_t CN, class FoldKernel, typename Arg1, typename Arg2>
	class colwise_fold_getter<CM, CN, FoldKernel, Arg1, Arg2>
	{
	public:
		typedef fold_policy<FoldKernel, matrix_shape<CM, 1>, Arg1, Arg2> pmap;
		typedef typename pmap::unit U;

		LMAT_ENSURE_INLINE
		colwise_fold_getter(const FoldKernel& kernel, const matrix_shape<CM, CN>& shape,
				const Arg1& arg1, const Arg2& arg2)
		: m_kernel(kernel)
		, m_coldim(shape.nrows())
		, m_rd1(make_multicol_accessor(U(), in_(arg1)))
		, m_rd2(make_multicol_accessor(U(), in_(arg2)))
		{ }

		LMAT_ENSURE_INLINE
		typename FoldKernel::accumulated_type operator[] (index_t j) const
		{
			return linear_fold_impl(m_coldim, U(), m_kernel, m_rd1.col(j), m_rd2.col(j));
		}

	private:
		const FoldKernel& m_kernel;
		dimension<CM> m_coldim;
		typename multicol_reader_map<Arg1, U>::type m_rd1;
		typename multicol_reader_map<Arg2, U>::type m_rd2;
	};


	template<index_t CM, index_t CN, class FoldKernel, typename Arg1, typename Arg2, typename Arg3>
	class colwise_fold_getter<CM, CN, FoldKernel, Arg1, Arg2, Arg3>
	{
	public:
		typedef fold_policy<FoldKernel, matrix_shape<CM, 1>, Arg1, Arg2, Arg3> pmap;
		typedef typename pmap::unit U;

		LMAT_ENSURE_INLINE
		colwise_fold_getter(const FoldKernel& kernel, const matrix_shape<CM, CN>& shape,
				const Arg1& arg1, const Arg2& arg2, const Arg3& arg3)
		: m_kernel(kernel)
		, m_coldim(shape.nrows())
		, m_rd1(make_multicol_accessor(U(), in_(arg1)))
		, m_rd2(make_multicol_accessor(U(), in_(arg2)))
		, m_rd3(make_multicol_accessor(U(), in_(arg3)))
		{ }

		LMAT_ENSURE_INLINE
		typename FoldKernel::accumulated_type operator[] (index_t j) const
		{
			return linear_fold_impl(m_coldim, U(), m_kernel, m_rd1.col(j), m_rd2.col(j), m_rd3.col(j));
		}

	private:
		const FoldKernel& m_kernel;
		dimension<CM> m_coldim;
		typename multicol_reader_map<Arg1, U>::type m_rd1;
		typename multicol_reader_map<Arg2, U>::type m_rd2;
		typename multicol_reader_map<Arg3, U>::type m_rd3;
	};


	template<index_t CM, index_t CN, class FoldKernel, typename Arg1, typename T1>
	LMAT_ENSURE_INLINE
	inline colwise_fold_getter<CM, CN, FoldKernel, Arg1>
	make_colwise_fold_getter(const FoldKernel& kernel, const matrix_shape<CM, CN>& shape,
			const IEWiseMatrix<Arg1, T1>& arg1)
	{
		typedef colwise_fold_getter<CM, CN, FoldKernel, Arg1> type;
		return type(kernel, shape, arg1.derived());
	}

	template<index_t CM, index_t CN, class FoldKernel, typename Arg1, typename T1, typename Arg2, typename T2>
	LMAT_ENSURE_INLINE
	inline colwise_fold_getter<CM, CN, FoldKernel, Arg1, Arg2>
	make_colwise_fold_getter(const FoldKernel& kernel, const matrix_shape<CM, CN>& shape,
			const IEWiseMatrix<Arg1, T1>& arg1,
			const IEWiseMatrix<Arg2, T2>& arg2)
	{
		typedef colwise_fold_getter<CM, CN, FoldKernel, Arg1, Arg2> type;
		return type(kernel, shape, arg1.derived(), arg2.derived());
	}

	template<index_t CM, index_t CN, class FoldKernel,
		typename Arg1, typename T1, typename Arg2, typename T2, typename Arg3, typename T3>
	LMAT_ENSURE_INLINE
	inline colwise_fold_getter<CM, CN, FoldKernel, Arg1, Arg2, Arg3>
	make_colwise_fold_getter(const FoldKernel& kernel, const matrix_shape<CM, CN>& shape,
			const IEWiseMatrix<Arg1, T1>& arg1,
			const IEWiseMatrix<Arg2, T2>& arg2,
			const IEWiseMatrix<Arg3, T2>& arg3)
	{
		typedef colwise_fold_getter<CM, CN, FoldKernel, Arg1, Arg2, Arg3> type;
		return type(kernel, shape, arg1.derived(), arg2.derived(), arg3.derived());
	}



	template<index_t CM, index_t CN, class FoldKernel, typename T, class DMat, class TExpr>
	inline void colwise_fold_impl(const matrix_shape<CM, CN>& shape,
			const FoldKernel& kernel, IRegularMatrix<DMat, T>& dmat, const IEWiseMatrix<TExpr, T>& texpr)
	{
		const index_t n = shape.ncolumns();
		LMAT_CHECK_DIMS( n == dmat.nelems() )

		auto g = make_colwise_fold_getter(kernel, shape, texpr);

		DMat& d_ = dmat.derived();
		for (index_t j = 0; j < n; ++j)
		{
			d_[j] = g[j];
		}
	}

	// row wise reduction

	template<index_t CM, index_t CN, class FoldKernel, typename T, class DMat, class TExpr>
	inline void rowwise_fold_impl(const matrix_shape<CM, CN>& shape,
			const FoldKernel& kernel, IRegularMatrix<DMat, T>& dmat, const IEWiseMatrix<TExpr, T>& texpr)
	{
		dimension<CM> col_dim(shape.nrows());
		const index_t n = shape.ncolumns();

		LMAT_CHECK_DIMS( col_dim.value() == dmat.nelems() )

		typedef preferred_macc_policy<matrix_shape<CM, 1>, FoldKernel, DMat, TExpr> pmap;
		// static_assert(pmap::use_simd, "should use SIMD here");

		typedef typename pmap::unit U;

		auto a = make_vec_accessor(U(), in_out_(dmat));
		auto rd = make_multicol_accessor(U(), in_(texpr));

		internal::_linear_ewise_eval(col_dim, U(), copy_kernel<T>(), rd.col(0), a);

		for (index_t j = 1; j < n; ++j)
		{
			internal::_linear_ewise_eval(col_dim, U(), kernel, a, rd.col(j));
		}
	}


} }

#endif 
