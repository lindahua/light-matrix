/**
 * @file dense_accessors.h
 *
 * Accessor classes for dense matrices
 *
 * @author Dahua Lin
 */

#ifdef _MSC_VER
#pragma once
#endif

#ifndef LIGHTMAT_DENSE_ACCESSORS_H_
#define LIGHTMAT_DENSE_ACCESSORS_H_

#include <light_mat/matexpr/matrix_access_eval.h>
#include <light_mat/matrix/dense_matrix.h>

namespace lmat
{
	// forward declarations

	template<class Mat, typename T> class dense_linear_scalar_accessor;
	template<class Mat, typename T> class dense_percol_scalar_accessor;
	template<class Mat, typename T> class grid_percol_scalar_accessor;
	template<class Mat, typename T> class cached_linear_scalar_accessor;
	template<class Mat, typename T> class cached_percol_scalar_accessor;


	/********************************************
	 *
	 *  accessor classes
	 *
	 ********************************************/

	// continuous linear accessing

	template<class Mat, typename T>
	class dense_linear_scalar_accessor
			: public ILinearMatrixScalarAccessor<dense_linear_scalar_accessor<Mat, T>, T>
	{
#ifdef LMAT_USE_STATIC_ASSERT
		static_assert(meta::supports_linear_index<Mat>::value,
				"Mat should support linear indexing");
#endif

	public:
		LMAT_ENSURE_INLINE
		dense_linear_scalar_accessor(const Mat& mat)
		: m_mat(mat)
		{ }

		LMAT_ENSURE_INLINE
		T get_scalar(const index_t i) const
		{
			return m_mat[i];
		}

	private:
		const Mat& m_mat;
	};


	// dense percol accessing

	template<class Mat, typename T>
	struct percol_macc_state_map<dense_percol_scalar_accessor<Mat, T> >
	{
		typedef const T* type;
	};

	template<class Mat, typename T>
	class dense_percol_scalar_accessor
			: public IPerColMatrixScalarAccessor<dense_percol_scalar_accessor<Mat, T>, T>
	{
#ifdef LMAT_USE_STATIC_ASSERT
		static_assert(meta::is_percol_continuous<Mat>::value,
				"Mat should be a compile-time percol-continuous matrix");
#endif

	public:
		LMAT_ENSURE_INLINE
		dense_percol_scalar_accessor(const Mat& mat)
		: m_data(mat.ptr_data())
		, m_stride(mat.col_stride())
		{ }

		LMAT_ENSURE_INLINE
		T get_scalar(const index_t i, const T* s) const
		{
			return s[i];
		}

		LMAT_ENSURE_INLINE
		const T* col_state(const index_t j) const
		{
			return m_data + j * m_stride;
		}

	private:
		const T *m_data;
		const index_t m_stride;
	};


	// grid percol accessing

	template<class Mat, typename T>
	struct percol_macc_state_map<grid_percol_scalar_accessor<Mat, T> >
	{
		typedef const T* type;
	};

	template<class Mat, typename T>
	class grid_percol_scalar_accessor
			: public IPerColMatrixScalarAccessor<grid_percol_scalar_accessor<Mat, T>, T>
	{
	public:
		LMAT_ENSURE_INLINE
		grid_percol_scalar_accessor(const Mat& mat)
		: m_data(mat.ptr_data())
		, m_step(mat.row_stride())
		, m_stride(mat.col_stride())
		{ }

		LMAT_ENSURE_INLINE
		T get_scalar(const index_t i, const T* s) const
		{
			return s[i * m_step];
		}

		LMAT_ENSURE_INLINE
		const T* col_state(const index_t j) const
		{
			return m_data + j * m_stride;
		}

	private:
		const T *m_data;
		const index_t m_step;
		const index_t m_stride;
	};



	// cache linear accessing

	template<class Mat, typename T>
	class cached_linear_scalar_accessor
			: public ILinearMatrixScalarAccessor<cached_linear_scalar_accessor<Mat, T>, T>
	{
	public:
		LMAT_ENSURE_INLINE
		cached_linear_scalar_accessor(const Mat& mat)
		: m_cache(mat)
		, m_data(m_cache.ptr_data())
		{ }

		LMAT_ENSURE_INLINE
		T get_scalar(const index_t i) const
		{
			return m_data[i];
		}

	private:
		dense_matrix<T, meta::nrows<Mat>::value, meta::ncols<Mat>::value> m_cache;
		const T *m_data;
	};


	// cache percol accessing

	template<class Mat, typename T>
	struct percol_macc_state_map<cached_percol_scalar_accessor<Mat, T> >
	{
		typedef const T* type;
	};

	template<class Mat, typename T>
	class cached_percol_scalar_accessor
			: public IPerColMatrixScalarAccessor<cached_percol_scalar_accessor<Mat, T>, T>
	{
	public:
		LMAT_ENSURE_INLINE
		cached_percol_scalar_accessor(const Mat& mat)
		: m_cache(mat)
		{ }

		LMAT_ENSURE_INLINE
		T get_scalar(const index_t i, const T* s) const
		{
			return s[i];
		}

		LMAT_ENSURE_INLINE
		const T* col_state(const index_t j) const
		{
			return m_cache.ptr_col(j);
		}

	private:
		dense_matrix<T, meta::nrows<Mat>::value, meta::ncols<Mat>::value> m_cache;
	};


	/********************************************
	 *
	 *  accessor maps
	 *
	 ********************************************/

	namespace internal
	{
		template<class Xpr, bool IsDense>
		struct macc_support_dense_linear;

		template<class Xpr>
		struct macc_support_dense_linear<Xpr, true>
		{
			static const bool value = meta::supports_linear_index<Xpr>::value;
		};

		template<class Xpr>
		struct macc_support_dense_linear<Xpr, false>
		{
			static const bool value = false;
		};

		template<class Xpr, bool IsDense>
		struct macc_support_dense_percol;

		template<class Xpr>
		struct macc_support_dense_percol<Xpr, true>
		{
			static const bool value = meta::is_percol_continuous<Xpr>::value;
		};

		template<class Xpr>
		struct macc_support_dense_percol<Xpr, false>
		{
			static const bool value = false;
		};
	}

	template<class SExpr, typename AccCate, typename KerCate>
	struct generic_macc_accessor_map;

	template<class SExpr>
	struct generic_macc_accessor_map<SExpr, linear_macc, scalar_ker>
	{
		typedef typename matrix_traits<SExpr>::value_type T;

		static const bool support_dense_linear =
				internal::macc_support_dense_linear<SExpr, meta::is_dense_mat<SExpr>::value>::value;

		typedef typename
				meta::if_c<support_dense_linear,
					dense_linear_scalar_accessor<SExpr, T>,
					cached_linear_scalar_accessor<SExpr, T>
				>::type type;
	};


	template<class SExpr>
	struct generic_macc_accessor_map<SExpr, percol_macc, scalar_ker>
	{
		typedef typename matrix_traits<SExpr>::value_type T;

		static const bool support_dense_percol =
				internal::macc_support_dense_percol<SExpr, meta::is_dense_mat<SExpr>::value>::value;

		typedef typename
				meta::if_c<meta::is_dense_mat<SExpr>::value,
					typename
					meta::if_c<support_dense_percol,
						dense_percol_scalar_accessor<SExpr, T>,
						grid_percol_scalar_accessor<SExpr, T>
					>::type,
					cached_percol_scalar_accessor<SExpr, T>
				>::type type;
	};


	template<class SExpr, typename AccCate, typename KerCate>
	struct macc_accessor_map
	{
		typedef typename generic_macc_accessor_map<SExpr, AccCate, KerCate>::type type;
	};


	/********************************************
	 *
	 *  cost model
	 *
	 ********************************************/

	template<class Xpr, typename AccCate, typename KerCate>
	struct generic_macc_cost;

	template<class Xpr, typename KerCate>
	struct generic_macc_cost<Xpr, linear_macc, KerCate>
	{
		static const bool support_dense_linear =
				internal::macc_support_dense_linear<Xpr, meta::is_dense_mat<Xpr>::value>::value;

		static const int value = support_dense_linear ? 0 : MACC_CACHE_COST;
	};

	template<class Xpr, typename KerCate>
	struct generic_macc_cost<Xpr, percol_macc, KerCate>
	{
		static const bool is_dense = meta::is_dense_mat<Xpr>::value;

		static const int value = is_dense ? 0 : MACC_CACHE_COST;
	};


	template<class Xpr, typename AccCate, typename KerCate>
	struct macc_cost
	{
		static const int value =
				generic_macc_cost<Xpr, AccCate, KerCate>::value;
	};

	template<typename T, typename AccCate, typename KerCate>
	struct macc_cost<scalar_expr<T>, AccCate, KerCate>
	{
		static const int value = 0;
	};


}

#endif
