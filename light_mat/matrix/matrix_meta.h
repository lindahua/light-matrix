/**
 * @file matrix_meta.h
 *
 * Meta programming tools for Matrix
 *
 * @author Dahua Lin
 */

#ifdef _MSC_VER
#pragma once
#endif

#ifndef LIGHTMAT_MATRIX_META_H_
#define LIGHTMAT_MATRIX_META_H_

#include <light_mat/matrix/matrix_fwd.h>

namespace lmat { namespace meta {

	/********************************************
	 *
	 *  Concept test tools
	 *
	 ********************************************/

	template<class T>
	struct is_supported_matrix_value_type
	{
		static const bool value = lmat::is_pod<T>::value;
	};

	template<class Derived, template<class D, typename T> class Interface>
	struct has_matrix_interface
	{
		typedef Interface<Derived, typename matrix_traits<Derived>::value_type> expect_base;
		static const bool value = is_base_of<expect_base, Derived>::value;
	};

	template<class Mat>
	struct is_mat_xpr
	{
		static const bool value = has_matrix_interface<Mat, IMatrixXpr>::value;
	};

	template<class Mat>
	struct is_dense_mat
	{
		static const bool value = has_matrix_interface<Mat, IDenseMatrix>::value;
	};


	/********************************************
	 *
	 *  Domain related tools
	 *
	 ********************************************/

	template<class Mat>
	struct domain_of
	{
		typedef typename matrix_traits<Mat>::domain type;
	};

	template<class LMat, class RMat>
	struct in_same_domain
	{
		static const bool value = is_same<
				typename matrix_traits<LMat>::domain,
				typename matrix_traits<RMat>::domain>::value;
	};


	template<class Mat>
	struct in_cpu_domain
	{
		static const bool value = is_same<
				typename matrix_traits<Mat>::domain,
				cpu_domain>::value;
	};

	template<class MatList>
	struct common_domain
	{
		typedef typename common_<MatList, domain_of>::type type;
	};

	/********************************************
	 *
	 *  Type related tools
	 *
	 ********************************************/

	template<class Mat>
	struct value_type_of
	{
		typedef typename matrix_traits<Mat>::value_type type;
	};

	template<class LMat, class RMat>
	struct have_same_value_type
	{
		static const bool value = is_same<
			typename matrix_traits<LMat>::value_type,
			typename matrix_traits<RMat>::value_type>::value;
	};

	template<class MatList>
	struct common_value_type
	{
		typedef typename common_<MatList, value_type_of>::type type;
	};

	/********************************************
	 *
	 *  Dimension related tools
	 *
	 ********************************************/

	template<class Mat>
	struct nrows
	{
		static const int value = matrix_traits<Mat>::ct_num_rows;
	};

	template<class Mat>
	struct ncols
	{
		static const int value = matrix_traits<Mat>::ct_num_cols;
	};

	template<class Mat>
	struct nelems
	{
		static const int value = nrows<Mat>::value * ncols<Mat>::value;
	};

	template<class Mat>
	struct shape
	{
		typedef matrix_shape<nrows<Mat>::value, ncols<Mat>::value> type;
	};

	template<class Mat>
	struct is_row
	{
		static const bool value = nrows<Mat>::value == 1;
	};

	template<class Mat>
	struct is_col
	{
		static const bool value = ncols<Mat>::value == 1;
	};

	template<class Mat>
	struct is_scalar
	{
		static const bool value = is_row<Mat>::value && is_col<Mat>::value;
	};

	template<class Mat>
	struct is_vector
	{
		static const bool value = is_row<Mat>::value || is_col<Mat>::value;
	};

	template<class Mat>
	struct has_static_nrows
	{
		static const bool value = nrows<Mat>::value > 0;
	};

	template<class Mat>
	struct has_static_ncols
	{
		static const bool value = ncols<Mat>::value > 0;
	};

	template<class Mat>
	struct has_dynamic_nrows
	{
		static const bool value = nrows<Mat>::value == 0;
	};

	template<class Mat>
	struct has_dynamic_ncols
	{
		static const bool value = ncols<Mat>::value == 0;
	};


	template<class Mat>
	struct has_static_size
	{
		static const bool value = has_static_nrows<Mat>::value && has_static_ncols<Mat>::value;
	};


	/********************************************
	 *
	 *  Operations on a list of matrices
	 *
	 ********************************************/

	// compatible ct_dim

	namespace internal
	{
		template<int M, int N>
		struct are_compatible_dims_c
		{
			static const bool value = M >= 0 && N >= 0 && (M == N || M == 0 || N == 0);
		};

		template<int M, int N>
		struct common_dim_c
		{
			static const int value =
					are_compatible_dims_c<M, N>::value ? (M == 0 ? N : M) : -1;
		};

		template<class DimList>
		struct common_dim_reduce
		{
			static const int value = fold_ints_<common_dim_c, DimList>::value;
		};
	}


	template<class DimList>
	struct are_compatble_dims
	{
		static const bool value = internal::common_dim_reduce<DimList>::value >= 0;
	};

	template<class DimList>
	struct common_dim
	{
		static const int _maybe_value = internal::common_dim_reduce<DimList>::value;
		static const int value = enable_int_if_c<(_maybe_value >= 0), _maybe_value>::value;
	};


	// common_ctdim


	// common rows

	template<class ArgList>
	struct common_nrows
	{
		static const int value = common_dim<
				typename types_to_ints_<nrows, ArgList>::type>::value;
	};

	template<class ArgList>
	struct common_ncols
	{
		static const int value = common_dim<
				typename types_to_ints_<ncols, ArgList>::type>::value;
	};

	template<class ArgList>
	struct common_nelems
	{
		static const int value = common_nrows<ArgList>::value * common_ncols<ArgList>::value;
	};

	template<class ArgList>
	struct common_shape
	{
		typedef matrix_shape<common_nrows<ArgList>::value, common_ncols<ArgList>::value> type;
	};



	// has compatible size

	template<class ArgList>
	struct have_compatible_nrows
	{
		static const bool value =
				internal::common_dim_reduce<
					typename types_to_ints_<nrows, ArgList>::type
				>::value >= 0;
	};


	template<class ArgList>
	struct have_compatible_ncols
	{
		static const bool value =
				internal::common_dim_reduce<
					typename types_to_ints_<ncols, ArgList>::type
				>::value >= 0;
	};

	template<class ArgList>
	struct have_compatible_shapes
	{
		static const bool value =
				have_compatible_nrows<ArgList>::value &&
				have_compatible_ncols<ArgList>::value;
	};


	/********************************************
	 *
	 *  Matrix access and manipulation
	 *
	 ********************************************/

	template<class Mat>
	struct is_readonly
	{
		static const bool value = matrix_traits<Mat>::is_readonly;
	};


	/********************************************
	 *
	 *  Layout attributes
	 *
	 ********************************************/

	template<class Mat>
	struct is_continuous
	{
		typedef typename matrix_traits<Mat>::layout_type layout_type;
		static const bool value = layout_traits<layout_type>::ct_is_continuous;
	};

	template<class Mat>
	struct is_percol_continuous
	{
		typedef typename matrix_traits<Mat>::layout_type layout_type;
		static const bool value = layout_traits<layout_type>::ct_is_percol_continuous;
	};

	template<class Mat>
	struct supports_linear_index
	{
		static const bool value = is_continuous<Mat>::value || is_vector<Mat>::value;
	};


	namespace internal
	{
		template<class Mat, bool IsDense>
		struct _is_continuous_ex_
		{
			static const bool value = is_continuous<Mat>::value;
		};

		template<class Mat, bool IsDense>
		struct _is_percol_continuous_ex_
		{
			static const bool value = is_percol_continuous<Mat>::value;
		};

		template<class Mat, bool IsDense>
		struct _supports_linear_index_ex_
		{
			static const bool value = supports_linear_index<Mat>::value;
		};

		template<class Mat>
		struct _is_continuous_ex_<Mat, false>
		{
			static const bool value = false;
		};

		template<class Mat>
		struct _is_percol_continuous_ex_<Mat, false>
		{
			static const bool value = false;
		};

		template<class Mat>
		struct _supports_linear_index_ex_<Mat, false>
		{
			static const bool value = false;
		};
	}

	template<class Mat>
	struct is_continuous_ex
	{
		static const bool value =
				internal::_is_continuous_ex_<Mat, is_dense_mat<Mat>::value>::value;
	};

	template<class Mat>
	struct is_percol_continuous_ex
	{
		static const bool value =
				internal::_is_percol_continuous_ex_<Mat, is_dense_mat<Mat>::value>::value;
	};

	template<class Mat>
	struct supports_linear_index_ex
	{
		static const bool value =
				internal::_supports_linear_index_ex_<Mat, is_dense_mat<Mat>::value>::value;
	};



	/********************************************
	 *
	 *  Evaluation
	 *
	 ********************************************/

	template<class SExpr, class DMat>
	struct is_mat_assignable
	{
		typedef LMAT_TYPELIST_2(SExpr, DMat) Lst;

		static const bool value =
				is_dense_mat<DMat>::value &&
				!is_readonly<DMat>::value &&
				is_mat_xpr<SExpr>::value &&
				in_same_domain<SExpr, DMat>::value &&
				have_same_value_type<SExpr, DMat>::value &&
				have_compatible_shapes<Lst>::value;
	};

} }

#endif /* MATRIX_META_H_ */
