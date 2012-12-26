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
#include <type_traits>

namespace lmat { namespace meta {

	/********************************************
	 *
	 *  Concept test tools
	 *
	 ********************************************/

	template<class T>
	struct is_supported_matrix_value_type
	{
		static const bool value = std::is_pod<T>::value;
	};

	template<class T>
	struct is_supported_matrix_value_type<mask_t<T> >
	{
		static const bool value = true;
	};

	template<class Derived, template<class D, typename T> class Interface>
	struct has_matrix_interface
	{
		typedef Interface<Derived, typename matrix_traits<Derived>::value_type> expect_base;
		static const bool value = std::is_base_of<expect_base, Derived>::value;
	};

	template<class Mat>
	struct is_mat_xpr
	{
		static const bool value = has_matrix_interface<Mat, IMatrixXpr>::value;
	};

	template<class Mat>
	struct is_regular_mat
	{
		static const bool value = has_matrix_interface<Mat, IRegularMatrix>::value;
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
		static const bool value = std::is_same<
				typename matrix_traits<LMat>::domain,
				typename matrix_traits<RMat>::domain>::value;
	};


	template<class Mat>
	struct in_cpu_domain
	{
		static const bool value = std::is_same<
				typename matrix_traits<Mat>::domain,
				cpu_domain>::value;
	};

	template<typename... Mat>
	struct common_domain
	{
		typedef typename common_type<typename domain_of<Mat>::type...>::type type;
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
		static const bool value = std::is_same<
			typename matrix_traits<LMat>::value_type,
			typename matrix_traits<RMat>::value_type>::value;
	};

	template<typename... Mat>
	struct common_value_type
	{
		typedef typename common_type<typename value_type_of<Mat>::type...>::type type;
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

		template<template<class X> class Fun, class... Mats>
		struct common_dim_reduce;

		template<template<class X> class Fun, class Mat>
		struct common_dim_reduce<Fun, Mat>
		{
			static const int value = Fun<Mat>::value;
		};

		template<template<class X> class Fun, class M0, class... Mats>
		struct common_dim_reduce<Fun, M0, Mats...>
		{
			static const int value = common_dim_c<
					Fun<M0>::value,
					common_dim_reduce<Fun, Mats...>::value>::value;
		};
	}

	template<int M, int N>
	struct common_dim
	{
		static const int _maybe_value = internal::common_dim_c<M, N>::value;
		static const int value = enable_int_if_c<(_maybe_value >= 0), _maybe_value>::value;
	};


	// common_ctdim


	// common rows

	template<typename... Mats>
	struct common_nrows
	{
		static const int _maybe_value = internal::common_dim_reduce<nrows, Mats...>::value;
		static const int value = enable_int_if_c<(_maybe_value >= 0), _maybe_value>::value;
	};

	template<typename... Mats>
	struct common_ncols
	{
		static const int _maybe_value = internal::common_dim_reduce<ncols, Mats...>::value;
		static const int value = enable_int_if_c<(_maybe_value >= 0), _maybe_value>::value;
	};

	template<typename... Mats>
	struct common_nelems
	{
		static const int value = common_nrows<Mats...>::value * common_ncols<Mats...>::value;
	};

	template<typename... Mats>
	struct common_shape
	{
		typedef matrix_shape<common_nrows<Mats...>::value, common_ncols<Mats...>::value> type;
	};



	// has compatible size

	template<typename... Mats>
	struct have_compatible_nrows
	{
		static const int _maybe_value = internal::common_dim_reduce<nrows, Mats...>::value;
		static const bool value = _maybe_value >= 0;
	};


	template<typename... Mats>
	struct have_compatible_ncols
	{
		static const int _maybe_value = internal::common_dim_reduce<ncols, Mats...>::value;
		static const bool value = _maybe_value >= 0;
	};

	template<typename... Mats>
	struct have_compatible_shapes
	{
		static const bool value =
				have_compatible_nrows<Mats...>::value &&
				have_compatible_ncols<Mats...>::value;
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
	struct continuous_level
	{
		typedef typename if_<is_continuous<Mat>, 		cont_level::whole,
				typename if_<is_percol_continuous<Mat>, cont_level::percol,
														cont_level::none
				>::type >::type type;
	};

	template<class LMat, class RMat>
	struct binary_continuous_level
	{
		static const bool is_cont = is_continuous<LMat>::value && is_continuous<RMat>::value;
		static const bool is_pc_cont =
				is_percol_continuous<LMat>::value && is_percol_continuous<RMat>::value;

		typedef typename if_c<is_cont, cont_level::whole,
				typename if_c<is_pc_cont, cont_level::percol, cont_level::none
				>::type >::type type;
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
				internal::_is_continuous_ex_<Mat, is_regular_mat<Mat>::value>::value;
	};

	template<class Mat>
	struct is_percol_continuous_ex
	{
		static const bool value =
				internal::_is_percol_continuous_ex_<Mat, is_regular_mat<Mat>::value>::value;
	};

	template<class Mat>
	struct supports_linear_index_ex
	{
		static const bool value =
				internal::_supports_linear_index_ex_<Mat, is_regular_mat<Mat>::value>::value;
	};



	/********************************************
	 *
	 *  Evaluation
	 *
	 ********************************************/

	template<class SExpr>
	struct supports_ewise_access
	{
		static const bool value = is_regular_mat<SExpr>::value;
	};


	template<class SExpr, class DMat>
	struct is_mat_assignable
	{
		static const bool value =
				is_regular_mat<DMat>::value &&
				!is_readonly<DMat>::value &&
				is_mat_xpr<SExpr>::value &&
				in_same_domain<SExpr, DMat>::value &&
				have_same_value_type<SExpr, DMat>::value &&
				have_compatible_shapes<SExpr, DMat>::value;
	};

} }

#endif /* MATRIX_META_H_ */
