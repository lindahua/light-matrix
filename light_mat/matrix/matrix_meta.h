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
	: public std::is_standard_layout<T>
	{ };

	template<class T>
	struct is_supported_matrix_value_type<mask_t<T> > : public true_ { };

	template<class Derived, template<class D, typename T> class Interface>
	struct has_matrix_interface
	{
		typedef Interface<Derived, typename matrix_traits<Derived>::value_type> expect_base;
		static const bool value = std::is_base_of<expect_base, Derived>::value;
	};

	template<class Mat>
	struct is_mat_xpr : public has_matrix_interface<Mat, IMatrixXpr> { };

	template<class Mat>
	struct is_ewise_mat : public has_matrix_interface<Mat, IEWiseMatrix> { };

	template<class Mat>
	struct is_regular_mat : public has_matrix_interface<Mat, IRegularMatrix> { };


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
	: public std::is_same<
	  	  typename matrix_traits<LMat>::domain,
	  	  typename matrix_traits<RMat>::domain> { };

	template<class Mat>
	struct in_cpu_domain
	: public std::is_same<
	  	  typename matrix_traits<Mat>::domain,
	  	  cpu_domain> { };

	template<typename... Mat>
	struct common_domain
	{
		typedef typename common_<typename domain_of<Mat>::type...>::type type;
	};

	/********************************************
	 *
	 *  Type related tools
	 *
	 ********************************************/

	template<class Mat>
	struct qualified_value_type_of
	{
		typedef typename matrix_traits<Mat>::qualified_value_type type;
	};

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
		typedef typename common_<typename value_type_of<Mat>::type...>::type type;
	};


	template<class Mat>
	struct pointer_of
	{
		typedef typename matrix_traits<Mat>::qualified_value_type* type;
	};

	template<class Mat>
	struct reference_of
	{
		typedef typename matrix_traits<Mat>::qualified_value_type& type;
	};

	template<class Mat>
	struct const_pointer_of
	{
		typedef const typename matrix_traits<Mat>::qualified_value_type* type;
	};

	template<class Mat>
	struct const_reference_of
	{
		typedef const typename matrix_traits<Mat>::qualified_value_type& type;
	};

	template<class Mat>
	struct is_readonly : public std::is_const<typename qualified_value_type_of<Mat>::type> { };


	/********************************************
	 *
	 *  Dimension related tools
	 *
	 ********************************************/

	template<class Mat>
	struct nrows : public int_<matrix_traits<Mat>::ct_num_rows> { };

	template<class Mat>
	struct ncols : public int_<matrix_traits<Mat>::ct_num_cols> { };

	template<class Mat>
	struct nelems : public mul_<nrows<Mat>, ncols<Mat> > { };

	template<class Mat>
	struct shape
	{
		typedef matrix_shape<nrows<Mat>::value, ncols<Mat>::value> type;
	};

	template<class Mat>
	struct is_row : public eq_<nrows<Mat>, int_<1> > { };

	template<class Mat>
	struct is_col : public eq_<ncols<Mat>, int_<1> > { };

	template<class Mat>
	struct is_scalar : public and_<is_row<Mat>, is_col<Mat> > { };

	template<class Mat>
	struct is_vector : public or_<is_row<Mat>, is_col<Mat> > { };

	template<class Mat>
	struct has_dynamic_nrows : public eq_<nrows<Mat>, int_<0> > { };

	template<class Mat>
	struct has_dynamic_ncols : public eq_<ncols<Mat>, int_<0> > { };

	template<class Mat>
	struct has_static_nrows : public not_<has_dynamic_nrows<Mat> > { };

	template<class Mat>
	struct has_static_ncols : public not_<has_dynamic_ncols<Mat> > { };

	template<class Mat>
	struct has_static_size
	: public and_<has_static_nrows<Mat>, has_static_ncols<Mat> > { };


	/********************************************
	 *
	 *  Operations on a list of matrices
	 *
	 ********************************************/

	namespace internal
	{
		template<int D1, int D2> struct _common_dim;

		template<int D>
		struct _common_dim<D, 0> : public int_<D> { };

		template<int D>
		struct _common_dim<0, D> : public int_<D> { };

		template<int D>
		struct _common_dim<D, D> : public int_<D> { };

		template<>
		struct _common_dim<0, 0> : public int_<0> { };
	}

	// dimension compatibility

	template<typename W1, typename... WR> struct are_compatible_dims;

	template<typename W1, typename W2>
	struct are_compatible_dims<W1, W2>
	: public bool_<(W1::value == W2::value || W1::value == 0 || W2::value == 0)> { };

	template<typename W1, typename W2, typename... WR>
	struct are_compatible_dims<W1, W2, WR...>
	: and_<are_compatible_dims<W2, WR...>, are_compatible_dims<W1, maximum_<W2, WR...> > > { };

	// common dimension

	template<typename... W> struct common_dim;

	template<typename W>
	struct common_dim<W> : public W { };

	template<typename W1, typename W2>
	struct common_dim<W1, W2> : public internal::_common_dim<W1::value, W2::value> { };

	template<typename W1, typename W2, typename... WR>
	struct common_dim<W1, W2, WR...>
	: public common_dim< common_dim<W1, W2>, common_dim<WR...> > { };


	template<typename... Mats>
	struct common_nrows : public common_dim<nrows<Mats>...> { };

	template<typename... Mats>
	struct common_ncols : public common_dim<ncols<Mats>...> { };

	template<typename... Mats>
	struct common_nelems : public mul_<common_nrows<Mats...>, common_ncols<Mats...> > { };

	template<typename... Mats>
	struct common_shape
	{
		typedef matrix_shape<common_nrows<Mats...>::value, common_ncols<Mats...>::value> type;
	};


	// has compatible size

	template<typename... Mats>
	struct have_compatible_nrows : public are_compatible_dims<nrows<Mats>...> { };

	template<typename... Mats>
	struct have_compatible_ncols : public are_compatible_dims<ncols<Mats>...> { };


	template<typename... Mats>
	struct have_compatible_shapes
	: public and_<
	  have_compatible_nrows<Mats...>,
	  have_compatible_ncols<Mats...> > { };


	template<class Mat>
	struct sq_dim
	: public common_dim<nrows<Mat>, ncols<Mat> > { };


	/********************************************
	 *
	 *  Layout attributes
	 *
	 ********************************************/

	template<class Mat>
	struct is_contiguous
	{
		typedef typename matrix_traits<Mat>::layout_type layout_type;
		static const bool value = layout_traits<layout_type>::ct_is_contiguous;
	};

	template<class Mat>
	struct is_percol_contiguous
	{
		typedef typename matrix_traits<Mat>::layout_type layout_type;
		static const bool value = layout_traits<layout_type>::ct_is_percol_contiguous;
	};


	template<typename... Mat> struct contiguousness;

	template<class A>
	struct contiguousness<A>
	: public select_<
	  	  is_contiguous<A>, cont_level::whole,
	  	  is_percol_contiguous<A>, cont_level::percol,
	  	  otherwise_, cont_level::none> { };


	template<class A, class B>
	struct contiguousness<A, B>
	: public select_<
	  	  and_<is_contiguous<A>, is_contiguous<B> >, cont_level::whole,
	  	  and_<is_percol_contiguous<A>, is_percol_contiguous<B> >, cont_level::percol,
	  	  otherwise_, cont_level::none> { };


	template<class Mat>
	struct supports_linear_index
	: public or_<is_contiguous<Mat>, is_vector<Mat> > { };


} }

#endif /* MATRIX_META_H_ */
