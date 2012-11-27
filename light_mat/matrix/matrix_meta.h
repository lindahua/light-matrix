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

namespace lmat
{

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
	struct is_mat_view
	{
		static const bool value = has_matrix_interface<Mat, IMatrixView>::value;
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

	template<class LMat, class RMat>
	struct has_same_domain
	{
		static const bool value = is_same<
				typename matrix_traits<LMat>::domain,
				typename matrix_traits<RMat>::domain>::value;
	};


	template<class Mat>
	struct has_cpu_domain
	{
		static const bool value = is_same<
				typename matrix_traits<Mat>::domain,
				cpu_domain>::value;
	};

	template<class LMat, class RMat>
	struct binary_domain
	{
#ifdef LMAT_USE_STATIC_ASSERT
		static_assert( has_same_domain<LMat, RMat>::value,
				"LMat and RMat have different domains." );
#endif

		typedef typename matrix_traits<LMat>::domain type;
	};


	/********************************************
	 *
	 *  Type related tools
	 *
	 ********************************************/

	template<class LMat, class RMat>
	struct has_same_value_type
	{
		static const bool value = is_same<
			typename matrix_traits<LMat>::value_type,
			typename matrix_traits<RMat>::value_type>::value;
	};


	template<class LMat, class RMat>
	struct binary_value_type
	{
#ifdef LMAT_USE_STATIC_ASSERT
		static_assert( has_same_value_type<LMat, RMat>::value,
				"LMat and RMat have different value_types." );
#endif

		typedef typename matrix_traits<LMat>::value_type type;
	};


	/********************************************
	 *
	 *  Dimension related tools
	 *
	 ********************************************/

	template<class Mat>
	struct ct_rows
	{
		static const int value = matrix_traits<Mat>::ct_num_rows;
	};

	template<class Mat>
	struct ct_cols
	{
		static const int value = matrix_traits<Mat>::ct_num_cols;
	};

	template<class Mat>
	struct ct_size
	{
		static const int value = ct_rows<Mat>::value * ct_cols<Mat>::value;
	};

	template<class Mat>
	struct matrix_shape_type
	{
		typedef matrix_shape<ct_rows<Mat>::value, ct_cols<Mat>::value> type;
	};

	template<int N1, int N2>
	struct is_compatible_ctdim
	{
		static const bool value = N1 >= 0 && N2 >= 0 && (N1 == 0 || N2 == 0 || N1 == N2);
	};

	template<int N1, int N2>
	struct binary_ctdim
	{
		static const int _maybe_value_ = N1 > N2 ? N1 : N2;
		static const int value = enable_int_if<is_compatible_ctdim<N1, N2>, _maybe_value_>::value;
	};


	template<class Mat1, class Mat2>
	struct binary_ct_rows
	{
		static const int value = binary_ctdim<ct_rows<Mat1>::value, ct_rows<Mat2>::value>::value;
	};

	template<class Mat1, class Mat2>
	struct binary_ct_cols
	{
		static const int value = binary_ctdim<ct_cols<Mat1>::value, ct_cols<Mat2>::value>::value;
	};

	template<class Mat1, class Mat2>
	struct binary_ct_size
	{
		static const int value = binary_ct_rows<Mat1, Mat2>::value * binary_ct_cols<Mat1, Mat2>::value;
	};

	template<class Mat1, class Mat2>
	struct binary_shape_type
	{
		static const int ct_nrows = binary_ct_rows<Mat1, Mat2>::value;
		static const int ct_ncols = binary_ct_cols<Mat1, Mat2>::value;
		typedef matrix_shape<ct_nrows, ct_ncols> type;
	};


	template<class LMat, class RMat>
	struct has_compatible_ct_size
	{
		static const bool value =
				is_compatible_ctdim<ct_rows<LMat>::value, ct_rows<RMat>::value>::value &&
				is_compatible_ctdim<ct_cols<LMat>::value, ct_cols<RMat>::value>::value;
	};

	template<class Mat>
	struct ct_is_row
	{
		static const bool value = ct_rows<Mat>::value == 1;
	};

	template<class Mat>
	struct ct_is_col
	{
		static const bool value = ct_cols<Mat>::value == 1;
	};

	template<class Mat>
	struct ct_is_vector
	{
		static const bool value = ct_is_row<Mat>::value || ct_is_col<Mat>::value;
	};

	template<class Mat>
	struct has_static_nrows
	{
		static const bool value = ct_rows<Mat>::value > 0;
	};

	template<class Mat>
	struct has_static_ncols
	{
		static const bool value = ct_cols<Mat>::value > 0;
	};

	template<class Mat>
	struct has_dynamic_nrows
	{
		static const bool value = ct_rows<Mat>::value == 0;
	};

	template<class Mat>
	struct has_dynamic_ncols
	{
		static const bool value = ct_cols<Mat>::value == 0;
	};


	template<class Mat>
	struct has_static_size
	{
		static const bool value = has_static_nrows<Mat>::value && has_static_ncols<Mat>::value;
	};

	template<class LMat, class RMat>
	struct ct_has_same_nrows
	{
		static const bool value =
				ct_rows<LMat>::value > 0 &&
				ct_rows<LMat>::value == ct_rows<RMat>::value;
	};

	template<class LMat, class RMat>
	struct ct_has_same_ncols
	{
		static const bool value =
				ct_cols<LMat>::value > 0 &&
				ct_cols<LMat>::value == ct_cols<RMat>::value;
	};


	/********************************************
	 *
	 *  Matrix access and manipulation
	 *
	 ********************************************/

	template<class Mat>
	struct is_readonly_mat
	{
		static const bool value = matrix_traits<Mat>::is_readonly;
	};

	template<class Mat>
	struct mat_access
	{
		typedef typename matrix_traits<Mat>::value_type value_type;
		typedef typename if_<is_readonly_mat<Mat>, const value_type, value_type>::type access_type;

		typedef access_type* pointer;
		typedef access_type& reference;
	};


	/********************************************
	 *
	 *  Layout attributes
	 *
	 ********************************************/

	template<class Mat>
	struct ct_is_continuous
	{
		typedef typename matrix_traits<Mat>::layout_type layout_type;
		static const bool value = layout_traits<layout_type>::ct_is_continuous;
	};

	template<class Mat>
	struct ct_is_percol_continuous
	{
		typedef typename matrix_traits<Mat>::layout_type layout_type;
		static const bool value = layout_traits<layout_type>::ct_is_percol_continuous;
	};

	template<class Mat>
	struct ct_supports_linear_index
	{
		static const bool value = ct_is_continuous<Mat>::value || ct_is_vector<Mat>::value;
	};


	/********************************************
	 *
	 *  Evaluation
	 *
	 ********************************************/

	template<class SExpr, class DMat>
	struct matrix_eval_verifier
	{
		static const bool value =
				is_dense_mat<DMat>::value &&
				!is_readonly_mat<DMat>::value &&
				is_mat_xpr<SExpr>::value &&
				has_same_domain<SExpr, DMat>::value &&
				has_same_value_type<SExpr, DMat>::value &&
				has_compatible_ct_size<DMat, SExpr>::value;
	};

}

#endif /* MATRIX_META_H_ */
