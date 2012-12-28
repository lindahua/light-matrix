/**
 * @file macc_policy.h
 *
 * @brief Matrix access policy
 *
 * @author Dahua Lin
 */

#ifndef LIGHTMAT_MACC_POLICY_H_
#define LIGHTMAT_MACC_POLICY_H_

#include <light_mat/mateval/mateval_fwd.h>
#include <light_mat/matrix/matrix_concepts.h>

namespace lmat
{
	// Policies

	template<typename ATag> struct linear_macc { };
	template<typename ATag> struct percol_macc { };

	// Policy resolution

	template<typename Mat>
	struct supports_linear_macc
	{
		static const bool value = meta::supports_linear_index<Mat>::value;
	};


	namespace internal
	{

		template<typename Mat, typename T, typename Kind, bool IsLinear>
		struct _supports_simd
		{
			static const bool value = false;
		};

		template<typename Mat, typename Kind>
		struct _supports_simd<Mat, float, Kind, true>
		{
			static const int pw = math::simd_traits<float, Kind>::pack_width;
			static const int mod_ = meta::nelems<Mat>::value % pw;

			static const bool value = meta::is_continuous<Mat>::value && (mod_ == 0);
		};

		template<typename Mat, typename Kind>
		struct _supports_simd<Mat, double, Kind, true>
		{
			static const int pw = math::simd_traits<float, Kind>::pack_width;
			static const int mod_ = meta::nelems<Mat>::value % pw;

			static const bool value = meta::is_continuous<Mat>::value && (mod_ == 0);
		};

		template<typename Mat, typename Kind>
		struct _supports_simd<Mat, float, Kind, false>
		{
			static const int pw = math::simd_traits<float, Kind>::pack_width;
			static const int mod_ = meta::nrows<Mat>::value % pw;

			static const bool value = meta::is_percol_continuous<Mat>::value && (mod_ == 0);
		};

		template<typename Mat, typename Kind>
		struct _supports_simd<Mat, double, Kind, false>
		{
			static const int pw = math::simd_traits<float, Kind>::pack_width;
			static const int mod_ = meta::nrows<Mat>::value % pw;

			static const bool value = meta::is_percol_continuous<Mat>::value && (mod_ == 0);
		};

		template<typename T1, typename T2>
		struct are_simd_compatible_types
		{
			static const bool value = false;
		};

		template<typename T>
		struct are_simd_compatible_types<T, T> { static const bool value = true; };

		template<typename T>
		struct are_simd_compatible_types<mask_t<T>, T> { static const bool value = true; };
	}

	template<typename Mat, typename T, typename Kind, bool IsLinear>
	struct supports_simd
	{
		typedef typename matrix_traits<Mat>::value_type VT;
		static const bool value = internal::are_simd_compatible_types<VT, T>::value
				&& internal::_supports_simd<Mat, T, Kind, IsLinear>::value;
	};


	template<class S>
	struct preferred_macc_policy
	{
		typedef typename matrix_traits<S>::value_type vtype;

		static const bool prefer_linear = supports_linear_macc<S>::value;

		static const bool prefer_simd =
				supports_simd<S, vtype, default_simd_kind, prefer_linear>::value;

		typedef typename meta::if_c<prefer_simd,
				atags::simd<default_simd_kind>,
				atags::scalar >::type atag;

		typedef typename meta::if_c<prefer_linear,
				linear_macc<atag>,
				percol_macc<atag> >::type type;
	};

	template<class S, class D>
	struct preferred_macc_policy_2
	{
		typedef typename meta::common_value_type<S, D>::type vtype;

		static const bool prefer_linear =
				supports_linear_macc<S>::value &&
				supports_linear_macc<D>::value;

		static const bool prefer_simd =
				supports_simd<S, vtype, default_simd_kind, prefer_linear>::value &&
				supports_simd<D, vtype, default_simd_kind, prefer_linear>::value;

		typedef typename meta::if_c<prefer_simd,
				atags::simd<default_simd_kind>,
				atags::scalar >::type atag;

		typedef typename meta::if_c<prefer_linear,
				linear_macc<atag>,
				percol_macc<atag> >::type type;
	};


	template<typename T, class A>
	LMAT_ENSURE_INLINE
	inline typename preferred_macc_policy<A>::type
	get_preferred_macc_policy(const IEWiseMatrix<A, T>& a)
	{
		typedef typename preferred_macc_policy<A>::type policy_t;
		return policy_t();
	}

	template<typename T, class A, class B>
	LMAT_ENSURE_INLINE
	inline typename preferred_macc_policy_2<A, B>::type
	get_preferred_macc_policy(const IEWiseMatrix<A, T>& a, const IEWiseMatrix<B, T>& b)
	{
		typedef typename preferred_macc_policy_2<A, B>::type policy_t;
		return policy_t();
	}

}

#endif /* MACC_POLICY_H_ */
