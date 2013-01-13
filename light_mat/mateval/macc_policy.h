/**
 * @file macc_policy.h
 *
 * @brief Matrix access policy
 *
 * @author Dahua Lin
 */

#ifdef _MSC_VER
#pragma once
#endif

#ifndef LIGHTMAT_MACC_POLICY_H_
#define LIGHTMAT_MACC_POLICY_H_

#include <light_mat/mateval/mateval_fwd.h>
#include <light_mat/matrix/matrix_concepts.h>

namespace lmat
{
	// Policies

	template<typename ATag> struct linear_macc { };
	template<typename ATag> struct percol_macc { };

	/********************************************
	 *
	 *  Linear MACC support
	 *
	 ********************************************/

	template<typename Mat>
	struct supports_linear_macc
	{
		static const bool value =
				meta::is_regular_mat<Mat>::value &&
				meta::supports_linear_index<Mat>::value;
	};


	/********************************************
	 *
	 *  SIMD support
	 *
	 ********************************************/

	template<typename T, typename Kind>
	struct supports_simd_access : public meta::false_ { };

	template<>
	struct supports_simd_access<float, sse_t> : public meta::true_ { };

	template<>
	struct supports_simd_access<double, sse_t> : public meta::true_ { };

	template<>
	struct supports_simd_access<float, avx_t> : public meta::true_ { };

	template<>
	struct supports_simd_access<double, avx_t> : public meta::true_ { };

	namespace internal
	{
		template<class Mat, typename Kind, bool IsLinear, bool ESimd>
		struct regular_mat_supports_simd : public meta::false_ { };

		template<class Mat, typename Kind>
		struct regular_mat_supports_simd<Mat, Kind, true, true>
		{
			typedef typename matrix_traits<Mat>::value_type T;
			static const unsigned int L = (unsigned int)meta::nelems<Mat>::value;
			static const unsigned int W = simd_traits<T, Kind>::pack_width;

			static const bool value = meta::is_contiguous<Mat>::value && (L % W == 0);
		};

		template<class Mat, typename Kind>
		struct regular_mat_supports_simd<Mat, Kind, false, true>
		{
			typedef typename matrix_traits<Mat>::value_type T;
			static const unsigned int L = (unsigned int)meta::nrows<Mat>::value;
			static const unsigned int W = simd_traits<T, Kind>::pack_width;

			static const bool value = meta::is_percol_contiguous<Mat>::value && (L % W == 0);
		};
	}

	template<typename Mat, typename Kind, bool IsLinear>
	struct supports_simd
	{
		static_assert(meta::is_regular_mat<Mat>::value, "Mat here should be a regular matrix.");

		typedef typename matrix_traits<Mat>::value_type T;
		static const bool value = internal::regular_mat_supports_simd<
				Mat, Kind, IsLinear, supports_simd_access<T, Kind>::value>::value;
	};


	/********************************************
	 *
	 *  preferred policy
	 *
	 ********************************************/

	template<class S>
	struct preferred_macc_policy
	{
		typedef typename matrix_traits<S>::value_type vtype;

		static const bool prefer_linear = supports_linear_macc<S>::value;

		static const bool prefer_simd =
				supports_simd<S, default_simd_kind, prefer_linear>::value;

		typedef typename std::conditional<prefer_simd,
				atags::simd<default_simd_kind>,
				atags::scalar >::type atag;

		typedef typename std::conditional<prefer_linear,
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
				supports_simd<S, default_simd_kind, prefer_linear>::value &&
				supports_simd<D, default_simd_kind, prefer_linear>::value;

		typedef typename std::conditional<prefer_simd,
				atags::simd<default_simd_kind>,
				atags::scalar >::type atag;

		typedef typename std::conditional<prefer_linear,
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
