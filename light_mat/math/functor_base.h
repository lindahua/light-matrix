/**
 * @file functor_base.h
 *
 * @brief The basis of evaluation functors
 *
 * @author Dahua Lin
 */

#ifdef _MSC_VER
#pragma once
#endif

#ifndef LIGHTMAT_FUNCTOR_BASE_H_
#define LIGHTMAT_FUNCTOR_BASE_H_

#include <light_mat/common/basic_defs.h>
#include <light_mat/math/simd_base.h>

namespace lmat {

	template<typename Fun, typename Kind>
	struct is_simdizable : public meta::false_ { };

	template<typename Fun, typename Kind>
	struct simdize_map;

	// generic result of predicates

	template<typename T>
	struct pred_result
	{
		typedef mask_t<T> type;
	};

	template<typename T>
	struct pred_result<mask_t<T> >
	{
		typedef mask_t<T> type;
	};

	template<>
	struct pred_result<bool>
	{
		typedef bool type;
	};

	template<typename T, typename Kind>
	struct pred_result<math::simd_pack<T, Kind> >
	{
		typedef math::simd_bpack<T, Kind> type;
	};

	template<typename T, typename Kind>
	struct pred_result<math::simd_bpack<T, Kind> >
	{
		typedef math::simd_bpack<T, Kind> type;
	};


}

#define LMAT_DECL_SIMDIZABLE_ON_REAL(FunT) \
		template<typename Kind> \
		struct is_simdizable<FunT<float>, Kind> : public std::true_type { }; \
		template<typename Kind> \
		struct is_simdizable<FunT<double>, Kind> : public std::true_type { };

#define LMAT_DEF_TRIVIAL_SIMDIZE_MAP(FunT) \
		template<typename Kind> \
		struct simdize_map<FunT<float>, Kind> { \
			typedef FunT<math::simd_pack<float, Kind> > type; \
			LMAT_ENSURE_INLINE \
			static type get(FunT<float> ) { return type(); } \
		}; \
		template<typename Kind> \
		struct simdize_map<FunT<double>, Kind> { \
			typedef FunT<math::simd_pack<double, Kind> > type; \
			LMAT_ENSURE_INLINE \
			static type get(FunT<double> ) { return type(); } \
		};

#define LMAT_DEF_SIMD_SUPPORT( FunT ) \
	LMAT_DECL_SIMDIZABLE_ON_REAL( FunT ) \
	LMAT_DEF_TRIVIAL_SIMDIZE_MAP( FunT )

#endif
