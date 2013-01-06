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

	template<typename Fun>
	struct is_simdizable
	{
		static const bool value = false;
	};

	template<typename Fun, typename Kind>
	struct simdize_map;

}

#define LMAT_DECL_SIMDIZABLE_ON_REAL(FunT) \
		template<> \
		struct is_simdizable<FunT<float> > \
		{ static const bool value = true; }; \
		template<> \
		struct is_simdizable<FunT<double> > \
		{ static const bool value = true; };

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
