/**
 * @file prim_types.h
 *
 * @brief Declaration of primitive types
 *
 * @author Dahua Lin
 */

#ifdef _MSC_VER
#pragma once
#endif

#ifndef LIGHTMAT_PRIM_TYPES_H_
#define LIGHTMAT_PRIM_TYPES_H_

#include <light_mat/config/config.h>

#include <cstddef>

#ifdef LMAT_USE_C11_STDLIB
#include <cstdint>
#else
#include <tr1/cstdint>
#endif

namespace lmat
{
	using LMAT_TR1::int8_t;
	using LMAT_TR1::int16_t;
	using LMAT_TR1::int32_t;
	using LMAT_TR1::int64_t;

	using LMAT_TR1::uint8_t;
	using LMAT_TR1::uint16_t;
	using LMAT_TR1::uint32_t;
	using LMAT_TR1::uint64_t;

	using std::ptrdiff_t;
	using std::size_t;

	typedef int32_t index_t;
}

#endif
