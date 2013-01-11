/**
 * @file math_constants.h
 *
 * Useful math constants
 *
 * @author Dahua Lin
 */

#ifdef _MSC_VER
#pragma once
#endif

#ifndef LIGHTMAT_MATH_CONSTANTS_H_
#define LIGHTMAT_MATH_CONSTANTS_H_

#include <light_mat/common/prim_types.h>

namespace lmat { namespace math {

	template<typename T> struct consts;

	template<> struct consts<double>
	{
		LMAT_ENSURE_INLINE
		static double pi() { return 3.141592653589793238462643383; }  	// Pi

		LMAT_ENSURE_INLINE
		static double two_pi() { return 6.283185307179586476925286767; }  	// 2 * Pi

		LMAT_ENSURE_INLINE
		static double half_pi() { return 1.570796326794896619231321692; }	// Pi / 2

		LMAT_ENSURE_INLINE
		static double quar_pi() { return 0.7853981633974483096156608458; }	// Pi / 4

		LMAT_ENSURE_INLINE
		static double rcp_pi() { return 0.3183098861837906715377675267; }	// 1 / Pi

		LMAT_ENSURE_INLINE
		static double two_rcp_pi() { return 0.6366197723675813430755350535; }	// 2 / Pi

		LMAT_ENSURE_INLINE
		static double rcp_two_rpi() { return 0.1591549430918953357688837634; }	// 1 / (2 * Pi)

		LMAT_ENSURE_INLINE
		static double e() { return 2.718281828459045235360287471; }	// E
	};


	template<> struct consts<float>
	{
		LMAT_ENSURE_INLINE
		static float pi() { return 3.14159265359f; }  	// Pi

		LMAT_ENSURE_INLINE
		static float two_pi() { return 6.28318530718f; }  	// 2 * Pi

		LMAT_ENSURE_INLINE
		static float half_pi() { return 1.57079632679f; }	// Pi / 2

		LMAT_ENSURE_INLINE
		static float quar_pi() { return 0.785398163397f; }	// Pi / 4

		LMAT_ENSURE_INLINE
		static float rcp_pi() { return 0.318309886184f; }	// 1 / Pi

		LMAT_ENSURE_INLINE
		static float two_rcp_pi() { return 0.636619772368f; }	// 2 / Pi

		LMAT_ENSURE_INLINE
		static float rcp_two_rpi() { return 0.159154943092f; }	// 1 / (2 * Pi)

		LMAT_ENSURE_INLINE
		static float e() { return 2.718281828459f; }	// E
	};

} }

#endif /* MATH_CONSTANTS_H_ */
