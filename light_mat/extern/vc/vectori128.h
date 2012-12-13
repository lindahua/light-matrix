/****************************  vectori128.h   *******************************
* Author:        Agner Fog
* Date created:  2012-05-30
* Last modified: 2012-08-01
* Version:       1.02 Beta
* Project:       vector classes
* Description:
* Header file defining integer vector classes as interface to intrinsic
* functions in x86 microprocessors with SSE2 and later instruction sets
* up to AVX.
*
* Instructions:
* Use Gnu, Intel or Microsoft C++ compiler. Compile for the desired
* instruction set, which must be at least SSE2. Specify the supported
* instruction set by a command line define, e.g. __SSE4_1__ if the
* compiler does not automatically do so.
*
* The following vector classes are defined here:
* Vec128b   Vector of 128  1-bit unsigned integers or Booleans
* Vec16c    Vector of  16  8-bit signed   integers
* Vec16uc   Vector of  16  8-bit unsigned integers
* Vec8s     Vector of   8 16-bit signed   integers
* Vec8us    Vector of   8 16-bit unsigned integers
* Vec4i     Vector of   4 32-bit signed   integers
* Vec4ui    Vector of   4 32-bit unsigned integers
* Vec2q     Vector of   2 64-bit signed   integers
* Vec2uq    Vector of   2 64-bit unsigned integers
*
* Each vector object is represented internally in the CPU as a 128-bit register.
* This header file defines operators and functions for these vectors.
*
* For example:
* Vec4i a(1,2,3,4), b(5,6,7,8), c;
* c = a + b;     // now c contains (6,8,10,12)
*
* For detailed instructions, see VectorClass.pdf
*
* (c) Copyright 2012 GNU General Public License http://www.gnu.org/licenses
*****************************************************************************/
#ifndef VECTORI128_H
#define VECTORI128_H

#include "instrset.h"  // Select supported instruction set


