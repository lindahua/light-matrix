
#ifdef _MSC_VER
#pragma once
#endif

#ifndef LIGHTMAT_VC_VEC_128B_H_
#define LIGHTMAT_VC_VEC_128B_H_

#include "vc_fwd.h"

namespace lmat { namespace vc {

	/*****************************************************************************
	*
	*          Vector of 128 1-bit unsigned integers or Booleans
	*
	*****************************************************************************/
	class Vec128b {
	protected:
	    __m128i xmm; // Integer vector
	public:
	    // Default constructor:
	    Vec128b() {
	    };
	    // Constructor to broadcast the same value into all elements:
	    Vec128b(int i) {
	        xmm = _mm_set1_epi32(-(i & 1));
	    };
	    // Constructor to convert from type __m128i used in intrinsics:
	    Vec128b(__m128i const & x) {
	        xmm = x;
	    };
	    // Assignment operator to convert from type __m128i used in intrinsics:
	    Vec128b & operator = (__m128i const & x) {
	        xmm = x;
	        return *this;
	    };
	    // Type cast operator to convert to __m128i used in intrinsics
	    operator __m128i() const {
	        return xmm;
	    }
	    // Member function to load from array (unaligned)
	    Vec128b & load(void const * p) {
	        xmm = _mm_loadu_si128((__m128i const*)p);
	        return *this;
	    }
	    // Member function to load from array, aligned by 16
	    // "load_a" is faster than "load" on older Intel processors (Pentium 4, Pentium M, Core 1,
	    // Merom, Wolfdale) and Atom, but not on other processors from Intel, AMD or VIA.
	    // You may use load_a instead of load if you are certain that p points to an address
	    // divisible by 16.
	    void load_a(void const * p) {
	        xmm = _mm_load_si128((__m128i const*)p);
	    }
	    // Member function to store into array (unaligned)
	    void store(void * p) const {
	        _mm_storeu_si128((__m128i*)p, xmm);
	    }
	    // Member function to store into array, aligned by 16
	    // "store_a" is faster than "store" on older Intel processors (Pentium 4, Pentium M, Core 1,
	    // Merom, Wolfdale) and Atom, but not on other processors from Intel, AMD or VIA.
	    // You may use store_a instead of store if you are certain that p points to an address
	    // divisible by 16.
	    void store_a(void * p) const {
	        _mm_store_si128((__m128i*)p, xmm);
	    }
	    // Member function to change a single bit
	    // Note: This function is inefficient. Use load function if changing more than one bit
	    Vec128b const & set_bit(uint32_t index, int value) {
	        static const union {
	            uint64_t i[4];
	            __m128i  x[2];
	        } u = {{1,0,0,1}};                 // 2 vectors with bit 0 and 64 set, respectively
	        int w = (index >> 6) & 1;          // qword index
	        int bi = index & 0x3F;             // bit index within qword w
	        __m128i mask = u.x[w];
	        mask = _mm_sll_epi64(mask,_mm_cvtsi32_si128(bi)); // mask with bit number b set
	        if (value & 1) {
	            xmm = _mm_or_si128(mask,xmm);
	        }
	        else {
	            xmm = _mm_andnot_si128(mask,xmm);
	        }
	        return *this;
	    }
	    // Member function to get a single bit
	    // Note: This function is inefficient. Use store function if reading more than one bit
	    int get_bit(uint32_t index) const {
	        union {
	            __m128i x;
	            uint8_t i[16];
	        } u;
	        u.x = xmm;
	        int w = (index >> 3) & 0xF;            // byte index
	        int bi = index & 7;                    // bit index within byte w
	        return (u.i[w] >> bi) & 1;
	    }
	    // Extract a single element. Use store function if extracting more than one element.
	    // Operator [] can only read an element, not write.
	    bool operator [] (uint32_t index) const {
	        return get_bit(index) != 0;
	    }
	};


	// Define operators for this class

	// vector operator & : bitwise and
	static inline Vec128b operator & (Vec128b const & a, Vec128b const & b) {
	    return _mm_and_si128(a, b);
	}
	static inline Vec128b operator && (Vec128b const & a, Vec128b const & b) {
	    return a & b;
	}

	// vector operator | : bitwise or
	static inline Vec128b operator | (Vec128b const & a, Vec128b const & b) {
	    return _mm_or_si128(a, b);
	}
	static inline Vec128b operator || (Vec128b const & a, Vec128b const & b) {
	    return a | b;
	}

	// vector operator ^ : bitwise xor
	static inline Vec128b operator ^ (Vec128b const & a, Vec128b const & b) {
	    return _mm_xor_si128(a, b);
	}

	// vector operator ~ : bitwise not
	static inline Vec128b operator ~ (Vec128b const & a) {
	    return _mm_xor_si128(a, _mm_set1_epi32(-1));
	}

	// vector operator &= : bitwise and
	static inline Vec128b & operator &= (Vec128b & a, Vec128b const & b) {
	    a = a & b;
	    return a;
	}

	// vector operator |= : bitwise or
	static inline Vec128b & operator |= (Vec128b & a, Vec128b const & b) {
	    a = a | b;
	    return a;
	}

	// vector operator ^= : bitwise xor
	static inline Vec128b & operator ^= (Vec128b & a, Vec128b const & b) {
	    a = a ^ b;
	    return a;
	}

	// Define functions for this class

	// function andnot: a & ~ b
	static inline Vec128b andnot (Vec128b const & a, Vec128b const & b) {
	    return _mm_andnot_si128(b, a);
	}


	/*****************************************************************************
	*
	*          Generate compile-time constant vector
	*
	*****************************************************************************/
	// Generate a constant vector of 4 integers stored in memory.
	// Can be converted to any integer vector type
	template <int i0, int i1, int i2, int i3>
	static inline __m128i constant4i() {
	    static const union {
	        int     i[4];
	        __m128i xmm;
	    } u = {{i0,i1,i2,i3}};
	    return u.xmm;
	}


	/*****************************************************************************
	*
	*          selectb function
	*
	*****************************************************************************/
	// Select between two sources, byte by byte. Used in various functions and operators
	// Corresponds to this pseudocode:
	// for (int i = 0; i < 16; i++) result[i] = s[i] ? a[i] : b[i];
	// Each byte in s must be either 0 (false) or 0xFF (true). No other values are allowed.
	// The implementation depends on the instruction set:
	// If SSE4.1 is supported then only bit 7 in each byte of s is checked,
	// otherwise all bits in s are used.
	static inline __m128i selectb (__m128i const & s, __m128i const & a, __m128i const & b) {
	#if INSTRSET >= 5   // SSE4.1 supported
	    return _mm_blendv_epi8 (b, a, s);
	#else
	    return _mm_or_si128(
	        _mm_and_si128(s,a),
	        _mm_andnot_si128(s,b));
	#endif
	}


	/*****************************************************************************
	*
	*          Horizontal Boolean functions
	*
	*****************************************************************************/

	// horizontal_and. Returns true if all bits are 1
	static inline bool horizontal_and (Vec128b const & a) {
	#if INSTRSET >= 5   // SSE4.1 supported. Use PTEST
	    return _mm_testc_si128(a,constant4i<-1,-1,-1,-1>()) != 0;
	#else
	    __m128i t1 = _mm_unpackhi_epi64(a,a);                  // get 64 bits down
	    __m128i t2 = _mm_and_si128(a,t1);                      // and 64 bits
	#ifdef __x86_64__
	    int64_t t5 = _mm_cvtsi128_si64(t2);                    // transfer 64 bits to integer
	    return  t5 == int64_t(-1);
	#else
	    __m128i t3 = _mm_srli_epi64(t2,32);                    // get 32 bits down
	    __m128i t4 = _mm_and_si128(t2,t3);                     // and 32 bits
	    int     t5 = _mm_cvtsi128_si32(t4);                    // transfer 32 bits to integer
	    return  t5 == -1;
	#endif  // __x86_64__
	#endif  // INSTRSET
	}

	// horizontal_or. Returns true if at least one bit is 1
	static inline bool horizontal_or (Vec128b const & a) {
	#if INSTRSET >= 5   // SSE4.1 supported. Use PTEST
	    return ! _mm_testz_si128(a,a);
	#else
	    __m128i t1 = _mm_unpackhi_epi64(a,a);                  // get 64 bits down
	    __m128i t2 = _mm_or_si128(a,t1);                       // and 64 bits
	#ifdef __x86_64__
	    int64_t t5 = _mm_cvtsi128_si64(t2);                    // transfer 64 bits to integer
	    return  t5 != int64_t(0);
	#else
	    __m128i t3 = _mm_srli_epi64(t2,32);                    // get 32 bits down
	    __m128i t4 = _mm_or_si128(t2,t3);                      // and 32 bits
	    int     t5 = _mm_cvtsi128_si32(t4);                    // transfer to integer
	    return  t5 != 0;
	#endif  // __x86_64__
	#endif  // INSTRSET
	}

} }

#endif 
