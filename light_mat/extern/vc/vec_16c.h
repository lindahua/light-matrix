
#ifdef _MSC_VER
#pragma once
#endif

#ifndef LIGHTMAT_VC_VEC_16C_H_
#define LIGHTMAT_VC_VEC_16C_H_

#include "vc_fwd.h"

namespace lmat { namespace vc {

	class Vec16c : public Vec128b {
	public:
	    // Default constructor:
	    Vec16c() {
	    }
	    // Constructor to broadcast the same value into all elements:
	    Vec16c(int i) {
	        xmm = _mm_set1_epi8(i);
	    }
	    // Constructor to build from all elements:
	    Vec16c(int8_t i0, int8_t i1, int8_t i2, int8_t i3, int8_t i4, int8_t i5, int8_t i6, int8_t i7,
	        int8_t i8, int8_t i9, int8_t i10, int8_t i11, int8_t i12, int8_t i13, int8_t i14, int8_t i15) {
	        xmm = _mm_setr_epi8(i0, i1, i2, i3, i4, i5, i6, i7, i8, i9, i10, i11, i12, i13, i14, i15);
	    }
	    // Constructor to convert from type __m128i used in intrinsics:
	    Vec16c(__m128i const & x) {
	        xmm = x;
	    }
	    // Assignment operator to convert from type __m128i used in intrinsics:
	    Vec16c & operator = (__m128i const & x) {
	        xmm = x;
	        return *this;
	    }
	    // Type cast operator to convert to __m128i used in intrinsics
	    operator __m128i() const {
	        return xmm;
	    }
	    // Member function to load from array (unaligned)
	    Vec16c & load(void const * p) {
	        xmm = _mm_loadu_si128((__m128i const*)p);
	        return *this;
	    }
	    // Member function to load from array (aligned)
	    Vec16c & load_a(void const * p) {
	        xmm = _mm_load_si128((__m128i const*)p);
	        return *this;
	    }
	    // Partial load. Load n elements and set the rest to 0
	    Vec16c & load_partial(int n, void const * p) {
	        if      (n >= 16) load(p);
	        else if (n <= 0)  *this = 0;
	        else if (((int)(intptr_t)p & 0xFFF) < 0xFF0) {
	            // p is at least 16 bytes from a page boundary. OK to read 16 bytes
	            load(p);
	        }
	        else {
	            // worst case. read 1 byte at a time and suffer store forwarding penalty
	            char x[16];
	            for (int i = 0; i < n; i++) x[i] = ((char *)p)[i];
	            load(x);
	        }
	        cutoff(n);
	        return *this;
	    }
	    // Partial store. Store n elements
	    void store_partial(int n, void * p) const {
	        if (n >= 16) {
	            store(p);
	            return;
	        }
	        if (n <= 0) return;
	        // we are not using _mm_maskmoveu_si128 because it is too slow on many processors
	        union {
	            int8_t  c[16];
	            int16_t s[8];
	            int32_t i[4];
	            int64_t q[2];
	        } u;
	        store(u.c);
	        int j = 0;
	        if (n & 8) {
	            *(int64_t*)p = u.q[0];
	            j += 8;
	        }
	        if (n & 4) {
	            ((int32_t*)p)[j/4] = u.i[j/4];
	            j += 4;
	        }
	        if (n & 2) {
	            ((int16_t*)p)[j/2] = u.s[j/2];
	            j += 2;
	        }
	        if (n & 1) {
	            ((int8_t*)p)[j]    = u.c[j];
	        }
	    }
	    // cut off vector to n elements. The last 16-n elements are set to zero
	    Vec16c & cutoff(int n) {
	        if (uint32_t(n) >= 16) return *this;
	        static const char mask[32] = {-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
	            0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
	        *this &= Vec16c().load(mask+16-n);
	        return *this;
	    }
	    // Member function to change a single element in vector
	    // Note: This function is inefficient. Use load function if changing more than one element
	    Vec16c const & insert(uint32_t index, int8_t value) {
	        static const int8_t maskl[32] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
	            -1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
	        __m128i broad = _mm_set1_epi8(value);  // broadcast value into all elements
	        __m128i mask  = _mm_loadu_si128((__m128i const*)(maskl+16-(index & 0x0F))); // mask with FF at index position
	        xmm = selectb(mask,broad,xmm);
	        return *this;
	    }
	    // Member function extract a single element from vector
	    int8_t extract(uint32_t index) const {
	        int8_t x[16];
	        store(x);
	        return x[index & 0x0F];
	    }
	    // Extract a single element. Use store function if extracting more than one element.
	    // Operator [] can only read an element, not write.
	    int8_t operator [] (uint32_t index) const {
	        return extract(index);
	    }
	};

	// Define operators for this class

	// vector operator + : add element by element
	static inline Vec16c operator + (Vec16c const & a, Vec16c const & b) {
	    return _mm_add_epi8(a, b);
	}

	// vector operator += : add
	static inline Vec16c & operator += (Vec16c & a, Vec16c const & b) {
	    a = a + b;
	    return a;
	}

	// postfix operator ++
	static inline Vec16c operator ++ (Vec16c & a, int) {
	    Vec16c a0 = a;
	    a = a + 1;
	    return a0;
	}

	// prefix operator ++
	static inline Vec16c & operator ++ (Vec16c & a) {
	    a = a + 1;
	    return a;
	}

	// vector operator - : subtract element by element
	static inline Vec16c operator - (Vec16c const & a, Vec16c const & b) {
	    return _mm_sub_epi8(a, b);
	}

	// vector operator - : unary minus
	static inline Vec16c operator - (Vec16c const & a) {
	    return _mm_sub_epi8(_mm_setzero_si128(), a);
	}

	// vector operator -= : add
	static inline Vec16c & operator -= (Vec16c & a, Vec16c const & b) {
	    a = a - b;
	    return a;
	}

	// postfix operator --
	static inline Vec16c operator -- (Vec16c & a, int) {
	    Vec16c a0 = a;
	    a = a - 1;
	    return a0;
	}

	// prefix operator --
	static inline Vec16c & operator -- (Vec16c & a) {
	    a = a - 1;
	    return a;
	}

	// vector operator * : multiply element by element
	static inline Vec16c operator * (Vec16c const & a, Vec16c const & b) {
	    // There is no 8-bit multiply in SSE2. Split into two 16-bit multiplies
	    __m128i aodd    = _mm_srli_epi16(a,8);                 // odd numbered elements of a
	    __m128i bodd    = _mm_srli_epi16(b,8);                 // odd numbered elements of b
	    __m128i muleven = _mm_mullo_epi16(a,b);                // product of even numbered elements
	    __m128i mulodd  = _mm_mullo_epi16(aodd,bodd);          // product of odd  numbered elements
	            mulodd  = _mm_slli_epi16(mulodd,8);            // put odd numbered elements back in place
	    __m128i mask    = _mm_set1_epi32(0x00FF00FF);          // mask for even positions
	    __m128i product = selectb(mask,muleven,mulodd);        // interleave even and odd
	    return product;
	}

	// vector operator *= : multiply
	static inline Vec16c & operator *= (Vec16c & a, Vec16c const & b) {
	    a = a * b;
	    return a;
	}

	// vector operator << : shift left all elements
	static inline Vec16c operator << (Vec16c const & a, int b) {
	    uint32_t mask = (uint32_t)0xFF >> (uint32_t)b;         // mask to remove bits that are shifted out
	    __m128i am    = _mm_and_si128(a,_mm_set1_epi8(mask));  // remove bits that will overflow
	    __m128i res   = _mm_sll_epi16(am,_mm_cvtsi32_si128(b));// 16-bit shifts
	    return res;
	}

	// vector operator <<= : shift left
	static inline Vec16c & operator <<= (Vec16c & a, int b) {
	    a = a << b;
	    return a;
	}

	// vector operator >> : shift right arithmetic all elements
	static inline Vec16c operator >> (Vec16c const & a, int b) {
	    __m128i aeven = _mm_slli_epi16(a,8);                   // even numbered elements of a. get sign bit in position
	            aeven = _mm_sra_epi16(aeven,_mm_cvtsi32_si128(b+8)); // shift arithmetic, back to position
	    __m128i aodd  = _mm_sra_epi16(a,_mm_cvtsi32_si128(b)); // shift odd numbered elements arithmetic
	    __m128i mask    = _mm_set1_epi32(0x00FF00FF);          // mask for even positions
	    __m128i res     = selectb(mask,aeven,aodd);            // interleave even and odd
	    return res;
	}

	// vector operator >>= : shift right arithmetic
	static inline Vec16c & operator >>= (Vec16c & a, int b) {
	    a = a >> b;
	    return a;
	}

	// vector operator == : returns true for elements for which a == b
	static inline Vec16c operator == (Vec16c const & a, Vec16c const & b) {
	    return _mm_cmpeq_epi8(a,b);
	}

	// vector operator != : returns true for elements for which a != b
	static inline Vec16c operator != (Vec16c const & a, Vec16c const & b) {
	#ifdef __XOP__  // AMD XOP instruction set
	    return _mm_comneq_epi8(a,b);
	#else  // SSE2 instruction set
	    return Vec16c(~(a == b));
	#endif
	}

	// vector operator > : returns true for elements for which a > b (signed)
	static inline Vec16c operator > (Vec16c const & a, Vec16c const & b) {
	    return _mm_cmpgt_epi8(a,b);
	}

	// vector operator < : returns true for elements for which a < b (signed)
	static inline Vec16c operator < (Vec16c const & a, Vec16c const & b) {
	    return b > a;
	}

	// vector operator >= : returns true for elements for which a >= b (signed)
	static inline Vec16c operator >= (Vec16c const & a, Vec16c const & b) {
	#ifdef __XOP__  // AMD XOP instruction set
	    return _mm_comge_epi8(a,b);
	#else  // SSE2 instruction set
	    return Vec16c(~(b > a));
	#endif
	}

	// vector operator <= : returns true for elements for which a <= b (signed)
	static inline Vec16c operator <= (Vec16c const & a, Vec16c const & b) {
	    return b >= a;
	}

	// vector operator & : bitwise and
	static inline Vec16c operator & (Vec16c const & a, Vec16c const & b) {
	    return Vec16c(Vec128b(a) & Vec128b(b));
	}
	static inline Vec16c operator && (Vec16c const & a, Vec16c const & b) {
	    return a & b;
	}

	// vector operator | : bitwise or
	static inline Vec16c operator | (Vec16c const & a, Vec16c const & b) {
	    return Vec16c(Vec128b(a) | Vec128b(b));
	}
	static inline Vec16c operator || (Vec16c const & a, Vec16c const & b) {
	    return a | b;
	}

	// vector operator ^ : bitwise xor
	static inline Vec16c operator ^ (Vec16c const & a, Vec16c const & b) {
	    return Vec16c(Vec128b(a) ^ Vec128b(b));
	}

	// vector operator ~ : bitwise not
	static inline Vec16c operator ~ (Vec16c const & a) {
	    return Vec16c( ~ Vec128b(a));
	}

	// vector operator ! : logical not, returns true for elements == 0
	static inline Vec16c operator ! (Vec16c const & a) {
	    return _mm_cmpeq_epi8(a,_mm_setzero_si128());
	}

	// Functions for this class

	// Select between two operands. Corresponds to this pseudocode:
	// for (int i = 0; i < 16; i++) result[i] = s[i] ? a[i] : b[i];
	// Each byte in s must be either 0 (false) or -1 (true). No other values are allowed.
	static inline Vec16c select (Vec16c const & s, Vec16c const & a, Vec16c const & b) {
	    return selectb(s,a,b);
	}

	// Horizontal add: Calculates the sum of all vector elements.
	// Overflow will wrap around
	static inline int32_t horizontal_add (Vec16c const & a) {
	    __m128i sum1 = _mm_sad_epu8(a,_mm_setzero_si128());
	    __m128i sum2 = _mm_shuffle_epi32(sum1,2);
	    __m128i sum3 = _mm_add_epi16(sum1,sum2);
	    int8_t  sum4 = _mm_cvtsi128_si32(sum3);      // truncate to 8 bits
	    return  sum4;                                // sign extend to 32 bits
	}

	// Horizontal add extended: Calculates the sum of all vector elements.
	// Each element is sign-extended before addition to avoid overflow
	static inline int32_t horizontal_add_x (Vec16c const & a) {
	#ifdef __XOP__       // AMD XOP instruction set
	    __m128i sum1  = _mm_haddq_epi8(a);
	    __m128i sum2  = _mm_shuffle_epi32(sum1,0x0E);          // high element
	    __m128i sum3  = _mm_add_epi32(sum1,sum2);              // sum
	    return          _mm_cvtsi128_si32(sum3);
	#elif  INSTRSET >= 4  // SSSE3
	    __m128i aeven = _mm_slli_epi16(a,8);                   // even numbered elements of a. get sign bit in position
	            aeven = _mm_srai_epi16(aeven,8);               // sign extend even numbered elements
	    __m128i aodd  = _mm_srai_epi16(a,8);                   // sign extend odd  numbered elements
	    __m128i sum1  = _mm_add_epi16(aeven,aodd);             // add even and odd elements
	    __m128i sum2  = _mm_hadd_epi16(sum1,sum1);             // horizontally add 8 elements in 3 steps
	    __m128i sum3  = _mm_hadd_epi16(sum2,sum2);
	    __m128i sum4  = _mm_hadd_epi16(sum3,sum3);
	    int16_t sum5  = _mm_cvtsi128_si32(sum4);               // 16 bit sum
	    return  sum5;                                          // sign extend to 32 bits
	#else                 // SSE2
	    __m128i aeven = _mm_slli_epi16(a,8);                   // even numbered elements of a. get sign bit in position
	            aeven = _mm_srai_epi16(aeven,8);               // sign extend even numbered elements
	    __m128i aodd  = _mm_srai_epi16(a,8);                   // sign extend odd  numbered elements
	    __m128i sum1  = _mm_add_epi16(aeven,aodd);             // add even and odd elements
	    __m128i sum2  = _mm_shuffle_epi32(sum1,0x0E);          // 4 high elements
	    __m128i sum3  = _mm_add_epi16(sum1,sum2);              // 4 sums
	    __m128i sum4  = _mm_shuffle_epi32(sum3,0x01);          // 2 high elements
	    __m128i sum5  = _mm_add_epi16(sum3,sum4);              // 2 sums
	    __m128i sum6  = _mm_shufflelo_epi16(sum5,0x01);        // 1 high element
	    __m128i sum7  = _mm_add_epi16(sum5,sum6);              // 1 sum
	    int16_t sum8  = _mm_cvtsi128_si32(sum7);               // 16 bit sum
	    return  sum8;                                          // sign extend to 32 bits
	#endif
	}


	// function add_saturated: add element by element, signed with saturation
	static inline Vec16c add_saturated(Vec16c const & a, Vec16c const & b) {
	    return _mm_adds_epi8(a, b);
	}

	// function sub_saturated: subtract element by element, signed with saturation
	static inline Vec16c sub_saturated(Vec16c const & a, Vec16c const & b) {
	    return _mm_subs_epi8(a, b);
	}

	// function max: a > b ? a : b
	static inline Vec16c max(Vec16c const & a, Vec16c const & b) {
	#if INSTRSET >= 5   // SSE4.1
	    return _mm_max_epi8(a,b);
	#else  // SSE2
	    __m128i signbit = _mm_set1_epi32(0x80808080);
	    __m128i a1      = _mm_xor_si128(a,signbit);            // add 0x80
	    __m128i b1      = _mm_xor_si128(b,signbit);            // add 0x80
	    __m128i m1      = _mm_max_epu8(a1,b1);                 // unsigned max
	    return  _mm_xor_si128(m1,signbit);                     // sub 0x80
	#endif
	}

	// function min: a < b ? a : b
	static inline Vec16c min(Vec16c const & a, Vec16c const & b) {
	#if INSTRSET >= 5   // SSE4.1
	    return _mm_min_epi8(a,b);
	#else  // SSE2
	    __m128i signbit = _mm_set1_epi32(0x80808080);
	    __m128i a1      = _mm_xor_si128(a,signbit);            // add 0x80
	    __m128i b1      = _mm_xor_si128(b,signbit);            // add 0x80
	    __m128i m1      = _mm_min_epu8(a1,b1);                 // unsigned min
	    return  _mm_xor_si128(m1,signbit);                     // sub 0x80
	#endif
	}

	// function abs: a >= 0 ? a : -a
	static inline Vec16c abs(Vec16c const & a) {
	#if INSTRSET >= 4     // SSSE3 supported
	    return _mm_sign_epi8(a,a);
	#else                 // SSE2
	    __m128i nega = _mm_sub_epi8(_mm_setzero_si128(), a);
	    return _mm_min_epu8(a, nega);   // unsigned min (the negative value is bigger when compared as unsigned)
	#endif
	}

	// function abs_saturated: same as abs, saturate if overflow
	static inline Vec16c abs_saturated(Vec16c const & a) {
	    __m128i absa   = abs(a);                               // abs(a)
	    __m128i overfl = _mm_cmpgt_epi8(_mm_setzero_si128(),absa);// 0 > a
	    return           _mm_add_epi8(absa,overfl);            // subtract 1 if 0x80
	}

	// function rotate_left: rotate each element left by b bits
	// Use negative count to rotate right
	static inline Vec16c rotate_left(Vec16c const & a, int b) {
	#ifdef __XOP__  // AMD XOP instruction set
	    return _mm_rot_epi8(a,_mm_set1_epi8(b));
	#else  // SSE2 instruction set
	    __m128i bb        = _mm_cvtsi32_si128(b & 7);          // b modulo 8
	    __m128i mbb       = _mm_cvtsi32_si128((8-b) & 7);      // 8-b modulo 8
	    __m128i maskeven  = _mm_set1_epi32(0x00FF00FF);        // mask for even numbered bytes
	    __m128i even      = _mm_and_si128(a,maskeven);         // even numbered bytes of a
	    __m128i odd       = _mm_andnot_si128(maskeven,a);      // odd numbered bytes of a
	    __m128i evenleft  = _mm_sll_epi16(even,bb);            // even bytes of a << b
	    __m128i oddleft   = _mm_sll_epi16(odd,bb);             // odd  bytes of a << b
	    __m128i evenright = _mm_srl_epi16(even,mbb);           // even bytes of a >> 8-b
	    __m128i oddright  = _mm_srl_epi16(odd,mbb);            // odd  bytes of a >> 8-b
	    __m128i evenrot   = _mm_or_si128(evenleft,evenright);  // even bytes of a rotated
	    __m128i oddrot    = _mm_or_si128(oddleft,oddright);    // odd  bytes of a rotated
	    __m128i allrot    = selectb(maskeven,evenrot,oddrot);  // all  bytes rotated
	    return  allrot;
	#endif
	}


	/*****************************************************************************
	*
	*          Vector of 16 8-bit unsigned integers
	*
	*****************************************************************************/

	class Vec16uc : public Vec16c {
	public:
	    // Default constructor:
	    Vec16uc() {
	    };
	    // Constructor to broadcast the same value into all elements:
	    Vec16uc(uint32_t i) {
	        xmm = _mm_set1_epi8(i);
	    };
	    // Constructor to build from all elements:
	    Vec16uc(uint8_t i0, uint8_t i1, uint8_t i2, uint8_t i3, uint8_t i4, uint8_t i5, uint8_t i6, uint8_t i7,
	        uint8_t i8, uint8_t i9, uint8_t i10, uint8_t i11, uint8_t i12, uint8_t i13, uint8_t i14, uint8_t i15) {
	        xmm = _mm_setr_epi8(i0, i1, i2, i3, i4, i5, i6, i7, i8, i9, i10, i11, i12, i13, i14, i15);
	    };
	    // Constructor to convert from type __m128i used in intrinsics:
	    Vec16uc(__m128i const & x) {
	        xmm = x;
	    };
	    // Assignment operator to convert from type __m128i used in intrinsics:
	    Vec16uc & operator = (__m128i const & x) {
	        xmm = x;
	        return *this;
	    };
	    // Member function to load from array (unaligned)
	    Vec16uc & load(void const * p) {
	        xmm = _mm_loadu_si128((__m128i const*)p);
	        return *this;
	    }
	    // Member function to load from array (aligned)
	    Vec16uc & load_a(void const * p) {
	        xmm = _mm_load_si128((__m128i const*)p);
	        return *this;
	    }
	    // Member function to change a single element in vector
	    // Note: This function is inefficient. Use load function if changing more than one element
	    Vec16uc const & insert(uint32_t index, uint8_t value) {
	        Vec16c::insert(index, value);
	        return *this;
	    }
	    // Member function extract a single element from vector
	    uint8_t extract(uint32_t index) const {
	        return Vec16c::extract(index);
	    }
	    // Extract a single element. Use store function if extracting more than one element.
	    // Operator [] can only read an element, not write.
	    uint8_t operator [] (uint32_t index) const {
	        return extract(index);
	    }
	};

	// Define operators for this class

	// vector operator << : shift left all elements
	static inline Vec16uc operator << (Vec16uc const & a, uint32_t b) {
	    uint32_t mask = (uint32_t)0xFF >> (uint32_t)b;         // mask to remove bits that are shifted out
	    __m128i am    = _mm_and_si128(a,_mm_set1_epi8(mask));  // remove bits that will overflow
	    __m128i res   = _mm_sll_epi16(am,_mm_cvtsi32_si128(b));// 16-bit shifts
	    return res;
	}

	// vector operator << : shift left all elements
	static inline Vec16uc operator << (Vec16uc const & a, int32_t b) {
	    return a << (uint32_t)b;
	}

	// vector operator >> : shift right logical all elements
	static inline Vec16uc operator >> (Vec16uc const & a, uint32_t b) {
	    uint32_t mask = (uint32_t)0xFF << (uint32_t)b;         // mask to remove bits that are shifted out
	    __m128i am    = _mm_and_si128(a,_mm_set1_epi8(mask));  // remove bits that will overflow
	    __m128i res   = _mm_srl_epi16(am,_mm_cvtsi32_si128(b));// 16-bit shifts
	    return res;
	}

	// vector operator >> : shift right logical all elements
	static inline Vec16uc operator >> (Vec16uc const & a, int32_t b) {
	    return a >> (uint32_t)b;
	}

	// vector operator >>= : shift right logical
	static inline Vec16uc & operator >>= (Vec16uc & a, int b) {
	    a = a >> b;
	    return a;
	}

	// vector operator >= : returns true for elements for which a >= b (unsigned)
	static inline Vec16c operator >= (Vec16uc const & a, Vec16uc const & b) {
	#ifdef __XOP__  // AMD XOP instruction set
	    return _mm_comge_epu8(a,b);
	#else  // SSE2 instruction set
	    return _mm_cmpeq_epi8(_mm_max_epu8(a,b),a); // a == max(a,b)
	#endif
	}

	// vector operator <= : returns true for elements for which a <= b (unsigned)
	static inline Vec16c operator <= (Vec16uc const & a, Vec16uc const & b) {
	    return b >= a;
	}

	// vector operator > : returns true for elements for which a > b (unsigned)
	static inline Vec16c operator > (Vec16uc const & a, Vec16uc const & b) {
	#ifdef __XOP__  // AMD XOP instruction set
	    return _mm_comgt_epu8(a,b);
	#else  // SSE2 instruction set
	    return Vec16c(~(b >= a));
	#endif
	}

	// vector operator < : returns true for elements for which a < b (unsigned)
	static inline Vec16c operator < (Vec16uc const & a, Vec16uc const & b) {
	    return b > a;
	}

	// vector operator + : add
	static inline Vec16uc operator + (Vec16uc const & a, Vec16uc const & b) {
	    return Vec16uc (Vec16c(a) + Vec16c(b));
	}

	// vector operator - : subtract
	static inline Vec16uc operator - (Vec16uc const & a, Vec16uc const & b) {
	    return Vec16uc (Vec16c(a) - Vec16c(b));
	}

	// vector operator * : multiply
	static inline Vec16uc operator * (Vec16uc const & a, Vec16uc const & b) {
	    return Vec16uc (Vec16c(a) * Vec16c(b));
	}

	// vector operator & : bitwise and
	static inline Vec16uc operator & (Vec16uc const & a, Vec16uc const & b) {
	    return Vec16uc(Vec128b(a) & Vec128b(b));
	}
	static inline Vec16uc operator && (Vec16uc const & a, Vec16uc const & b) {
	    return a & b;
	}

	// vector operator | : bitwise or
	static inline Vec16uc operator | (Vec16uc const & a, Vec16uc const & b) {
	    return Vec16uc(Vec128b(a) | Vec128b(b));
	}
	static inline Vec16uc operator || (Vec16uc const & a, Vec16uc const & b) {
	    return a | b;
	}

	// vector operator ^ : bitwise xor
	static inline Vec16uc operator ^ (Vec16uc const & a, Vec16uc const & b) {
	    return Vec16uc(Vec128b(a) ^ Vec128b(b));
	}

	// vector operator ~ : bitwise not
	static inline Vec16uc operator ~ (Vec16uc const & a) {
	    return Vec16uc( ~ Vec128b(a));
	}

	// Functions for this class

	// Select between two operands. Corresponds to this pseudocode:
	// for (int i = 0; i < 16; i++) result[i] = s[i] ? a[i] : b[i];
	// Each byte in s must be either 0 (false) or -1 (true). No other values are allowed.
	// (s is signed)
	static inline Vec16uc select (Vec16c const & s, Vec16uc const & a, Vec16uc const & b) {
	    return selectb(s,a,b);
	}

	// Horizontal add: Calculates the sum of all vector elements.
	// Overflow will wrap around
	// (Note: horizontal_add_x(Vec16uc) is slightly faster)
	static inline uint32_t horizontal_add (Vec16uc const & a) {
	    __m128i sum1 = _mm_sad_epu8(a,_mm_setzero_si128());
	    __m128i sum2 = _mm_shuffle_epi32(sum1,2);
	    __m128i sum3 = _mm_add_epi16(sum1,sum2);
	    uint8_t sum4 = _mm_cvtsi128_si32(sum3);      // truncate to 8 bits
	    return  sum4;
	}

	// Horizontal add extended: Calculates the sum of all vector elements.
	// Each element is zero-extended before addition to avoid overflow
	static inline uint32_t horizontal_add_x (Vec16uc const & a) {
	    __m128i sum1 = _mm_sad_epu8(a,_mm_setzero_si128());
	    __m128i sum2 = _mm_shuffle_epi32(sum1,2);
	    __m128i sum3 = _mm_add_epi16(sum1,sum2);
	    return _mm_cvtsi128_si32(sum3);
	}

	// function add_saturated: add element by element, unsigned with saturation
	static inline Vec16uc add_saturated(Vec16uc const & a, Vec16uc const & b) {
	    return _mm_adds_epu8(a, b);
	}

	// function sub_saturated: subtract element by element, unsigned with saturation
	static inline Vec16uc sub_saturated(Vec16uc const & a, Vec16uc const & b) {
	    return _mm_subs_epu8(a, b);
	}

	// function max: a > b ? a : b
	static inline Vec16uc max(Vec16uc const & a, Vec16uc const & b) {
	    return _mm_max_epu8(a,b);
	}

	// function min: a < b ? a : b
	static inline Vec16uc min(Vec16uc const & a, Vec16uc const & b) {
	    return _mm_min_epu8(a,b);
	}

} }

#endif 
