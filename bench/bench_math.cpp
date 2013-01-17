/**
 * @file bench_math.cpp
 *
 * @brief Benchmark the evaluation of math functions
 *
 * @author Dahua Lin
 */


#include "bench_base.h"
#include <light_mat/matexpr/mat_arith.h>
#include <light_mat/matexpr/mat_emath.h>
#include <light_mat/matexpr/mat_special.h>

using namespace lmat;
using namespace ltest;
using namespace lmat::bench;


template<typename T>
struct bench_math1_base
{
	cref_matrix<T> a;
	mutable ref_matrix<T> dst;

	bench_math1_base(index_t m, index_t n, const T *pa, T *pd)
	: a(pa, m, n)
	, dst(pd, m, n) { }

	size_t size() const
	{
		return (size_t)(a.nelems());
	}
};


template<typename T>
struct bench_math2_base
{
	cref_matrix<T> a;
	cref_matrix<T> b;
	mutable ref_matrix<T> dst;

	bench_math2_base(index_t m, index_t n, const T *pa, const T* pb, T *pd)
	: a(pa, m, n)
	, b(pb, m, n)
	, dst(pd, m, n) { }

	size_t size() const
	{
		return (size_t)(a.nelems());
	}
};


#define DEF_BENCH_MATH1(Fun) \
		template<typename T> \
		struct bench_##Fun##_scalar : public bench_math1_base<T> { \
			bench_##Fun##_scalar(const bench_math1_base<T>& base) \
			: bench_math1_base<T>(base) { } \
			const char *name() const { return #Fun "-scalar"; } \
			LMAT_ENSURE_INLINE \
			void operator() () const { \
				const cref_matrix<T>& a = this->a; \
				ref_matrix<T>& dst = this->dst; \
				evaluate(Fun(a), dst, linear_macc<atags::scalar>()); } \
		}; \
		template<typename T> \
		struct bench_##Fun##_simd : public bench_math1_base<T> { \
			bench_##Fun##_simd(const bench_math1_base<T>& base) \
			: bench_math1_base<T>(base) { } \
			const char *name() const { return #Fun "-simd"; } \
			LMAT_ENSURE_INLINE \
			void operator() () const { \
				const cref_matrix<T>& a = this->a; \
				ref_matrix<T>& dst = this->dst; \
				dst = Fun(a); } \
		};


#define DEF_BENCH_MATH2(Fun) \
		template<typename T> \
		struct bench_##Fun##_scalar : public bench_math2_base<T> { \
			bench_##Fun##_scalar(const bench_math2_base<T>& base) \
			: bench_math2_base<T>(base) { } \
			const char *name() const { return #Fun "-scalar"; } \
			LMAT_ENSURE_INLINE \
			void operator() () const { \
				const cref_matrix<T>& a = this->a; \
				const cref_matrix<T>& b = this->b; \
				ref_matrix<T>& dst = this->dst; \
				evaluate(Fun(a, b), dst, linear_macc<atags::scalar>()); } \
		}; \
		template<typename T> \
		struct bench_##Fun##_simd : public bench_math2_base<T> { \
			bench_##Fun##_simd(const bench_math2_base<T>& base) \
			: bench_math2_base<T>(base) { } \
			const char *name() const { return #Fun "-simd"; } \
			LMAT_ENSURE_INLINE \
			void operator() () const { \
				const cref_matrix<T>& a = this->a; \
				const cref_matrix<T>& b = this->b; \
				ref_matrix<T>& dst = this->dst; \
				dst = Fun(a, b); } \
		};


#define ADD_MBENCH1(Fun) \
		run_benchmark(bench_##Fun##_scalar<T>(base1), mon, opt); \
		run_benchmark(bench_##Fun##_simd<T>(base1), mon, opt);

#define ADD_MBENCH2(Fun) \
		run_benchmark(bench_##Fun##_scalar<T>(base2), mon, opt); \
		run_benchmark(bench_##Fun##_simd<T>(base2), mon, opt);


DEF_BENCH_MATH1( sqrt )
DEF_BENCH_MATH1( cbrt )
DEF_BENCH_MATH2( pow )
DEF_BENCH_MATH2( hypot )

DEF_BENCH_MATH1( exp )
DEF_BENCH_MATH1( log )
DEF_BENCH_MATH1( log10 )
DEF_BENCH_MATH1( xlogx )
DEF_BENCH_MATH2( xlogy )
DEF_BENCH_MATH1( exp2 )
DEF_BENCH_MATH1( log2 )
DEF_BENCH_MATH1( expm1 )
DEF_BENCH_MATH1( log1p )

DEF_BENCH_MATH1( sin )
DEF_BENCH_MATH1( cos )
DEF_BENCH_MATH1( tan )
DEF_BENCH_MATH1( asin )
DEF_BENCH_MATH1( acos )
DEF_BENCH_MATH1( atan )
DEF_BENCH_MATH2( atan2 )

DEF_BENCH_MATH1( erf )
DEF_BENCH_MATH1( erfc )
DEF_BENCH_MATH1( norminv )


const index_t m = 128;
const index_t n = 128;
const size_t pbsiz = 500;

#define MATH_BENCH_TEMPL "{{jobname : %-16s}}:  {{times: %10lu}}  | {{secs: %10.4f}} s  | {{mps: %10.2f}} MPS\n"

template<typename T>
void run_bench()
{
	dense_matrix<T> a(n, n);
	dense_matrix<T> b(n, n);
	dense_matrix<T> dst(n, n, zero());
	const T *pa = a.ptr_data();
	const T *pb = b.ptr_data();
	T *pd = dst.ptr_data();
	fill_rand(a);
	fill_rand(b);

	std_bench_monitor mon(MATH_BENCH_TEMPL);
	benchmark_option opt(pbsiz);

	bench_math1_base<T> base1(m, n, pa, pd);
	bench_math2_base<T> base2(m, n, pa, pb, pd);

	std::cout << "power functions" << std::endl;
	std::cout << "---------------------" << std::endl;

	ADD_MBENCH1( sqrt )
	ADD_MBENCH1( cbrt )
	ADD_MBENCH2( pow )
	ADD_MBENCH2( hypot )

	std::cout << "\nexp & log functions" << std::endl;
	std::cout << "---------------------" << std::endl;

	ADD_MBENCH1( exp )
	ADD_MBENCH1( log )
	ADD_MBENCH1( log10 )
	ADD_MBENCH1( xlogx )
	ADD_MBENCH2( xlogy )
	ADD_MBENCH1( exp2 )
	ADD_MBENCH1( log2 )
	ADD_MBENCH1( expm1 )
	ADD_MBENCH1( log1p )

	std::cout << "\ntrigonometry functions" << std::endl;
	std::cout << "---------------------" << std::endl;

	ADD_MBENCH1( sin )
	ADD_MBENCH1( cos )
	ADD_MBENCH1( tan )
	ADD_MBENCH1( asin )
	ADD_MBENCH1( acos )
	ADD_MBENCH1( atan )
	ADD_MBENCH2( atan2 )

	std::cout << "\nspecial functions" << std::endl;
	std::cout << "---------------------" << std::endl;

	ADD_MBENCH1( erf )
	ADD_MBENCH1( erfc )
	ADD_MBENCH1( norminv )

	std::cout << std::endl;
}


int main(int argc, char *argv[])
{
	std::printf("On float (size = %d x %d)\n", (int)m, (int)n);
	std::printf("**************************************\n");
	run_bench<float>();

	std::printf("\n");

	std::printf("On double (size = %d x %d)\n", (int)m, (int)n);
	std::printf("**************************************\n");
	run_bench<double>();

	std::printf("\n");
}





