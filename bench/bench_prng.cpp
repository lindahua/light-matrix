/**
 * @file bench_prng.cpp
 *
 * Benchmark of Pseudo Random Number Generators
 * 
 * @author Dahua Lin 
 */

#include "bench_base.h"
#include <light_mat/random/distributions.h>

using namespace lmat;
using namespace ltest;
using namespace lmat::random;
using namespace lmat::bench;

template<typename T>
struct bench_prng_base
{
	const char *_name;
	index_t _n;
	T *dst;

	bench_prng_base(const index_t n, T *d)
	: _name(0), _n(n), dst(d)
	{ }

	const char *name() const
	{
		return _name;
	}

	size_t size() const
	{
		return (size_t)(_n);
	}
};


default_rand_stream rstream;

template<class Distr>
struct bench_prng_scalar
: public bench_prng_base<typename Distr::result_type>
{
	typedef bench_prng_base<typename Distr::result_type> base_t;
	Distr distr;

	bench_prng_scalar(const Distr& d, const char *name, const base_t& base)
	: base_t(base), distr(d) { this->_name = name; }

	LMAT_ENSURE_INLINE
	void operator() () const
	{
		typename Distr::result_type *pd = this->dst;
		const index_t n = this->_n;

		for (index_t i = 0; i < n; ++i)
		{
			pd[i] = distr(rstream);
		}
	}
};


template<class Distr>
struct bench_prng_simd
: public bench_prng_base<typename Distr::result_type>
{
	typedef bench_prng_base<typename Distr::result_type> base_t;
	typedef typename simdize_map<Distr, default_simd_kind>::type simd_distr_t;

	simd_distr_t sdistr;

	bench_prng_simd(const Distr& d, const char *name, const base_t& base)
	: base_t(base)
	, sdistr( simdize_map<Distr, default_simd_kind>::get(d) )
	{ this->_name = name; }

	LMAT_ENSURE_INLINE
	void operator() () const
	{
		typename Distr::result_type *pd = this->dst;
		const index_t n = this->_n;

		typedef typename Distr::result_type T;
		typedef simd_pack<T, default_simd_kind> pack_t;
		const index_t pw = pack_t::pack_width;

		for (index_t i = 0; i < n; i += pw)
		{
			pack_t pk = sdistr(rstream);
			pk.store_u(pd + i);
		}
	}
};



const index_t Length = 1024;  // Length must be multiples of 16

#define ADD_DISTR_BENCH_P0( Name ) \
	run_benchmark(bench_prng_scalar<Name##_distr<T> >(Name##_distr<T>( ), #Name "-scalar", base), mon, opt)

#define ADD_DISTR_BENCH_P1( Name, P1 ) \
	run_benchmark(bench_prng_scalar<Name##_distr<T> >(Name##_distr<T>( P1 ), #Name "-scalar", base), mon, opt)

#define ADD_DISTR_BENCH_P2( Name, P1, P2 ) \
	run_benchmark(bench_prng_scalar<Name##_distr<T> >(Name##_distr<T>( P1, P2 ), #Name "-scalar", base), mon, opt)

#define ADD_DISTR_SIMD_BENCH_P0( Name ) \
	run_benchmark(bench_prng_simd<Name##_distr<T> >(Name##_distr<T>( ), #Name "-simd", base), mon, opt)

#define ADD_DISTR_SIMD_BENCH_P1( Name, P1 ) \
	run_benchmark(bench_prng_simd<Name##_distr<T> >(Name##_distr<T>( P1 ), #Name "-simd", base), mon, opt)

#define ADD_DISTR_SIMD_BENCH_P2( Name, P1, P2 ) \
	run_benchmark(bench_prng_simd<Name##_distr<T> >(Name##_distr<T>( P1, P2 ), #Name "-simd", base), mon, opt)


void bench_discrete_distrs()
{
	typedef uint32_t T;
	dense_col<T> dst(Length);

	std_bench_monitor mon;
	size_t pbsiz = 1024000 / Length;
	benchmark_option opt(pbsiz);

	bench_prng_base<T> base(dst.nelems(), dst.ptr_data());

	ADD_DISTR_BENCH_P1( std_uniform_int, 100 );
	ADD_DISTR_BENCH_P2( uniform_int, 12, 36 );
	ADD_DISTR_BENCH_P2( binomial, 5, 0.4 );
	ADD_DISTR_BENCH_P1( geometric, 0.4 );
}

template<typename T>
void bench_real_distrs()
{
	dense_col<T> dst(Length);

	std_bench_monitor mon;
	size_t pbsiz = 1024000 / Length;
	benchmark_option opt(pbsiz);

	bench_prng_base<T> base(dst.nelems(), dst.ptr_data());

	std::cout << "uniform:\n";
	std::cout << "---------------------\n";

	ADD_DISTR_BENCH_P0( std_uniform_real );
	ADD_DISTR_SIMD_BENCH_P0( std_uniform_real );

	ADD_DISTR_BENCH_P2( uniform_real, T(1.2), T(4.8) );
	ADD_DISTR_SIMD_BENCH_P2( uniform_real, T(1.2), T(4.8) );

	std::cout << "\n";
	std::cout << "exponential:\n";
	std::cout << "---------------------\n";

	ADD_DISTR_BENCH_P0( std_exponential );
	ADD_DISTR_SIMD_BENCH_P0( std_exponential );

	ADD_DISTR_BENCH_P1( exponential, T(2.5) );
	ADD_DISTR_SIMD_BENCH_P1( exponential, T(2.5) );

	std::cout << "\n";
	std::cout << "normal:\n";
	std::cout << "---------------------\n";

	ADD_DISTR_BENCH_P0( std_normal );
	ADD_DISTR_SIMD_BENCH_P0( std_normal );

	ADD_DISTR_BENCH_P2( normal, T(1.6), T(2.5) );
	ADD_DISTR_SIMD_BENCH_P2( normal, T(1.6), T(2.5) );

	std::cout << "\n";
	std::cout << "gamma:\n";
	std::cout << "---------------------\n";

	ADD_DISTR_BENCH_P1( std_gamma, T(2.0) );
	ADD_DISTR_BENCH_P2( gamma, T(2.0), T(2.5) );
}


int main(int argc, char *argv[])
{
	std::cout << "Discrete distributions [uint32_t]\n";
	std::cout << "**************************************\n";
	bench_discrete_distrs();
	std::cout << "\n";

	std::cout << "Continuous distributions [float]\n";
	std::cout << "**************************************\n";
	bench_real_distrs<float>();
	std::cout << "\n";

	std::cout << "Continuous distributions [double]\n";
	std::cout << "**************************************\n";
	bench_real_distrs<double>();
	std::cout << "\n";
}

