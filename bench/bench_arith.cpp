/**
 * @file bench_arith.cpp
 *
 * Benchmark of different ways of doing arithmetics
 *
 * reference formula
 *
 *   cube(x + y) + sqrt(sqr(x) + sqr(y) + x * y);
 * 
 * @author Dahua Lin 
 */


#include "bench_base.h"
#include <light_mat/matexpr/mat_arith.h>

using namespace lmat;
using namespace ltest;
using namespace lmat::bench;

template<typename T>
LMAT_ENSURE_INLINE
inline void calc(const T& x, const T& y, T& r)
{
	r = math::cube(x + y) + math::sqrt(math::sqr(x) + math::sqr(y)) + x * y;
}

template<typename T>
struct calc_kernel
{
	typedef T value_type;

	LMAT_ENSURE_INLINE
	void operator()(const T& x, const T& y, T& r) const
	{
		calc(x, y, r);
	}
};

namespace lmat
{
	LMAT_DEF_SIMD_SUPPORT(calc_kernel)
}

template<typename T>
struct bench_arith_base
{
	cref_matrix<T> a;
	cref_matrix<T> b;
	mutable ref_matrix<T> dst;

	bench_arith_base(index_t m, index_t n, const T *pa, const T* pb, T *pd)
	: a(pa, m, n)
	, b(pb, m, n)
	, dst(pd, m, n) { }

	size_t size() const
	{
		return (size_t)(a.nelems());
	}
};

template<typename T>
struct bench_rawloop : public bench_arith_base<T>
{
	bench_rawloop(const bench_arith_base<T>& base)
	: bench_arith_base<T>(base) { }

	const char *name() const { return "arith-rawloop"; }

	LMAT_ENSURE_INLINE
	void operator() () const
	{
		const T* pa = this->a.ptr_data();
		const T *pb = this->b.ptr_data();
		T *pd = this->dst.ptr_data();

		const index_t N = this->a.nelems();

		for (index_t i = 0; i < N; ++i)
			calc(pa[i], pb[i], pd[i]);
	}
};


template<typename T>
struct bench_rawsimd : public bench_arith_base<T>
{
	bench_rawsimd(const bench_arith_base<T>& base)
	: bench_arith_base<T>(base) { }

	const char *name() const { return "arith-rawsimd"; }

	LMAT_ENSURE_INLINE
	void operator() () const
	{
		const T* pa = this->a.ptr_data();
		const T *pb = this->b.ptr_data();
		T *pd = this->dst.ptr_data();

		typedef simd_pack<T, default_simd_kind> pack_t;
		pack_t x, y, z;

		const index_t N = this->a.nelems();
		const index_t W = (index_t)pack_t::pack_width;

		for (index_t i = 0; i < N; ++i, i += W)
		{
			x.load_u(pa + i);
			y.load_u(pb + i);
			calc(x, y, z);
			z.store_u(pd + i);
		}
	}
};


template<typename T>
struct bench_kernel_linear_scalar : public bench_arith_base<T>
{
	bench_kernel_linear_scalar(const bench_arith_base<T>& base)
	: bench_arith_base<T>(base) { }

	const char *name() const { return "arith-kernel-linear-scalar"; }

	LMAT_ENSURE_INLINE
	void operator() () const
	{
		const index_t N = this->a.nelems();
		ewise(calc_kernel<T>(), atags::scalar())(N,
				in_(this->a), in_(this->b), out_(this->dst));
	}
};

template<typename T>
struct bench_kernel_linear_simd : public bench_arith_base<T>
{
	bench_kernel_linear_simd(const bench_arith_base<T>& base)
	: bench_arith_base<T>(base) { }

	const char *name() const { return "arith-kernel-linear-simd"; }

	LMAT_ENSURE_INLINE
	void operator() () const
	{
		const index_t N = this->a.nelems();
		ewise(calc_kernel<T>())(N,
				in_(this->a), in_(this->b), out_(this->dst));
	}
};


template<typename T>
struct bench_kernel_percol_scalar : public bench_arith_base<T>
{
	bench_kernel_percol_scalar(const bench_arith_base<T>& base)
	: bench_arith_base<T>(base) { }

	const char *name() const { return "arith-kernel-percol-scalar"; }

	LMAT_ENSURE_INLINE
	void operator() () const
	{
		const index_t m = this->a.nrows();
		const index_t n = this->a.ncolumns();
		percol(ewise(calc_kernel<T>(), atags::scalar()), m, n,
				in_(this->a), in_(this->b), out_(this->dst));
	}
};

template<typename T>
struct bench_kernel_percol_simd : public bench_arith_base<T>
{
	bench_kernel_percol_simd(const bench_arith_base<T>& base)
	: bench_arith_base<T>(base) { }

	const char *name() const { return "arith-kernel-percol-simd"; }

	LMAT_ENSURE_INLINE
	void operator() () const
	{
		const index_t m = this->a.nrows();
		const index_t n = this->a.ncolumns();
		percol(ewise(calc_kernel<T>()), m, n,
				in_(this->a), in_(this->b), out_(this->dst));
	}
};

// r = cube(x + y) + sqrt(sqr(x) + sqr(y)) + x * y;

template<typename T>
struct bench_expr_linear_scalar : public bench_arith_base<T>
{
	bench_expr_linear_scalar(const bench_arith_base<T>& base)
	: bench_arith_base<T>(base) { }

	const char *name() const { return "arith-expr-linear-scalar"; }

	LMAT_ENSURE_INLINE
	void operator() () const
	{
		const cref_matrix<T>& a = this->a;
		const cref_matrix<T>& b = this->b;
		ref_matrix<T>& dst = this->dst;

		linear_macc<atags::scalar> policy;
		evaluate(cube(a + b) + sqrt(sqr(a) + sqr(b)) + a * b, dst, policy);
	}
};

template<typename T>
struct bench_expr_linear_simd : public bench_arith_base<T>
{
	bench_expr_linear_simd(const bench_arith_base<T>& base)
	: bench_arith_base<T>(base) { }

	const char *name() const { return "arith-expr-linear-simd"; }

	LMAT_ENSURE_INLINE
	void operator() () const
	{
		const cref_matrix<T>& a = this->a;
		const cref_matrix<T>& b = this->b;
		ref_matrix<T>& dst = this->dst;

		linear_macc<atags::simd<default_simd_kind> > policy;
		evaluate(cube(a + b) + sqrt(sqr(a) + sqr(b)) + a * b, dst, policy);
	}
};


template<typename T>
struct bench_expr_percol_scalar : public bench_arith_base<T>
{
	bench_expr_percol_scalar(const bench_arith_base<T>& base)
	: bench_arith_base<T>(base) { }

	const char *name() const { return "arith-expr-percol-scalar"; }

	LMAT_ENSURE_INLINE
	void operator() () const
	{
		const cref_matrix<T>& a = this->a;
		const cref_matrix<T>& b = this->b;
		ref_matrix<T>& dst = this->dst;

		percol_macc<atags::scalar> policy;
		evaluate(cube(a + b) + sqrt(sqr(a) + sqr(b)) + a * b, dst, policy);
	}
};

template<typename T>
struct bench_expr_percol_simd : public bench_arith_base<T>
{
	bench_expr_percol_simd(const bench_arith_base<T>& base)
	: bench_arith_base<T>(base) { }

	const char *name() const { return "arith-expr-percol-simd"; }

	LMAT_ENSURE_INLINE
	void operator() () const
	{
		const cref_matrix<T>& a = this->a;
		const cref_matrix<T>& b = this->b;
		ref_matrix<T>& dst = this->dst;

		percol_macc<atags::simd<default_simd_kind> > policy;
		evaluate(cube(a + b) + sqrt(sqr(a) + sqr(b)) + a * b, dst, policy);
	}
};


template<typename T>
struct bench_expr_default : public bench_arith_base<T>
{
	bench_expr_default(const bench_arith_base<T>& base)
	: bench_arith_base<T>(base) { }

	const char *name() const { return "arith-expr-default"; }

	LMAT_ENSURE_INLINE
	void operator() () const
	{
		const cref_matrix<T>& a = this->a;
		const cref_matrix<T>& b = this->b;
		ref_matrix<T>& dst = this->dst;

		dst = cube(a + b) + sqrt(sqr(a) + sqr(b)) + a * b;
	}
};



index_t max_size = 2048;
index_t sizes[] = {4, 8, 16, 32, 64, 128, 256, 512, 1024 };
const size_t nsizes = sizeof(sizes) / sizeof(index_t);


template<typename T>
void run_bench()
{
	dense_matrix<T> a(max_size, max_size);
	dense_matrix<T> b(max_size, max_size);
	dense_matrix<T> dst(max_size, max_size, zero());
	const T *pa = a.ptr_data();
	const T *pb = b.ptr_data();
	T *pd = dst.ptr_data();
	fill_rand(a);
	fill_rand(b);

	std_bench_monitor mon;

	for (size_t k = 0; k < nsizes; ++k)
	{
		index_t siz = sizes[k];
		index_t m = siz;
		index_t n = siz;
		size_t pbsiz = 20000000 / size_t(m * n);

		benchmark_option opt(pbsiz);

		std::cout << "size = " << m << " x " << n << "\n";
		std::cout << "=======================================\n";

		bench_arith_base<T> base(m, n, pa, pb, pd);

		run_benchmark(bench_rawloop<T>(base), mon, opt);
		run_benchmark(bench_rawsimd<T>(base), mon, opt);

		run_benchmark(bench_kernel_linear_scalar<T>(base), mon, opt);
		run_benchmark(bench_kernel_linear_simd<T>  (base), mon, opt);
		run_benchmark(bench_kernel_percol_scalar<T>(base), mon, opt);
		run_benchmark(bench_kernel_percol_simd<T>  (base), mon, opt);

		run_benchmark(bench_expr_linear_scalar<T>(base), mon, opt);
		run_benchmark(bench_expr_linear_simd<T>(base), mon, opt);
		run_benchmark(bench_expr_percol_scalar<T>(base), mon, opt);
		run_benchmark(bench_expr_percol_simd<T>(base), mon, opt);
		run_benchmark(bench_expr_default<T>(base), mon, opt);

		std::cout << "\n";
	}
}


int main(int argc, char *argv[])
{
	std::printf("On float\n");
	std::printf("**************************************\n");
	run_bench<float>();

	std::printf("\n");

	std::printf("On double\n");
	std::printf("**************************************\n");
	run_bench<double>();

	std::printf("\n");
}

