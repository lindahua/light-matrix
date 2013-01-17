/**
 * @file bench_copy.cpp
 *
 * @brief Comparison of the performance of different copying way
 *
 * @author Dahua Lin
 */

#include "bench_base.h"
#include <light_mat/mateval/ewise_eval.h>

using namespace lmat;
using namespace ltest;
using namespace lmat::bench;

// Jobs

template<typename T>
struct direct_copy
{
	const index_t nrows;
	const index_t ncols;
	const T *src;
	T *dst;

	direct_copy(index_t m, index_t n, const T* s, T *d)
	: nrows(m), ncols(n), src(s), dst(d) { }

	const char *name() const
	{
		return "direct_copy";
	}

	size_t size() const
	{
		return (size_t)(nrows * ncols);
	}

	void operator() () const
	{
		copy_vec(nrows * ncols, src, dst);
	}
};


template<typename T>
struct forloop_copy
{
	const index_t nrows;
	const index_t ncols;
	const T *src;
	T *dst;

	forloop_copy(index_t m, index_t n, const T* s, T *d)
	: nrows(m), ncols(n), src(s), dst(d) { }

	const char *name() const
	{
		return "forloop_copy";
	}

	size_t size() const
	{
		return (size_t)(nrows * ncols);
	}

	void operator() () const
	{
		index_t N = nrows * ncols;
		for (index_t i = 0; i < N; ++i) dst[i] = src[i];
	}
};


template<typename T>
struct matrix_copy
{
	cref_matrix<T> src;
	mutable ref_matrix<T> dst;

	matrix_copy(index_t m, index_t n, const T* s, T *d)
	: src(s, m, n), dst(d, m, n) { }

	const char *name() const
	{
		return "matrix_copy";
	}

	size_t size() const
	{
		return (size_t)src.nelems();
	}

	void operator() () const
	{
		copy(src, dst);
	}
};


template<typename T>
struct percol_copy
{
	cref_matrix<T> src;
	mutable ref_matrix<T> dst;

	percol_copy(index_t m, index_t n, const T* s, T *d)
	: src(s, m, n), dst(d, m, n) { }

	const char *name() const
	{
		return "percol_copy";
	}

	size_t size() const
	{
		return (size_t)src.nelems();
	}

	void operator() () const
	{
		const index_t n = src.ncolumns();
		for (index_t j = 0; j < n; ++j)
		{
			ref_col<T> dcol = dst.column(j);
			copy(src.column(j), dcol);
		}
	}
};


template<typename T>
struct percol_assign
{
	cref_matrix<T> src;
	mutable ref_matrix<T> dst;

	percol_assign(index_t m, index_t n, const T* s, T *d)
	: src(s, m, n), dst(d, m, n) { }

	const char *name() const
	{
		return "percol_assign";
	}

	size_t size() const
	{
		return (size_t)src.nelems();
	}

	void operator() () const
	{
		const index_t n = src.ncolumns();
		for (index_t j = 0; j < n; ++j)
		{
			dst.column(j) = src.column(j);
		}
	}
};


template<typename T>
struct perelem_assign
{
	cref_matrix<T> src;
	mutable ref_matrix<T> dst;

	perelem_assign(index_t m, index_t n, const T* s, T *d)
	: src(s, m, n), dst(d, m, n) { }

	const char *name() const
	{
		return "perelem_assign";
	}

	size_t size() const
	{
		return (size_t)src.nelems();
	}

	void operator() () const
	{
		const index_t m = src.nrows();
		const index_t n = src.ncolumns();

		for (index_t j = 0; j < n; ++j)
		{
			for (index_t i = 0; i < m; ++i)
				dst(i, j) = src(i, j);
		}
	}
};



template<typename T>
struct linearscalar_copy
{
	cref_matrix<T> src;
	mutable ref_matrix<T> dst;

	linearscalar_copy(index_t m, index_t n, const T* s, T *d)
	: src(s, m, n), dst(d, m, n) { }

	const char *name() const
	{
		return "linearscalar_copy";
	}

	size_t size() const
	{
		return (size_t)src.nelems();
	}

	void operator() () const
	{
		typedef atags::scalar tag;
		ewise(copy_kernel<T>(), tag())(src.nelems(), in_(src), out_(dst));
	}
};


template<typename T>
struct linearsimd_copy
{
	cref_matrix<T> src;
	mutable ref_matrix<T> dst;

	linearsimd_copy(index_t m, index_t n, const T* s, T *d)
	: src(s, m, n), dst(d, m, n) { }

	const char *name() const
	{
		return "linearsimd_copy";
	}

	size_t size() const
	{
		return (size_t)src.nelems();
	}

	void operator() () const
	{
		typedef atags::simd<default_simd_kind> tag;
		ewise(copy_kernel<T>(), tag())(src.nelems(), in_(src), out_(dst));
	}
};


template<typename T>
struct percolscalar_copy
{
	cref_matrix<T> src;
	mutable ref_matrix<T> dst;

	percolscalar_copy(index_t m, index_t n, const T* s, T *d)
	: src(s, m, n), dst(d, m, n) { }

	const char *name() const
	{
		return "percolscalar_copy";
	}

	size_t size() const
	{
		return (size_t)src.nelems();
	}

	void operator() () const
	{
		typedef atags::scalar tag;
		percol(ewise(copy_kernel<T>(), tag()), src.shape(), in_(src), out_(dst));
	}
};


template<typename T>
struct percolsimd_copy
{
	cref_matrix<T> src;
	mutable ref_matrix<T> dst;

	percolsimd_copy(index_t m, index_t n, const T* s, T *d)
	: src(s, m, n), dst(d, m, n) { }

	const char *name() const
	{
		return "percolsimd_copy";
	}

	size_t size() const
	{
		return (size_t)src.nelems();
	}

	void operator() () const
	{
		typedef atags::simd<default_simd_kind> tag;
		percol(ewise(copy_kernel<T>(), tag()), src.shape(), in_(src), out_(dst));
	}
};


index_t sizes[] = {2, 5, 10, 20, 50, 100, 200, 500, 1000, 2000 };
const size_t nsizes = sizeof(sizes) / sizeof(index_t);


template<typename T>
void run_bench()
{
	dense_matrix<T> src(2000, 2000);
	dense_matrix<T> dst(2000, 2000, zero());
	const T *ps = src.ptr_data();
	T *pd = dst.ptr_data();
	fill_rand(src);

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

		direct_copy<T> job_direct_copy(m, n, ps, pd);
		run_benchmark(job_direct_copy, mon, opt);

		forloop_copy<T> job_forloop_copy(m, n, ps, pd);
		run_benchmark(job_forloop_copy, mon, opt);

		matrix_copy<T> job_matrix_copy(m, n, ps, pd);
		run_benchmark(job_matrix_copy, mon, opt);

		percol_copy<T> job_percol_copy(m, n, ps, pd);
		run_benchmark(job_percol_copy, mon, opt);

		percol_assign<T> job_percol_assign(m, n, ps, pd);
		run_benchmark(job_percol_assign, mon, opt);

		perelem_assign<T> job_perelem_assign(m, n, ps, pd);
		run_benchmark(job_perelem_assign, mon, opt);

		linearscalar_copy<T> job_linearscalar_copy(m, n, ps, pd);
		run_benchmark(job_linearscalar_copy, mon, opt);

		linearsimd_copy<T> job_linearsimd_copy(m, n, ps, pd);
		run_benchmark(job_linearsimd_copy, mon, opt);

		percolscalar_copy<T> job_percolscalar_copy(m, n, ps, pd);
		run_benchmark(job_percolscalar_copy, mon, opt);

		percolsimd_copy<T> job_percolsimd_copy(m, n, ps, pd);
		run_benchmark(job_percolsimd_copy, mon, opt);

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


