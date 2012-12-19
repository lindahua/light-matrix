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

	unsigned int size() const
	{
		return (unsigned int)(nrows * ncols);
	}

	void run()
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

	unsigned int size() const
	{
		return (unsigned int)(nrows * ncols);
	}

	void run()
	{
		index_t N = nrows * ncols;
		for (index_t i = 0; i < N; ++i) dst[i] = src[i];
	}
};


template<typename T>
struct matrix_copy
{
	ref_matrix<T> src;
	ref_matrix<T> dst;

	matrix_copy(index_t m, index_t n, const T* s, T *d)
	: src(s, m, n), dst(d, m, n) { }

	const char *name() const
	{
		return "matrix_copy";
	}

	unsigned int size() const
	{
		return (unsigned int)src.nelems();
	}

	void run()
	{
		copy(src, dst);
	}
};


template<typename T>
struct percol_copy
{
	ref_matrix<T> src;
	ref_matrix<T> dst;

	percol_copy(index_t m, index_t n, const T* s, T *d)
	: src(s, m, n), dst(d, m, n) { }

	const char *name() const
	{
		return "matrix_copy";
	}

	unsigned int size() const
	{
		return (unsigned int)src.nelems();
	}

	void run()
	{
		const index_t n = src.ncolumns();
		for (index_t j = 0; j < n; ++j)
		{
			copy(src.column(j), dst.column(j));
		}
	}
};


template<typename T>
struct linearscalar_copy
{
	ref_matrix<T> src;
	ref_matrix<T> dst;

	linearscalar_copy(index_t m, index_t n, const T* s, T *d)
	: src(s, m, n), dst(d, m, n) { }

	const char *name() const
	{
		return "lineareval_copy";
	}

	unsigned int size() const
	{
		return (unsigned int)src.nelems();
	}

	void run()
	{
		typedef atags::scalar tag;
		linear_ewise(tag(), src.nelems()).copy(in_(src), out_(dst));
	}
};


template<typename T>
struct linearsimd_copy
{
	ref_matrix<T> src;
	ref_matrix<T> dst;

	linearsimd_copy(index_t m, index_t n, const T* s, T *d)
	: src(s, m, n), dst(d, m, n) { }

	const char *name() const
	{
		return "lineareval_copy";
	}

	unsigned int size() const
	{
		return (unsigned int)src.nelems();
	}

	void run()
	{
		typedef atags::simd<T, default_simd_kind> tag;
		linear_ewise(tag(), src.nelems()).copy(in_(src), out_(dst));
	}
};


template<typename T>
struct percolscalar_copy
{
	ref_matrix<T> src;
	ref_matrix<T> dst;

	percolscalar_copy(index_t m, index_t n, const T* s, T *d)
	: src(s, m, n), dst(d, m, n) { }

	const char *name() const
	{
		return "lineareval_copy";
	}

	unsigned int size() const
	{
		return (unsigned int)src.nelems();
	}

	void run()
	{
		typedef atags::scalar tag;
		percol_ewise(tag(), src.nelems()).copy(in_(src), out_(dst));
	}
};


template<typename T>
struct percolsimd_copy
{
	ref_matrix<T> src;
	ref_matrix<T> dst;

	percolsimd_copy(index_t m, index_t n, const T* s, T *d)
	: src(s, m, n), dst(d, m, n) { }

	const char *name() const
	{
		return "lineareval_copy";
	}

	unsigned int size() const
	{
		return (unsigned int)src.nelems();
	}

	void run()
	{
		typedef atags::simd<T, default_simd_kind> tag;
		percol_ewise(tag(), src.nelems()).copy(in_(src), out_(dst));
	}
};


index_t sizes[] = {2, 5, 10, 20, 50, 100, 200, 500, 1000, 2000 };
unsigned int nruns[] = {
		100000000,
		 40000000,
		 10000000,
		  2500000,
		   400000,
		   100000,
		    25000,
		     4000,
		     1000,
		      250};

const unsigned int nsizes = sizeof(sizes) / sizeof(index_t);


template<typename T>
void run_bench()
{
	dense_matrix<T> src(2000, 2000);
	dense_matrix<T> dst(2000, 2000, zero());
	const T *ps = src.ptr_data();
	T *pd = dst.ptr_data();

	fill_rand(src);

	for (unsigned int k = 0; k < nsizes; ++k)
	{
		index_t siz = sizes[k];
		unsigned int nrun = nruns[k];
		index_t m = siz;
		index_t n = siz;

		std::printf("size = %d x %d:\n", m, n);

		direct_copy<T> job_direct_copy(m, n, ps, pd);
		time_job(job_direct_copy, 1, nrun, PUNIT_MPS);
	}
}


int main(int argc, char *argv[])
{
	run_bench<double>();
}


