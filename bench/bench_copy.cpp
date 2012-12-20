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
	cref_matrix<T> src;
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
	cref_matrix<T> src;
	ref_matrix<T> dst;

	percol_copy(index_t m, index_t n, const T* s, T *d)
	: src(s, m, n), dst(d, m, n) { }

	const char *name() const
	{
		return "percol_copy";
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
			ref_col<T> dcol = dst.column(j);
			copy(src.column(j), dcol);
		}
	}
};


template<typename T>
struct percol_assign
{
	cref_matrix<T> src;
	ref_matrix<T> dst;

	percol_assign(index_t m, index_t n, const T* s, T *d)
	: src(s, m, n), dst(d, m, n) { }

	const char *name() const
	{
		return "percol_assign";
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
			dst.column(j) = src.column(j);
		}
	}
};


template<typename T>
struct perelem_assign
{
	cref_matrix<T> src;
	ref_matrix<T> dst;

	perelem_assign(index_t m, index_t n, const T* s, T *d)
	: src(s, m, n), dst(d, m, n) { }

	const char *name() const
	{
		return "perelem_assign";
	}

	unsigned int size() const
	{
		return (unsigned int)src.nelems();
	}

	void run()
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
	ref_matrix<T> dst;

	linearscalar_copy(index_t m, index_t n, const T* s, T *d)
	: src(s, m, n), dst(d, m, n) { }

	const char *name() const
	{
		return "linearscalar_copy";
	}

	unsigned int size() const
	{
		return (unsigned int)src.nelems();
	}

	void run()
	{
		typedef atags::scalar tag;
		ewise(copy_kernel<T>(), tag())(src.nelems(), in_(src), out_(dst));
	}
};


template<typename T>
struct linearsimd_copy
{
	cref_matrix<T> src;
	ref_matrix<T> dst;

	linearsimd_copy(index_t m, index_t n, const T* s, T *d)
	: src(s, m, n), dst(d, m, n) { }

	const char *name() const
	{
		return "linearsimd_copy";
	}

	unsigned int size() const
	{
		return (unsigned int)src.nelems();
	}

	void run()
	{
		typedef atags::simd<default_simd_kind> tag;
		ewise(copy_kernel<T>(), tag())(src.nelems(), in_(src), out_(dst));
	}
};


template<typename T>
struct percolscalar_copy
{
	cref_matrix<T> src;
	ref_matrix<T> dst;

	percolscalar_copy(index_t m, index_t n, const T* s, T *d)
	: src(s, m, n), dst(d, m, n) { }

	const char *name() const
	{
		return "percolscalar_copy";
	}

	unsigned int size() const
	{
		return (unsigned int)src.nelems();
	}

	void run()
	{
		typedef atags::scalar tag;
		percol(ewise(copy_kernel<T>(), tag()), src.shape(), in_(src), out_(dst));
	}
};


template<typename T>
struct percolsimd_copy
{
	cref_matrix<T> src;
	ref_matrix<T> dst;

	percolsimd_copy(index_t m, index_t n, const T* s, T *d)
	: src(s, m, n), dst(d, m, n) { }

	const char *name() const
	{
		return "percolsimd_copy";
	}

	unsigned int size() const
	{
		return (unsigned int)src.nelems();
	}

	void run()
	{
		typedef atags::simd<default_simd_kind> tag;
		percol(ewise(copy_kernel<T>(), tag()), src.shape(), in_(src), out_(dst));
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
		      500,
		      100};

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
		time_job(job_direct_copy, 1, nrun, PUNIT_GPS);

		forloop_copy<T> job_forloop_copy(m, n, ps, pd);
		time_job(job_forloop_copy, 1, nrun, PUNIT_GPS);

		matrix_copy<T> job_matrix_copy(m, n, ps, pd);
		time_job(job_matrix_copy, 1, nrun, PUNIT_GPS);

		percol_copy<T> job_percol_copy(m, n, ps, pd);
		time_job(job_percol_copy, 1, nrun, PUNIT_GPS);

		percol_assign<T> job_percol_assign(m, n, ps, pd);
		time_job(job_percol_assign, 1, nrun, PUNIT_GPS);

		perelem_assign<T> job_perelem_assign(m, n, ps, pd);
		time_job(job_perelem_assign, 1, nrun, PUNIT_GPS);

		linearscalar_copy<T> job_linearscalar_copy(m, n, ps, pd);
		time_job(job_linearscalar_copy, 1, nrun, PUNIT_GPS);

		linearsimd_copy<T> job_linearsimd_copy(m, n, ps, pd);
		time_job(job_linearsimd_copy, 1, nrun, PUNIT_GPS);

		percolscalar_copy<T> job_percolscalar_copy(m, n, ps, pd);
		time_job(job_percolscalar_copy, 1, nrun, PUNIT_GPS);

		percolsimd_copy<T> job_percolsimd_copy(m, n, ps, pd);
		time_job(job_percolsimd_copy, 1, nrun, PUNIT_GPS);
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


