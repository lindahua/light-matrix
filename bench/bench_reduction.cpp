/**
 * @file bench_reduction.cpp
 *
 * @brief Benchmarking of reduction
 *
 * @author Dahua Lin
 */


#include "bench_base.h"
#include <light_mat/mateval/mat_reduce.h>

using namespace lmat;
using namespace ltest;
using namespace lmat::bench;


template<typename T>
struct bench_fullreduc
{
	const char *_name;
	cref_matrix<T> a;
	cref_matrix<T> b;
	mutable volatile T res;

	bench_fullreduc(index_t m, index_t n, const T *pa, const T *pb)
	: _name(0), a(pa, m, n), b(pb, m, n), res(0) { }

	void set_name(const char *name) { _name = name; }

	const char *name() const
	{
		return _name;
	}

	size_t size() const
	{
		return (size_t)(a.nelems());
	}

	LMAT_ENSURE_INLINE
	void force_res(T r) const
	{
		res = r;
	}
};


template<typename T>
struct bench_cwreduc
{
	const char *_name;
	cref_matrix<T> a;
	cref_matrix<T> b;
	mutable ref_matrix<T> r;
	mutable volatile T res;

	bench_cwreduc(index_t m, index_t n, const T *pa, const T *pb, T *pr)
	: _name(0), a(pa, m, n), b(pb, m, n), r(pr, 1, n) { }

	void set_name(const char *name) { _name = name; }

	const char *name() const
	{
		return _name;
	}

	size_t size() const
	{
		return (size_t)(a.nelems());
	}

	LMAT_ENSURE_INLINE
	void force_res() const
	{
		res = r[0];
	}
};

template<typename T>
struct bench_rwreduc
{
	const char *_name;
	cref_matrix<T> a;
	cref_matrix<T> b;
	mutable ref_matrix<T> r;
	mutable volatile T res;

	bench_rwreduc(index_t m, index_t n, const T *pa, const T *pb, T *pr)
	: _name(0), a(pa, m, n), b(pb, m, n), r(pr, m, 1) { }

	void set_name(const char *name) { _name = name; }

	const char *name() const
	{
		return _name;
	}

	size_t size() const
	{
		return (size_t)(a.nelems());
	}

	LMAT_ENSURE_INLINE
	void force_res() const
	{
		res = r[0];
	}
};


template<typename T>
struct bench_full_sum_rawloop : public bench_fullreduc<T>
{
	bench_full_sum_rawloop( const bench_fullreduc<T>& base )
	: bench_fullreduc<T>(base)
	{ this->set_name("full-sum-rawloop");  }

	LMAT_ENSURE_INLINE
	void operator() () const
	{
		index_t N = this->a.nelems();
		const T *pa = this->a.ptr_data();
		T s(0);
		for (index_t i = 0; i < N; ++i) s += pa[i];
		this->force_res(s);
	}
};

template<typename T>
struct bench_full_sum_eval : public bench_fullreduc<T>
{
	bench_full_sum_eval( const bench_fullreduc<T>& base )
	: bench_fullreduc<T>(base)
	{ this->set_name("full-sum-eval");  }

	LMAT_ENSURE_INLINE
	void operator() () const
	{
		T s = sum(this->a);
		this->force_res(s);
	}
};

template<typename T>
struct bench_full_dot_rawloop : public bench_fullreduc<T>
{
	bench_full_dot_rawloop( const bench_fullreduc<T>& base )
	: bench_fullreduc<T>(base)
	{ this->set_name("full-dot-rawloop");  }

	LMAT_ENSURE_INLINE
	void operator() () const
	{
		index_t N = this->a.nelems();
		const T *pa = this->a.ptr_data();
		const T *pb = this->b.ptr_data();
		T s(0);
		for (index_t i = 0; i < N; ++i) s += pa[i] * pb[i];
		this->force_res(s);
	}
};

template<typename T>
struct bench_full_dot_eval : public bench_fullreduc<T>
{
	bench_full_dot_eval( const bench_fullreduc<T>& base )
	: bench_fullreduc<T>(base)
	{ this->set_name("full-dot-eval");  }

	LMAT_ENSURE_INLINE
	void operator() () const
	{
		T s = dot(this->a, this->b);
		this->force_res(s);
	}
};


template<typename T>
struct bench_colwise_sum_rawloop : public bench_cwreduc<T>
{
	bench_colwise_sum_rawloop( const bench_cwreduc<T>& base )
	: bench_cwreduc<T>(base)
	{ this->set_name("colwise-sum-rawloop");  }

	LMAT_ENSURE_INLINE
	void operator() () const
	{
		index_t m = this->a.nrows();
		index_t n = this->a.ncolumns();
		T *pr = this->r.ptr_data();

		for (index_t j = 0; j < n; ++j)
		{
			const T *pa = this->a.ptr_col(j);
			T s(0);
			for (index_t i = 0; i < m; ++i)
				s += pa[i];
			pr[j] = s;
		}
		this->force_res();
	}
};

template<typename T>
struct bench_colwise_sum_eval : public bench_cwreduc<T>
{
	bench_colwise_sum_eval( const bench_cwreduc<T>& base )
	: bench_cwreduc<T>(base)
	{ this->set_name("colwise-sum-eval");  }

	LMAT_ENSURE_INLINE
	void operator() () const
	{
		colwise_sum(this->a, this->r);
		this->force_res();
	}
};


template<typename T>
struct bench_colwise_dot_rawloop : public bench_cwreduc<T>
{
	bench_colwise_dot_rawloop( const bench_cwreduc<T>& base )
	: bench_cwreduc<T>(base)
	{ this->set_name("colwise-dot-rawloop");  }

	LMAT_ENSURE_INLINE
	void operator() () const
	{
		index_t m = this->a.nrows();
		index_t n = this->a.ncolumns();
		T *pr = this->r.ptr_data();

		for (index_t j = 0; j < n; ++j)
		{
			const T *pa = this->a.ptr_col(j);
			const T *pb = this->b.ptr_col(j);
			T s(0);
			for (index_t i = 0; i < m; ++i)
				s += pa[i] * pb[i];
			pr[j] = s;
		}
		this->force_res();
	}
};

template<typename T>
struct bench_colwise_dot_eval : public bench_cwreduc<T>
{
	bench_colwise_dot_eval( const bench_cwreduc<T>& base )
	: bench_cwreduc<T>(base)
	{ this->set_name("colwise-dot-eval");  }

	LMAT_ENSURE_INLINE
	void operator() () const
	{
		colwise_dot(this->a, this->b, this->r);
		this->force_res();
	}
};


template<typename T>
struct bench_rowwise_sum_naiveloop : public bench_rwreduc<T>
{
	bench_rowwise_sum_naiveloop( const bench_rwreduc<T>& base )
	: bench_rwreduc<T>(base)
	{ this->set_name("rowwise-sum-naiveloop");  }

	LMAT_ENSURE_INLINE
	void operator() () const
	{
		index_t m = this->a.nrows();
		index_t n = this->a.ncolumns();
		const cref_matrix<T>& a_ = this->a;
		T *pr = this->r.ptr_data();

		for (index_t i = 0; i < m; ++i)
		{
			T s(0);
			for (index_t j = 0; j < n; ++j)
				s += a_(i, j);
			pr[i] = s;
		}
		this->force_res();
	}
};


template<typename T>
struct bench_rowwise_sum_rawloop : public bench_rwreduc<T>
{
	bench_rowwise_sum_rawloop( const bench_rwreduc<T>& base )
	: bench_rwreduc<T>(base)
	{ this->set_name("rowwise-sum-rawloop");  }

	LMAT_ENSURE_INLINE
	void operator() () const
	{
		index_t m = this->a.nrows();
		index_t n = this->a.ncolumns();
		const T *pa = this->a.ptr_data();
		T *pr = this->r.ptr_data();

		for (index_t i = 0; i < m; ++i)
			pr[i] = pa[i];

		for (index_t j = 1; j < n; ++j)
		{
			for (index_t i = 0; i < m; ++i)
			{
				pr[i] += pa[i];
			}
		}
		this->force_res();
	}
};


template<typename T>
struct bench_rowwise_sum_eval : public bench_rwreduc<T>
{
	bench_rowwise_sum_eval( const bench_rwreduc<T>& base )
	: bench_rwreduc<T>(base)
	{ this->set_name("rowwise-sum-eval");  }

	LMAT_ENSURE_INLINE
	void operator() () const
	{
		rowwise_sum(this->a, this->r);
		this->force_res();
	}
};


template<typename T>
struct bench_rowwise_dot_naiveloop : public bench_rwreduc<T>
{
	bench_rowwise_dot_naiveloop( const bench_rwreduc<T>& base )
	: bench_rwreduc<T>(base)
	{ this->set_name("rowwise-dot-naiveloop");  }

	LMAT_ENSURE_INLINE
	void operator() () const
	{
		index_t m = this->a.nrows();
		index_t n = this->a.ncolumns();
		const cref_matrix<T>& a_ = this->a;
		const cref_matrix<T>& b_ = this->b;
		T *pr = this->r.ptr_data();

		for (index_t i = 0; i < m; ++i)
		{
			T s(0);
			for (index_t j = 0; j < n; ++j)
				s += a_(i, j) * b_(i, j);
			pr[i] = s;
		}
		this->force_res();
	}
};


template<typename T>
struct bench_rowwise_dot_rawloop : public bench_rwreduc<T>
{
	bench_rowwise_dot_rawloop( const bench_rwreduc<T>& base )
	: bench_rwreduc<T>(base)
	{ this->set_name("rowwise-dot-rawloop");  }

	LMAT_ENSURE_INLINE
	void operator() () const
	{
		index_t m = this->a.nrows();
		index_t n = this->a.ncolumns();
		const T *pa = this->a.ptr_data();
		const T *pb = this->b.ptr_data();
		T *pr = this->r.ptr_data();

		for (index_t i = 0; i < m; ++i)
			pr[i] = pa[i] * pb[i];

		for (index_t j = 1; j < n; ++j)
		{
			for (index_t i = 0; i < m; ++i)
			{
				pr[i] += pa[i] * pb[i];
			}
		}
		this->force_res();
	}
};


template<typename T>
struct bench_rowwise_dot_eval : public bench_rwreduc<T>
{
	bench_rowwise_dot_eval( const bench_rwreduc<T>& base )
	: bench_rwreduc<T>(base)
	{ this->set_name("rowwise-dot-eval");  }

	LMAT_ENSURE_INLINE
	void operator() () const
	{
		rowwise_dot(this->a, this->b, this->r);
		this->force_res();
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
	dense_matrix<T> dst(max_size, 1, zero());
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
		size_t pbsiz = 10000000 / size_t(m * n);

		benchmark_option opt(pbsiz);

		std::cout << "size = " << m << " x " << n << "\n";
		std::cout << "=======================================\n";

		bench_fullreduc<T> full_base(m, n, pa, pb);
		bench_cwreduc<T> cw_base(m, n, pa, pb, pd);
		bench_rwreduc<T> rw_base(m, n, pa, pb, pd);

		run_benchmark(bench_full_sum_rawloop<T>(full_base), mon, opt);
		run_benchmark(bench_full_sum_eval<T>(full_base),    mon, opt);
		run_benchmark(bench_full_dot_rawloop<T>(full_base), mon, opt);
		run_benchmark(bench_full_dot_eval<T>(full_base),    mon, opt);

		run_benchmark(bench_colwise_sum_rawloop<T>(cw_base), mon, opt);
		run_benchmark(bench_colwise_sum_eval<T>(cw_base),    mon, opt);
		run_benchmark(bench_colwise_dot_rawloop<T>(cw_base), mon, opt);
		run_benchmark(bench_colwise_dot_eval<T>(cw_base),    mon, opt);

		run_benchmark(bench_rowwise_sum_naiveloop<T>(rw_base), mon, opt);
		run_benchmark(bench_rowwise_sum_rawloop<T>(rw_base),   mon, opt);
		run_benchmark(bench_rowwise_sum_eval<T>(rw_base),      mon, opt);
		run_benchmark(bench_rowwise_dot_naiveloop<T>(rw_base), mon, opt);
		run_benchmark(bench_rowwise_dot_rawloop<T>(rw_base),   mon, opt);
		run_benchmark(bench_rowwise_dot_eval<T>(rw_base),      mon, opt);

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





