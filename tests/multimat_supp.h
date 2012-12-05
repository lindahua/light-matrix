/**
 * @file multimat_supp.h
 *
 * @brief Support the testing of multiple matrix type
 *
 * @author Dahua Lin
 */

#ifndef MULTIMAT_SUPP_H_
#define MULTIMAT_SUPP_H_

#include <light_mat/matrix/matrix_classes.h>
#include <cstdlib>

using namespace lmat;

#ifndef DEFAULT_M_VALUE
#define DEFAULT_M_VALUE 7
#endif

#ifndef DEFAULT_N_VALUE
#define DEFAULT_N_VALUE 8
#endif

const index_t DM = DEFAULT_M_VALUE;
const index_t DN = DEFAULT_N_VALUE;
const index_t rs = 2;
const index_t cs = (DM > DN ? DM : DN) * rs + 1;
const index_t LDim = cs;


struct fill_lin { };
struct fill_rand { };

template<typename T>
void do_fill_lin(T *p, const index_t n, T a, T b)
{
	for (index_t i = 0; i < n; ++i)
	{
		p[i] = a * T(i) + b;
	}
}

template<typename T>
void do_fill_lin(T *p, const index_t n)
{
	for (index_t i = 0; i < n; ++i)
	{
		p[i] = T(i + 1);
	}
}

template<typename T>
void do_fill_rand(T *p, const index_t n)
{
	for (index_t i = 0; i < n; ++i)
	{
		p[i] = (T)std::rand() / RAND_MAX;
	}
}

template<typename T>
void do_fill_rand(T *p, const index_t n, T lb, T ub)
{
	for (index_t i = 0; i < n; ++i)
	{
		T u = (T)std::rand() / RAND_MAX;
		p[i] = lb + u * (ub - lb);
	}
}


template<typename VT>
class mat_host_base
{
public:
	mat_host_base(index_t m, index_t n, index_t max_siz)
	: m_nrows(m), m_ncols(n), m_blk(max_siz, zero()) { }

	const VT *pdata() const { return m_blk.ptr_data(); }

	VT *pdata() { return m_blk.ptr_data(); }

	index_t nrows() const { return m_nrows; }

	index_t ncolumns() const { return m_ncols; }

	index_t capacity() const { return m_blk.nelems(); }

public:
	void fill_lin()
	{
		do_fill_lin(pdata(), capacity());
	}

	void fill_lin(VT a, VT b)
	{
		do_fill_lin(pdata(), capacity(), a, b);
	}

	void fill_rand()
	{
		do_fill_rand(pdata(), capacity());
	}

	void fill_rand(VT lb, VT ub)
	{
		do_fill_rand(pdata(), capacity(), lb, ub);
	}
private:
	index_t m_nrows;
	index_t m_ncols;
	dblock<VT> m_blk;
};


struct cont {};
struct bloc {};
struct grid {};

template<typename Tag, typename VT, int M, int N>
class mat_host;

template<typename VT, int M, int N>
class mat_host<cont, VT, M, N> : public mat_host_base<VT>
{
public:
	typedef cref_matrix<double, M, N> cmat_t;
	typedef ref_matrix<double, M, N> mat_t;

	mat_host(index_t m, index_t n)
	: mat_host_base<VT>(m, n, m * n)
	{ }

	cmat_t get_cmat() const
	{
		return cmat_t(this->pdata(), this->nrows(), this->ncolumns());
	}

	mat_t get_mat()
	{
		return mat_t(this->pdata(), this->nrows(), this->ncolumns());
	}

private:
	dblock<VT> m_blk;
};

template<typename VT, int M, int N>
class mat_host<bloc, VT, M, N> : public mat_host_base<VT>
{
public:
	typedef cref_block<double, M, N> cmat_t;
	typedef ref_block<double, M, N> mat_t;

	mat_host(index_t m, index_t n)
	: mat_host_base<VT>(m, n, LDim * n)
	{ }

	cmat_t get_cmat() const
	{
		return cmat_t(this->pdata(), this->nrows(), this->ncolumns(), cs);
	}

	mat_t get_mat()
	{
		return mat_t(this->pdata(), this->nrows(), this->ncolumns(), cs);
	}

private:
	dblock<VT> m_blk;
};


template<typename VT, int M, int N>
class mat_host<grid, VT, M, N> : public mat_host_base<VT>
{
public:
	typedef cref_grid<double, M, N> cmat_t;
	typedef ref_grid<double, M, N> mat_t;

	mat_host(index_t m, index_t n)
	: mat_host_base<VT>(m, n, LDim * n)
	{ }

	cmat_t get_cmat() const
	{
		return cmat_t(this->pdata(), this->nrows(), this->ncolumns(), rs, cs);
	}

	mat_t get_mat()
	{
		return mat_t(this->pdata(), this->nrows(), this->ncolumns(), rs, cs);
	}

private:
	dblock<VT> m_blk;
};


#endif /* MULTIMAT_SUPP_H_ */
