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

template<template<typename T, int R, int C> class ClassT, int M, int N>
struct mat_maker;

template<int M, int N>
struct mat_maker<ref_matrix, M, N>
{
	typedef cref_matrix<double, M, N> cmat_t;
	typedef ref_matrix<double, M, N> mat_t;

	static index_t max_size(index_t m, index_t n)
	{
		return m * n;
	}

	static cmat_t get_cmat(const double *p, index_t m, index_t n)
	{
		return cmat_t(p, m, n);
	}

	static mat_t get_mat(double *p, index_t m, index_t n)
	{
		return mat_t(p, m, n);
	}
};

template<int M, int N>
struct mat_maker<ref_block, M, N>
{
	typedef cref_block<double, M, N> cmat_t;
	typedef ref_block<double, M, N> mat_t;

	static index_t max_size(index_t m, index_t n)
	{
		return LDim * n;
	}

	static cmat_t get_cmat(const double *p, index_t m, index_t n)
	{
		return cmat_t(p, m, n, cs);
	}

	static mat_t get_mat(double *p, index_t m, index_t n)
	{
		return mat_t(p, m, n, cs);
	}
};

template<int M, int N>
struct mat_maker<ref_grid, M, N>
{
	typedef cref_grid<double, M, N> cmat_t;
	typedef ref_grid<double, M, N> mat_t;

	static index_t max_size(index_t m, index_t n)
	{
		return LDim * n;
	}

	static cmat_t get_cmat(const double *p, index_t m, index_t n)
	{
		return cmat_t(p, m, n, rs, cs);
	}

	static mat_t get_mat(double *p, index_t m, index_t n)
	{
		return mat_t(p, m, n, rs, cs);
	}
};


#endif /* MULTIMAT_SUPP_H_ */
