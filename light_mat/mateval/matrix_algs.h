/**
 * @file matrix_algs.h
 *
 * @brief Generic algorithms on matrices
 *
 * @author Dahua Lin
 */

#ifndef LIGHTMAT_MATRIX_ALGS_H_
#define LIGHTMAT_MATRIX_ALGS_H_

#include <light_mat/matrix/matrix_classes.h>
#include "internal/matrix_algs_internal.h"

#include <functional>
#include <algorithm>


namespace lmat
{

	/********************************************
	 *
	 *  counting
	 *
	 ********************************************/

	template<class A, typename T>
	inline size_t count(const IEWiseMatrix<A, T>& a)
	{
		return internal::count_impl<supports_linear_macc<A>::value>::run(a.derived());
	}

	template<class A, typename T, class D, typename TD>
	inline void colwise_count(const IEWiseMatrix<A, T>& a, IRegularMatrix<D, TD>& dmat)
	{
		D& d = dmat.derived();
		const index_t m = a.nrows();
		const index_t n = a.ncolumns();

		LMAT_CHECK_DIMS( dmat.nelems() == n )
		auto rd = make_multicol_accessor(atags::scalar(), in_(a.derived()));

		for (index_t j = 0; j < n; ++j)
		{
			dmat[j] = internal::count_by_reader(m, rd.col(j));
		}
	}


	/********************************************
	 *
	 *  zip & unzip value and index
	 *
	 ********************************************/

	template<class A, typename T>
	dense_matrix<std::pair<T, index_t> > zip_with_idx(const IRegularMatrix<A, T>& a)
	{
		typedef std::pair<T, index_t> RT;
		dense_matrix<RT> r(a.nrows(), a.ncolumns());

		auto p = begin(a);
		const RT* r_ = r.ptr_data();
		const index_t n = r.nelems();

		for (index_t i = 0; i < n; ++i, ++p)
			r_[i] = std::make_pair(*p, i);
	}

	template<class A, typename T>
	dense_matrix<std::pair<T, std::pair<index_t, index_t> > >
	zip_with_ij(const IRegularMatrix<A, T>& a)
	{
		typedef std::pair<T, std::pair<index_t, index_t> > RT;

		const index_t m = a.nrows();
		const index_t n = a.ncolumns();
		dense_matrix<RT> r(m, n);

		for (index_t j = 0; j < n; ++j)
		{
			const RT* r_ = r.ptr_col(j);
			auto p = a.col_begin(j);

			for (index_t i = 0; i < m; ++i)
				r_[i] = std::make_pair(*p, std::make_pair(i, j));
		}
	}


	/********************************************
	 *
	 *  sorting
	 *
	 ********************************************/

	enum sort_order
	{
		SORD_ASC,
		SORD_DESC
	};

	// test is_sorted

	template<class A, typename T, typename Compare>
	inline bool
	is_sorted(const IRegularMatrix<A, T>& a, const Compare& comp)
	{
		return std::is_sorted(begin(a), end(a), comp);
	}

	template<class A, typename T>
	inline bool
	is_sorted(const IRegularMatrix<A, T>& a, sort_order ord=SORD_ASC)
	{
		return ord == SORD_ASC ?
				is_sorted(a, std::less_equal<T>()) :
				is_sorted(a, std::greater_equal<T>());
	}

	template<class A, typename T, class D, typename Compare>
	inline void
	colwise_is_sorted(const IRegularMatrix<A, T>& a, IRegularMatrix<D, bool>& r, const Compare& comp)
	{
		const index_t n = a.ncolumns();
		LMAT_CHECK_DIMS( r.nelems() == n )

		for (index_t j = 0; j < n; ++j)
			r[j] = std::is_sorted(a.col_begin(j), a.col_end(j), comp);
	}

	template<class A, typename T, class D>
	inline void
	colwise_is_sorted(const IRegularMatrix<A, T>& a, IRegularMatrix<D, bool>& r, sort_order ord=SORD_ASC)
	{
		if (ord == SORD_ASC)
			colwise_is_sorted(a, r, std::less_equal<T>());
		else
			colwise_is_sorted(a, r, std::greater_equal<T>());
	}


	// sorting algorithms

	struct std_sort
	{
		template<typename Iterator, typename Compare>
		LMAT_ENSURE_INLINE
		void sort(Iterator first, Iterator last, Compare comp)
		{
			std::sort(first, last, comp);
		}
	};

	struct stable_sort
	{
		template<typename Iterator, typename Compare>
		LMAT_ENSURE_INLINE
		void sort(Iterator first, Iterator last, Compare comp)
		{
			std::stable_sort(first, last, comp);
		}
	};


	struct partial_sort
	{
		const index_t k;

		partial_sort(index_t k_) : k(k_) { }

		template<typename Iterator, typename Compare>
		LMAT_ENSURE_INLINE
		void sort(Iterator first, Iterator last, Compare comp)
		{
			std::partial_sort(first, first + k, last);
		}
	};


	// Inplace sorting

	template<class A, typename T, typename Alg, typename Compare>
	inline typename meta::enable_if<meta::supports_random_access<A>,
	void>::type
	gsort(IRegularMatrix<A, T>& a, const Alg& alg, const Compare& comp)
	{
		alg.sort(begin(a), end(a), comp);
	}


	template<class A, typename T, typename Alg>
	inline typename meta::enable_if<meta::supports_random_access<A>,
	void>::type
	gsort(IRegularMatrix<A, T>& a, const Alg& alg, sort_order ord=SORD_ASC)
	{
		if (ord == SORD_ASC)
			gsort(a, alg, std::less<T>());
		else
			gsort(a, alg, std::greater<T>());
	}

	template<class A, typename T, typename Alg, typename Compare>
	inline void
	colwise_gsort(IRegularMatrix<A, T>& a, const Alg& alg, const Compare& comp)
	{
		const index_t n = a.ncolumns();
		for (index_t j = 0; j < n; ++j)
			alg.sort(a.col_begin(j), a.col_end(j), comp);
	}


	template<class A, typename T, typename Alg>
	inline void
	colwise_gsort(IRegularMatrix<A, T>& a, const Alg& alg, sort_order ord=SORD_ASC)
	{
		if (ord == SORD_ASC)
			colwise_gsort(a, alg, std::less<T>());
		else
			colwise_gsort(a, alg, std::greater<T>());
	}


	template<class A, typename T, typename Compare>
	inline void sort(IRegularMatrix<A, T>& a, const Compare& comp)
	{
		gsort(a, std_sort(), comp);
	}


	template<class A, typename T>
	inline void sort(IRegularMatrix<A, T>& a, sort_order ord=SORD_ASC)
	{
		gsort(a, std_sort(), ord);
	}

	template<class A, typename T, typename Compare>
	inline void colwise_sort(IRegularMatrix<A, T>& a, const Compare& comp)
	{
		colwise_gsort(a, std_sort(), comp);
	}

	template<class A, typename T>
	inline void colwise_gsort(IRegularMatrix<A, T>& a, sort_order ord=SORD_ASC)
	{
		colwise_gsort(a, std_sort(), ord);
	}




	template<class A, typename T, class I, typename IDX, typename Compare>
	inline typename meta::enable_if<
		meta::and_<meta::supports_random_access<A>, meta::supports_random_access<I> >,
	void>::type
	sort_idx(const IRegularMatrix<A, T>& a, IRegularMatrix<I, IDX>& idx, const Compare& comp)
	{
		idx.require_size(a.nrows(), a.ncolumns());
		internal::_fill_inds(idx.derived());

		auto a_ = begin(a);
		std::sort(begin(idx), end(idx),
				[](IDX i, IDX j) { return comp(a_[i], a_[j]); } );
	}

	template<class A, typename T, class I, typename IDX, typename Compare>
	inline typename meta::enable_if<
		meta::and_<meta::supports_random_access<A>, meta::supports_random_access<I> >,
	void>::type
	sort_idx(const IRegularMatrix<A, T>& a, IRegularMatrix<I, IDX>& idx, sort_order ord=SORD_ASC)
	{
		if (ord == SORD_ASC)
			sort_idx(a, idx, std::less<T>());
		else
			sort_idx(a, idx, std::greater<T>());
	}

	template<class A, typename T, class I, typename IDX, typename Compare>
	inline void
	colwise_sort_idx(const IRegularMatrix<A, T>& a, IRegularMatrix<I, IDX>& idx, const Compare& comp)
	{
		idx.require_size(a.nrows(), a.ncolumns());
		internal::_fill_subs_i(idx.derived());

		const index_t n = a.ncolumns();
		for (index_t j = 0; j < n; ++j)
		{
			auto a_ = a.col_begin(j);
			std::sort(idx.col_begin(j), idx.col_end(j),
					[](IDX u, IDX v) { return comp(a_[u], a_[v]); } );
		}
	}

	template<class A, typename T, class I, typename IDX, typename Compare>
	inline void
	colwise_sort_idx(const IRegularMatrix<A, T>& a, IRegularMatrix<I, IDX>& idx, sort_order ord=SORD_ASC)
	{
		if (ord == SORD_ASC)
			colwise_sort_idx(a, idx, std::less<T>());
		else
			colwise_sort_idx(a, idx, std::greater<T>());
	}






}

#endif /* MATRIX_ALGS_H_ */
