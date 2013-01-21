/**
 * @file matrix_sort.h
 *
 * Matrix sorting algorithms
 * 
 * @author Dahua Lin 
 */

#ifdef _MSC_VER
#pragma once
#endif

#ifndef LIGHTMAT_MATRIX_SORT_H_
#define LIGHTMAT_MATRIX_SORT_H_

#include <light_mat/matrix/matrix_classes.h>
#include <light_mat/matexpr/subs_expr.h>
#include <light_mat/matexpr/mat_zip.h>

#include <functional>
#include <algorithm>

namespace lmat
{

	/********************************************
	 *
	 *  definitions of common types
	 *
	 ********************************************/

	struct asc_ { };
	struct desc_ { };

	// forward declaration of expressions

	template<class Arg, class Alg, class Compare> class sort_expr;
	template<class Arg, class Alg, class Compare> class colwise_sort_expr;

	template<class Arg, class Alg, class Compare> class sort_idx_expr;
	template<class Arg, class Alg, class Compare> class colwise_sort_idx_expr;

	template<class A>
	struct sortx_value
	{
		typedef std::pair<typename matrix_traits<A>::value_type, index_t> type;
	};

	template<class Arg, class Alg, class Compare> class sortx_expr;
	template<class Arg, class Alg, class Compare> class colwise_sortx_expr;


	// sorting algorithms

	struct std_sort
	{
		template<typename Iterator, typename Compare>
		LMAT_ENSURE_INLINE
		void sort(Iterator first, Iterator last, Compare comp) const
		{
			std::sort(first, last, comp);
		}
	};

	struct stable_sort
	{
		template<typename Iterator, typename Compare>
		LMAT_ENSURE_INLINE
		void sort(Iterator first, Iterator last, Compare comp) const
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
		void sort(Iterator first, Iterator last, Compare comp) const
		{
			std::partial_sort(first, first + k, last, comp);
		}
	};

	typedef std_sort default_sort_alg;


	/********************************************
	 *
	 *  test is_sorted
	 *
	 ********************************************/

	template<class A, typename T, typename Compare>
	inline bool is_sorted(const IRegularMatrix<A, T>& a, const Compare& comp)
	{
		return std::is_sorted(begin(a), end(a), comp);
	}

	template<class A, typename T>
	inline bool is_sorted(const IRegularMatrix<A, T>& a, asc_)
	{
		return is_sorted(a, std::less<T>());
	}

	template<class A, typename T>
	inline bool is_sorted(const IRegularMatrix<A, T>& a, desc_)
	{
		return is_sorted(a, std::greater<T>());
	}

	template<class A, typename T>
	inline bool is_sorted(const IRegularMatrix<A, T>& a)
	{
		return is_sorted(a, asc_());
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
	inline void colwise_is_sorted(const IRegularMatrix<A, T>& a, IRegularMatrix<D, bool>& r, asc_)
	{
		colwise_is_sorted(a, r, std::less_equal<T>());
	}

	template<class A, typename T, class D>
	inline void colwise_is_sorted(const IRegularMatrix<A, T>& a, IRegularMatrix<D, bool>& r, desc_)
	{
		colwise_is_sorted(a, r, std::greater_equal<T>());
	}

	template<class A, typename T, class D>
	inline void colwise_is_sorted(const IRegularMatrix<A, T>& a, IRegularMatrix<D, bool>& r)
	{
		colwise_is_sorted(a, r, asc_());
	}


	/********************************************
	 *
	 *  inplace sorting (core)
	 *
	 ********************************************/

	template<class A, typename T, typename Alg, typename Compare>
	inline typename meta::enable_if_<meta::supports_random_access<A>,
	void>::type
	gsort(IRegularMatrix<A, T>& a, const Alg& alg, const Compare& comp)
	{
		alg.sort(begin(a), end(a), comp);
	}

	template<class A, typename T, typename Alg>
	void gsort(IRegularMatrix<A, T>& a, const Alg& alg, asc_)
	{
		gsort(a, alg, std::less<T>());
	}

	template<class A, typename T, typename Alg>
	void gsort(IRegularMatrix<A, T>& a, const Alg& alg, desc_)
	{
		gsort(a, alg, std::greater<T>());
	}

	template<class A, typename T, typename Spec>
	inline void sort(IRegularMatrix<A, T>& a, const Spec& s)
	{
		gsort(a, default_sort_alg(), s);
	}

	template<class A, typename T>
	inline void sort(IRegularMatrix<A, T>& a)
	{
		gsort(a, default_sort_alg(), asc_());
	}


	// colwise

	template<class A, typename T, typename Alg, typename Compare>
	inline void
	colwise_gsort(IRegularMatrix<A, T>& a, const Alg& alg, const Compare& comp)
	{
		const index_t n = a.ncolumns();
		for (index_t j = 0; j < n; ++j)
			alg.sort(a.col_begin(j), a.col_end(j), comp);
	}

	template<class A, typename T, typename Alg>
	inline void colwise_gsort(IRegularMatrix<A, T>& a, const Alg& alg, asc_)
	{
		colwise_gsort(a, alg, std::less<T>());
	}

	template<class A, typename T, typename Alg>
	inline void colwise_gsort(IRegularMatrix<A, T>& a, const Alg& alg, desc_)
	{
		colwise_gsort(a, alg, std::greater<T>());
	}

	template<class A, typename T, typename Spec>
	inline void colwise_sort(IRegularMatrix<A, T>& a, const Spec& s)
	{
		colwise_gsort(a, default_sort_alg(), s);
	}

	template<class A, typename T>
	inline void colwise_sort(IRegularMatrix<A, T>& a)
	{
		colwise_gsort(a, default_sort_alg(), asc_());
	}


	/********************************************
	 *
	 *  sort expression
	 *
	 ********************************************/

	template<class Arg, class Alg, class Compare>
	struct matrix_traits<sort_expr<Arg, Alg, Compare> >
	{
		static const int num_dimensions = 2;
		static const int ct_num_rows = meta::nrows<Arg>::value;
		static const int ct_num_cols = meta::ncols<Arg>::value;

		static const bool is_readonly = true;

		typedef matrix_shape<ct_num_rows, ct_num_cols> shape_type;
		typedef typename matrix_traits<Arg>::value_type value_type;
		typedef typename matrix_traits<Arg>::domain domain;
	};

	template<class Arg, class Alg, class Compare>
	struct matrix_traits<colwise_sort_expr<Arg, Alg, Compare> >
	{
		static const int num_dimensions = 2;
		static const int ct_num_rows = meta::nrows<Arg>::value;
		static const int ct_num_cols = meta::ncols<Arg>::value;

		static const bool is_readonly = true;

		typedef matrix_shape<ct_num_rows, ct_num_cols> shape_type;
		typedef typename matrix_traits<Arg>::value_type value_type;
		typedef typename matrix_traits<Arg>::domain domain;
	};


	template<class Arg, class Alg, class Compare>
	class sort_expr
	: public IMatrixXpr<sort_expr<Arg, Alg, Compare>, typename matrix_traits<Arg>::value_type>
	{
	public:
		sort_expr(const Arg& arg, const Alg& alg, const Compare& cmp)
		: m_arg(arg), m_alg(alg), m_cmp(cmp) { }

		LMAT_ENSURE_INLINE const Arg& arg() const
		{
			return m_arg;
		}

		LMAT_ENSURE_INLINE const Alg& algorithm() const
		{
			return m_alg;
		}

		LMAT_ENSURE_INLINE Compare comparer() const
		{
			return m_cmp;
		}

		LMAT_ENSURE_INLINE index_t nrows() const
		{
			return m_arg.nrows();
		}

		LMAT_ENSURE_INLINE index_t ncolumns() const
		{
			return m_arg.ncolumns();
		}

		LMAT_ENSURE_INLINE index_t nelems() const
		{
			return m_arg.nelems();
		}

		LMAT_ENSURE_INLINE
		matrix_shape<meta::nrows<Arg>::value, meta::ncols<Arg>::value>
		shape() const
		{
			return m_arg.shape();
		}

	private:
		const Arg& m_arg;
		Alg m_alg;
		Compare m_cmp;
	};

	template<class Arg, class Alg, class Compare>
	class colwise_sort_expr
	: public IMatrixXpr<colwise_sort_expr<Arg, Alg, Compare>, typename matrix_traits<Arg>::value_type>
	{
	public:
		colwise_sort_expr(const Arg& arg, const Alg& alg, const Compare& cmp)
		: m_arg(arg), m_alg(alg), m_cmp(cmp) { }

		LMAT_ENSURE_INLINE const Arg& arg() const
		{
			return m_arg;
		}

		LMAT_ENSURE_INLINE const Alg& algorithm() const
		{
			return m_alg;
		}

		LMAT_ENSURE_INLINE Compare comparer() const
		{
			return m_cmp;
		}

		LMAT_ENSURE_INLINE index_t nrows() const
		{
			return m_arg.nrows();
		}

		LMAT_ENSURE_INLINE index_t ncolumns() const
		{
			return m_arg.ncolumns();
		}

		LMAT_ENSURE_INLINE index_t nelems() const
		{
			return m_arg.nelems();
		}

		LMAT_ENSURE_INLINE
		matrix_shape<meta::nrows<Arg>::value, meta::ncols<Arg>::value>
		shape() const
		{
			return m_arg.shape();
		}

	private:
		const Arg& m_arg;
		Alg m_alg;
		Compare m_cmp;
	};

	template<class Arg, class Alg, class Compare, class DMat>
	LMAT_ENSURE_INLINE
	inline void evaluate(const sort_expr<Arg, Alg, Compare>& expr,
			IRegularMatrix<DMat, typename matrix_traits<Arg>::value_type>& dmat)
	{
		dmat.derived() = expr.arg();
		gsort(dmat.derived(), expr.algorithm(), expr.comparer());
	}

	template<class Arg, class Alg, class Compare, class DMat>
	LMAT_ENSURE_INLINE
	inline void evaluate(const colwise_sort_expr<Arg, Alg, Compare>& expr,
			IRegularMatrix<DMat, typename matrix_traits<Arg>::value_type>& dmat)
	{
		dmat.derived() = expr.arg();
		colwise_gsort(dmat.derived(), expr.algorithm(), expr.comparer());
	}


	// expression construction

	template<typename T, class Mat, class Alg, class Compare>
	LMAT_ENSURE_INLINE
	inline sort_expr<Mat, Alg, Compare> gsorted(const IMatrixXpr<Mat, T>& a, const Alg& alg, const Compare& cmp)
	{
		return sort_expr<Mat, Alg, Compare>(a.derived(), alg, cmp);
	}

	template<typename T, class Mat, class Alg>
	LMAT_ENSURE_INLINE
	inline sort_expr<Mat, Alg, std::less<T> > gsorted(const IMatrixXpr<Mat, T>& a, const Alg& alg, asc_)
	{
		return sort_expr<Mat, Alg, std::less<T> >(a.derived(), alg, std::less<T>());
	}

	template<typename T, class Mat, class Alg>
	LMAT_ENSURE_INLINE
	inline sort_expr<Mat, Alg, std::greater<T> > gsorted(const IMatrixXpr<Mat, T>& a, const Alg& alg, desc_)
	{
		return sort_expr<Mat, Alg, std::greater<T> >(a.derived(), alg, std::greater<T>());
	}

	template<typename T, class Mat, class Compare>
	LMAT_ENSURE_INLINE
	inline sort_expr<Mat, default_sort_alg, Compare> sorted(const IMatrixXpr<Mat, T>& a, const Compare& cmp)
	{
		return gsorted(a, default_sort_alg(), cmp);
	}

	template<typename T, class Mat>
	LMAT_ENSURE_INLINE
	inline sort_expr<Mat, default_sort_alg, std::less<T> > sorted(const IMatrixXpr<Mat, T>& a, asc_)
	{
		return gsorted(a, default_sort_alg(), asc_());
	}

	template<typename T, class Mat>
	LMAT_ENSURE_INLINE
	inline sort_expr<Mat, default_sort_alg, std::greater<T> > sorted(const IMatrixXpr<Mat, T>& a, desc_)
	{
		return gsorted(a, default_sort_alg(), desc_());
	}

	template<typename T, class Mat>
	LMAT_ENSURE_INLINE
	inline sort_expr<Mat, default_sort_alg, std::less<T> > sorted(const IMatrixXpr<Mat, T>& a)
	{
		return gsorted(a, default_sort_alg(), asc_());
	}


	template<typename T, class Mat, class Alg, class Compare>
	LMAT_ENSURE_INLINE
	inline colwise_sort_expr<Mat, Alg, Compare> colwise_gsorted(const IMatrixXpr<Mat, T>& a, const Alg& alg, const Compare& cmp)
	{
		return colwise_sort_expr<Mat, Alg, Compare>(a.derived(), alg, cmp);
	}

	template<typename T, class Mat, class Alg>
	LMAT_ENSURE_INLINE
	inline colwise_sort_expr<Mat, Alg, std::less<T> > colwise_gsorted(const IMatrixXpr<Mat, T>& a, const Alg& alg, asc_)
	{
		return colwise_sort_expr<Mat, Alg, std::less<T> >(a.derived(), alg, std::less<T>());
	}

	template<typename T, class Mat, class Alg>
	LMAT_ENSURE_INLINE
	inline colwise_sort_expr<Mat, Alg, std::greater<T> > colwise_gsorted(const IMatrixXpr<Mat, T>& a, const Alg& alg, desc_)
	{
		return colwise_sort_expr<Mat, Alg, std::greater<T> >(a.derived(), alg, std::greater<T>());
	}

	template<typename T, class Mat, class Compare>
	LMAT_ENSURE_INLINE
	inline colwise_sort_expr<Mat, default_sort_alg, Compare> colwise_sorted(const IMatrixXpr<Mat, T>& a, const Compare& cmp)
	{
		return colwise_gsorted(a, default_sort_alg(), cmp);
	}

	template<typename T, class Mat>
	LMAT_ENSURE_INLINE
	inline colwise_sort_expr<Mat, default_sort_alg, std::less<T> > colwise_sorted(const IMatrixXpr<Mat, T>& a, asc_)
	{
		return colwise_gsorted(a, default_sort_alg(), asc_());
	}

	template<typename T, class Mat>
	LMAT_ENSURE_INLINE
	inline colwise_sort_expr<Mat, default_sort_alg, std::greater<T> > colwise_sorted(const IMatrixXpr<Mat, T>& a, desc_)
	{
		return colwise_gsorted(a, default_sort_alg(), desc_());
	}

	template<typename T, class Mat>
	LMAT_ENSURE_INLINE
	inline colwise_sort_expr<Mat, default_sort_alg, std::less<T> > colwise_sorted(const IMatrixXpr<Mat, T>& a)
	{
		return colwise_gsorted(a, default_sort_alg(), asc_());
	}


	/********************************************
	 *
	 *  sort indices
	 *
	 ********************************************/

	template<class Arg, class Alg, class Compare>
	struct matrix_traits<sort_idx_expr<Arg, Alg, Compare> >
	{
		static const int num_dimensions = 2;
		static const int ct_num_rows = meta::nrows<Arg>::value;
		static const int ct_num_cols = meta::ncols<Arg>::value;

		static const bool is_readonly = true;

		typedef matrix_shape<ct_num_rows, ct_num_cols> shape_type;
		typedef index_t value_type;
		typedef typename matrix_traits<Arg>::domain domain;
	};

	template<class Arg, class Alg, class Compare>
	struct matrix_traits<colwise_sort_idx_expr<Arg, Alg, Compare> >
	{
		static const int num_dimensions = 2;
		static const int ct_num_rows = meta::nrows<Arg>::value;
		static const int ct_num_cols = meta::ncols<Arg>::value;

		static const bool is_readonly = true;

		typedef matrix_shape<ct_num_rows, ct_num_cols> shape_type;
		typedef index_t value_type;
		typedef typename matrix_traits<Arg>::domain domain;
	};


	template<class Arg, class Alg, class Compare>
	class sort_idx_expr
	: public IMatrixXpr<sort_idx_expr<Arg, Alg, Compare>, index_t>
	{
		static_assert(meta::supports_random_access<Arg>::value,
				"Arg for sort_idx_expr should support random access iteration.");

	public:
		sort_idx_expr(const Arg& arg, const Alg& alg, const Compare& cmp)
		: m_arg(arg), m_alg(alg), m_cmp(cmp) { }

		LMAT_ENSURE_INLINE const Arg& arg() const
		{
			return m_arg;
		}

		LMAT_ENSURE_INLINE const Alg& algorithm() const
		{
			return m_alg;
		}

		LMAT_ENSURE_INLINE Compare comparer() const
		{
			return m_cmp;
		}

		LMAT_ENSURE_INLINE index_t nrows() const
		{
			return m_arg.nrows();
		}

		LMAT_ENSURE_INLINE index_t ncolumns() const
		{
			return m_arg.ncolumns();
		}

		LMAT_ENSURE_INLINE index_t nelems() const
		{
			return m_arg.nelems();
		}

		LMAT_ENSURE_INLINE
		matrix_shape<meta::nrows<Arg>::value, meta::ncols<Arg>::value>
		shape() const
		{
			return m_arg.shape();
		}

	private:
		const Arg& m_arg;
		Alg m_alg;
		Compare m_cmp;
	};

	template<class Arg, class Alg, class Compare>
	class colwise_sort_idx_expr
	: public IMatrixXpr<colwise_sort_idx_expr<Arg, Alg, Compare>, index_t>
	{
	public:
		colwise_sort_idx_expr(const Arg& arg, const Alg& alg, const Compare& cmp)
		: m_arg(arg), m_alg(alg), m_cmp(cmp) { }

		LMAT_ENSURE_INLINE const Arg& arg() const
		{
			return m_arg;
		}

		LMAT_ENSURE_INLINE const Alg& algorithm() const
		{
			return m_alg;
		}

		LMAT_ENSURE_INLINE Compare comparer() const
		{
			return m_cmp;
		}

		LMAT_ENSURE_INLINE index_t nrows() const
		{
			return m_arg.nrows();
		}

		LMAT_ENSURE_INLINE index_t ncolumns() const
		{
			return m_arg.ncolumns();
		}

		LMAT_ENSURE_INLINE index_t nelems() const
		{
			return m_arg.nelems();
		}

		LMAT_ENSURE_INLINE
		matrix_shape<meta::nrows<Arg>::value, meta::ncols<Arg>::value>
		shape() const
		{
			return m_arg.shape();
		}

	private:
		const Arg& m_arg;
		Alg m_alg;
		Compare m_cmp;
	};

	template<class Arg, class Alg, class Compare, class DMat>
	LMAT_ENSURE_INLINE
	inline void evaluate(const sort_idx_expr<Arg, Alg, Compare>& expr, IRegularMatrix<DMat, index_t>& dmat)
	{
		auto a = begin(expr.arg());
		Compare cmp = expr.comparer();

		DMat& dm = dmat.derived();
		dm = inds(expr.shape());

		gsort(dm, expr.algorithm(),
				[&](const index_t& u, const index_t& v)
				{ return cmp(a[u], a[v]); }
		);
	}

	template<class Arg, class Alg, class Compare, class DMat>
	inline void evaluate(const colwise_sort_idx_expr<Arg, Alg, Compare>& expr, IRegularMatrix<DMat, index_t>& dmat)
	{
		Compare cmp = expr.comparer();
		const index_t n = expr.ncolumns();

		DMat& dm = dmat.derived();
		dm = subs_i(expr.shape());

		for (index_t j = 0; j < n; ++j)
		{
			auto a = expr.arg().col_begin(j);
			expr.algorithm().sort(dm.col_begin(j), dm.col_end(j),
					[&](const index_t& u, const index_t& v)
					{ return cmp(a[u], a[v]); }
			);
		}
	}


	// expression construction

	template<typename T, class Mat, class Alg, class Compare>
	LMAT_ENSURE_INLINE
	inline sort_idx_expr<Mat, Alg, Compare> gsorted_idx(const IRegularMatrix<Mat, T>& a, const Alg& alg, const Compare& cmp)
	{
		return sort_idx_expr<Mat, Alg, Compare>(a.derived(), alg, cmp);
	}

	template<typename T, class Mat, class Alg>
	LMAT_ENSURE_INLINE
	inline sort_idx_expr<Mat, Alg, std::less<T> > gsorted_idx(const IRegularMatrix<Mat, T>& a, const Alg& alg, asc_)
	{
		return sort_idx_expr<Mat, Alg, std::less<T> >(a.derived(), alg, std::less<T>());
	}

	template<typename T, class Mat, class Alg>
	LMAT_ENSURE_INLINE
	inline sort_idx_expr<Mat, Alg, std::greater<T> > gsorted_idx(const IRegularMatrix<Mat, T>& a, const Alg& alg, desc_)
	{
		return sort_idx_expr<Mat, Alg, std::greater<T> >(a.derived(), alg, std::greater<T>());
	}

	template<typename T, class Mat, class Compare>
	LMAT_ENSURE_INLINE
	inline sort_idx_expr<Mat, default_sort_alg, Compare> sorted_idx(const IRegularMatrix<Mat, T>& a, const Compare& cmp)
	{
		return gsorted_idx(a, default_sort_alg(), cmp);
	}

	template<typename T, class Mat>
	LMAT_ENSURE_INLINE
	inline sort_idx_expr<Mat, default_sort_alg, std::less<T> > sorted_idx(const IRegularMatrix<Mat, T>& a, asc_)
	{
		return gsorted_idx(a, default_sort_alg(), asc_());
	}

	template<typename T, class Mat>
	LMAT_ENSURE_INLINE
	inline sort_idx_expr<Mat, default_sort_alg, std::greater<T> > sorted_idx(const IRegularMatrix<Mat, T>& a, desc_)
	{
		return gsorted_idx(a, default_sort_alg(), desc_());
	}

	template<typename T, class Mat>
	LMAT_ENSURE_INLINE
	inline sort_idx_expr<Mat, default_sort_alg, std::less<T> > sorted_idx(const IRegularMatrix<Mat, T>& a)
	{
		return gsorted_idx(a, default_sort_alg(), asc_());
	}


	template<typename T, class Mat, class Alg, class Compare>
	LMAT_ENSURE_INLINE
	inline colwise_sort_idx_expr<Mat, Alg, Compare> colwise_gsorted_idx(const IRegularMatrix<Mat, T>& a, const Alg& alg, const Compare& cmp)
	{
		return colwise_sort_idx_expr<Mat, Alg, Compare>(a.derived(), alg, cmp);
	}

	template<typename T, class Mat, class Alg>
	LMAT_ENSURE_INLINE
	inline colwise_sort_idx_expr<Mat, Alg, std::less<T> > colwise_gsorted_idx(const IRegularMatrix<Mat, T>& a, const Alg& alg, asc_)
	{
		return colwise_sort_idx_expr<Mat, Alg, std::less<T> >(a.derived(), alg, std::less<T>());
	}

	template<typename T, class Mat, class Alg>
	LMAT_ENSURE_INLINE
	inline colwise_sort_idx_expr<Mat, Alg, std::greater<T> > colwise_gsorted_idx(const IRegularMatrix<Mat, T>& a, const Alg& alg, desc_)
	{
		return colwise_sort_idx_expr<Mat, Alg, std::greater<T> >(a.derived(), alg, std::greater<T>());
	}

	template<typename T, class Mat, class Compare>
	LMAT_ENSURE_INLINE
	inline colwise_sort_idx_expr<Mat, default_sort_alg, Compare> colwise_sorted_idx(const IRegularMatrix<Mat, T>& a, const Compare& cmp)
	{
		return colwise_gsorted_idx(a, default_sort_alg(), cmp);
	}

	template<typename T, class Mat>
	LMAT_ENSURE_INLINE
	inline colwise_sort_idx_expr<Mat, default_sort_alg, std::less<T> > colwise_sorted_idx(const IRegularMatrix<Mat, T>& a, asc_)
	{
		return colwise_gsorted_idx(a, default_sort_alg(), asc_());
	}

	template<typename T, class Mat>
	LMAT_ENSURE_INLINE
	inline colwise_sort_idx_expr<Mat, default_sort_alg, std::greater<T> > colwise_sorted_idx(const IRegularMatrix<Mat, T>& a, desc_)
	{
		return colwise_gsorted_idx(a, default_sort_alg(), desc_());
	}

	template<typename T, class Mat>
	LMAT_ENSURE_INLINE
	inline colwise_sort_idx_expr<Mat, default_sort_alg, std::less<T> > colwise_sorted_idx(const IRegularMatrix<Mat, T>& a)
	{
		return colwise_gsorted_idx(a, default_sort_alg(), asc_());
	}


	/********************************************
	 *
	 *  Extended (sort value + index)
	 *
	 ********************************************/

	template<class Arg, class Alg, class Compare>
	struct matrix_traits<sortx_expr<Arg, Alg, Compare> >
	{
		static const int num_dimensions = 2;
		static const int ct_num_rows = meta::nrows<Arg>::value;
		static const int ct_num_cols = meta::ncols<Arg>::value;

		static const bool is_readonly = true;

		typedef matrix_shape<ct_num_rows, ct_num_cols> shape_type;
		typedef typename sortx_value<Arg>::type value_type;
		typedef typename matrix_traits<Arg>::domain domain;
	};

	template<class Arg, class Alg, class Compare>
	struct matrix_traits<colwise_sortx_expr<Arg, Alg, Compare> >
	{
		static const int num_dimensions = 2;
		static const int ct_num_rows = meta::nrows<Arg>::value;
		static const int ct_num_cols = meta::ncols<Arg>::value;

		static const bool is_readonly = true;

		typedef matrix_shape<ct_num_rows, ct_num_cols> shape_type;
		typedef typename sortx_value<Arg>::type value_type;
		typedef typename matrix_traits<Arg>::domain domain;
	};

	template<class Arg, class Alg, class Compare>
	class sortx_expr
	: public IMatrixXpr<sortx_expr<Arg, Alg, Compare>, typename sortx_value<Arg>::type>
	{
	public:
		sortx_expr(const Arg& arg, const Alg& alg, const Compare& cmp)
		: m_arg(arg), m_alg(alg), m_cmp(cmp) { }

		LMAT_ENSURE_INLINE const Arg& arg() const
		{
			return m_arg;
		}

		LMAT_ENSURE_INLINE const Alg& algorithm() const
		{
			return m_alg;
		}

		LMAT_ENSURE_INLINE Compare comparer() const
		{
			return m_cmp;
		}

		LMAT_ENSURE_INLINE index_t nrows() const
		{
			return m_arg.nrows();
		}

		LMAT_ENSURE_INLINE index_t ncolumns() const
		{
			return m_arg.ncolumns();
		}

		LMAT_ENSURE_INLINE index_t nelems() const
		{
			return m_arg.nelems();
		}

		LMAT_ENSURE_INLINE
		matrix_shape<meta::nrows<Arg>::value, meta::ncols<Arg>::value>
		shape() const
		{
			return m_arg.shape();
		}

	private:
		const Arg& m_arg;
		Alg m_alg;
		Compare m_cmp;
	};

	template<class Arg, class Alg, class Compare>
	class colwise_sortx_expr
	: public IMatrixXpr<colwise_sortx_expr<Arg, Alg, Compare>, typename sortx_value<Arg>::type>
	{
	public:
		colwise_sortx_expr(const Arg& arg, const Alg& alg, const Compare& cmp)
		: m_arg(arg), m_alg(alg), m_cmp(cmp) { }

		LMAT_ENSURE_INLINE const Arg& arg() const
		{
			return m_arg;
		}

		LMAT_ENSURE_INLINE const Alg& algorithm() const
		{
			return m_alg;
		}

		LMAT_ENSURE_INLINE Compare comparer() const
		{
			return m_cmp;
		}

		LMAT_ENSURE_INLINE index_t nrows() const
		{
			return m_arg.nrows();
		}

		LMAT_ENSURE_INLINE index_t ncolumns() const
		{
			return m_arg.ncolumns();
		}

		LMAT_ENSURE_INLINE index_t nelems() const
		{
			return m_arg.nelems();
		}

		LMAT_ENSURE_INLINE
		matrix_shape<meta::nrows<Arg>::value, meta::ncols<Arg>::value>
		shape() const
		{
			return m_arg.shape();
		}

	private:
		const Arg& m_arg;
		Alg m_alg;
		Compare m_cmp;
	};


	template<class Arg, class Alg, class Compare, class DMat>
	LMAT_ENSURE_INLINE
	inline void evaluate(const sortx_expr<Arg, Alg, Compare>& expr,
			IRegularMatrix<DMat, typename sortx_value<Arg>::type>& dmat)
	{
		typedef typename sortx_value<Arg>::type T;

		dmat.derived() = zip_pair(expr.arg(), inds(expr.shape()));
		Compare cmp = expr.comparer();
		gsort(dmat.derived(), expr.algorithm(),
				[&](const T& u, const T& v) { return cmp(u.first, v.first); } );
	}

	template<class Arg, class Alg, class Compare, class DMat>
	inline void evaluate(const colwise_sortx_expr<Arg, Alg, Compare>& expr,
			IRegularMatrix<DMat, typename sortx_value<Arg>::type>& dmat)
	{
		typedef typename sortx_value<Arg>::type T;

		dmat.derived() = zip_pair(expr.arg(), subs_i(expr.shape()));
		Compare cmp = expr.comparer();
		colwise_gsort(dmat.derived(), expr.algorithm(),
				[&](const T& u, const T& v) { return cmp(u.first, v.first); } );
	}


	// expression construction

	template<typename T, class Mat, class Alg, class Compare>
	LMAT_ENSURE_INLINE
	inline sortx_expr<Mat, Alg, Compare> gsorted_ex(const IEWiseMatrix<Mat, T>& a, const Alg& alg, const Compare& cmp)
	{
		return sortx_expr<Mat, Alg, Compare>(a.derived(), alg, cmp);
	}

	template<typename T, class Mat, class Alg>
	LMAT_ENSURE_INLINE
	inline sortx_expr<Mat, Alg, std::less<T> > gsorted_ex(const IEWiseMatrix<Mat, T>& a, const Alg& alg, asc_)
	{
		return sortx_expr<Mat, Alg, std::less<T> >(a.derived(), alg, std::less<T>());
	}

	template<typename T, class Mat, class Alg>
	LMAT_ENSURE_INLINE
	inline sortx_expr<Mat, Alg, std::greater<T> > gsorted_ex(const IEWiseMatrix<Mat, T>& a, const Alg& alg, desc_)
	{
		return sortx_expr<Mat, Alg, std::greater<T> >(a.derived(), alg, std::greater<T>());
	}

	template<typename T, class Mat, class Compare>
	LMAT_ENSURE_INLINE
	inline sortx_expr<Mat, default_sort_alg, Compare> sorted_ex(const IEWiseMatrix<Mat, T>& a, const Compare& cmp)
	{
		return gsorted_ex(a, default_sort_alg(), cmp);
	}

	template<typename T, class Mat>
	LMAT_ENSURE_INLINE
	inline sortx_expr<Mat, default_sort_alg, std::less<T> > sorted_ex(const IEWiseMatrix<Mat, T>& a, asc_)
	{
		return gsorted_ex(a, default_sort_alg(), asc_());
	}

	template<typename T, class Mat>
	LMAT_ENSURE_INLINE
	inline sortx_expr<Mat, default_sort_alg, std::greater<T> > sorted_ex(const IEWiseMatrix<Mat, T>& a, desc_)
	{
		return gsorted_ex(a, default_sort_alg(), desc_());
	}

	template<typename T, class Mat>
	LMAT_ENSURE_INLINE
	inline sortx_expr<Mat, default_sort_alg, std::less<T> > sorted_ex(const IEWiseMatrix<Mat, T>& a)
	{
		return gsorted_ex(a, default_sort_alg(), asc_());
	}


	template<typename T, class Mat, class Alg, class Compare>
	LMAT_ENSURE_INLINE
	inline colwise_sortx_expr<Mat, Alg, Compare> colwise_gsorted_ex(const IEWiseMatrix<Mat, T>& a, const Alg& alg, const Compare& cmp)
	{
		return colwise_sortx_expr<Mat, Alg, Compare>(a.derived(), alg, cmp);
	}

	template<typename T, class Mat, class Alg>
	LMAT_ENSURE_INLINE
	inline colwise_sortx_expr<Mat, Alg, std::less<T> > colwise_gsorted_ex(const IEWiseMatrix<Mat, T>& a, const Alg& alg, asc_)
	{
		return colwise_sortx_expr<Mat, Alg, std::less<T> >(a.derived(), alg, std::less<T>());
	}

	template<typename T, class Mat, class Alg>
	LMAT_ENSURE_INLINE
	inline colwise_sortx_expr<Mat, Alg, std::greater<T> > colwise_gsorted_ex(const IEWiseMatrix<Mat, T>& a, const Alg& alg, desc_)
	{
		return colwise_sortx_expr<Mat, Alg, std::greater<T> >(a.derived(), alg, std::greater<T>());
	}

	template<typename T, class Mat, class Compare>
	LMAT_ENSURE_INLINE
	inline colwise_sortx_expr<Mat, default_sort_alg, Compare> colwise_sorted_ex(const IEWiseMatrix<Mat, T>& a, const Compare& cmp)
	{
		return colwise_gsorted_ex(a, default_sort_alg(), cmp);
	}

	template<typename T, class Mat>
	LMAT_ENSURE_INLINE
	inline colwise_sortx_expr<Mat, default_sort_alg, std::less<T> > colwise_sorted_ex(const IEWiseMatrix<Mat, T>& a, asc_)
	{
		return colwise_gsorted_ex(a, default_sort_alg(), asc_());
	}

	template<typename T, class Mat>
	LMAT_ENSURE_INLINE
	inline colwise_sortx_expr<Mat, default_sort_alg, std::greater<T> > colwise_sorted_ex(const IEWiseMatrix<Mat, T>& a, desc_)
	{
		return colwise_gsorted_ex(a, default_sort_alg(), desc_());
	}

	template<typename T, class Mat>
	LMAT_ENSURE_INLINE
	inline colwise_sortx_expr<Mat, default_sort_alg, std::less<T> > colwise_sorted_ex(const IEWiseMatrix<Mat, T>& a)
	{
		return colwise_gsorted_ex(a, default_sort_alg(), asc_());
	}

}

#endif 
