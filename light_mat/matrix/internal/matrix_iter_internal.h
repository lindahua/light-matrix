/**
 * @file matrix_iter_internal.h
 *
 * @brief Internal implementation of matrix iteration
 *
 * @author Dahua Lin
 */

#ifdef _MSC_VER
#pragma once
#endif

#ifndef LIGHTMAT_MATRIX_ITER_INTERNAL_H_
#define LIGHTMAT_MATRIX_ITER_INTERNAL_H_

#include <light_mat/matrix/matrix_meta.h>

namespace lmat { namespace internal {

	// tags

	struct cont_iter_tag { };
	struct step_col_iter_tag { };
	struct step_row_iter_tag { };
	struct pcc_iter_tag { };
	struct gen_iter_tag { };

	struct cont_coliter_tag { };
	struct step_coliter_tag { };

	template<class Mat>
	struct iter_tag
	: public meta::select_<
	  meta::is_contiguous<Mat>, cont_iter_tag,
	  meta::is_col<Mat>, step_col_iter_tag,
	  meta::is_row<Mat>, step_row_iter_tag,
	  meta::is_percol_contiguous<Mat>, pcc_iter_tag,
	  meta::otherwise_, gen_iter_tag> { };

	template<class Mat>
	struct coliter_tag
	: public meta::if_<
	  meta::is_percol_contiguous<Mat>,
	  cont_coliter_tag,
	  step_coliter_tag> { };


	// forward

	template<class Mat, typename Tag> struct matrix_iter_helper;
	template<class Mat, typename Tag> struct matrix_coliter_helper;


	/********************************************
	 *
	 *  iteration
	 *
	 ********************************************/

	template<class Mat>
	struct matrix_iter_helper<Mat, cont_iter_tag>
	{
		typedef typename matrix_traits<Mat>::value_type T;
		typedef const T* const_iterator;
		typedef typename meta::if_<meta::is_readonly<Mat>,
				const T*, T*>::type iterator;

		static const bool supp_random_access = true;

		LMAT_ENSURE_INLINE
		static const_iterator begin(const Mat& a)
		{
			return a.ptr_data();
		}

		LMAT_ENSURE_INLINE
		static const_iterator end(const Mat& a)
		{
			return a.ptr_data() + a.nelems();
		}

		LMAT_ENSURE_INLINE
		static iterator begin(Mat& a)
		{
			return a.ptr_data();
		}

		LMAT_ENSURE_INLINE
		static iterator end(Mat& a)
		{
			return a.ptr_data() + a.nelems();
		}
	};


	template<class Mat>
	struct matrix_iter_helper<Mat, step_col_iter_tag>
	{
		typedef typename matrix_traits<Mat>::value_type T;
		typedef step_ptr_t<const T> const_iterator;
		typedef typename meta::if_<meta::is_readonly<Mat>,
				step_ptr_t<const T>, step_ptr_t<T> >::type iterator;

		static const bool supp_random_access = true;

		LMAT_ENSURE_INLINE
		static const_iterator begin(const Mat& a)
		{
			return const_iterator(a.ptr_data(), a.row_stride());
		}

		LMAT_ENSURE_INLINE
		static const_iterator end(const Mat& a)
		{
			return const_iterator(
					a.ptr_data() + a.nrows() * a.row_stride(),
					a.row_stride());
		}

		LMAT_ENSURE_INLINE
		static iterator begin(Mat& a)
		{
			return iterator(a.ptr_data(), a.row_stride());
		}

		LMAT_ENSURE_INLINE
		static iterator end(Mat& a)
		{
			return iterator(
					a.ptr_data() + a.nrows() * a.row_stride(),
					a.row_stride());
		}
	};


	template<class Mat>
	struct matrix_iter_helper<Mat, step_row_iter_tag>
	{
		typedef typename matrix_traits<Mat>::value_type T;
		typedef step_ptr_t<const T> const_iterator;
		typedef typename meta::if_<meta::is_readonly<Mat>,
				step_ptr_t<const T>, step_ptr_t<T> >::type iterator;

		static const bool supp_random_access = true;

		LMAT_ENSURE_INLINE
		static const_iterator begin(const Mat& a)
		{
			return const_iterator(a.ptr_data(), a.col_stride());
		}

		LMAT_ENSURE_INLINE
		static const_iterator end(const Mat& a)
		{
			return const_iterator(
					a.ptr_data() + a.ncolumns() * a.col_stride(),
					a.col_stride());
		}

		LMAT_ENSURE_INLINE
		static iterator begin(Mat& a)
		{
			return iterator(a.ptr_data(), a.col_stride());
		}

		LMAT_ENSURE_INLINE
		static iterator end(Mat& a)
		{
			return iterator(
					a.ptr_data() + a.ncolumns() * a.col_stride(),
					a.col_stride());
		}
	};


	template<typename T>
	class pcc_iterator_t
	{
	public:
		typedef T value_type;
		typedef ptrdiff_t difference_type;
		typedef T* pointer;
		typedef T& reference;
		typedef std::forward_iterator_tag iterator_category;

		LMAT_ENSURE_INLINE
		pcc_iterator_t()
		: m_pcol(nullptr), m_i(0), m_m(0), m_cs(0)
		{ }

		LMAT_ENSURE_INLINE
		pcc_iterator_t(T *p, index_t m, index_t cs)
		: m_pcol(p), m_i(0), m_m(m), m_cs(cs)
		{ }

		LMAT_ENSURE_INLINE
		T& operator * () const
		{
			return m_pcol[m_i];
		}

		LMAT_ENSURE_INLINE
		T* operator -> () const
		{
			return m_pcol + m_i;
		}

		LMAT_ENSURE_INLINE
		pcc_iterator_t& operator++ ()
		{
			next();
			return *this;
		}

		LMAT_ENSURE_INLINE
		pcc_iterator_t operator++ (int)
		{
			pcc_iterator_t old(*this);
			next();
			return old;
		}

		LMAT_ENSURE_INLINE
		bool operator == (const pcc_iterator_t& r) const
		{
			return m_pcol == r.m_pcol && m_i == r.m_i;
		}

		LMAT_ENSURE_INLINE
		bool operator != (const pcc_iterator_t& r) const
		{
			return !(operator == (r));
		}

	private:
		LMAT_ENSURE_INLINE
		void next()
		{
			if (++ m_i == m_m)
			{
				m_pcol += m_cs;
				m_i = 0;
			}
		}

		T* m_pcol;
		index_t m_i;
		index_t m_m;
		index_t m_cs;
	};


	template<class Mat>
	struct matrix_iter_helper<Mat, pcc_iter_tag>
	{
		typedef typename matrix_traits<Mat>::value_type T;
		typedef pcc_iterator_t<const T> const_iterator;
		typedef typename meta::if_<meta::is_readonly<Mat>,
				pcc_iterator_t<const T>, pcc_iterator_t<T> >::type iterator;

		static const bool supp_random_access = false;

		LMAT_ENSURE_INLINE
		static const_iterator begin(const Mat& a)
		{
			return const_iterator(a.ptr_data(), a.nrows(), a.col_stride());
		}

		LMAT_ENSURE_INLINE
		static const_iterator end(const Mat& a)
		{
			return const_iterator(a.ptr_data() + a.col_stride() * a.ncolumns(),
					a.nrows(), a.col_stride());
		}

		LMAT_ENSURE_INLINE
		static iterator begin(Mat& a)
		{
			return iterator(a.ptr_data(), a.nrows(), a.col_stride());
		}

		LMAT_ENSURE_INLINE
		static iterator end(Mat& a)
		{
			return iterator(a.ptr_data() + a.col_stride() * a.ncolumns(),
					a.nrows(), a.col_stride());
		}
	};


	template<typename T>
	class gen_iterator_t
	{
	public:
		typedef T value_type;
		typedef ptrdiff_t difference_type;
		typedef T* pointer;
		typedef T& reference;
		typedef std::forward_iterator_tag iterator_category;

		LMAT_ENSURE_INLINE
		gen_iterator_t()
		: m_pcol(nullptr), m_i(0), m_m(0), m_rs(0), m_cs(0)
		{ }

		LMAT_ENSURE_INLINE
		gen_iterator_t(T *p, index_t m, index_t rs, index_t cs)
		: m_pcol(p), m_i(0), m_m(m), m_rs(rs), m_cs(cs)
		{ }

		LMAT_ENSURE_INLINE
		T& operator * () const
		{
			return m_pcol[m_i * m_rs];
		}

		LMAT_ENSURE_INLINE
		T* operator -> () const
		{
			return m_pcol + m_i * m_rs;
		}

		LMAT_ENSURE_INLINE
		gen_iterator_t& operator++ ()
		{
			next();
			return *this;
		}

		LMAT_ENSURE_INLINE
		gen_iterator_t operator++ (int)
		{
			gen_iterator_t old(*this);
			next();
			return old;
		}

		LMAT_ENSURE_INLINE
		bool operator == (const gen_iterator_t& r) const
		{
			return m_pcol == r.m_pcol && m_i == r.m_i;
		}

		LMAT_ENSURE_INLINE
		bool operator != (const gen_iterator_t& r) const
		{
			return !(operator == (r));
		}

	private:
		LMAT_ENSURE_INLINE
		void next()
		{
			if (++ m_i == m_m)
			{
				m_pcol += m_cs;
				m_i = 0;
			}
		}

		T* m_pcol;
		index_t m_i;
		index_t m_m;
		index_t m_rs;
		index_t m_cs;
	};


	template<class Mat>
	struct matrix_iter_helper<Mat, gen_iter_tag>
	{
		typedef typename matrix_traits<Mat>::value_type T;
		typedef gen_iterator_t<const T> const_iterator;
		typedef typename meta::if_<meta::is_readonly<Mat>,
				gen_iterator_t<const T>, gen_iterator_t<T> >::type iterator;

		static const bool supp_random_access = false;

		LMAT_ENSURE_INLINE
		static const_iterator begin(const Mat& a)
		{
			return const_iterator(a.ptr_data(), a.nrows(),
					a.row_stride(), a.col_stride());
		}

		LMAT_ENSURE_INLINE
		static const_iterator end(const Mat& a)
		{
			return const_iterator(a.ptr_data() + a.col_stride() * a.ncolumns(),
					a.nrows(), a.row_stride(), a.col_stride());
		}

		LMAT_ENSURE_INLINE
		static iterator begin(Mat& a)
		{
			return iterator(a.ptr_data(), a.nrows(), a.row_stride(), a.col_stride());
		}

		LMAT_ENSURE_INLINE
		static iterator end(Mat& a)
		{
			return iterator(a.ptr_data() + a.col_stride() * a.ncolumns(),
					a.nrows(), a.row_stride(), a.col_stride());
		}
	};




	/********************************************
	 *
	 *  per-column iteration
	 *
	 ********************************************/

	template<class Mat>
	struct matrix_coliter_helper<Mat, cont_coliter_tag>
	{
		typedef typename matrix_traits<Mat>::value_type T;
		typedef const T* const_iterator;
		typedef typename meta::if_<meta::is_readonly<Mat>,
				const T*, T*>::type iterator;

		LMAT_ENSURE_INLINE
		static const_iterator begin(const Mat& a, index_t j)
		{
			return a.ptr_col(j);
		}

		LMAT_ENSURE_INLINE
		static const_iterator end(const Mat& a, index_t j)
		{
			return a.ptr_col(j) + a.nrows();
		}

		LMAT_ENSURE_INLINE
		static iterator begin(Mat& a, index_t j)
		{
			return a.ptr_col(j);
		}

		LMAT_ENSURE_INLINE
		static iterator end(Mat& a, index_t j)
		{
			return a.ptr_col(j) + a.nrows();
		}
	};


	template<class Mat>
	struct matrix_coliter_helper<Mat, step_coliter_tag>
	{
		typedef typename matrix_traits<Mat>::value_type T;
		typedef step_ptr_t<const T> const_iterator;
		typedef typename meta::if_<meta::is_readonly<Mat>,
				step_ptr_t<const T>, step_ptr_t<T> >::type iterator;

		LMAT_ENSURE_INLINE
		static const_iterator begin(const Mat& a, index_t j)
		{
			return const_iterator(a.ptr_col(j), a.row_stride());
		}

		LMAT_ENSURE_INLINE
		static const_iterator end(const Mat& a, index_t j)
		{
			return const_iterator(a.ptr_col(j) + a.row_stride() * a.nrows(), a.row_stride());
		}

		LMAT_ENSURE_INLINE
		static iterator begin(Mat& a, index_t j)
		{
			return iterator(a.ptr_col(j), a.row_stride());
		}

		LMAT_ENSURE_INLINE
		static iterator end(Mat& a, index_t j)
		{
			return iterator(a.ptr_col(j) + a.row_stride() * a.nrows(), a.row_stride());
		}
	};


} }

#endif /* MATRIX_ITER_INTERNAL_H_ */
