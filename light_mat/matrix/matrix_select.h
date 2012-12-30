/**
 * @file matrix_select.h
 *
 * @brief Selection of subset of elements
 *
 * @author Dahua Lin
 */

#ifndef LIGHTMAT_MATRIX_SELECT_H_
#define LIGHTMAT_MATRIX_SELECT_H_

#include <light_mat/matrix/matrix_properties.h>

namespace lmat
{
	// forward

	template<class Mat, class L> class selectl_expr;
	template<class Mat, class I, class J> class selectl2_expr;
	template<class Mat, class I, class J> class select_expr;

	template<class Mat, class L>
	struct matrix_traits<selectl_expr<Mat, L> >
	{
		static const int num_dimensions = 2;
		static const int ct_num_rows = meta::nrows<L>::value;
		static const int ct_num_cols = meta::ncols<L>::value;

		static const bool is_readonly = true;

		typedef matrix_shape<ct_num_rows, ct_num_cols> shape_type;
		typedef typename matrix_traits<Mat>::value_type value_type;
		typedef typename meta::common_domain<Mat, L>::type domain;
	};

	template<class Mat, class I, class J>
	struct matrix_traits<selectl2_expr<Mat, I, J> >
	{
		static const int num_dimensions = 2;
		static const int ct_num_rows = meta::common_nrows<I, J>::value;
		static const int ct_num_cols = meta::common_ncols<I, J>::value;

		static const bool is_readonly = true;

		typedef matrix_shape<ct_num_rows, ct_num_cols> shape_type;
		typedef typename matrix_traits<Mat>::value_type value_type;
		typedef typename meta::common_domain<Mat, I, J>::type domain;
	};

	template<class Mat, class I, class J>
	struct matrix_traits<select_expr<Mat, I, J> >
	{
		static const int num_dimensions = 2;
		static const int ct_num_rows = meta::nelems<I>::value;
		static const int ct_num_cols = meta::nelems<J>::value;

		static const bool is_readonly = true;

		typedef matrix_shape<ct_num_rows, ct_num_cols> shape_type;
		typedef typename matrix_traits<Mat>::value_type value_type;
		typedef typename meta::common_domain<Mat, I, J>::type domain;
	};


	/********************************************
	 *
	 *  Linear selection
	 *
	 ********************************************/

	template<class Mat, class L>
	class selectl_expr
	: public IMatrixXpr<selectl_expr<Mat, L>, typename matrix_traits<Mat>::value_type>
	{
		static_assert( meta::supports_linear_index<Mat>::value,
				"Mat should support linear indexing" );

	public:
		LMAT_ENSURE_INLINE
		selectl_expr(const Mat& src, const L& idx)
		: m_src(src), m_inds(idx) { }

		LMAT_ENSURE_INLINE const Mat& source() const
		{
			return m_src;
		}

		LMAT_ENSURE_INLINE const L& indices() const
		{
			return m_inds;
		}

		LMAT_ENSURE_INLINE index_t nrows() const
		{
			return m_inds.nrows();
		}

		LMAT_ENSURE_INLINE index_t ncolumns() const
		{
			return m_inds.ncolumns();
		}

		LMAT_ENSURE_INLINE index_t nelems() const
		{
			return m_inds.nelems();
		}

		LMAT_ENSURE_INLINE
		matrix_shape<meta::nrows<L>::value, meta::ncols<L>::value>
		shape() const
		{
			return m_inds.shape();
		}

	private:
		const Mat& m_src;
		const L& m_inds;
	};


	template<class Mat, class I, class J>
	class selectl2_expr
	: public IMatrixXpr<selectl2_expr<Mat, I, J>, typename matrix_traits<Mat>::value_type>
	{
		typedef matrix_shape<
				meta::common_nrows<I, J>::value,
				meta::common_ncols<I, J>::value> shape_type;

	public:
		LMAT_ENSURE_INLINE
		selectl2_expr(const Mat& src, const I& si, const J& sj)
		: m_shape(common_nrows(si, sj), common_ncols(si, sj))
		, m_src(src), m_si(si), m_sj(sj) { }

		LMAT_ENSURE_INLINE const Mat& source() const
		{
			return m_src;
		}

		LMAT_ENSURE_INLINE const I& i_subs() const
		{
			return m_si;
		}

		LMAT_ENSURE_INLINE const J& j_subs() const
		{
			return m_sj;
		}

		LMAT_ENSURE_INLINE index_t nrows() const
		{
			return m_shape.nrows();
		}

		LMAT_ENSURE_INLINE index_t ncolumns() const
		{
			return m_shape.ncolumns();
		}

		LMAT_ENSURE_INLINE index_t nelems() const
		{
			return m_shape.nelems();
		}

		LMAT_ENSURE_INLINE
		shape_type shape() const
		{
			return m_shape;
		}

	private:
		shape_type m_shape;
		const Mat& m_src;
		const I& m_si;
		const J& m_sj;
	};


	// construction

	template<class S, typename T, class L, typename TI>
	LMAT_ENSURE_INLINE
	inline typename meta::enable_if<meta::supports_linear_index<S>,
	selectl_expr<S, L> >::type
	selectl(const IRegularMatrix<S, T>& s, const IRegularMatrix<L, TI>& idx)
	{
		return selectl_expr<S, L>(s.derived(), idx.derived());
	}

	template<class S, typename T, class I, typename TI, class J, typename TJ>
	LMAT_ENSURE_INLINE
	inline typename meta::enable_if<meta::supports_linear_index<S>,
	selectl2_expr<S, I, J> >::type
	selectl(const IRegularMatrix<S, T>& s, const IRegularMatrix<I, TI>& si, const IRegularMatrix<J, TJ>& sj)
	{
		return selectl2_expr<S, I, J>(s.derived(), si.derived(), sj.derived());
	}


	// evaluation

	namespace internal
	{
		template<class S, class L, class D, bool UseLinear>
		struct selectl_eval;

		template<class S, class L, class D>
		struct selectl_eval<S, L, D, true>
		{
			LMAT_ENSURE_INLINE
			static void run(const S& s, const L& idx, D& dst)
			{
				const index_t n = idx.nelems();
				for (index_t i = 0; i < n; ++i)
					dst[i] = s[(index_t)idx[i]];
			}
		};

		template<class S, class L, class D>
		struct selectl_eval<S, L, D, false>
		{
			LMAT_ENSURE_INLINE
			static void run(const S& s, const L& idx, D& dst)
			{
				const index_t m = idx.nrows();
				const index_t n = idx.ncolumns();
				for (index_t j = 0; j < n; ++j)
				{
					for (index_t i = 0; i < m; ++i)
						dst.elem(i, j) = s[(index_t)idx.elem(i, j)];
				}
			}
		};


		template<class S, class I, class J, class D, bool UseLinear>
		struct selectl2_eval;

		template<class S, class I, class J, class D>
		struct selectl2_eval<S, I, J, D, true>
		{
			LMAT_ENSURE_INLINE
			static void run(const S& s, const I& si, const J& sj, D& dst)
			{
				const index_t n = dst.nelems();
				for (index_t i = 0; i < n; ++i)
					dst[i] = s((index_t)si[i], (index_t)sj[i]);
			}
		};

		template<class S, class I, class J, class D>
		struct selectl2_eval<S, I, J, D, false>
		{
			LMAT_ENSURE_INLINE
			static void run(const S& s, const I& si, const J& sj, D& dst)
			{
				const index_t m = dst.nrows();
				const index_t n = dst.ncolumns();
				for (index_t j = 0; j < n; ++j)
				{
					for (index_t i = 0; i < m; ++i)
						dst.elem(i, j) = s((index_t)si.elem(i, j), (index_t)sj.elem(i, j));
				}
			}
		};
	}


	template<typename T, class S, class L, class D>
	LMAT_ENSURE_INLINE
	inline void evaluate(const selectl_expr<S, L>& expr, IRegularMatrix<D, T>& dst)
	{
		const bool use_linear =
				meta::supports_linear_index<L>::value &&
				meta::supports_linear_index<D>::value;

		internal::selectl_eval<S, L, D, use_linear>::run(
				expr.source(), expr.indices(), dst.derived());
	}

	template<typename T, class S, class I, class J, class D>
	LMAT_ENSURE_INLINE
	inline void evaluate(const selectl2_expr<S, I, J>& expr, IRegularMatrix<D, T>& dst)
	{
		const bool use_linear =
				meta::supports_linear_index<I>::value &&
				meta::supports_linear_index<J>::value &&
				meta::supports_linear_index<D>::value;

		internal::selectl2_eval<S, I, J, D, use_linear>::run(
				expr.source(), expr.i_subs(), expr.j_subs(), dst.derived());
	}



	/********************************************
	 *
	 *  subscript selection
	 *
	 ********************************************/

	template<class Mat, class I, class J>
	class select_expr
	: public IMatrixXpr<select_expr<Mat, I, J>, typename matrix_traits<Mat>::value_type>
	{
		typedef matrix_shape<meta::nelems<I>::value, meta::nelems<J>::value> shape_type;

		static_assert( meta::supports_linear_index<I>::value, "I should support linear indexing" );
		static_assert( meta::supports_linear_index<J>::value, "J should support linear indexing" );

	public:
		LMAT_ENSURE_INLINE
		select_expr(const Mat& src, const I& si, const J& sj)
		: m_src(src), m_si(si), m_sj(sj) { }

		LMAT_ENSURE_INLINE const Mat& source() const
		{
			return m_src;
		}

		LMAT_ENSURE_INLINE const Mat& i_subs() const
		{
			return m_si;
		}

		LMAT_ENSURE_INLINE const Mat& j_subs() const
		{
			return m_sj;
		}

		LMAT_ENSURE_INLINE index_t nrows() const
		{
			return m_si.nelems();
		}

		LMAT_ENSURE_INLINE index_t ncolumns() const
		{
			return m_sj.nelems();
		}

		LMAT_ENSURE_INLINE index_t nelems() const
		{
			return nrows() * ncolumns();
		}

		LMAT_ENSURE_INLINE
		shape_type shape() const
		{
			return shape_type(nrows(), ncolumns());
		}

	private:
		const Mat& m_src;
		const I& m_si;
		const J& m_sj;
	};


	template<class S, typename T, class I, typename TI, class J, typename TJ>
	LMAT_ENSURE_INLINE
	inline typename meta::enable_if<
		meta::and_<meta::supports_linear_index<I>, meta::supports_linear_index<J> >,
	select_expr<S, I, J> >::type
	select(const IRegularMatrix<S, T>& s, const IRegularMatrix<I, TI>& si, const IRegularMatrix<J, TJ>& sj)
	{
		return select_expr<S, I, J>(s, si, sj);
	}


	template<typename T, class S, class I, class J, class D>
	LMAT_ENSURE_INLINE
	inline void evaluate(const select_expr<S, I, J>& expr, IRegularMatrix<D, T>& dst)
	{
		const S& s = expr.source();
		const I& si = expr.i_subs();
		const J& sj = expr.j_subs();
		D& d = dst.derived();

		const index_t m = dst.nrows();
		const index_t n = dst.ncolumns();

		for (index_t j = 0; j < n; ++j)
		{
			for (index_t i = 0; i < m; ++i)
			{
				d.elem(i, j) = s.elem((index_t)si[i], (index_t)sj[j]);
			}
		}
	}



}

#endif /* MATRIX_SELECT_H_ */
