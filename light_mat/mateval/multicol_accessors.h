/**
 * @file multicol_accessors.h
 *
 * @brief Multi-column accessor classes
 *
 * @author Dahua Lin
 */

#ifndef LIGHTMAT_MULTICOL_ACCESSORS_H_
#define LIGHTMAT_MULTICOL_ACCESSORS_H_

#include <light_mat/mateval/vec_accessors.h>

namespace lmat
{

	// forward declarations

	template<typename T, typename U> class multi_contcol_reader;
	template<typename T, typename U> class multi_stepcol_reader;

	template<typename T, typename U> class repeat_contcol_reader;
	template<typename T, typename U> class repeat_stepcol_reader;
	template<typename T, typename U> class repeat_controw_reader;
	template<typename T, typename U> class repeat_steprow_reader;

	template<typename T, typename U> class multi_contcol_writer;
	template<typename T, typename U> class multi_stepcol_writer;


	class multicol_accessor_base
	{
	public:
		LMAT_ENSURE_INLINE
		nil_t finalize() const { return nil_t(); }
	};


	/********************************************
	 *
	 *  reader classes
	 *
	 ********************************************/

	template<typename T, typename U>
	class multi_contcol_reader : public multicol_accessor_base
	{
	public:
		template<class Mat>
		LMAT_ENSURE_INLINE
		explicit multi_contcol_reader(const Mat& mat)
		: m_pbase(mat.ptr_data()), m_colstride(mat.col_stride())
		{ }

		LMAT_ENSURE_INLINE
		contvec_reader<T, U> col(index_t j) const
		{
			return contvec_reader<T, U>(m_pbase + m_colstride * j);
		}

	private:
		const T *m_pbase;
		index_t m_colstride;
	};


	template<typename T, typename U>
	class multi_stepcol_reader : public multicol_accessor_base
	{
	public:
		template<class Mat>
		LMAT_ENSURE_INLINE
		explicit multi_stepcol_reader(const Mat& mat)
		: m_pbase(mat.ptr_data())
		, m_rowstride(mat.row_stride())
		, m_colstride(mat.col_stride())
		{ }

		LMAT_ENSURE_INLINE
		stepvec_reader<T, U> col(index_t j) const
		{
			return stepvec_reader<T, U>(m_pbase + m_colstride * j, m_rowstride);
		}

	private:
		const T *m_pbase;
		index_t m_rowstride;
		index_t m_colstride;
	};


	template<typename T, typename U>
	class repeat_contcol_reader : public multicol_accessor_base
	{
	public:
		template<class Mat>
		LMAT_ENSURE_INLINE
		explicit repeat_contcol_reader(const Mat& col)
		: m_pbase(col.ptr_data())
		{ }

		LMAT_ENSURE_INLINE
		contvec_reader<T, U> col(index_t ) const
		{
			return contvec_reader<T, U>(m_pbase);
		}

	private:
		const T *m_pbase;
	};


	template<typename T, typename U>
	class repeat_stepcol_reader : public multicol_accessor_base
	{
	public:
		template<class Mat>
		LMAT_ENSURE_INLINE
		explicit repeat_stepcol_reader(const Mat& col)
		: m_pbase(col.ptr_data())
		, m_step(col.row_stride())
		{ }

		LMAT_ENSURE_INLINE
		stepvec_reader<T, U> col(index_t ) const
		{
			return stepvec_reader<T, U>(m_pbase, m_step);
		}

	private:
		const T *m_pbase;
		index_t m_step;
	};


	template<typename T, typename U>
	class repeat_controw_reader : public multicol_accessor_base
	{
	public:
		template<class Mat>
		LMAT_ENSURE_INLINE
		explicit repeat_controw_reader(const Mat& row)
		: m_pbase(row.ptr_data())
		{ }

		LMAT_ENSURE_INLINE
		single_reader<T, U> col(index_t j) const
		{
			return single_reader<T, U>(m_pbase[j]);
		}

	private:
		const T *m_pbase;
	};


	template<typename T, typename U>
	class repeat_steprow_reader : public multicol_accessor_base
	{
	public:
		template<class Mat>
		LMAT_ENSURE_INLINE
		explicit repeat_steprow_reader(const Mat& row)
		: m_pbase(row.ptr_data())
		, m_step(row.col_stride())
		{ }

		LMAT_ENSURE_INLINE
		single_reader<T, U> col(index_t j) const
		{
			return single_reader<T, U>(m_pbase[j * m_step]);
		}

	private:
		const T *m_pbase;
		index_t m_step;
	};



	/********************************************
	 *
	 *  reader maps
	 *
	 ********************************************/

	namespace internal
	{
		template<class Mat, typename U>
		struct multicol_reader_map
		{
			typedef typename matrix_traits<Mat>::value_type T;

			typedef typename meta::if_<
					meta::is_percol_continuous<Mat>,
					multi_contcol_reader<T, U>,
					multi_stepcol_reader<T, U>
			>::type type;
		};

		template<class Mat, typename U>
		struct repeatcol_reader_map
		{
			typedef typename matrix_traits<Mat>::value_type T;

			typedef typename meta::if_<
					meta::is_percol_continuous<Mat>,
					repeat_contcol_reader<T, U>,
					repeat_stepcol_reader<T, U>
			>::type type;
		};

		template<class Mat, typename U>
		struct repeatrow_reader_map
		{
			typedef typename matrix_traits<Mat>::value_type T;

			typedef typename meta::if_<
					meta::is_continuous<Mat>,
					repeat_controw_reader<T, U>,
					repeat_steprow_reader<T, U>
			>::type type;
		};
	}

	template<class Mat, typename U>
	LMAT_ENSURE_INLINE
	inline typename internal::multicol_reader_map<Mat, U>::type
	make_multicol_accessor(U, const in_wrap<Mat, atags::normal>& wrap)
	{
		typedef typename internal::multicol_reader_map<Mat, U>::type type;
		return type(wrap.arg());
	}

	template<typename T, typename U>
	LMAT_ENSURE_INLINE
	inline single_reader<T, U>
	make_multicol_accessor(U, const in_wrap<T, atags::single>& wrap)
	{
		return single_reader<T, U>(wrap.arg());
	}

	template<class Mat, typename U>
	LMAT_ENSURE_INLINE
	inline typename internal::repeatcol_reader_map<Mat, U>::type
	make_multicol_accessor(U, const in_wrap<Mat, atags::repcol>& wrap)
	{
		typedef typename internal::repeatcol_reader_map<Mat, U>::type type;
		return type(wrap.arg());
	}

	template<class Mat, typename U>
	LMAT_ENSURE_INLINE
	inline typename internal::repeatrow_reader_map<Mat, U>::type
	make_multicol_accessor(U, const in_wrap<Mat, atags::reprow>& wrap)
	{
		typedef typename internal::repeatrow_reader_map<Mat, U>::type type;
		return type(wrap.arg());
	}


	/********************************************
	 *
	 *  writer classes
	 *
	 ********************************************/

	template<typename T, typename U>
	class multi_contcol_writer : public multicol_accessor_base
	{
	public:
		template<class Mat>
		LMAT_ENSURE_INLINE
		explicit multi_contcol_writer(Mat& mat)
		: m_pbase(mat.ptr_data()), m_colstride(mat.col_stride())
		{ }

		LMAT_ENSURE_INLINE
		contvec_writer<T, U> col(index_t j) const
		{
			return contvec_writer<T, U>(m_pbase + m_colstride * j);
		}

	private:
		T *m_pbase;
		index_t m_colstride;
	};


	template<typename T, typename U>
	class multi_stepcol_writer : public multicol_accessor_base
	{
	public:
		template<class Mat>
		LMAT_ENSURE_INLINE
		explicit multi_stepcol_writer(Mat& mat)
		: m_pbase(mat.ptr_data())
		, m_rowstride(mat.row_stride())
		, m_colstride(mat.col_stride())
		{ }

		LMAT_ENSURE_INLINE
		stepvec_writer<T, U> col(index_t j) const
		{
			return stepvec_writer<T, U>(m_pbase + m_colstride * j, m_rowstride);
		}

	private:
		T *m_pbase;
		index_t m_rowstride;
		index_t m_colstride;
	};


	/********************************************
	 *
	 *  writer maps
	 *
	 ********************************************/

	namespace internal
	{
		template<class Mat, typename U>
		struct multicol_writer_map
		{
			typedef typename matrix_traits<Mat>::value_type T;

			typedef typename meta::if_<
					meta::is_percol_continuous<Mat>,
					multi_contcol_writer<T, U>,
					multi_stepcol_writer<T, U>
			>::type type;
		};
	}

	template<class Mat, typename U>
	LMAT_ENSURE_INLINE
	inline typename internal::multicol_writer_map<Mat, U>::type
	make_multicol_accessor(U, const out_wrap<Mat, atags::normal>& wrap)
	{
		typedef typename internal::multicol_writer_map<Mat, U>::type type;
		return type(wrap.arg());
	}


	/********************************************
	 *
	 *  updater classes
	 *
	 ********************************************/

	template<typename T, typename U>
	class multi_contcol_updater : public multicol_accessor_base
	{
	public:
		template<class Mat>
		LMAT_ENSURE_INLINE
		explicit multi_contcol_updater(Mat& mat)
		: m_pbase(mat.ptr_data()), m_colstride(mat.col_stride())
		{ }

		LMAT_ENSURE_INLINE
		contvec_updater<T, U> col(index_t j) const
		{
			return contvec_updater<T, U>(m_pbase + m_colstride * j);
		}

	private:
		T *m_pbase;
		index_t m_colstride;
	};


	template<typename T, typename U>
	class multi_stepcol_updater : public multicol_accessor_base
	{
	public:
		template<class Mat>
		LMAT_ENSURE_INLINE
		explicit multi_stepcol_updater(Mat& mat)
		: m_pbase(mat.ptr_data())
		, m_rowstride(mat.row_stride())
		, m_colstride(mat.col_stride())
		{ }

		LMAT_ENSURE_INLINE
		stepvec_updater<T, U> col(index_t j) const
		{
			return stepvec_updater<T, U>(m_pbase + m_colstride * j, m_rowstride);
		}

	private:
		T *m_pbase;
		index_t m_rowstride;
		index_t m_colstride;
	};


	/********************************************
	 *
	 *  updater maps
	 *
	 ********************************************/

	namespace internal
	{
		template<class Mat, typename U>
		struct multicol_updater_map
		{
			typedef typename matrix_traits<Mat>::value_type T;

			typedef typename meta::if_<
					meta::is_percol_continuous<Mat>,
					multi_contcol_updater<T, U>,
					multi_stepcol_updater<T, U>
			>::type type;
		};
	}

	template<class Mat, typename U>
	LMAT_ENSURE_INLINE
	inline typename internal::multicol_updater_map<Mat, U>::type
	make_multicol_accessor(U, const in_out_wrap<Mat, atags::normal>& wrap)
	{
		typedef typename internal::multicol_updater_map<Mat, U>::type type;
		return type(wrap.arg());
	}


}


#endif /* MULTICOL_ACCESSORS_H_ */
