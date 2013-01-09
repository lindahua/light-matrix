/**
 * @file multicol_accessors.h
 *
 * @brief Multi-column accessor classes
 *
 * @author Dahua Lin
 */

#ifdef _MSC_VER
#pragma once
#endif

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
	class multicol_single_reader : public multicol_accessor_base
	{
	public:
		typedef single_reader<T, U> col_accessor_type;

		LMAT_ENSURE_INLINE
		explicit multicol_single_reader(const T& v)
		: m_intern(v) { }

		LMAT_ENSURE_INLINE
		col_accessor_type col(index_t j) const
		{
			return m_intern;
		}

	private:
		col_accessor_type m_intern;
	};

	template<typename T, typename U>
	class multi_contcol_reader : public multicol_accessor_base
	{
	public:
		typedef contvec_reader<T, U> col_accessor_type;

		template<class Mat>
		LMAT_ENSURE_INLINE
		explicit multi_contcol_reader(const Mat& mat)
		: m_pbase(mat.ptr_data()), m_colstride(mat.col_stride())
		{ }

		LMAT_ENSURE_INLINE
		col_accessor_type col(index_t j) const
		{
			return col_accessor_type(m_pbase + m_colstride * j);
		}

	private:
		const T *m_pbase;
		index_t m_colstride;
	};


	template<typename T, typename U>
	class multi_stepcol_reader : public multicol_accessor_base
	{
	public:
		typedef stepvec_reader<T, U> col_accessor_type;

		template<class Mat>
		LMAT_ENSURE_INLINE
		explicit multi_stepcol_reader(const Mat& mat)
		: m_pbase(mat.ptr_data())
		, m_rowstride(mat.row_stride())
		, m_colstride(mat.col_stride())
		{ }

		LMAT_ENSURE_INLINE
		col_accessor_type col(index_t j) const
		{
			return col_accessor_type(m_pbase + m_colstride * j, m_rowstride);
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
		typedef contvec_reader<T, U> col_accessor_type;

		template<class Mat>
		LMAT_ENSURE_INLINE
		explicit repeat_contcol_reader(const Mat& col)
		: m_pbase(col.ptr_data())
		{ }

		LMAT_ENSURE_INLINE
		col_accessor_type col(index_t ) const
		{
			return col_accessor_type(m_pbase);
		}

	private:
		const T *m_pbase;
	};


	template<typename T, typename U>
	class repeat_stepcol_reader : public multicol_accessor_base
	{
	public:
		typedef stepvec_reader<T, U> col_accessor_type;

		template<class Mat>
		LMAT_ENSURE_INLINE
		explicit repeat_stepcol_reader(const Mat& col)
		: m_pbase(col.ptr_data())
		, m_step(col.row_stride())
		{ }

		LMAT_ENSURE_INLINE
		col_accessor_type col(index_t ) const
		{
			return col_accessor_type(m_pbase, m_step);
		}

	private:
		const T *m_pbase;
		index_t m_step;
	};


	template<typename T, typename U>
	class repeat_controw_reader : public multicol_accessor_base
	{
	public:
		typedef single_reader<T, U> col_accessor_type;

		template<class Mat>
		LMAT_ENSURE_INLINE
		explicit repeat_controw_reader(const Mat& row)
		: m_pbase(row.ptr_data())
		{ }

		LMAT_ENSURE_INLINE
		col_accessor_type col(index_t j) const
		{
			return col_accessor_type(m_pbase[j]);
		}

	private:
		const T *m_pbase;
	};


	template<typename T, typename U>
	class repeat_steprow_reader : public multicol_accessor_base
	{
	public:
		typedef single_reader<T, U> col_accessor_type;

		template<class Mat>
		LMAT_ENSURE_INLINE
		explicit repeat_steprow_reader(const Mat& row)
		: m_pbase(row.ptr_data())
		, m_step(row.col_stride())
		{ }

		LMAT_ENSURE_INLINE
		col_accessor_type col(index_t j) const
		{
			return col_accessor_type(m_pbase[j * m_step]);
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
		struct invalid_multicol_reader { };

		template<class Mat, typename U>
		struct multicol_reader_map
		{
			typedef typename matrix_traits<Mat>::value_type T;

			typedef typename meta::if_<
					meta::is_regular_mat<Mat>,
					typename meta::if_<meta::is_percol_contiguous<Mat>,
						multi_contcol_reader<T, U>,
						multi_stepcol_reader<T, U> >::type,
					invalid_multicol_reader
			>::type type;

			LMAT_ENSURE_INLINE
			static type get(const Mat& mat)
			{
				return type(mat);
			}
		};

		template<class Mat, typename U>
		struct repeatcol_reader_map
		{
			typedef typename matrix_traits<Mat>::value_type T;

			typedef typename meta::if_<
					meta::is_regular_mat<Mat>,
						typename meta::if_<meta::is_percol_contiguous<Mat>,
						repeat_contcol_reader<T, U>,
						repeat_stepcol_reader<T, U> >::type,
					invalid_multicol_reader
			>::type type;
		};

		template<class Mat, typename U>
		struct repeatrow_reader_map
		{
			typedef typename matrix_traits<Mat>::value_type T;

			typedef typename meta::if_<
					meta::is_regular_mat<Mat>,
						typename meta::if_<meta::is_contiguous<Mat>,
						repeat_controw_reader<T, U>,
						repeat_steprow_reader<T, U> >::type,
					invalid_multicol_reader
			>::type type;
		};
	}

	template<class Mat, typename U>
	LMAT_ENSURE_INLINE
	inline typename internal::multicol_reader_map<Mat, U>::type
	make_multicol_accessor(U, const in_wrap<Mat, atags::normal>& wrap)
	{
		return internal::multicol_reader_map<Mat, U>::get(wrap.arg());
	}

	template<typename T, typename U>
	LMAT_ENSURE_INLINE
	inline multicol_single_reader<T, U>
	make_multicol_accessor(U, const in_wrap<T, atags::single>& wrap)
	{
		return multicol_single_reader<T, U>(wrap.arg());
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
		typedef contvec_writer<T, U> col_accessor_type;

		template<class Mat>
		LMAT_ENSURE_INLINE
		explicit multi_contcol_writer(Mat& mat)
		: m_pbase(mat.ptr_data()), m_colstride(mat.col_stride())
		{ }

		LMAT_ENSURE_INLINE
		col_accessor_type col(index_t j) const
		{
			return col_accessor_type(m_pbase + m_colstride * j);
		}

	private:
		T *m_pbase;
		index_t m_colstride;
	};


	template<typename T, typename U>
	class multi_stepcol_writer : public multicol_accessor_base
	{
	public:
		typedef stepvec_writer<T, U> col_accessor_type;

		template<class Mat>
		LMAT_ENSURE_INLINE
		explicit multi_stepcol_writer(Mat& mat)
		: m_pbase(mat.ptr_data())
		, m_rowstride(mat.row_stride())
		, m_colstride(mat.col_stride())
		{ }

		LMAT_ENSURE_INLINE
		col_accessor_type col(index_t j) const
		{
			return col_accessor_type(m_pbase + m_colstride * j, m_rowstride);
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
					meta::is_percol_contiguous<Mat>,
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
		typedef contvec_updater<T, U> col_accessor_type;

		template<class Mat>
		LMAT_ENSURE_INLINE
		explicit multi_contcol_updater(Mat& mat)
		: m_pbase(mat.ptr_data()), m_colstride(mat.col_stride())
		{ }

		LMAT_ENSURE_INLINE
		col_accessor_type col(index_t j) const
		{
			return col_accessor_type(m_pbase + m_colstride * j);
		}

	private:
		T *m_pbase;
		index_t m_colstride;
	};


	template<typename T, typename U>
	class multi_stepcol_updater : public multicol_accessor_base
	{
	public:
		typedef stepvec_updater<T, U> col_accessor_type;

		template<class Mat>
		LMAT_ENSURE_INLINE
		explicit multi_stepcol_updater(Mat& mat)
		: m_pbase(mat.ptr_data())
		, m_rowstride(mat.row_stride())
		, m_colstride(mat.col_stride())
		{ }

		LMAT_ENSURE_INLINE
		col_accessor_type col(index_t j) const
		{
			return col_accessor_type(m_pbase + m_colstride * j, m_rowstride);
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
					meta::is_percol_contiguous<Mat>,
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


	/********************************************
	 *
	 *  accumulator classes
	 *
	 ********************************************/

	// full accumulator

	template<typename T, typename U>
	class multicol_sum_accumulator : public multicol_accessor_base
	{
	public:
		typedef sum_accumulator<T, U> col_accessor_type;

		LMAT_ENSURE_INLINE
		explicit multicol_sum_accumulator(T& s)
		: m_p(&s) { }

		LMAT_ENSURE_INLINE
		col_accessor_type col(index_t) const
		{
			return col_accessor_type(*m_p);
		}

	private:
		T *m_p;
	};

	template<typename T, typename U>
	class multicol_max_accumulator : public multicol_accessor_base
	{
	public:
		typedef max_accumulator<T, U> col_accessor_type;

		LMAT_ENSURE_INLINE
		explicit multicol_max_accumulator(T& s)
		: m_p(&s) { }

		LMAT_ENSURE_INLINE
		col_accessor_type col(index_t) const
		{
			return col_accessor_type(*m_p);
		}

	private:
		T *m_p;
	};

	template<typename T, typename U>
	class multicol_min_accumulator : public multicol_accessor_base
	{
	public:
		typedef min_accumulator<T, U> col_accessor_type;

		LMAT_ENSURE_INLINE
		explicit multicol_min_accumulator(T& s)
		: m_p(&s) { }

		LMAT_ENSURE_INLINE
		col_accessor_type col(index_t) const
		{
			return col_accessor_type(*m_p);
		}

	private:
		T *m_p;
	};


	// colwise accumulator

	template<typename T, typename U>
	class colwise_sum_accumulator : public multicol_accessor_base
	{
	public:
		typedef sum_accumulator<T, U> col_accessor_type;

		template<class Mat>
		LMAT_ENSURE_INLINE
		explicit colwise_sum_accumulator(Mat& row)
		: m_pbase(row.ptr_data()) { }

		LMAT_ENSURE_INLINE
		col_accessor_type col(index_t j) const
		{
			return col_accessor_type(m_pbase[j]);
		}

	private:
		T *m_pbase;
	};

	template<typename T, typename U>
	class colwise_sum_accumulator_x : public multicol_accessor_base
	{
	public:
		typedef sum_accumulator<T, U> col_accessor_type;

		template<class Mat>
		LMAT_ENSURE_INLINE
		explicit colwise_sum_accumulator_x(Mat& row)
		: m_pbase(row.ptr_data()), m_step(row.col_stride()) { }

		LMAT_ENSURE_INLINE
		col_accessor_type col(index_t j) const
		{
			return col_accessor_type(m_pbase[j * m_step]);
		}

	private:
		T *m_pbase;
		index_t m_step;
	};


	template<typename T, typename U>
	class colwise_max_accumulator : public multicol_accessor_base
	{
	public:
		typedef max_accumulator<T, U> col_accessor_type;

		template<class Mat>
		LMAT_ENSURE_INLINE
		explicit colwise_max_accumulator(Mat& row)
		: m_pbase(row.ptr_data()) { }

		LMAT_ENSURE_INLINE
		col_accessor_type col(index_t j) const
		{
			return col_accessor_type(m_pbase[j]);
		}

	private:
		T *m_pbase;
	};

	template<typename T, typename U>
	class colwise_max_accumulator_x : public multicol_accessor_base
	{
	public:
		typedef max_accumulator<T, U> col_accessor_type;

		template<class Mat>
		LMAT_ENSURE_INLINE
		explicit colwise_max_accumulator_x(Mat& row)
		: m_pbase(row.ptr_data()), m_step(row.col_stride()) { }

		LMAT_ENSURE_INLINE
		col_accessor_type col(index_t j) const
		{
			return col_accessor_type(m_pbase[j * m_step]);
		}

	private:
		T *m_pbase;
		index_t m_step;
	};


	template<typename T, typename U>
	class colwise_min_accumulator : public multicol_accessor_base
	{
	public:
		typedef min_accumulator<T, U> col_accessor_type;

		template<class Mat>
		LMAT_ENSURE_INLINE
		explicit colwise_min_accumulator(Mat& row)
		: m_pbase(row.ptr_data()) { }

		LMAT_ENSURE_INLINE
		col_accessor_type col(index_t j) const
		{
			return col_accessor_type(m_pbase[j]);
		}

	private:
		T *m_pbase;
	};

	template<typename T, typename U>
	class colwise_min_accumulator_x : public multicol_accessor_base
	{
	public:
		typedef min_accumulator<T, U> col_accessor_type;

		template<class Mat>
		LMAT_ENSURE_INLINE
		explicit colwise_min_accumulator_x(Mat& row)
		: m_pbase(row.ptr_data()), m_step(row.col_stride()) { }

		LMAT_ENSURE_INLINE
		col_accessor_type col(index_t j) const
		{
			return col_accessor_type(m_pbase[j * m_step]);
		}

	private:
		T *m_pbase;
		index_t m_step;
	};


	// rowwise accumulator

	template<class Col, typename U>
	class rowwise_accumulator : public multicol_accessor_base
	{
	public:
		typedef typename internal::vec_updater_map<Col, U>::type col_accessor_type;

		LMAT_ENSURE_INLINE
		explicit rowwise_accumulator(Col& col)
		: m_col(col) { }

		LMAT_ENSURE_INLINE
		col_accessor_type col(index_t ) const
		{
			return internal::vec_updater_map<Col, U>::get(m_col);
		}

	private:
		Col& m_col;
	};



	/********************************************
	 *
	 *  accumulator maps
	 *
	 ********************************************/

	namespace internal
	{
		template<class Mat, typename U>
		struct multicol_colwise_accumulator_map
		{
			typedef typename matrix_traits<Mat>::value_type T;

			typedef typename meta::if_<
					meta::is_contiguous<Mat>,
					colwise_sum_accumulator<T, U>,
					colwise_sum_accumulator_x<T, U>
			>::type sum_type;

			typedef typename meta::if_<
					meta::is_contiguous<Mat>,
					colwise_max_accumulator<T, U>,
					colwise_max_accumulator_x<T, U>
			>::type max_type;

			typedef typename meta::if_<
					meta::is_contiguous<Mat>,
					colwise_min_accumulator<T, U>,
					colwise_min_accumulator_x<T, U>
			>::type min_type;
		};
	}

	template<typename T, typename U>
	LMAT_ENSURE_INLINE
	inline multicol_sum_accumulator<T, U>
	make_multicol_accessor(U, const in_out_wrap<T, atags::sum>& wrap)
	{
		return multicol_sum_accumulator<T, U>(wrap.arg());
	}

	template<typename T, typename U>
	LMAT_ENSURE_INLINE
	inline multicol_max_accumulator<T, U>
	make_multicol_accessor(U, const in_out_wrap<T, atags::max>& wrap)
	{
		return multicol_max_accumulator<T, U>(wrap.arg());
	}

	template<typename T, typename U>
	LMAT_ENSURE_INLINE
	inline multicol_min_accumulator<T, U>
	make_multicol_accessor(U, const in_out_wrap<T, atags::min>& wrap)
	{
		return multicol_min_accumulator<T, U>(wrap.arg());
	}


	template<class Mat, typename U>
	LMAT_ENSURE_INLINE
	inline typename internal::multicol_colwise_accumulator_map<Mat, U>::sum_type
	make_multicol_accessor(U, const in_out_wrap<Mat, atags::colwise_sum>& wrap)
	{
		typedef typename internal::multicol_colwise_accumulator_map<Mat, U>::sum_type type;
		return type(wrap.arg());
	}

	template<class Mat, typename U>
	LMAT_ENSURE_INLINE
	inline typename internal::multicol_colwise_accumulator_map<Mat, U>::max_type
	make_multicol_accessor(U, const in_out_wrap<Mat, atags::colwise_max>& wrap)
	{
		typedef typename internal::multicol_colwise_accumulator_map<Mat, U>::max_type type;
		return type(wrap.arg());
	}

	template<class Mat, typename U>
	LMAT_ENSURE_INLINE
	inline typename internal::multicol_colwise_accumulator_map<Mat, U>::min_type
	make_multicol_accessor(U, const in_out_wrap<Mat, atags::colwise_min>& wrap)
	{
		typedef typename internal::multicol_colwise_accumulator_map<Mat, U>::min_type type;
		return type(wrap.arg());
	}


	template<class Mat, typename U>
	LMAT_ENSURE_INLINE
	inline rowwise_accumulator<Mat, U>
	make_multicol_accessor(U, const in_out_wrap<Mat, atags::rowwise_sum>& wrap)
	{
		return rowwise_accumulator<Mat, U>(wrap.arg());
	}

	template<class Mat, typename U>
	LMAT_ENSURE_INLINE
	inline rowwise_accumulator<Mat, U>
	make_multicol_accessor(U, const in_out_wrap<Mat, atags::rowwise_max>& wrap)
	{
		return rowwise_accumulator<Mat, U>(wrap.arg());
	}

	template<class Mat, typename U>
	LMAT_ENSURE_INLINE
	inline rowwise_accumulator<Mat, U>
	make_multicol_accessor(U, const in_out_wrap<Mat, atags::rowwise_min>& wrap)
	{
		return rowwise_accumulator<Mat, U>(wrap.arg());
	}

}


#endif /* MULTICOL_ACCESSORS_H_ */
