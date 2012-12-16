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

			LMAT_ENSURE_INLINE
			static type get(const Mat& mat)
			{
				return type(mat);
			}
		};
	}

	template<class Mat, typename U>
	LMAT_ENSURE_INLINE
	inline typename internal::multicol_reader_map<Mat, U>::type
	make_multicol_accessor(U, const in_wrap<Mat, atags::normal>& wrap)
	{
		return internal::multicol_reader_map<Mat, U>::get(wrap.arg());
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

			LMAT_ENSURE_INLINE
			static type get(Mat& mat)
			{
				return type(mat);
			}
		};
	}

	template<class Mat, typename U>
	LMAT_ENSURE_INLINE
	inline typename internal::multicol_writer_map<Mat, U>::type
	make_multicol_accessor(U, const out_wrap<Mat, atags::normal>& wrap)
	{
		return internal::multicol_writer_map<Mat, U>::get(wrap.arg());
	}


}


#endif /* MULTICOL_ACCESSORS_H_ */
