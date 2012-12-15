/**
 * @file vec_accessors.h
 *
 * @brief single-vector accessors
 *
 * @author Dahua Lin
 */

#ifndef LIGHTMAT_VEC_ACCESSORS_H_
#define LIGHTMAT_VEC_ACCESSORS_H_

#include <light_mat/mateval/mateval_fwd.h>
#include <light_mat/matrix/matrix_concepts.h>

#define LMAT_TRIVIAL_DONE_( Nam ) \
	LMAT_ENSURE_INLINE nil_t done_##Nam(index_t ) const { return nil_t(); }

namespace lmat
{

	// forward declarations

	template<typename T, typename U> class contvec_reader;
	template<typename T, typename U> class stepvec_reader;

	template<typename T, typename U> class contvec_writer;
	template<typename T, typename U> class stepvec_writer;

	/********************************************
	 *
	 *  reader classes
	 *
	 ********************************************/

	// contvec_reader

	template<typename T>
	class contvec_reader<T, atags::scalar>
	{
	public:
		LMAT_ENSURE_INLINE
		contvec_reader(const T* p) : m_pdata(p) { }

		T scalar(index_t i) const
		{
			return m_pdata[i];
		}

		LMAT_TRIVIAL_DONE_(scalar)

	private:
		const T* m_pdata;
	};


	// stepvec_reader

	template<typename T>
	class stepvec_reader<T, atags::scalar>
	{
	public:
		LMAT_ENSURE_INLINE
		stepvec_reader(const T* p, index_t step)
		: m_pdata(p), m_step(step) { }

		T scalar(index_t i) const
		{
			return m_pdata[i * m_step];
		}

		LMAT_TRIVIAL_DONE_(scalar)

	private:
		const T* m_pdata;
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
		struct contvec_reader_map
		{
			typedef contvec_reader<Mat, U> type;

			LMAT_ENSURE_INLINE
			static type get(const Mat& mat)
			{
				return type(mat.ptr_data());
			}
		};

		template<class Mat, typename U>
		struct stepcol_reader_map
		{
			typedef stepvec_reader<Mat, U> type;

			LMAT_ENSURE_INLINE
			static type get(const Mat& mat)
			{
				return type(mat.ptr_data(), mat.row_stride());
			}
		};

		template<class Mat, typename U>
		struct steprow_reader_map
		{
			typedef stepvec_reader<Mat, U> type;

			LMAT_ENSURE_INLINE
			static type get(const Mat& mat)
			{
				return type(mat.ptr_data(), mat.col_stride());
			}
		};

		template<class Mat, typename U>
		struct invalid_reader_map
		{
			typedef nil_t type;

			LMAT_ENSURE_INLINE
			static type get(const Mat& mat)
			{
				throw invalid_operation("The input matrix is not eligible for vector-reading.");
			}
		};

		template<class Mat, typename U>
		struct vec_reader_map
		{
			typedef
			typename meta::if_<meta::is_continuous<Mat>,	contvec_reader_map<Mat, U>,
			typename meta::if_<meta::is_col<Mat>, 			stepcol_reader_map<Mat, U>,
			typename meta::if_<meta::is_row<Mat>, 			steprow_reader_map<Mat, U>,
															invalid_reader_map<Mat, U>
			>::type >::type >::type internal_map;

			typedef typename internal_map::type type;

			LMAT_ENSURE_INLINE
			static type get(const Mat& mat)
			{
				return internal_map::get(mat);
			}
		};
	};


	template<class Mat, typename U>
	LMAT_ENSURE_INLINE
	inline typename internal::vec_reader_map<Mat, U>::type
	make_vec_accessor(U, const in_wrap<Mat, atags::normal>& wrap)
	{
		return internal::vec_reader_map<Mat, U>(wrap.arg());
	}



	/********************************************
	 *
	 *  writer classes
	 *
	 ********************************************/

	// contvec_writer

	template<typename T>
	class contvec_writer<T, atags::scalar>
	{
	public:
		LMAT_ENSURE_INLINE
		contvec_writer(const T* p) : m_pdata(p) { }

		T& scalar(index_t i) const
		{
			return m_pdata[i];
		}

		LMAT_TRIVIAL_DONE_(scalar)

	private:
		const T* m_pdata;
	};


	// stepvec_writer

	template<typename T>
	class stepvec_writer<T, atags::scalar>
	{
	public:
		LMAT_ENSURE_INLINE
		stepvec_writer(T* p, index_t step)
		: m_pdata(p), m_step(step) { }

		T& scalar(index_t i) const
		{
			return m_pdata[i * m_step];
		}

		LMAT_TRIVIAL_DONE_(scalar)

	private:
		T* m_pdata;
		index_t m_step;
	};



	/********************************************
	 *
	 *  writer maps
	 *
	 ********************************************/

	namespace internal
	{
		template<class Mat, typename U>
		struct contvec_writer_map
		{
			typedef contvec_writer<Mat, U> type;

			LMAT_ENSURE_INLINE
			static type get(Mat& mat)
			{
				return type(mat.ptr_data());
			}
		};

		template<class Mat, typename U>
		struct stepcol_writer_map
		{
			typedef stepvec_writer<Mat, U> type;

			LMAT_ENSURE_INLINE
			static type get(Mat& mat)
			{
				return type(mat.ptr_data(), mat.row_stride());
			}
		};

		template<class Mat, typename U>
		struct steprow_writer_map
		{
			typedef stepvec_writer<Mat, U> type;

			LMAT_ENSURE_INLINE
			static type get(Mat& mat)
			{
				return type(mat.ptr_data(), mat.col_stride());
			}
		};

		template<class Mat, typename U>
		struct invalid_writer_map
		{
			typedef nil_t type;

			LMAT_ENSURE_INLINE
			static type get(Mat& mat)
			{
				throw invalid_operation("The input matrix is not eligible for vector-writing.");
			}
		};

		template<class Mat, typename U>
		struct vec_writer_map
		{
			typedef
			typename meta::if_<meta::is_continuous<Mat>,	contvec_writer_map<Mat, U>,
			typename meta::if_<meta::is_col<Mat>, 			stepcol_writer_map<Mat, U>,
			typename meta::if_<meta::is_row<Mat>, 			steprow_writer_map<Mat, U>,
															invalid_writer_map<Mat, U>
			>::type >::type >::type internal_map;

			typedef typename internal_map::type type;

			LMAT_ENSURE_INLINE
			static type get(const Mat& mat)
			{
				return internal_map::get(mat);
			}
		};
	};


	template<class Mat, typename U>
	LMAT_ENSURE_INLINE
	inline typename internal::vec_writer_map<Mat, U>::type
	make_vec_accessor(U, const out_wrap<Mat, atags::normal>& wrap)
	{
		return internal::vec_writer_map<Mat, U>(wrap.arg());
	}




}

#endif /* VEC_ACCESSORS_H_ */
