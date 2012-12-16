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
#include <light_mat/math/simd_packs.h>

namespace lmat
{

	// forward declarations

	template<typename T, typename U> class contvec_reader;
	template<typename T, typename U> class stepvec_reader;

	template<typename T, typename U> class contvec_writer;
	template<typename T, typename U> class stepvec_writer;

	class scalar_vec_accessor_base
	{
	public:
		LMAT_ENSURE_INLINE
		nil_t done_scalar(index_t ) const { return nil_t(); }

		LMAT_ENSURE_INLINE
		nil_t finalize() const { return nil_t(); }
	};

	class simd_vec_accessor_base
	{
	public:
		LMAT_ENSURE_INLINE
		nil_t done_scalar(index_t ) const { return nil_t(); }

		LMAT_ENSURE_INLINE
		nil_t done_pack(index_t ) const { return nil_t(); }

		LMAT_ENSURE_INLINE
		nil_t finalize() const { return nil_t(); }
	};


	/********************************************
	 *
	 *  reader classes
	 *
	 ********************************************/

	// contvec_reader

	template<typename T>
	class contvec_reader<T, atags::scalar> : public scalar_vec_accessor_base
	{
	public:
		LMAT_ENSURE_INLINE
		explicit contvec_reader(const T* p) : m_pdata(p) { }

		LMAT_ENSURE_INLINE
		T scalar(index_t i) const
		{
			return m_pdata[i];
		}

	private:
		const T* m_pdata;
	};


	template<typename T, typename Kind>
	class contvec_reader<T, atags::simd<T, Kind> > : public simd_vec_accessor_base
	{
	public:
		typedef math::simd_pack<T, Kind> pack_type;

		LMAT_ENSURE_INLINE
		explicit contvec_reader(const T* p) : m_pdata(p) { }

		LMAT_ENSURE_INLINE
		T scalar(index_t i) const
		{
			return m_pdata[i];
		}

		LMAT_ENSURE_INLINE
		pack_type pack(index_t i) const
		{
			return pack_type(m_pdata + i);
		}

	private:
		const T* m_pdata;
	};



	// stepvec_reader

	template<typename T>
	class stepvec_reader<T, atags::scalar> : public scalar_vec_accessor_base
	{
	public:
		LMAT_ENSURE_INLINE
		explicit stepvec_reader(const T* p, index_t step)
		: m_pdata(p), m_step(step) { }

		T scalar(index_t i) const
		{
			return m_pdata[i * m_step];
		}

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
			typedef typename matrix_traits<Mat>::value_type T;
			typedef contvec_reader<T, U> type;

			LMAT_ENSURE_INLINE
			static type get(const Mat& mat)
			{
				return type(mat.ptr_data());
			}
		};

		template<class Mat, typename U>
		struct stepcol_reader_map
		{
			typedef typename matrix_traits<Mat>::value_type T;
			typedef stepvec_reader<T, U> type;

			LMAT_ENSURE_INLINE
			static type get(const Mat& mat)
			{
				return type(mat.ptr_data(), mat.row_stride());
			}
		};

		template<class Mat, typename U>
		struct steprow_reader_map
		{
			typedef typename matrix_traits<Mat>::value_type T;
			typedef stepvec_reader<T, U> type;

			LMAT_ENSURE_INLINE
			static type get(const Mat& mat)
			{
				return type(mat.ptr_data(), mat.col_stride());
			}
		};


		struct invalid_reader_map { };

		template<class Mat, typename U>
		struct vec_reader_map
		{
			typedef
			typename meta::if_<meta::is_continuous<Mat>,	contvec_reader_map<Mat, U>,
			typename meta::if_<meta::is_col<Mat>, 			stepcol_reader_map<Mat, U>,
			typename meta::if_<meta::is_row<Mat>, 			steprow_reader_map<Mat, U>,
															invalid_reader_map
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
		return internal::vec_reader_map<Mat, U>::get(wrap.arg());
	}



	/********************************************
	 *
	 *  writer classes
	 *
	 ********************************************/

	// contvec_writer

	template<typename T>
	class contvec_writer<T, atags::scalar> : public scalar_vec_accessor_base
	{
	public:
		LMAT_ENSURE_INLINE
		explicit contvec_writer(T* p) : m_pdata(p) { }

		T& scalar(index_t ) const
		{
			return m_temp;
		}

		nil_t done_scalar(index_t i) const
		{
			m_pdata[i] = m_temp;
			return nil_t();
		}

	private:
		mutable T m_temp;
		T* m_pdata;
	};


	template<typename T, typename Kind>
	class contvec_writer<T, atags::simd<T, Kind> > : public simd_vec_accessor_base
	{
	public:
		typedef math::simd_pack<T, Kind> pack_type;

		LMAT_ENSURE_INLINE
		explicit contvec_writer(const T* p) : m_pdata(p) { }

		LMAT_ENSURE_INLINE
		T& scalar(index_t) const
		{
			return m_stemp;
		}

		LMAT_ENSURE_INLINE
		pack_type& pack(index_t) const
		{
			return m_ptemp;
		}

		LMAT_ENSURE_INLINE
		nil_t done_scalar(index_t i) const
		{
			m_pdata[i] = m_stemp;
			return nil_t();
		}

		LMAT_ENSURE_INLINE
		nil_t done_pack(index_t i) const
		{
			m_ptemp.store(m_pdata + i);
			return nil_t();
		}

	private:
		mutable T m_stemp;
		mutable pack_type m_ptemp;
		const T* m_pdata;
	};

	// stepvec_writer

	template<typename T>
	class stepvec_writer<T, atags::scalar> : public scalar_vec_accessor_base
	{
	public:
		LMAT_ENSURE_INLINE
		explicit stepvec_writer(T* p, index_t step)
		: m_pdata(p), m_step(step) { }

		T& scalar(index_t ) const
		{
			return m_temp;
		}

		nil_t done_scalar(index_t i) const
		{
			m_pdata[i * m_step] = m_temp;
			return nil_t();
		}

	private:
		mutable T m_temp;
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
			typedef typename matrix_traits<Mat>::value_type T;
			typedef contvec_writer<T, U> type;

			LMAT_ENSURE_INLINE
			static type get(Mat& mat)
			{
				return type(mat.ptr_data());
			}
		};

		template<class Mat, typename U>
		struct stepcol_writer_map
		{
			typedef typename matrix_traits<Mat>::value_type T;
			typedef stepvec_writer<T, U> type;

			LMAT_ENSURE_INLINE
			static type get(Mat& mat)
			{
				return type(mat.ptr_data(), mat.row_stride());
			}
		};

		template<class Mat, typename U>
		struct steprow_writer_map
		{
			typedef typename matrix_traits<Mat>::value_type T;
			typedef stepvec_writer<T, U> type;

			LMAT_ENSURE_INLINE
			static type get(Mat& mat)
			{
				return type(mat.ptr_data(), mat.col_stride());
			}
		};

		struct invalid_writer_map { };

		template<class Mat, typename U>
		struct vec_writer_map
		{
			typedef
			typename meta::if_<meta::is_continuous<Mat>,	contvec_writer_map<Mat, U>,
			typename meta::if_<meta::is_col<Mat>, 			stepcol_writer_map<Mat, U>,
			typename meta::if_<meta::is_row<Mat>, 			steprow_writer_map<Mat, U>,
															invalid_writer_map
			>::type >::type >::type internal_map;

			typedef typename internal_map::type type;

			LMAT_ENSURE_INLINE
			static type get(Mat& mat)
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
		return internal::vec_writer_map<Mat, U>::get(wrap.arg());
	}
}

#endif /* VEC_ACCESSORS_H_ */
