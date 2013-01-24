/**
 * @file vec_accessors.h
 *
 * @brief single-vector accessors
 *
 * @author Dahua Lin
 */

#ifdef _MSC_VER
#pragma once
#endif

#ifndef LIGHTMAT_VEC_ACCESSORS_H_
#define LIGHTMAT_VEC_ACCESSORS_H_

#include <light_mat/mateval/mateval_fwd.h>
#include <light_mat/matrix/matrix_concepts.h>

#include <light_mat/math/math_base.h>

namespace lmat
{

	// forward declarations

	template<typename T, typename U> class contvec_reader;
	template<typename T, typename U> class stepvec_reader;
	template<typename T, typename U> class single_reader;

	template<typename T, typename U> class contvec_writer;
	template<typename T, typename U> class stepvec_writer;

	template<typename T, typename U> class contvec_updater;
	template<typename T, typename U> class stepvec_updater;

	template<typename T, typename U> class sum_accumulator;
	template<typename T, typename U> class max_accumulator;
	template<typename T, typename U> class min_accumulator;


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
		nil_t begin_packs() const { return nil_t(); }

		LMAT_ENSURE_INLINE
		nil_t end_packs() const { return nil_t(); }

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
	class contvec_reader<T, scalar_> : public scalar_vec_accessor_base
	{
	public:
		typedef T scalar_type;

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
	class contvec_reader<T, simd_<Kind> > : public simd_vec_accessor_base
	{
	public:
		typedef T scalar_type;
		typedef Kind simd_kind;
		typedef simd_pack<T, Kind> pack_type;

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
	class stepvec_reader<T, scalar_> : public scalar_vec_accessor_base
	{
	public:
		typedef T scalar_type;

		LMAT_ENSURE_INLINE
		explicit stepvec_reader(const T* p, index_t step)
		: m_pdata(p), m_step(step) { }

		LMAT_ENSURE_INLINE
		T scalar(index_t i) const
		{
			return m_pdata[i * m_step];
		}

	private:
		const T* m_pdata;
		index_t m_step;
	};


	// single_reader

	template<typename T>
	class single_reader<T, scalar_> : public scalar_vec_accessor_base
	{
	public:
		typedef T scalar_type;

		LMAT_ENSURE_INLINE
		explicit single_reader(const T& v)
		: m_val(v) { }

		LMAT_ENSURE_INLINE
		T scalar(index_t ) const
		{
			return m_val;
		}

	private:
		const T m_val;
	};


	template<typename T, typename Kind>
	class single_reader<T, simd_<Kind> > : public simd_vec_accessor_base
	{
	public:
		typedef T scalar_type;
		typedef Kind simd_kind;
		typedef simd_pack<T, Kind> pack_type;

		LMAT_ENSURE_INLINE
		explicit single_reader(const T& v)
		: m_pack(v), m_val(v) { }

		LMAT_ENSURE_INLINE
		T scalar(index_t ) const
		{
			return m_val;
		}

		LMAT_ENSURE_INLINE
		pack_type pack(index_t ) const
		{
			return m_pack;
		}

	private:
		pack_type m_pack;
		T m_val;
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
			typedef typename meta::select_<
					meta::is_contiguous<Mat>, contvec_reader_map<Mat, U>,
					meta::is_col<Mat>, stepcol_reader_map<Mat, U>,
					meta::is_row<Mat>, steprow_reader_map<Mat, U>,
					meta::otherwise_, invalid_reader_map>::type internal_map_;

			typedef typename meta::if_<meta::is_regular_mat<Mat>,
					internal_map_, invalid_reader_map>::type internal_map;

			typedef typename internal_map::type type;

			LMAT_ENSURE_INLINE
			static type get(const Mat& mat)
			{
				return internal_map::get(mat);
			}
		};
	}


	template<class Mat, typename U>
	LMAT_ENSURE_INLINE
	inline typename internal::vec_reader_map<Mat, U>::type
	make_vec_accessor(U, const arg_wrap<Mat, atags::in>& wrap)
	{
		return internal::vec_reader_map<Mat, U>::get(wrap.arg());
	}


	template<typename T, typename U>
	LMAT_ENSURE_INLINE
	inline single_reader<T, U>
	make_vec_accessor(U, const arg_wrap<T, atags::single>& wrap)
	{
		return single_reader<T, U>(wrap.arg());
	}


	/********************************************
	 *
	 *  writer classes
	 *
	 ********************************************/

	// contvec_writer

	template<typename T>
	class contvec_writer<T, scalar_> : public scalar_vec_accessor_base
	{
	public:
		LMAT_ENSURE_INLINE
		explicit contvec_writer(T* p) : m_pdata(p) { }

		LMAT_ENSURE_INLINE
		T& scalar(index_t ) const
		{
			return m_temp;
		}

		LMAT_ENSURE_INLINE
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
	class contvec_writer<T, simd_<Kind> > : public simd_vec_accessor_base
	{
	public:
		typedef T scalar_type;
		typedef Kind simd_kind;
		typedef simd_pack<T, Kind> pack_type;

		LMAT_ENSURE_INLINE
		explicit contvec_writer(T* p) : m_pdata(p) { }

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
			m_ptemp.store_u(m_pdata + i);
			return nil_t();
		}

	private:
		mutable pack_type m_ptemp;
		mutable T m_stemp;
		T* m_pdata;
	};

	// stepvec_writer

	template<typename T>
	class stepvec_writer<T, scalar_> : public scalar_vec_accessor_base
	{
	public:
		typedef T scalar_type;

		LMAT_ENSURE_INLINE
		explicit stepvec_writer(T* p, index_t step)
		: m_pdata(p), m_step(step) { }

		LMAT_ENSURE_INLINE
		T& scalar(index_t ) const
		{
			return m_temp;
		}

		LMAT_ENSURE_INLINE
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
			typedef typename meta::select_<
					meta::is_contiguous<Mat>, contvec_writer_map<Mat, U>,
					meta::is_col<Mat>, stepcol_writer_map<Mat, U>,
					meta::is_row<Mat>, steprow_writer_map<Mat, U>,
					meta::otherwise_, invalid_writer_map>::type internal_map;

			typedef typename internal_map::type type;

			LMAT_ENSURE_INLINE
			static type get(Mat& mat)
			{
				return internal_map::get(mat);
			}
		};
	}


	template<class Mat, typename U>
	LMAT_ENSURE_INLINE
	inline typename internal::vec_writer_map<Mat, U>::type
	make_vec_accessor(U, const arg_wrap<Mat, atags::out>& wrap)
	{
		return internal::vec_writer_map<Mat, U>::get(wrap.arg());
	}


	/********************************************
	 *
	 *  updater classes
	 *
	 ********************************************/

	// contvec_updater

	template<typename T>
	class contvec_updater<T, scalar_> : public scalar_vec_accessor_base
	{
	public:
		typedef T scalar_type;

		LMAT_ENSURE_INLINE
		explicit contvec_updater(T* p) : m_pdata(p) { }

		LMAT_ENSURE_INLINE
		T& scalar(index_t i) const
		{
			return m_temp = m_pdata[i];
		}

		LMAT_ENSURE_INLINE
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
	class contvec_updater<T, simd_<Kind> > : public simd_vec_accessor_base
	{
	public:
		typedef T scalar_type;
		typedef Kind simd_kind;
		typedef simd_pack<T, Kind> pack_type;

		LMAT_ENSURE_INLINE
		explicit contvec_updater(T* p) : m_pdata(p) { }

		LMAT_ENSURE_INLINE
		T& scalar(index_t i) const
		{
			return m_stemp = m_pdata[i];
		}

		LMAT_ENSURE_INLINE
		pack_type& pack(index_t i) const
		{
			m_ptemp.load_u(m_pdata + i);
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
			m_ptemp.store_u(m_pdata + i);
			return nil_t();
		}

	private:
		mutable pack_type m_ptemp;
		mutable T m_stemp;
		T* m_pdata;
	};

	// stepvec_updater

	template<typename T>
	class stepvec_updater<T, scalar_> : public scalar_vec_accessor_base
	{
	public:
		typedef T scalar_type;

		LMAT_ENSURE_INLINE
		explicit stepvec_updater(T* p, index_t step)
		: m_pdata(p), m_step(step) { }

		LMAT_ENSURE_INLINE
		T& scalar(index_t i) const
		{
			return m_temp = m_pdata[i * m_step];
		}

		LMAT_ENSURE_INLINE
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
	 *  updater maps
	 *
	 ********************************************/

	namespace internal
	{
		template<class Mat, typename U>
		struct contvec_updater_map
		{
			typedef typename matrix_traits<Mat>::value_type T;
			typedef contvec_updater<T, U> type;

			LMAT_ENSURE_INLINE
			static type get(Mat& mat)
			{
				return type(mat.ptr_data());
			}
		};

		template<class Mat, typename U>
		struct stepcol_updater_map
		{
			typedef typename matrix_traits<Mat>::value_type T;
			typedef stepvec_updater<T, U> type;

			LMAT_ENSURE_INLINE
			static type get(Mat& mat)
			{
				return type(mat.ptr_data(), mat.row_stride());
			}
		};

		template<class Mat, typename U>
		struct steprow_updater_map
		{
			typedef typename matrix_traits<Mat>::value_type T;
			typedef stepvec_updater<T, U> type;

			LMAT_ENSURE_INLINE
			static type get(Mat& mat)
			{
				return type(mat.ptr_data(), mat.col_stride());
			}
		};

		struct invalid_updater_map { };

		template<class Mat, typename U>
		struct vec_updater_map
		{
			typedef typename meta::select_<
					meta::is_contiguous<Mat>, contvec_updater_map<Mat, U>,
					meta::is_col<Mat>, stepcol_updater_map<Mat, U>,
					meta::is_row<Mat>, steprow_updater_map<Mat, U>,
					meta::otherwise_, invalid_updater_map>::type internal_map;

			typedef typename internal_map::type type;

			LMAT_ENSURE_INLINE
			static type get(Mat& mat)
			{
				return internal_map::get(mat);
			}
		};
	}


	template<class Mat, typename U>
	LMAT_ENSURE_INLINE
	inline typename internal::vec_updater_map<Mat, U>::type
	make_vec_accessor(U, const arg_wrap<Mat, atags::in_out>& wrap)
	{
		return internal::vec_updater_map<Mat, U>::get(wrap.arg());
	}


	/********************************************
	 *
	 *  accumulators
	 *
	 ********************************************/

	template<typename T>
	class sum_accumulator<T, scalar_> : public scalar_vec_accessor_base
	{
	public:
		LMAT_ENSURE_INLINE
		sum_accumulator(T& a) : m_val(0), m_p(&a) { }

		LMAT_ENSURE_INLINE
		T& scalar(index_t ) const
		{
			return m_val;
		}

		LMAT_ENSURE_INLINE
		nil_t finalize() const
		{
			*m_p += m_val;
			return nil_t();
		}

	private:
		mutable T m_val;
		T *m_p;
	};

	template<typename T, typename Kind>
	class sum_accumulator<T, simd_<Kind> > : public simd_vec_accessor_base
	{
	public:
		typedef simd_pack<T, Kind> pack_type;

		LMAT_ENSURE_INLINE
		sum_accumulator(T& a) : m_val(0), m_p(&a) { }

		LMAT_ENSURE_INLINE
		T& scalar(index_t ) const
		{
			return m_val;
		}

		LMAT_ENSURE_INLINE
		pack_type& pack(index_t ) const
		{
			return m_pack;
		}

		LMAT_ENSURE_INLINE
		nil_t begin_packs() const
		{
			m_pack = pack_type::zeros();
			return nil_t();
		}

		LMAT_ENSURE_INLINE
		nil_t end_packs() const
		{
			m_val += sum(m_pack);
			return nil_t();
		}

		LMAT_ENSURE_INLINE
		nil_t finalize() const
		{
			*m_p += m_val;
			return nil_t();
		}

	private:
		mutable pack_type m_pack;
		mutable T m_val;
		T *m_p;
	};


	template<typename T>
	class max_accumulator<T, scalar_> : public scalar_vec_accessor_base
	{
	public:
		LMAT_ENSURE_INLINE
		max_accumulator(T& a)
		: m_val(-std::numeric_limits<T>::infinity()), m_p(&a) { }

		LMAT_ENSURE_INLINE
		T& scalar(index_t ) const
		{
			return m_val;
		}

		LMAT_ENSURE_INLINE
		nil_t finalize() const
		{
			*m_p = math::max(*m_p, m_val);
			return nil_t();
		}

	private:
		mutable T m_val;
		T *m_p;
	};

	template<typename T, typename Kind>
	class max_accumulator<T, simd_<Kind> > : public simd_vec_accessor_base
	{
	public:
		typedef simd_pack<T, Kind> pack_type;

		LMAT_ENSURE_INLINE
		max_accumulator(T& a)
		: m_val(-std::numeric_limits<T>::infinity()), m_p(&a) { }

		LMAT_ENSURE_INLINE
		T& scalar(index_t ) const
		{
			return m_val;
		}

		LMAT_ENSURE_INLINE
		pack_type& pack(index_t ) const
		{
			return m_pack;
		}

		LMAT_ENSURE_INLINE
		nil_t begin_packs() const
		{
			m_pack = pack_type::neg_inf();
			return nil_t();
		}

		LMAT_ENSURE_INLINE
		nil_t end_packs() const
		{
			m_val = math::max(m_val, maximum(m_pack));
			return nil_t();
		}

		LMAT_ENSURE_INLINE
		nil_t finalize() const
		{
			*m_p = math::max(*m_p, m_val);
			return nil_t();
		}

	private:
		mutable pack_type m_pack;
		mutable T m_val;
		T *m_p;
	};


	template<typename T>
	class min_accumulator<T, scalar_> : public scalar_vec_accessor_base
	{
	public:
		LMAT_ENSURE_INLINE
		min_accumulator(T& a)
		: m_val(std::numeric_limits<T>::infinity()), m_p(&a) { }

		LMAT_ENSURE_INLINE
		T& scalar(index_t ) const
		{
			return m_val;
		}

		LMAT_ENSURE_INLINE
		nil_t finalize() const
		{
			*m_p = math::min(*m_p, m_val);
			return nil_t();
		}

	private:
		mutable T m_val;
		T *m_p;
	};

	template<typename T, typename Kind>
	class min_accumulator<T, simd_<Kind> > : public simd_vec_accessor_base
	{
	public:
		typedef simd_pack<T, Kind> pack_type;

		LMAT_ENSURE_INLINE
		min_accumulator(T& a)
		: m_val(std::numeric_limits<T>::infinity()), m_p(&a) { }

		LMAT_ENSURE_INLINE
		T& scalar(index_t ) const
		{
			return m_val;
		}

		LMAT_ENSURE_INLINE
		pack_type& pack(index_t ) const
		{
			return m_pack;
		}

		LMAT_ENSURE_INLINE
		nil_t begin_packs() const
		{
			m_pack = pack_type::inf();
			return nil_t();
		}

		LMAT_ENSURE_INLINE
		nil_t end_packs() const
		{
			m_val = math::min(m_val, minimum(m_pack));
			return nil_t();
		}

		LMAT_ENSURE_INLINE
		nil_t finalize() const
		{
			*m_p = math::min(*m_p, m_val);
			return nil_t();
		}

	private:
		mutable pack_type m_pack;
		mutable T m_val;
		T *m_p;
	};


	/********************************************
	 *
	 *  accumulator maps
	 *
	 ********************************************/

	template<typename T, typename U>
	LMAT_ENSURE_INLINE
	inline sum_accumulator<T, U>
	make_vec_accessor(U, const arg_wrap<T, atags::sum>& wrap)
	{
		return sum_accumulator<T, U>(wrap.arg());
	}

	template<typename T, typename U>
	LMAT_ENSURE_INLINE
	inline max_accumulator<T, U>
	make_vec_accessor(U, const arg_wrap<T, atags::max>& wrap)
	{
		return max_accumulator<T, U>(wrap.arg());
	}

	template<typename T, typename U>
	LMAT_ENSURE_INLINE
	inline min_accumulator<T, U>
	make_vec_accessor(U, const arg_wrap<T, atags::min>& wrap)
	{
		return min_accumulator<T, U>(wrap.arg());
	}


}

#endif /* VEC_ACCESSORS_H_ */
