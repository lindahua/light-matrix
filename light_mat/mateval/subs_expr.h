/**
 * @file subs_expr.h
 *
 * Subscript matrix expression
 * 
 * @author Dahua Lin 
 */

#ifdef _MSC_VER
#pragma once
#endif

#ifndef LIGHTMAT_SUBS_EXPR_H_
#define LIGHTMAT_SUBS_EXPR_H_

#include <light_mat/matrix/matrix_base.h>
#include <light_mat/mateval/ewise_eval.h>

namespace lmat
{
	// forward declaration

	template<typename T, int M, int N> class inds_expr;
	template<typename T, int M, int N> class subs_i_expr;
	template<typename T, int M, int N> class subs_j_expr;

	/********************************************
	 *
	 *  Expression classes
	 *
	 ********************************************/

	template<typename T, int M, int N>
	struct matrix_traits<inds_expr<T, M, N> >
	{
		static const int num_dimensions = 2;
		static const int ct_num_rows = M;
		static const int ct_num_cols = N;

		static const bool is_readonly = true;

		typedef matrix_shape<M, N> shape_type;
		typedef T value_type;
		typedef cpu_domain domain;
	};

	template<typename T, int M, int N>
	struct matrix_traits<subs_i_expr<T, M, N> >
	{
		static const int num_dimensions = 2;
		static const int ct_num_rows = M;
		static const int ct_num_cols = N;

		static const bool is_readonly = true;

		typedef matrix_shape<M, N> shape_type;
		typedef T value_type;
		typedef cpu_domain domain;
	};

	template<typename T, int M, int N>
	struct matrix_traits<subs_j_expr<T, M, N> >
	{
		static const int num_dimensions = 2;
		static const int ct_num_rows = M;
		static const int ct_num_cols = N;

		static const bool is_readonly = true;

		typedef matrix_shape<M, N> shape_type;
		typedef T value_type;
		typedef cpu_domain domain;
	};

	// classes

	template<typename T, int M, int N>
	class inds_expr : public IEWiseMatrix<inds_expr<T, M, N>, T>
	{
	public:
		typedef matrix_shape<M, N> shape_type;

		LMAT_ENSURE_INLINE inds_expr(index_t m, index_t n)
		: m_shape(m, n) { }

		LMAT_ENSURE_INLINE explicit inds_expr(const shape_type& shape)
		: m_shape(shape) { }

		LMAT_ENSURE_INLINE index_t nelems() const
		{
			return m_shape.nelems();
		}

		LMAT_ENSURE_INLINE index_t nrows() const
		{
			return m_shape.nrows();
		}

		LMAT_ENSURE_INLINE index_t ncolumns() const
		{
			return m_shape.ncolumns();
		}

		LMAT_ENSURE_INLINE shape_type shape() const
		{
			return m_shape;
		}

	private:
		shape_type m_shape;
	};

	template<typename T, int M, int N>
	class subs_i_expr : public IEWiseMatrix<subs_i_expr<T, M, N>, T>
	{
	public:
		typedef matrix_shape<M, N> shape_type;

		LMAT_ENSURE_INLINE subs_i_expr(index_t m, index_t n)
		: m_shape(m, n) { }

		LMAT_ENSURE_INLINE explicit subs_i_expr(const shape_type& shape)
		: m_shape(shape) { }

		LMAT_ENSURE_INLINE index_t nelems() const
		{
			return m_shape.nelems();
		}

		LMAT_ENSURE_INLINE index_t nrows() const
		{
			return m_shape.nrows();
		}

		LMAT_ENSURE_INLINE index_t ncolumns() const
		{
			return m_shape.ncolumns();
		}

		LMAT_ENSURE_INLINE shape_type shape() const
		{
			return m_shape;
		}

	private:
		shape_type m_shape;
	};

	template<typename T, int M, int N>
	class subs_j_expr : public IEWiseMatrix<subs_j_expr<T, M, N>, T>
	{
	public:
		typedef matrix_shape<M, N> shape_type;

		LMAT_ENSURE_INLINE subs_j_expr(index_t m, index_t n)
		: m_shape(m, n) { }

		LMAT_ENSURE_INLINE explicit subs_j_expr(const shape_type& shape)
		: m_shape(shape) { }

		LMAT_ENSURE_INLINE index_t nelems() const
		{
			return m_shape.nelems();
		}

		LMAT_ENSURE_INLINE index_t nrows() const
		{
			return m_shape.nrows();
		}

		LMAT_ENSURE_INLINE index_t ncolumns() const
		{
			return m_shape.ncolumns();
		}

		LMAT_ENSURE_INLINE shape_type shape() const
		{
			return m_shape;
		}

	private:
		shape_type m_shape;
	};


	/********************************************
	 *
	 *  Expression construction
	 *
	 ********************************************/

	LMAT_ENSURE_INLINE
	inline inds_expr<index_t, 0, 0>
	inds(index_t m, index_t n)
	{
		return inds_expr<index_t, 0, 0>(m, n);
	}

	template<typename T>
	LMAT_ENSURE_INLINE
	inline inds_expr<T, 0, 0>
	inds(type_<T>, index_t m, index_t n)
	{
		return inds_expr<T, 0, 0>(m, n);
	}

	template<int M, int N>
	LMAT_ENSURE_INLINE
	inline inds_expr<index_t, M, N>
	inds(const matrix_shape<M, N>& shape)
	{
		return inds_expr<index_t, M, N>(shape);
	}

	template<typename T, int M, int N>
	LMAT_ENSURE_INLINE
	inline inds_expr<T, M, N>
	inds(type_<T>, const matrix_shape<M, N>& shape)
	{
		return inds_expr<T, M, N>(shape);
	}

	LMAT_ENSURE_INLINE
	inline subs_i_expr<index_t, 0, 0>
	subs_i(index_t m, index_t n)
	{
		return subs_i_expr<index_t, 0, 0>(m, n);
	}

	template<typename T>
	LMAT_ENSURE_INLINE
	inline subs_i_expr<T, 0, 0>
	subs_i(type_<T>, index_t m, index_t n)
	{
		return subs_i_expr<T, 0, 0>(m, n);
	}

	template<int M, int N>
	LMAT_ENSURE_INLINE
	inline subs_i_expr<index_t, M, N>
	subs_i(const matrix_shape<M, N>& shape)
	{
		return subs_i_expr<index_t, M, N>(shape);
	}

	template<typename T, int M, int N>
	LMAT_ENSURE_INLINE
	inline subs_i_expr<T, M, N>
	subs_i(type_<T>, const matrix_shape<M, N>& shape)
	{
		return subs_i_expr<T, M, N>(shape);
	}

	LMAT_ENSURE_INLINE
	inline subs_j_expr<index_t, 0, 0>
	subs_j(index_t m, index_t n)
	{
		return subs_j_expr<index_t, 0, 0>(m, n);
	}

	template<typename T>
	LMAT_ENSURE_INLINE
	inline subs_j_expr<T, 0, 0>
	subs_j(type_<T>, index_t m, index_t n)
	{
		return subs_j_expr<T, 0, 0>(m, n);
	}

	template<int M, int N>
	LMAT_ENSURE_INLINE
	inline subs_j_expr<index_t, M, N>
	subs_j(const matrix_shape<M, N>& shape)
	{
		return subs_j_expr<index_t, M, N>(shape);
	}

	template<typename T, int M, int N>
	LMAT_ENSURE_INLINE
	inline subs_j_expr<T, M, N>
	subs_j(type_<T>, const matrix_shape<M, N>& shape)
	{
		return subs_j_expr<T, M, N>(shape);
	}


	/********************************************
	 *
	 *  Accessor classes
	 *
	 ********************************************/

	namespace internal
	{
		template<typename T, typename U>
		class iota_vec_reader;

		template<typename T, typename U> class inds_multicol_reader;
		template<typename T, typename U> class subs_i_multicol_reader;
		template<typename T, typename U> class subs_j_multicol_reader;

		template<typename T>
		class iota_vec_reader<T, atags::scalar> : public scalar_vec_accessor_base
		{
		public:
			typedef T scalar_type;

			LMAT_ENSURE_INLINE
			explicit iota_vec_reader() : m_base(0) { }

			LMAT_ENSURE_INLINE
			explicit iota_vec_reader(T b) : m_base(b) { }

			LMAT_ENSURE_INLINE
			T scalar(index_t i) const
			{
				return m_base + static_cast<T>(i);
			}

		private:
			const T m_base;
		};


		template<typename T, typename U>
		class inds_multicol_reader : public multicol_accessor_base
		{
		public:
			typedef T scalar_type;
			typedef iota_vec_reader<T, U> col_accessor_type;

			LMAT_ENSURE_INLINE
			explicit inds_multicol_reader(index_t m) : m_nrows(m) { }

			LMAT_ENSURE_INLINE
			col_accessor_type col(index_t j) const
			{ return col_accessor_type((T)(m_nrows * j)); }

		private:
			const index_t m_nrows;
		};

		template<typename T, typename U>
		class subs_i_multicol_reader : public multicol_accessor_base
		{
		public:
			typedef T scalar_type;
			typedef iota_vec_reader<T, U> col_accessor_type;

			LMAT_ENSURE_INLINE
			explicit subs_i_multicol_reader() : m_colacc(0) { }

			LMAT_ENSURE_INLINE
			col_accessor_type col(index_t ) const { return m_colacc; }

		private:
			col_accessor_type m_colacc;
		};

		template<typename T, typename U>
		class subs_j_multicol_reader : public multicol_accessor_base
		{
		public:
			typedef T scalar_type;
			typedef single_reader<T, U> col_accessor_type;

			LMAT_ENSURE_INLINE
			explicit subs_j_multicol_reader() { }

			LMAT_ENSURE_INLINE
			col_accessor_type col(index_t j) const
			{ return single_reader<T, U>((T)j); }
		};


		template<typename T, int M, int N, typename U>
		struct vec_reader_map<inds_expr<T, M, N>, U>
		{
			typedef iota_vec_reader<T, U> type;

			LMAT_ENSURE_INLINE
			static type get(const inds_expr<T, M, N>&)
			{
				return iota_vec_reader<T, U>();
			}
		};

		template<typename T, int M, int N, typename U>
		struct multicol_reader_map<inds_expr<T, M, N>, U>
		{
			typedef inds_multicol_reader<T, U> type;

			LMAT_ENSURE_INLINE
			static type get(const inds_expr<T, M, N>& expr)
			{
				return type(expr.nrows());
			}
		};

		template<typename T, int M, int N, typename U>
		struct multicol_reader_map<subs_i_expr<T, M, N>, U>
		{
			typedef subs_i_multicol_reader<T, U> type;

			LMAT_ENSURE_INLINE
			static type get(const subs_i_expr<T, M, N>&)
			{
				return type();
			}
		};

		template<typename T, int M, int N, typename U>
		struct multicol_reader_map<subs_j_expr<T, M, N>, U>
		{
			typedef subs_j_multicol_reader<T, U> type;

			LMAT_ENSURE_INLINE
			static type get(const subs_j_expr<T, M, N>&)
			{
				return type();
			}
		};
	}


	/********************************************
	 *
	 *  Evaluation
	 *
	 ********************************************/

	template<typename T, int M, int N>
	struct supports_linear_macc<inds_expr<T, M, N> >
	{
		static const bool value = true;
	};

	template<typename T, int M, int N>
	struct supports_linear_macc<subs_i_expr<T, M, N> >
	{
		static const bool value = false;
	};

	template<typename T, int M, int N>
	struct supports_linear_macc<subs_j_expr<T, M, N> >
	{
		static const bool value = false;
	};

	template<typename VT, int M, int N, typename T, typename Kind, bool IsLinear>
	struct supports_simd<inds_expr<VT, M, N>, T, Kind, IsLinear>
	{
		static const bool value = false;
	};

	template<typename VT, int M, int N, typename T, typename Kind, bool IsLinear>
	struct supports_simd<subs_i_expr<VT, M, N>, T, Kind, IsLinear>
	{
		static const bool value = false;
	};

	template<typename VT, int M, int N, typename T, typename Kind, bool IsLinear>
	struct supports_simd<subs_j_expr<VT, M, N>, T, Kind, IsLinear>
	{
		static const bool value = false;
	};


	template<typename T, int M, int N, class DMat>
	LMAT_ENSURE_INLINE
	inline void evaluate(const inds_expr<T, M, N>& sexpr, IRegularMatrix<DMat, T>& dmat)
	{
		evaluate_by_map(sexpr, dmat);
	}

	template<typename T, int M, int N, class DMat>
	LMAT_ENSURE_INLINE
	inline void evaluate(const subs_i_expr<T, M, N>& sexpr, IRegularMatrix<DMat, T>& dmat)
	{
		evaluate_by_map(sexpr, dmat);
	}

	template<typename T, int M, int N, class DMat>
	LMAT_ENSURE_INLINE
	inline void evaluate(const subs_j_expr<T, M, N>& sexpr, IRegularMatrix<DMat, T>& dmat)
	{
		evaluate_by_map(sexpr, dmat);
	}


}

#endif 
