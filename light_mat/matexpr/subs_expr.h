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

	template<typename T, int CM, int CN> class inds_expr;
	template<typename T, int CM, int CN> class subs_i_expr;
	template<typename T, int CM, int CN> class subs_j_expr;

	/********************************************
	 *
	 *  Expression classes
	 *
	 ********************************************/

	template<typename T, int CM, int CN>
	struct matrix_traits<inds_expr<T, CM, CN> >
	: public matrix_xpr_traits_base<T, CM, CN, cpu_domain> { };


	template<typename T, int CM, int CN>
	struct matrix_traits<subs_i_expr<T, CM, CN> >
	: public matrix_xpr_traits_base<T, CM, CN, cpu_domain> { };

	template<typename T, int CM, int CN>
	struct matrix_traits<subs_j_expr<T, CM, CN> >
	: public matrix_xpr_traits_base<T, CM, CN, cpu_domain> { };


	// classes

	template<typename T, int CM, int CN>
	class inds_expr
	: public ewise_matrix_base<inds_expr<T, CM, CN> >
	{
		typedef ewise_matrix_base<inds_expr<T, CM, CN> > base_t;
		using typename base_t::shape_type;

	public:
		LMAT_ENSURE_INLINE inds_expr(index_t m, index_t n)
		: base_t(m, n) { }

		LMAT_ENSURE_INLINE explicit inds_expr(const shape_type& shape)
		: base_t(shape) { }
	};

	template<typename T, int CM, int CN>
	class subs_i_expr
	: public ewise_matrix_base<subs_i_expr<T, CM, CN> >
	{
		typedef ewise_matrix_base<subs_i_expr<T, CM, CN> > base_t;
		using typename base_t::shape_type;

	public:
		LMAT_ENSURE_INLINE subs_i_expr(index_t m, index_t n)
		: base_t(m, n) { }

		LMAT_ENSURE_INLINE explicit subs_i_expr(const shape_type& shape)
		: base_t(shape) { }
	};

	template<typename T, int CM, int CN>
	class subs_j_expr
	: public ewise_matrix_base<subs_j_expr<T, CM, CN> >
	{
		typedef ewise_matrix_base<subs_j_expr<T, CM, CN> > base_t;
		using typename base_t::shape_type;

	public:
		LMAT_ENSURE_INLINE subs_j_expr(index_t m, index_t n)
		: base_t(m, n) { }

		LMAT_ENSURE_INLINE explicit subs_j_expr(const shape_type& shape)
		: base_t(shape) { }
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

	template<int CM, int CN>
	LMAT_ENSURE_INLINE
	inline inds_expr<index_t, CM, CN>
	inds(const matrix_shape<CM, CN>& shape)
	{
		return inds_expr<index_t, CM, CN>(shape);
	}

	template<typename T, int CM, int CN>
	LMAT_ENSURE_INLINE
	inline inds_expr<T, CM, CN>
	inds(type_<T>, const matrix_shape<CM, CN>& shape)
	{
		return inds_expr<T, CM, CN>(shape);
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

	template<int CM, int CN>
	LMAT_ENSURE_INLINE
	inline subs_i_expr<index_t, CM, CN>
	subs_i(const matrix_shape<CM, CN>& shape)
	{
		return subs_i_expr<index_t, CM, CN>(shape);
	}

	template<typename T, int CM, int CN>
	LMAT_ENSURE_INLINE
	inline subs_i_expr<T, CM, CN>
	subs_i(type_<T>, const matrix_shape<CM, CN>& shape)
	{
		return subs_i_expr<T, CM, CN>(shape);
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

	template<int CM, int CN>
	LMAT_ENSURE_INLINE
	inline subs_j_expr<index_t, CM, CN>
	subs_j(const matrix_shape<CM, CN>& shape)
	{
		return subs_j_expr<index_t, CM, CN>(shape);
	}

	template<typename T, int CM, int CN>
	LMAT_ENSURE_INLINE
	inline subs_j_expr<T, CM, CN>
	subs_j(type_<T>, const matrix_shape<CM, CN>& shape)
	{
		return subs_j_expr<T, CM, CN>(shape);
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


		template<typename T, int CM, int CN, typename U>
		struct vec_reader_map<inds_expr<T, CM, CN>, U>
		{
			typedef iota_vec_reader<T, U> type;

			LMAT_ENSURE_INLINE
			static type get(const inds_expr<T, CM, CN>&)
			{
				return iota_vec_reader<T, U>();
			}
		};

		template<typename T, int CM, int CN, typename U>
		struct multicol_reader_map<inds_expr<T, CM, CN>, U>
		{
			typedef inds_multicol_reader<T, U> type;

			LMAT_ENSURE_INLINE
			static type get(const inds_expr<T, CM, CN>& expr)
			{
				return type(expr.nrows());
			}
		};

		template<typename T, int CM, int CN, typename U>
		struct multicol_reader_map<subs_i_expr<T, CM, CN>, U>
		{
			typedef subs_i_multicol_reader<T, U> type;

			LMAT_ENSURE_INLINE
			static type get(const subs_i_expr<T, CM, CN>&)
			{
				return type();
			}
		};

		template<typename T, int CM, int CN, typename U>
		struct multicol_reader_map<subs_j_expr<T, CM, CN>, U>
		{
			typedef subs_j_multicol_reader<T, U> type;

			LMAT_ENSURE_INLINE
			static type get(const subs_j_expr<T, CM, CN>&)
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

	template<typename T, int CM, int CN>
	struct supports_linear_macc<inds_expr<T, CM, CN> > : public meta::true_ { };

	template<typename T, int CM, int CN>
	struct supports_linear_macc<subs_i_expr<T, CM, CN> > : public meta::false_ { };

	template<typename T, int CM, int CN>
	struct supports_linear_macc<subs_j_expr<T, CM, CN> > : public meta::false_ { };

	template<typename VT, int CM, int CN, typename Kind, bool IsLinear>
	struct supports_simd<inds_expr<VT, CM, CN>, Kind, IsLinear> : public meta::false_ { };

	template<typename VT, int CM, int CN, typename Kind, bool IsLinear>
	struct supports_simd<subs_i_expr<VT, CM, CN>, Kind, IsLinear> : public meta::false_ { };

	template<typename VT, int CM, int CN, typename Kind, bool IsLinear>
	struct supports_simd<subs_j_expr<VT, CM, CN>, Kind, IsLinear> : public meta::false_ { };


	template<typename T, int CM, int CN, class DMat>
	LMAT_ENSURE_INLINE
	inline void evaluate(const inds_expr<T, CM, CN>& sexpr, IRegularMatrix<DMat, T>& dmat)
	{
		evaluate_by_map(sexpr, dmat);
	}

	template<typename T, int CM, int CN, class DMat>
	LMAT_ENSURE_INLINE
	inline void evaluate(const subs_i_expr<T, CM, CN>& sexpr, IRegularMatrix<DMat, T>& dmat)
	{
		evaluate_by_map(sexpr, dmat);
	}

	template<typename T, int CM, int CN, class DMat>
	LMAT_ENSURE_INLINE
	inline void evaluate(const subs_j_expr<T, CM, CN>& sexpr, IRegularMatrix<DMat, T>& dmat)
	{
		evaluate_by_map(sexpr, dmat);
	}


}

#endif 
