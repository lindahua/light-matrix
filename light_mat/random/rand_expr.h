/**
 * @file rand_expr.h
 *
 * @brief Random matrix expression
 *
 * @author Dahua Lin
 */

#ifdef _MSC_VER
#pragma once
#endif

#ifndef LIGHTMAT_RAND_EXPR_H_
#define LIGHTMAT_RAND_EXPR_H_

#include <light_mat/matrix/matrix_classes.h>
#include <light_mat/simd/simd.h>
#include <light_mat/mateval/ewise_eval.h>
#include <light_mat/random/rand_accessors.h>

#include <light_mat/random/uniform_int_distr.h>
#include <light_mat/random/uniform_real_distr.h>
#include <light_mat/random/exponential_distr.h>
#include <light_mat/random/normal_distr.h>
#include <light_mat/random/gamma_distr.h>

namespace lmat
{
	// forward declaration

	template<class Distr, class RStream, int CM=0, int CN=0> class rand_expr;

	/********************************************
	 *
	 *  random expression classes
	 *
	 ********************************************/

	template<class Distr, class RStream, int CM, int CN>
	struct matrix_traits<rand_expr<Distr, RStream, CM, CN> >
	: public matrix_xpr_traits_base<
	  typename Distr::result_type, CM, CN, cpu_domain> { };

	template<class Distr, class RStream, int CM, int CN>
	class rand_expr
	: public ewise_matrix_base<rand_expr<Distr, RStream, CM, CN> >
	{
		typedef ewise_matrix_base<rand_expr<Distr, RStream, CM, CN> > base_t;
		using typename base_t::shape_type;

	public:
		typedef Distr distribution_type;

		LMAT_ENSURE_INLINE
		rand_expr(const Distr& distr, RStream& rs, index_t m, index_t n)
		: base_t(m, n), m_rstream(rs), m_distr(distr) { }

		LMAT_ENSURE_INLINE
		rand_expr(const Distr& distr, RStream& rs, const shape_type& shape)
		: base_t(shape), m_rstream(rs), m_distr(distr) { }

	public:
		LMAT_ENSURE_INLINE
		RStream& stream() const
		{
			return m_rstream;
		}

		LMAT_ENSURE_INLINE
		const Distr& distr() const
		{
			return m_distr;
		}

	private:
		RStream& m_rstream;
		Distr m_distr;
	};


	/********************************************
	 *
	 *  Expression construction functions
	 *
	 ********************************************/

	template<class Distr, class RStream>
	LMAT_ENSURE_INLINE
	inline rand_expr<Distr, RStream> rand_mat(const Distr& distr, RStream& rs, index_t m, index_t n)
	{
		return rand_expr<Distr, RStream>(distr, rs, m, n);
	}

	template<class Distr, class RStream, int M, int N>
	LMAT_ENSURE_INLINE
	inline rand_expr<Distr, RStream, M, N> rand_mat(const Distr& distr, RStream& rs, const matrix_shape<M, N>& shape)
	{
		return rand_expr<Distr, RStream, M, N>(distr, rs, shape);
	}

	// randu

	template<class RStream>
	LMAT_ENSURE_INLINE
	inline rand_expr<random::std_uniform_real_distr<double>, RStream>
	randu(RStream& rs, index_t m, index_t n)
	{
		return rand_mat(random::std_uniform_real_distr<double>(), rs, m, n);
	}

	template<class RStream>
	LMAT_ENSURE_INLINE
	inline rand_expr<random::std_uniform_real_distr<float>, RStream>
	randuf(RStream& rs, index_t m, index_t n)
	{
		return rand_mat(random::std_uniform_real_distr<float>(), rs, m, n);
	}

	template<typename T, class RStream>
	LMAT_ENSURE_INLINE
	inline typename std::enable_if<std::is_floating_point<T>::value,
	rand_expr<random::uniform_real_distr<T>, RStream> >::type
	randu(RStream& rs, index_t m, index_t n, T a, T b)
	{
		return rand_mat(random::uniform_real_distr<T>(a, b), rs, m, n);
	}

	// randn

	template<class RStream>
	LMAT_ENSURE_INLINE
	inline rand_expr<random::std_normal_distr<double>, RStream>
	randn(RStream& rs, index_t m, index_t n)
	{
		return rand_mat(random::std_normal_distr<double>(), rs, m, n);
	}

	template<class RStream>
	LMAT_ENSURE_INLINE
	inline rand_expr<random::std_normal_distr<float>, RStream>
	randnf(RStream& rs, index_t m, index_t n)
	{
		return rand_mat(random::std_normal_distr<float>(), rs, m, n);
	}

	template<typename T, class RStream>
	LMAT_ENSURE_INLINE
	inline typename std::enable_if<std::is_floating_point<T>::value,
	rand_expr<random::normal_distr<T>, RStream> >::type
	randn(RStream& rs, index_t m, index_t n, T mu, T sigma)
	{
		return rand_mat(random::normal_distr<T>(mu, sigma), rs, m, n);
	}


	// rande

	template<class RStream>
	LMAT_ENSURE_INLINE
	inline rand_expr<random::std_exponential_distr<double>, RStream>
	rande(RStream& rs, index_t m, index_t n)
	{
		return rand_mat(random::std_exponential_distr<double>(), rs, m, n);
	}

	template<class RStream>
	LMAT_ENSURE_INLINE
	inline rand_expr<random::std_exponential_distr<float>, RStream>
	randef(RStream& rs, index_t m, index_t n)
	{
		return rand_mat(random::std_exponential_distr<float>(), rs, m, n);
	}

	template<typename T, class RStream>
	LMAT_ENSURE_INLINE
	inline typename std::enable_if<std::is_floating_point<T>::value,
	rand_expr<random::exponential_distr<T>, RStream> >::type
	rande(RStream& rs, index_t m, index_t n, T lambda)
	{
		return rand_mat(random::exponential_distr<T>(lambda), rs, m, n);
	}


	// randg

	template<typename T, class RStream>
	LMAT_ENSURE_INLINE
	inline typename std::enable_if<std::is_floating_point<T>::value,
	rand_expr<random::std_gamma_distr<T>, RStream> >::type
	randg(RStream& rs, index_t m, index_t n, T alpha)
	{
		return rand_mat(random::std_gamma_distr<T>(alpha), rs, m, n);
	}

	template<typename T, class RStream>
	LMAT_ENSURE_INLINE
	inline typename std::enable_if<std::is_floating_point<T>::value,
	rand_expr<random::gamma_distr<T>, RStream> >::type
	randg(RStream& rs, index_t m, index_t n, T alpha, T beta)
	{
		return rand_mat(random::gamma_distr<T>(alpha, beta), rs, m, n);
	}


	/********************************************
	 *
	 *  Accessor classes
	 *
	 ********************************************/

	namespace internal
	{
		template<class Distr, class RStream, int CM, int CN, typename U>
		struct vec_reader_map<rand_expr<Distr, RStream, CM, CN>, U>
		{
			typedef rand_expr<Distr, RStream, CM, CN> expr_type;
			typedef rand_vec_reader<RStream, Distr, U> type;

			LMAT_ENSURE_INLINE
			static type get(const expr_type& expr)
			{
				return type(expr.stream(), expr.distr());
			}
		};

		template<class Distr, class RStream, int CM, int CN, typename U>
		struct multicol_reader_map<rand_expr<Distr, RStream, CM, CN>, U>
		{
			typedef rand_expr<Distr, RStream, CM, CN> expr_type;
			typedef rand_multicol_reader<RStream, Distr, U> type;

			LMAT_ENSURE_INLINE
			static type get(const expr_type& expr)
			{
				return type(expr.stream(), expr.distr());
			}
		};
	}


	/********************************************
	 *
	 *  Evaluation
	 *
	 ********************************************/

	template<class Distr, class RStream, int CM, int CN>
	struct supports_linear_macc<rand_expr<Distr, RStream, CM, CN> >
	{
		static const bool value = true;
	};

	template<class Distr, class RStream, int CM, int CN, typename Kind, bool IsLinear>
	struct supports_simd<rand_expr<Distr, RStream, CM, CN>, Kind, IsLinear>
	{
		static const bool value = is_simdizable<Distr, Kind>::value;
	};

	template<class Distr, class RStream, int CM, int CN, class DMat>
	LMAT_ENSURE_INLINE
	inline void evaluate(const rand_expr<Distr, RStream, CM, CN>& sexpr,
			IRegularMatrix<DMat, typename Distr::result_type>& dmat)
	{
		evaluate_by_map(sexpr, dmat);
	}


}

#endif
