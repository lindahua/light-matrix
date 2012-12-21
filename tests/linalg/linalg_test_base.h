/**
 * @file linalg_test_base.h
 *
 * The testing facilities for linear algebra unit testing
 *
 * @author Dahua Lin
 */

#ifdef _MSC_VER
#pragma once
#endif

#ifndef LINALG_TEST_BASE_H_
#define LINALG_TEST_BASE_H_

#include "test_base.h"
#include <light_mat/linalg/linalg_base.h>
#include <cstdlib>

namespace lmat { namespace test {

	template<typename T> struct value_type_name;

	template<> struct value_type_name<float>
	{
		static const char *get() { return "f32"; }
	};

	template<> struct value_type_name<double>
	{
		static const char *get() { return "f64"; }
	};


	template<typename T, int N>
	class TN_case : public ltest::test_case
	{
		std::string m_name;

	public:
		TN_case(const char *nam)
		{
			std::stringstream ss;
			ss << nam
				<< " [" << N << " : " << value_type_name<T>::get() << "]";
			m_name = ss.str();
		}

		virtual ~TN_case() { }

		const char *name() const
		{
			return m_name.c_str();
		}
	};

	template<typename T, int M, int N>
	class TMN_case : public ltest::test_case
	{
		std::string m_name;

	public:
		TMN_case(const char *nam)
		{
			std::stringstream ss;
			ss << nam
				<< " [" << M << " x " << N << " : " << value_type_name<T>::get() << "]";
			m_name = ss.str();
		}

		virtual ~TMN_case() { }

		const char *name() const
		{
			return m_name.c_str();
		}
	};


	template<typename T, class Mat>
	void fill_rand(IDenseMatrix<Mat, T>& A, const T a, const T b)
	{
		const index_t m = A.nrows();
		const index_t n = A.ncolumns();

		for (index_t j = 0; j < n; ++j)
		{
			for (index_t i = 0; i < m; ++i)
			{
				double u = double(std::rand()) / double(RAND_MAX);
				A(i, j) = a + T(u) * (b - a);
			}
		}
	}


	template<typename T, class MatA, class MatB, class MatC>
	void naive_mtimes(
			const IDenseMatrix<MatA, T>& A,
			const IDenseMatrix<MatB, T>& B,
			IDenseMatrix<MatC, T>& C)
	{
		const index_t m = A.nrows();
		const index_t k = A.ncolumns();
		const index_t n = B.ncolumns();

		for (index_t j = 0; j < n; ++j)
		{
			for (index_t i = 0; i < m; ++i)
			{
				double s = 0;
				for (index_t l = 0; l < k; ++l)
				{
					s += double( A(i, l) * B(l, j) );
				}
				C(i, j) = T(s);
			}
		}
	}


} }


#define TN_CASE( pname, tname ) \
	template<typename T, int N> \
	class TCASE_CLASS(pname, tname) : public lmat::test::TN_case<T, N> { \
	public: \
		TCASE_CLASS(pname, tname)() : lmat::test::TN_case<T, N>( #tname ) { } \
		virtual ~TCASE_CLASS(pname, tname)() { } \
		virtual void run(); \
	}; \
	template<typename T, int N> \
	void TCASE_CLASS(pname, tname)<T, N>::run()

#define ADD_TN_CASE( pname, tname, ty, n ) \
		tpack->add( new TCASE_CLASS( pname, tname )<ty, n>()  );


#define TMN_CASE( pname, tname ) \
	template<typename T, int M, int N> \
	class TCASE_CLASS(pname, tname) : public lmat::test::TMN_case<T, M, N> { \
	public: \
		TCASE_CLASS(pname, tname)() : lmat::test::TMN_case<T, M, N>( #tname ) { } \
		virtual ~TCASE_CLASS(pname, tname)() { } \
		virtual void run(); \
	}; \
	template<typename T, int M, int N> \
	void TCASE_CLASS(pname, tname)<T, M, N>::run()

#define ADD_TMN_CASE( pname, tname, ty, m, n ) \
		tpack->add( new TCASE_CLASS( pname, tname )<ty, m,n>()  );

#define ADD_TMN_CASE_3X3( pname, tname, ty, m, n ) \
		ADD_TMN_CASE( pname, tname, ty, 0, 0 ) \
		ADD_TMN_CASE( pname, tname, ty, 0, 1 ) \
		ADD_TMN_CASE( pname, tname, ty, 0, n ) \
		ADD_TMN_CASE( pname, tname, ty, 1, 0 ) \
		ADD_TMN_CASE( pname, tname, ty, 1, 1 ) \
		ADD_TMN_CASE( pname, tname, ty, 1, n ) \
		ADD_TMN_CASE( pname, tname, ty, m, 0 ) \
		ADD_TMN_CASE( pname, tname, ty, m, 1 ) \
		ADD_TMN_CASE( pname, tname, ty, m, n )

#endif /* LINALG_TEST_BASE_H_ */



