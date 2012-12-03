/**
 * @file test_base.h
 *
 * @brief The basis for testing
 *
 * @author Dahua Lin
 */

#ifdef _MSC_VER
#pragma once
#endif

#ifndef TEST_BASE_H_
#define TEST_BASE_H_

#include <light_test/tests.h>
#include <string>
#include <sstream>

namespace lmat { namespace test {

	class simple_case : public ltest::test_case
	{
		const char *m_name;

	public:
		simple_case(const char *nam)
		: m_name(nam)
		{
		}

		virtual ~simple_case() { }

		const char *name() const
		{
			return m_name;
		}
	};


	template<int N>
	class N_case : public ltest::test_case
	{
		std::string m_name;

	public:
		N_case(const char *nam)
		{
			std::stringstream ss;
			ss << nam << " [" << N << "]";
			m_name = ss.str();
		}

		virtual ~N_case() { }

		const char *name() const
		{
			return m_name.c_str();
		}
	};

	template<int M, int N>
	class MN_case : public ltest::test_case
	{
		std::string m_name;

	public:
		MN_case(const char *nam)
		{
			std::stringstream ss;
			ss << nam << " [" << M << " x " << N << "]";
			m_name = ss.str();
		}

		virtual ~MN_case() { }

		const char *name() const
		{
			return m_name.c_str();
		}
	};

} }


#define ASSERT_CT_VALUE( T, V ) \
	if (!((T::value) == (V))) throw ::ltest::assertion_failure(__FILE__, __LINE__, #T "::value == " #V)

#define ASSERT_CT_TYPE( T, R ) \
	if (!std::is_same<T, R>::value) throw ::ltest::assertion_failure(__FILE__, __LINE__, #T " is " #R)


#define BEGIN_MAIN_SUITE \
	ltest::test_suite lmat_main_suite( "Main" ); \
	void lmat_add_test_packs() {

#define END_MAIN_SUITE }

#define ADD_TPACK( pname ) lmat_main_suite.add( create_tpack_##pname() );


#define TCASE_CLASS( pname, tname ) pname##_##tname##_tests

#define SIMPLE_CASE( pname, tname ) \
	class TCASE_CLASS(pname, tname) : public lmat::test::simple_case { \
	public: \
		TCASE_CLASS(pname, tname)() : lmat::test::simple_case( #tname ) { } \
		virtual ~TCASE_CLASS(pname, tname)() { } \
		virtual void run(); \
	}; \
	void TCASE_CLASS(pname, tname)::run()

#define ADD_SIMPLE_CASE( pname, tname ) \
		tpack->add( new TCASE_CLASS( pname, tname )()  );

#define N_CASE( pname, tname ) \
	template<int N> \
	class TCASE_CLASS(pname, tname) : public lmat::test::N_case<N> { \
	public: \
		TCASE_CLASS(pname, tname)() : lmat::test::N_case<N>( #tname ) { } \
		virtual ~TCASE_CLASS(pname, tname)() { } \
		virtual void run(); \
	}; \
	template<int N> \
	void TCASE_CLASS(pname, tname)<N>::run()

#define ADD_N_CASE( pname, tname, n ) \
		tpack->add( new TCASE_CLASS( pname, tname )<n>()  );


#define MN_CASE( pname, tname ) \
	template<int M, int N> \
	class TCASE_CLASS(pname, tname) : public lmat::test::MN_case<M, N> { \
	public: \
		TCASE_CLASS(pname, tname)() : lmat::test::MN_case<M, N>( #tname ) { } \
		virtual ~TCASE_CLASS(pname, tname)() { } \
		virtual void run(); \
	}; \
	template<int M, int N> \
	void TCASE_CLASS(pname, tname)<M, N>::run()

#define ADD_MN_CASE( pname, tname, m, n ) \
		tpack->add( new TCASE_CLASS( pname, tname )<m,n>()  );

#define ADD_MN_CASE_3X3( pname, tname, m, n ) \
		ADD_MN_CASE( pname, tname, 0, 0 ) \
		ADD_MN_CASE( pname, tname, 0, 1 ) \
		ADD_MN_CASE( pname, tname, 0, n ) \
		ADD_MN_CASE( pname, tname, 1, 0 ) \
		ADD_MN_CASE( pname, tname, 1, 1 ) \
		ADD_MN_CASE( pname, tname, 1, n ) \
		ADD_MN_CASE( pname, tname, m, 0 ) \
		ADD_MN_CASE( pname, tname, m, 1 ) \
		ADD_MN_CASE( pname, tname, m, n )


#define BEGIN_TPACK( pname ) \
	ltest::test_pack* create_tpack_##pname() { \
		ltest::test_pack *tpack = new ltest::test_pack( #pname );

#define END_TPACK return tpack; }




#endif /* TEST_BASE_H_ */
