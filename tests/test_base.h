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

	/********************************************
	 *
	 *  type names
	 *
	 ********************************************/

	template<typename T>
	struct type_name;

	template<> struct type_name<float>
	{ static const char *get() { return "f32"; } };

	template<> struct type_name<double>
	{ static const char *get() { return "f64"; } };

	template<> struct type_name<int64_t>
	{ static const char *get() { return "i64"; } };

	template<> struct type_name<uint64_t>
	{ static const char *get() { return "u64"; } };

	template<> struct type_name<int32_t>
	{ static const char *get() { return "i32"; } };

	template<> struct type_name<uint32_t>
	{ static const char *get() { return "u32"; } };

	template<> struct type_name<int16_t>
	{ static const char *get() { return "i16"; } };

	template<> struct type_name<uint16_t>
	{ static const char *get() { return "u16"; } };

	template<> struct type_name<int8_t>
	{ static const char *get() { return "i8"; } };

	template<> struct type_name<uint8_t>
	{ static const char *get() { return "u8"; } };


	/********************************************
	 *
	 *  case classes
	 *
	 ********************************************/

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


	template<typename T>
	class T_case : public ltest::test_case
	{
		std::string m_name;

	public:
		T_case(const char *nam)
		{
			std::stringstream ss;
			ss << nam << " [" << type_name<T>::get() << "]";
			m_name = ss.str();
		}

		virtual ~T_case() { }

		const char *name() const
		{
			return m_name.c_str();
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


	template<typename T, int N>
	class TN_case : public ltest::test_case
	{
		std::string m_name;

	public:
		TN_case(const char *nam)
		{
			std::stringstream ss;
			ss << nam << " [" << type_name<T>::get() << " x " << N << "]";
			m_name = ss.str();
		}

		virtual ~TN_case() { }

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


	template<typename T, int M, int N>
	class TMN_case : public ltest::test_case
	{
		std::string m_name;

	public:
		TMN_case(const char *nam)
		{
			std::stringstream ss;
			ss << nam << " [" << type_name<T>::get() << " x " << M << " x " << N << "]";
			m_name = ss.str();
		}

		virtual ~TMN_case() { }

		const char *name() const
		{
			return m_name.c_str();
		}
	};

} }


/********************************************
 *
 *  general constructs
 *
 ********************************************/

#define BEGIN_MAIN_SUITE \
	ltest::test_suite lmat_main_suite( "Main" ); \
	void lmat_add_test_packs() {

#define END_MAIN_SUITE }

#define ADD_TPACK( pname ) lmat_main_suite.add( create_tpack_##pname() );

#define TCASE_CLASS( pname, tname ) pname##_##tname##_tests

#define BEGIN_TPACK( pname ) \
	ltest::test_pack* create_tpack_##pname() { \
		ltest::test_pack *tpack = new ltest::test_pack( #pname );

#define END_TPACK return tpack; }

/********************************************
 *
 *  specific constructs
 *
 ********************************************/

// simple cases

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

// T cases

#define T_CASE( pname, tname ) \
	template<typename T> \
	class TCASE_CLASS(pname, tname) : public lmat::test::T_case<T> { \
	public: \
		TCASE_CLASS(pname, tname)() : lmat::test::T_case<T>( #tname ) { } \
		virtual ~TCASE_CLASS(pname, tname)() { } \
		virtual void run(); \
	}; \
	template<typename T> \
	void TCASE_CLASS(pname, tname)<T>::run()

#define ADD_T_CASE( pname, tname, ty ) \
	tpack->add( new TCASE_CLASS( pname, tname )<ty>() );

// N cases

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


// TN cases

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
	tpack->add( new TCASE_CLASS( pname, tname )<ty, n>() );


// MN cases

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


// TMN cases

#define TMN_CASE( pname, tname ) \
	template<typename T, int M, int N> \
	class TCASE_CLASS(pname, tname) : public lmat::test::TN_case<T, M, N> { \
	public: \
		TCASE_CLASS(pname, tname)() : lmat::test::TN_case<T, M, N>( #tname ) { } \
		virtual ~TCASE_CLASS(pname, tname)() { } \
		virtual void run(); \
	}; \
	template<typename T, int M, int N> \
	void TCASE_CLASS(pname, tname)<T, M, N>::run()

#define ADD_TMN_CASE( pname, tname, ty, m, n ) \
	tpack->add( new TCASE_CLASS( pname, tname )<ty, m, n>() );


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



#endif /* TEST_BASE_H_ */



