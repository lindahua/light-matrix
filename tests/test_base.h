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
 *  specific constructs
 *
 ********************************************/

// simple cases

#define SIMPLE_CASE( Name ) \
	class Name: public lmat::test::simple_case { \
	public: \
		Name() : lmat::test::simple_case( #Name ) { } \
		virtual ~Name() { } \
		virtual void run(); \
	}; \
	void Name::run()

#define ADD_SIMPLE_CASE( Name ) this->add( new Name() );

// T cases

#define T_CASE( Name ) \
	template<typename T> \
	class Name : public lmat::test::T_case<T> { \
	public: \
		Name() : lmat::test::T_case<T>( #Name ) { } \
		virtual ~Name() { } \
		virtual void run(); \
	}; \
	template<typename T> \
	void Name<T>::run()

#define ADD_T_CASE( Name, ty ) this->add( new Name<ty>());

#define ADD_T_CASE_FP( Name ) \
		ADD_T_CASE( Name, double ) \
		ADD_T_CASE( Name, float )

// N cases

#define N_CASE( Name ) \
	template<int N> \
	class Name : public lmat::test::N_case<N> { \
	public: \
		Name() : lmat::test::N_case<N>( #Name ) { } \
		virtual ~Name() { } \
		virtual void run(); \
	}; \
	template<int N> \
	void Name<N>::run()

#define ADD_N_CASE( Name, n ) this->add( new Name<n>() );

#define ADD_N_CASE_3( Name, n ) \
		ADD_N_CASE( Name, 0 ) \
		ADD_N_CASE( Name, 1 ) \
		ADD_N_CASE( Name, n )


// TN cases

#define TN_CASE( Name ) \
	template<typename T, int N> \
	class Name : public lmat::test::TN_case<T, N> { \
	public: \
		Name() : lmat::test::TN_case<T, N>( #Name ) { } \
		virtual ~Name() { } \
		virtual void run(); \
	}; \
	template<typename T, int N> \
	void Name<T, N>::run()

#define ADD_TN_CASE( Name, ty, n ) this->add( new Name<ty, n>() );

#define ADD_TN_CASE_3( Name, ty, n ) \
		ADD_TN_CASE( Name, ty, 0 ) \
		ADD_TN_CASE( Name, ty, 1 ) \
		ADD_TN_CASE( Name, ty, n )


// MN cases

#define MN_CASE( Name ) \
	template<int M, int N> \
	class Name : public lmat::test::MN_case<M, N> { \
	public: \
		Name() : lmat::test::MN_case<M, N>( #Name ) { } \
		virtual ~Name() { } \
		virtual void run(); \
	}; \
	template<int M, int N> \
	void Name<M, N>::run()

#define ADD_MN_CASE( Name, m, n ) this->add( new Name<m,n>()  );

#define ADD_MN_CASE_3X3( Name, m, n ) \
		ADD_MN_CASE( Name, 0, 0 ) \
		ADD_MN_CASE( Name, 0, 1 ) \
		ADD_MN_CASE( Name, 0, n ) \
		ADD_MN_CASE( Name, 1, 0 ) \
		ADD_MN_CASE( Name, 1, 1 ) \
		ADD_MN_CASE( Name, 1, n ) \
		ADD_MN_CASE( Name, m, 0 ) \
		ADD_MN_CASE( Name, m, 1 ) \
		ADD_MN_CASE( Name, m, n )


// TMN cases

#define TMN_CASE( Name ) \
	template<typename T, int M, int N> \
	class Name : public lmat::test::TMN_case<T, M, N> { \
	public: \
		Name() : lmat::test::TMN_case<T, M, N>( #Name ) { } \
		virtual ~Name() { } \
		virtual void run(); \
	}; \
	template<typename T, int M, int N> \
	void Name<T, M, N>::run()

#define ADD_TMN_CASE( Name, ty, m, n ) this->add( new Name<ty, m, n>() );

#define ADD_TMN_CASE_3X3( Name, ty, m, n ) \
		ADD_TMN_CASE( Name, ty, 0, 0 ) \
		ADD_TMN_CASE( Name, ty, 0, 1 ) \
		ADD_TMN_CASE( Name, ty, 0, n ) \
		ADD_TMN_CASE( Name, ty, 1, 0 ) \
		ADD_TMN_CASE( Name, ty, 1, 1 ) \
		ADD_TMN_CASE( Name, ty, 1, n ) \
		ADD_TMN_CASE( Name, ty, m, 0 ) \
		ADD_TMN_CASE( Name, ty, m, 1 ) \
		ADD_TMN_CASE( Name, ty, m, n )



#endif /* TEST_BASE_H_ */



