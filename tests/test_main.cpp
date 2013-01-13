/**
 * @file test_main.cpp
 *
 * The main file for a testing program
 *
 * @author Dahua Lin
 */

#include <light_test/auto_suite.h>
#include <light_test/std_test_mon.h>

#ifdef _MSC_VER
#pragma warning(disable:4100)
#endif

using namespace ltest;

int main(int argc, char *argv[])
{
	if (std_test_main(*auto_test_suite::main_suite()))
	{
		return 0;
	}
	else
	{
		return -1;
	}
}




