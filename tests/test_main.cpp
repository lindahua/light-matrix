/**
 * @file test_main.cpp
 *
 * The main file for a testing program
 *
 * @author Dahua Lin
 */

#include <light_test/std_test_mon.h>

#ifdef _MSC_VER
#pragma warning(disable:4100)
#endif

using namespace ltest;

extern void lmat_add_test_packs();
extern ltest::test_suite lmat_main_suite;

int main(int argc, char *argv[])
{
	::lmat_add_test_packs();

	if (std_test_main(lmat_main_suite))
	{
		return 0;
	}
	else
	{
		return -1;
	}
}




