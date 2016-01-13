#define BOOST_TEST_MODULE MyTest
#include <boost/test/unit_test.hpp>

int add(int i, int j){return i + j;}

BOOST_AUTO_TEST_CASE(my_test)
{
	BOOST_CHECK(add(2, 2) == 4);
}

