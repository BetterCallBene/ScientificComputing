#define BOOST_TEST_MODULE MultiShooting Dynamic 
#include <boost/test/included/unit_test.hpp>

int add(int i, int j){return i + j;}

BOOST_AUTO_TEST_CASE(Constructor)
{
	BOOST_CHECK(add(2, 2) == 4);
}

