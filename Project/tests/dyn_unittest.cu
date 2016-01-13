#define BOOST_TEST_MODULE MultiShooting Dynamic 
#include <boost/lexical_cast.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/ini_parser.hpp>

#include <boost/test/included/unit_test.hpp>

#include "config.hpp"
#include "dyn.hpp"
int add(int i, int j){return i + j;}


using namespace std;
using namespace boost::lexical_cast;
using namespace boost::property_tree;


double string_to_double(string value)
{
	return lexical_cast<double>(value);
}

void getInitialValues(value_type* g, value_type* m, value_type* kT, 
				 value_type* kQ, value_type* d, value_type* IM,
				 value_type* Iges)
{
	ptree pt;
	ini_parser::read_ini("rtopt.ini", pt);
	cout << pt.get<string>("environment.Iges");
}

BOOST_AUTO_TEST_CASE(Constructor)
{
	
}

