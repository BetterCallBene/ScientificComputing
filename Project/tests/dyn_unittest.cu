#define BOOST_TEST_MODULE MultiShooting Dynamic 
#include <boost/lexical_cast.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/ini_parser.hpp>
#include <boost/algorithm/string.hpp>

#include <boost/test/included/unit_test.hpp>

#include "config.hpp"
#include "dyn.hpp"
int add(int i, int j){return i + j;}


using namespace std;
using namespace boost::property_tree;


double string_to_double(string value)
{
	return boost::lexical_cast<double>(value);
}

void getInitialValues(value_type* g, value_type* m, value_type* kT, 
				 value_type* kQ, value_type* d, value_type* IM,
				 value_type* Iges)
{
	ptree pt;
	ini_parser::read_ini("rtopt.ini", pt);
  
	string strIges =  pt.get<string>("environment.Iges");
	vector<string> strs;
        boost::split(strs,strIges,boost::is_any_of("(,"));
	cout << strs[0];
}

BOOST_AUTO_TEST_CASE(Constructor)
{
	value_type g, m, kT, kQ, d, IM;
	value_type Iges[3] = {0, 0, 0};
	g = 0; m = 0; kT = 0; kQ = 0; d = 0; IM = 0;
	getInitialValues(&g, &m, &kT, &kQ, &d, &IM, Iges);	
}

