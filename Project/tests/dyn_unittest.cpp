#define BOOST_TEST_MODULE MultiShooting Dynamic 
#include <boost/lexical_cast.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/ini_parser.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/numeric/ublas/io.hpp>

#include <boost/test/included/unit_test.hpp>

#include "config.hpp"

#include "dyn.hpp"

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
	std::vector<string> strs;
        boost::split(strs,strIges,boost::is_any_of("(,"));
	
	*g = string_to_double(pt.get<string>("environment.g"));
	*m = string_to_double(pt.get<string>("environment.m"));
	*kT = string_to_double(pt.get<string>("environment.kT"));
	*kQ = string_to_double(pt.get<string>("environment.kQ"));
	*d = string_to_double(pt.get<string>("environment.d"));
	*IM = string_to_double(pt.get<string>("environment.IM"));

	Iges[0] =  string_to_double(strs[1]);
	Iges[1] =  string_to_double(strs[2]);
	Iges[2] =  string_to_double(strs[3].substr(0, strs[3].length()-1));
}
/*
BOOST_AUTO_TEST_CASE(TestFex)
{
	size_t n_intervals = 0;
	value_type g, m, kT, kQ, d, IM;
	value_type Iges[3] = {0, 0, 0};
	g = 0; m = 0; kT = 0; kQ = 0; d = 0; IM = 0;
	getInitialValues(&g, &m, &kT, &kQ, &d, &IM, Iges);
	cout << "---		--" << endl;
	cout << "Test: TestFex" << endl;
	cout << "Gravitation  g=" << g << endl;
	cout << "Masse        m=" << m << endl;
	cout << "            kT=" << kT <<endl;
	cout << "            kQ=" << kQ << endl;
	cout << "             d=" << d <<endl;
	cout << "            IM=" << IM <<endl;
	cout << "           Iges={" << Iges[0] << "," << Iges[1] << "," << Iges[2] <<"}"  << endl;


	const int mN = 13;
	const int mM = 17;

	value_type pd[mN * mM];

	const value_type r[3] = {0, 0, 0};
	const value_type q[4] = {0, 0, 0, 1};
	const value_type v[3] = {0, 0, 0};
	const value_type omega[3] ={0, 0, 0};

	dyn::TestFuncFunctionality(n_intervals, g, m, kT, kQ, d, IM, Iges);
}
*/

BOOST_AUTO_TEST_CASE(TestFex)
{
	size_t n_intervals = 0;
	value_type g, m, kT, kQ, d, IM;
	value_type Iges[3] = {0, 0, 0};
	g = 0; m = 0; kT = 0; kQ = 0; d = 0; IM = 0;
	getInitialValues(&g, &m, &kT, &kQ, &d, &IM, Iges);
	cout << "---		--" << endl;
	cout << "Test: TestFex" << endl;
	cout << "Gravitation  g=" << g << endl;
	cout << "Masse        m=" << m << endl;
	cout << "            kT=" << kT <<endl;
	cout << "            kQ=" << kQ << endl;
	cout << "             d=" << d <<endl;
	cout << "            IM=" << IM <<endl;
	cout << "           Iges={" << Iges[0] << "," << Iges[1] << "," << Iges[2] <<"}"  << endl;


	dyn::TestFuncFunctionality(n_intervals, g, m, kT, kQ, d, IM, Iges);
}


BOOST_AUTO_TEST_CASE(TestJacobi)
{
	size_t n_intervals = 20;
	value_type g, m, kT, kQ, d, IM;
	value_type Iges[3] = {0, 0, 0};
	g = 0; m = 0; kT = 0; kQ = 0; d = 0; IM = 0;
	getInitialValues(&g, &m, &kT, &kQ, &d, &IM, Iges);
	cout << "---		--" << endl;
	cout << "Test: Jacobi" << endl;
	cout << "Gravitation  g=" << g << endl;
	cout << "Masse        m=" << m << endl;
	cout << "            kT=" << kT <<endl;
	cout << "            kQ=" << kQ << endl;
	cout << "             d=" << d <<endl;
	cout << "            IM=" << IM <<endl;
	cout << "           Iges={" << Iges[0] << "," << Iges[1] << "," << Iges[2] <<"}"  << endl;

	dyn::TestJacobi(n_intervals, g, m, kT, kQ, d, IM, Iges);
}