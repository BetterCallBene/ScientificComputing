#include <boost/numeric/ublas/vector.hpp>
#include "config.hpp"
#include "drand48_generator.hpp"
#pragma once
typedef boost::numeric::ublas::vector< value_type > vector_type;
class TestEnv
{
protected:
	static void generator(int nvec, vector_type& vec)
	{
		
		srand48_generator s48;
		for(int i = 0; i < nvec; i++) vec[i] = s48();
		//nrmQ(timepoints, vec);
	}
};