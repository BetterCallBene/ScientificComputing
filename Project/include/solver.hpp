
#include <boost/numeric/odeint.hpp>
#include "dyn.hpp"
#include "TestEnv.hpp"

#pragma once 

using namespace std;
using namespace boost::numeric::odeint;

class solver : public TestEnv
{

	struct quad_system : dyn
	{
		quad_system(size_t N, const vector_type &u, value_type g, value_type m, value_type kT, 
	 			 value_type kQ, value_type d, value_type IM,
	 			 value_type* Iges) : dyn(N, u, g, m, kT, kQ, d, IM, Iges)
	 	{

	 	}
	 	void operator()( const vector_type &x , vector_type &dxdt , double /* t */ )
 		{
 			fex(x, dxdt);
 		}

		
	};

	struct quad_system_jacobi : dyn
	{
		quad_system_jacobi(size_t N, const vector_type &u, value_type g, value_type m, value_type kT, 
				 value_type kQ, value_type d, value_type IM,
				 value_type* Iges) : dyn(N, u, g, m, kT, kQ, d, IM, Iges)
		{

		}

		void operator()( const vector_type &x, matrix_type &J , const double & /* t */ , vector_type &dfdt )
     	{
      		for(int i = 0; i < mN; i++)
      			dfdt[i] = 0;
      		jac(x, J);
      	} 
	};
protected:
	static void nrmQ(int ntimepoints, vector_type& vec)
	{
		for(int i =0; i < ntimepoints; i++ )
		{
			vector_type s 
				= subrange(vec, i * NEQ + nqFirst, i * NEQ + nqLast + 1); 
			value_type n2 = norm_2(s);
			for(int j = 0; j < nquad; j++)
			{
				int indx = i * NEQ + nqFirst + j;
				vec[indx] = vec[indx] / n2;
			}
		}
	}
public:
	void integrate(vector_type state0)
	{
		
	}
public:
	static void TestIntegratorRungeKutta(int nintervall, value_type g,
	value_type m,
	value_type kT,
	value_type kQ,
	value_type d,
	value_type IM,
	value_type* Iges
	)
	{
		size_t N; 
		cout << "Test: Integrator" << endl;

		int ntimepoints = (nintervall + 1);
		int nstate0 = NEQ * ntimepoints;
		int nallcontr = ncontr * ntimepoints;
		
		vector_type state0(nstate0);
		vector_type u0(nallcontr);

		generator(nstate0, state0);
		generator(nallcontr, u0);
		nrmQ(ntimepoints, state0);

		quad_system quad_sys(ntimepoints, u0, g, m, kT, kQ, d, IM, Iges);
		quad_system_jacobi quad_sys_jac(ntimepoints, u0, g, m, kT, kQ, d, IM, Iges);

		//typedef runge_kutta_dopri5< vector_type , vector_type , vector_type , vector_type > stepper_type;
		typedef runge_kutta_dopri5< vector_type , value_type , vector_type , value_type > stepper_type;
		 
		integrate_adaptive( make_controlled( 1.0e-6 , 1.0e-4 , stepper_type()) , quad_sys , state0 , 0.0 , 1.0 , 0.1 );
		cout << "YES" << endl;
	}

	static void TestIntegratorRosenbrock(int nintervall, value_type g,
	value_type m,
	value_type kT,
	value_type kQ,
	value_type d,
	value_type IM,
	value_type* Iges
	)
	{
		size_t N; 
		cout << "Test: Integrator" << endl;

		int ntimepoints = (nintervall + 1);
		int nstate0 = NEQ * ntimepoints;
		int nallcontr = ncontr * ntimepoints;
		
		vector_type state0(nstate0);
		vector_type u0(nallcontr);

		generator(nstate0, state0);
		generator(nallcontr, u0);

		nrmQ(ntimepoints, state0);
		
		quad_system quad_sys(ntimepoints, u0, g, m, kT, kQ, d, IM, Iges);
		quad_system_jacobi quad_sys_jac(ntimepoints, u0, g, m, kT, kQ, d, IM, Iges);

		//typedef runge_kutta_dopri5< vector_type , vector_type , vector_type , vector_type > stepper_type;
		typedef rosenbrock4< value_type > stepper_type;
		 
		integrate_adaptive( make_controlled( 1.0e-6 , 1.0e-4 , stepper_type()), make_pair(quad_sys, quad_sys_jac), state0 , 0.0 , 1.0 , 0.1 );
		cout << "YES" << endl;
	}
};