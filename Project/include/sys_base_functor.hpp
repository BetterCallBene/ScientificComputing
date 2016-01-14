
#include <cmath>
#include "config.hpp"
#pragma once

class sys_base_functor
{
//parameter (Iges = (/0.0093886,0.0093886,0.018406/))
//            parameter (g = 9.81, m = 1.022,  kT  = 1.5e1)
//            parameter (kQ = 3e-01, d = 0.22, IM = 4.4466e-06)
	private:
		value_type g;
		value_type m;
		value_type kT;
		value_type kQ;
		value_type d;
		value_type IM;
		value_type Iges[3];
		value_type u[4];
	public:
		sys_base_functor(value_type g, value_type m, value_type kT, 
				 value_type kQ, value_type d, value_type IM,
				 value_type* Iges, value_type* u)
		{
			this->g = g;
			this->m = m;
			this->kT = kT;
			this->kQ = kQ;
			this->d = d;
			this->IM = IM;
			this->Iges[0] = Iges[0];
			this->Iges[1] = Iges[1];
			this->Iges[2] = Iges[2];
			this->u[0] = u[0];
			this->u[1] = u[1];
			this->u[2] = u[2];
			this->u[3] = u[3];
		}


	protected:
		void func(const value_type* r, const value_type* q, const value_type* v, 
			  const value_type* omega, value_type* rhs)
		{
	
		}

		void jacobi(const value_type* r, const value_type* q, const value_type* v, 
		        const value_type* omega, value_type** pd)
		{
		}

		void hesse(const value_type* r, const value_type* q, const value_type* v, 
			const value_type* omega, value_type*** pdd)
		{
		}

};
