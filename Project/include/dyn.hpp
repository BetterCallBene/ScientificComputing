#include <thrust/device_vector.h>
#include <boost/numeric/odeint/external/thrust/thrust.hpp>

#include "config.hpp"
#include "sys_base_functor.hpp"

#pragma once
using namespace std;

class dyn
{
	value_type m_g;
	value_type m_m;
	value_type m_kT;
	value_type m_kQ;
	value_type m_d;
	value_type m_IM;
	value_type m_Iges[3];
	value_type m_u[3];
public:

	struct sys_functor : sys_base_functor
	{
		
		sys_functor(value_type g, value_type m, value_type kT,
			value_type kQ, value_type d, value_type IM, 
			value_type* Iges, value_type *u) : sys_base_functor(g, m, kT, kQ, d, IM, Iges, u)
		{
			
		}
		template <class Tuple>
		__host__ __device__ 
		void operator()(Tuple t)
		{
			const value_type r[3] = {thrust::get<0>(t), thrust::get<1>(t), thrust::get<2>(t)};
			const value_type q[4] = {thrust::get<3>(t), thrust::get<4>(t), thrust::get<5>(t), thrust::get<6>(t)};
			const value_type v[3] = {thrust::get<7>(t), thrust::get<8>(t), thrust::get<9>(t)};
			const value_type omega[3] = {thrust::get<10>(t), thrust::get<11>(t), thrust::get<12>(t)};
		}
	};
	dyn(value_type g, value_type m, value_type kT, 
				 value_type kQ, value_type d, value_type IM,
				 value_type* Iges, value_type* u) : m_g(g), m_m(m), m_kT(kT),
				 m_kQ(kQ), m_d(d), m_IM(IM)
	{
		this->m_Iges[0] = Iges[0];
		this->m_Iges[1] = Iges[1];
		this->m_Iges[2] = Iges[2];

		this->m_u[0] = u[0];
		this->m_u[1] = u[1];
		this->m_u[2] = u[2];
		this->m_u[3] = u[3];
	}



};
