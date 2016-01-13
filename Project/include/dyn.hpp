#include "config.hpp"
#include "dyn_base.h"

#pragma once

class dyn : protected dyn_base
{
public:
	dyn_base(value_type g, value_type m, value_type kT, 
				 value_type kQ, value_type d, value_type IM,
				 value_type* Iges) : dyn_base(g, m, kT, 
				 kQ, d, IM,	Iges)
	{

	}

};