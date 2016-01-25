
#include <stdlib.h>
#pragma once
struct srand48_generator
{
    typedef float result_type;
    result_type operator()( void ) const { return (result_type) drand48(); }
    result_type min( void ) const { return 0.0; }
    result_type max( void ) const { return 1.0; }
};