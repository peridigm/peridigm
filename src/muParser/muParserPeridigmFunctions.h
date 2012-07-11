#pragma GCC system_header

#ifndef MU_PARSER_PERIDIGM_FUNCTIONS_H
#define MU_PARSER_PERIDIGM_FUNCTIONS_H

#include "muParser.h"
#include "Peridigm_RandomNumber.hpp"

namespace mu {
  static value_type Rnd(value_type magnitude) {
    value_type tmp;
    value_type maxval = std::numeric_limits<value_type>::max();
    tmp = (2.0*( (value_type)PeridigmNS::rand_num() / (value_type)maxval )) - 1.0;
    return ( magnitude*tmp );
  }
} 

#endif


