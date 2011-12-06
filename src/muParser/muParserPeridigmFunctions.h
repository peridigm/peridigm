#ifndef MU_PARSER_PERIDIGM_FUNCTIONS_H
#define MU_PARSER_PERIDIGM_FUNCTIONS_H

#include "muParser.h"

namespace mu {
  static value_type Rnd(value_type magnitude) {
    value_type tmp;
    tmp = (2.0*( (value_type)rand() / (value_type)RAND_MAX )) - 1.0;
    return ( magnitude*tmp );
  }
} 

#endif


