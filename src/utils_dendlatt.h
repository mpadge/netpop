// utils.h

#ifndef UTILSDL_H
#define UTILSDL_H

#include "utils.h"

#include "boost/multi_array.hpp"

const double minqc = 0.001; // NOTE: Current results used 0.02
typedef boost::multi_array<bool, 3> barr3; 

std::vector <int> randseq (int n, base_generator_type * generator);
void pseq (dvec * pvals, bool q);
ivec sort (ivec sortvec, int veclen);
void timeout (double tseconds);
void dumpseed ();

#endif
