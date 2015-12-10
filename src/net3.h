/*
 * diff_sim3.h
 */

#include "utils.h"

struct Parameters 
{
    int nTrials;
    double k0, k0sd, ksd, alpha0, alphasd, r;
};

struct Results 
{
    double connectivity, nmn_node [4], nsd_node [4], nmn_network, nsd_network, 
           cov [3] [3];
};

Results runPop (Parameters pars, base_generator_type * generator);

const int runin = 100;
const double minqc = 0.001;
// minqc is important to prevent dividing by sometimes *really* small numbers
