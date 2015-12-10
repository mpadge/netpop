/*
 * diff_sim4.h
 *
 * Adapated from diff_sim4_basic.cc to include
 * randomly varying habitat qualities and 
 * connectivities.
 * 
 * Simulates the four basic network types:
 * 
 * A =	*---*---*---*		B =	*---*---*
 * 							    |
 * 							    *
 * 					          *
 * C =	*---*		         /|
 *		|   |		D = *---* |
 *		*---*		         \|
 *					          *
 *
 * Connectivities are readily analytically derivable,
 * and the sole point of this is to estimate average
 * nodal variance for each type. Calculations are
 * done in all cases over the full range of values of
 * connectivity (= alpha).
 */

#include "utils.h"

struct Parameters 
{
    int nTrials;
    double k0, k0sd, ksd, alpha0, alphasd, r;
};

struct Results 
{
    double connectivity, nmn_node [5], nsd_node [5], nmn_network, nsd_network, 
           cov [4] [4];
};

Results runPop (Parameters pars, int network_type, 
        base_generator_type * generator);

const int runin = 100;
const double minqc = 0.001;
// minqc is important to prevent dividing by sometimes *really* small numbers
