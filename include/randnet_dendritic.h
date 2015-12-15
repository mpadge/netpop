// randnet_dendritic.h

#include "utils.h"
#include "utils_dendlatt.h"

// There used to be 8 radii, with the last plist = 0.0, but that generates the
// same network every time, so that last one has been ditched.
const int nradii = 7, rlist[nradii] = {5,6,7,8,9,10,11};
const double plist[nradii] = {0.8, 0.4, 0.24, 0.155, 0.1, 0.06, 0.026};
/*
 * mincq stops q values from being below this.
 * plist are the branching probabilities that give the corresponding average radii.
 */

struct Parameters
{
    /* directed orders connectivites so cdown < cup, simulating directed flow
	 * sigma is environmental stochasticity.
     * N0 is the Allee effect parameter (= 0), and k0 is the mean value of k. */
	bool directed;
	int nnodes, nTrials, maxSims;
	double pbranch, r, N0, k0, alpha0, pscale, sigma, alphasd;
};

struct NetStructure
{
	int nextus1, nextus2, nextds, xpos, ypos;
	double q, cdown, cup1, cup2, N;
};

struct NetResults
{
	/* NetResults are stored for each of the 101 values of alpha0.
	 * edgedist holds counts of nodes with [1..2] edges.
	 * edgelength is total length of connected, unbranching
	 * edges extending from terminal nodes.*/
	int edgedist [3], radius, edgelength;
	double connectivity [101], diff_q [101], N_mn [101], N_sd [101];
};

struct NetResults_OneNet
{
	int edgedist [3], radius, edgelength;
};

NetResults do1trial (Parameters params, base_generator_type * generator, 
        barr3 * varcheck);
void randnet (NetStructure *NetPt, Parameters params, 
        base_generator_type * generator);
void fillqc (NetStructure* NetPt, Parameters params, 
        base_generator_type * generator, bool unif);
NetResults_OneNet netstats (NetStructure *NetPt, Parameters params);
dmat populateNet (dmat * pmat, NetStructure *NetPt, Parameters params, 
        base_generator_type * generator);
dmat populateNet_eq (dmat * pmat, Parameters params);
void makepmat (const NetStructure *NetPt, dmat * pmat);
imat makedmat (const NetStructure *NetPt, Parameters params);
