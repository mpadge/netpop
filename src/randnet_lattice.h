// randnet_lattice.h

#include "utils.h"
#include "utils_dendlatt.h"

const int nradii = 6, rlist[nradii] = {3,4,5,6,7,8};

struct Parameters
{
	/* N0 is the Allee effect parameter (= 0),
	 * and k0 is the mean value of k. */
	int nnodes[3], nTrials, maxSims;
	double r, N0, k0, alpha0, pscale, sigma, alphasd;
};

struct NetStructure
{
	int xpos, ypos, nextnode[8];
	double q, conn[8], N;
};

struct NetResults
{
	/* edgedist holds counts of nodes with [1..8] edges;
	 * edgelength is total length of connected, unbranching
	 * edges extending from terminal nodes.*/
	int edgedist [8], radius, edgelength, nloops;
	double connectivity [101], diff_q [101], N_mn [101], N_sd [101];
};

struct NetResults_OneNet 
{
	int edgedist [8], radius, edgelength, nloops;
};


NetResults do1trial (Parameters params, base_generator_type * generator, 
        barr3 * varcheck);
void randnet (NetStructure *NetPt, Parameters params, 
        base_generator_type * generator);
void fillqc (NetStructure* NetPt, Parameters params, 
        base_generator_type * generator, bool unif);
NetResults_OneNet netstats (NetStructure *NetPt, Parameters params);
int countloops (const NetStructure* NetPt, Parameters params);
dmat populateNet (dmat * pmat, NetStructure *NetPt, Parameters params, 
        base_generator_type * generator);
dmat populateNet_eq (dmat * pmat, Parameters params);
void makepmat (const NetStructure *NetPt, dmat * pmat);
imat makedmat (const NetStructure *NetPt, Parameters params);
