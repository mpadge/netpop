/***************************************************************************
 *  Project:    netpop
 *  File:       randnet_lattice.h
 *  Language:   C++
 *
 *  netpop is free software: you can redistribute it and/or modify it under the
 *  terms of the GNU General Public License as published by the Free Software
 *  Foundation, either version 3 of the License, or (at your option) any later
 *  version.
 *
 *  netpop is distributed in the hope that it will be useful, but WITHOUT ANY
 *  WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 *  FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 *  details.
 *
 *  You should have received a copy of the GNU General Public License along with
 *  NeutralClusters.  If not, see <http://www.gnu.org/licenses/>.
 *
 *  Copyright   Mark Padgham December 2015
 *  Author:     Mark Padgham
 *  E-Mail:     mark.padgham@email.com
 *
 *  Description:    Simulates populations reproducing according to cubic
 *                  logistic dynamic and moving throughout weighted and directed
 *                  networks of three and four nodes, and more complex dendritic
 *                  and lattice networks of 25 nodes.
 *
 *  Project Structure:  
 *      Routines are divided between the four main programs:
 *      1. net3:        Triangular network of 3 nodes
 *      2. net4:        Four kinds of networks of 4 nodes (see net.c++ for
 *                      diagrams).
 *      3. net_dend:    Dendritic network of 25 nodes
 *      4. net_latt:    Square lattice network of 25 nodes
 *
 *  Limitations:
 *
 *  Dependencies:       libboost
 *
 *  Compiler Options:   -std=c++11 -lboost_program_options
 ***************************************************************************/


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
