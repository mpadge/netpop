/***************************************************************************
 *  Project:    netpop
 *  File:       randnet_dendritic.h
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
