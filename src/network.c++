/***************************************************************************
 *  Project:    netpop
 *  File:       network.c++
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

#include "network.h"


/************************************************************************
 ************************************************************************
 **                                                                    **
 **                         GET_FILENAME                               **
 **                                                                    **
 ************************************************************************
 ************************************************************************/

void Network::get_filename (int nnodes)
{

    std::stringstream ss;
    ss.str ("");
    ss << nnodes;
    filename = "net" + ss.str () + "_results_r";
    ss.str ("");

    if (pars.r < 0.1)
        filename += "0";
    if (pars.r < 1.0)
        filename += "0";
    ss.str (""); 
    //ss << floor (pars.r * 10.0); // presumes 0.1 <= r < 1
    ss << round (pars.r * 100.0); 
    filename += ss.str () + "_ksd";
    if (pars.ksd < 0.01) 
        filename += "0";
    if (pars.ksd < 0.1) 
        filename += "0";
    ss.str (""); ss << round (1000.0 * pars.ksd);
    filename += ss.str () + "_k0sd";
    if (pars.k0sd < 0.01) 
        filename += "0";
    if (pars.k0sd < 0.1) 
        filename += "0";
    ss.str (""); ss << round (1000.0 * pars.k0sd);
    filename += ss.str() + "_alphasd";
    if (pars.alphasd < 0.01) 
        filename += "0";
    if (pars.alphasd < 0.1) 
        filename += "0";
    ss.str (""); ss << round (1000.0 * pars.alphasd);
    filename += ss.str() + ".txt";
    ss.str ("");
}


/************************************************************************
 ************************************************************************
 **                                                                    **
 **                         ITERATE_POPULATION                         **
 **                                                                    **
 ************************************************************************
 ************************************************************************/

void Network::iterate_population (base_generator_type * generator, int nnodes)
{
    int tempi, count = 0;
    bool flag = true, bigflag = false;
    double tempd, n [nnodes], nold [nnodes], kvals [nnodes];

    boost::normal_distribution<> norm_dist (0.0, 1.0);
    boost::variate_generator<base_generator_type&,
        boost::normal_distribution<> > rnorm((*generator), norm_dist);
    // Burn generator in
    for (int i=0; i<20; i++) 
        tempd = rnorm();

    results1.nmn_node.resize (nnodes + 1); 
    results1.nsd_node.resize (nnodes + 1);
    results1.cov.resize (nnodes, nnodes);

    while (flag) 
    {
        for (int j=0; j<nnodes; j++) 
            nold [j] = k0 [j];
        for (int j=0; j<(nnodes + 1); j++) 
        {
            results1.nmn_node (j) = 0.0;
            results1.nsd_node (j) = 0.0;
        }
        results1.nmn_network = 0.0;
        results1.nsd_network = 0.0;
        flag = false;
        for (int i=0; i<(nnodes - 1); i++)
            for (int j=(i+1); j<nnodes; j++)
                results1.cov (i, j) = 0.0;

        for (int i=0; i<(runin + pars.nTrials); i++) 
        {
            // Movement through network
            for (int j=0; j<nnodes; j++) 
            {
                n [j] = 0.0;
                for (int k=0; k<nnodes; k++) 
                    n [j] += pmat (k, j) * nold [k];
            }
            // Change carrying capacities
            for (int j=0; j<nnodes; j++) 
            {
                kvals [j] = k0 [j] + pars.ksd * rnorm ();
                while (kvals [j] < minqc || 
                        kvals [j] > (2.0 * k0 [j] - minqc))
                    kvals [j] = k0 [j] + pars.ksd * rnorm ();
            }
            // Population dynamic
            tempi = 0;
            for (int j=0; j<nnodes; j++) 
            {
                nold [j] = n [j] + pars.r * n [j] * n [j] / kvals [j] - 
                    pars.r * n [j] * n [j] * n [j] / (kvals [j] * kvals [j]);
                if (nold [j] < 0.0) 
                    tempi++;
            } // end for j
            if (tempi == nnodes) 
            {
                flag = true;
                break;
            } else if (i > runin) {
                tempd = 0.0;
                for (int j=0; j<nnodes; j++) 
                {
                    tempd += nold [j];
                    results1.nmn_node (j) += nold [j];
                    results1.nsd_node (j) += nold [j] * nold [j];
                }
                results1.nmn_network += tempd;
                results1.nsd_network += tempd * tempd;
                for (int j=0; j<(nnodes - 1); j++)
                    for (int k=(j+1); k<nnodes; k++)
                        results1.cov (j, k) += nold [j] * nold [k];
            }
            // Then set all negative nodal abundances to zero
            for (int j=0; j<nnodes; j++)
                if (nold [j] < 0.0)
                    nold [j] = 0.0;
        } // end for i over nTrials
        if (!flag) 
        {
            for (int j=0; j<nnodes; j++) 
            {
                results1.nmn_node (j) = results1.nmn_node (j) 
                    / (double) pars.nTrials;
                results1.nsd_node (j) = results1.nsd_node (j) / 
                    (double) pars.nTrials -
                    results1.nmn_node (j) * results1.nmn_node (j);
            }
            for (int j=0; j<(nnodes - 1); j++)
                for (int k=(j+1); k<nnodes; k++)
                    results1.cov (j, k) = results1.cov (j, k) / 
                        (double) pars.nTrials -
                        results1.nmn_node (j) * results1.nmn_node (k);
        } else { 
            count++;
        }
        if (count >= pars.nTrials) 
        {
            flag = false;
            bigflag = true;	
        }
    } // end while flag
    if (bigflag) 
    {
        for (int i=0; i<(nnodes + 1); i++) 
        {
            results1.nmn_node (i) = DOUBLE_MIN;
            results1.nsd_node (i) = DOUBLE_MIN;
        }
        results1.nmn_network = DOUBLE_MIN;
        results1.nsd_network = DOUBLE_MIN;
    } else {
        results1.nmn_node (nnodes) = 0.0;
        results1.nsd_node (nnodes) = 0.0;
        for (int i=0; i<(nnodes + 1); i++) 
        {
            results1.nmn_node (nnodes) += results1.nmn_node [i];
            results1.nsd_node (nnodes) += results1.nsd_node [i];
        }
        results1.nmn_node (nnodes) = results1.nmn_node (nnodes) / (double) nnodes;
        results1.nsd_node (nnodes) = results1.nsd_node (nnodes) / (double) nnodes;
        results1.nmn_network = results1.nmn_network  / (double) pars.nTrials;
        results1.nsd_network = results1.nsd_network / (double) pars.nTrials -
            results1.nmn_network * results1.nmn_network;
        results1.nmn_network = results1.nmn_network / (double) nnodes;
        // So it's on the same scale as nodal abundance.
    }
}
