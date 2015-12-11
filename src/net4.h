/***************************************************************************
 *  Project:    netpop
 *  File:       net4.h
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
 *                  networks of three and four nodes.
 *
 *  Project Structure:  
 *      Routines are divided between the two main programs:
 *      1. net3:    Triangular network of 3 nodes
 *      2. net4:    Four kinds of networks of 4 nodes (see net.c++ for
 *                  diagrams).
 *
 *  Limitations:
 *
 *  Dependencies:       libboost
 *
 *  Compiler Options:   -std=c++11 -lboost_program_options
 ***************************************************************************/

#include "utils.h"

class Network
{
    private:
        static constexpr int nnodes = 4, maxTrials = 100, runin = 100;
        static constexpr double minqc = 0.001;
        // minqc is important to prevent dividing by sometimes *really* small
        // numbers
    public:
        int network_type;
        double k0 [nnodes], alpha [nnodes] [nnodes], pmat [nnodes] [nnodes];
        std::string filename = "aaasim3_results_r"; 
        // filename is extended in "get_filename"

        struct Parameters 
        {
            int nTrials, nRepeats;
            double k0, k0sd, ksd, alpha0, alphasd, r;
        };
        Parameters pars;

        struct Results1Net
        {
            double connectivity, nmn_node [nnodes + 1], nsd_node [nnodes + 1], 
                   nmn_network, nsd_network, cov [nnodes] [nnodes];
        };
        Results1Net results1;

        struct ResultsAll
        {
            double conn_mn, nmn_node [nnodes + 1], nsd_node [nnodes + 1], 
                   nmn_net, nsd_net;
        };
        ResultsAll results;

        Network ()
        {
        }
        ~Network ()
        {
        }

        int get_nnodes () { return nnodes;  }

        void get_filename ();
        void fill_alpha (base_generator_type * generator);
        void make_pmat (base_generator_type * generator);
        void iterate_population (base_generator_type * generator);
};

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
