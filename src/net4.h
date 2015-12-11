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
#include "network.h"

class Net4 : public Network
{
    private:
        static constexpr int nnodes = 4;
    
    public:
        int network_type;
        struct ResultsAll
        {
            double conn_mn, nmn_node [nnodes + 1], nsd_node [nnodes + 1], 
                   nmn_net, nsd_net;
        };
        ResultsAll results;

        Net4 ()
        {
            alpha.resize (nnodes, nnodes);
            pmat.resize (nnodes, nnodes);
        }
        ~Net4 ()
        {
        }

        int get_nnodes () { return nnodes;  }

        void fill_alpha (base_generator_type * generator);
        void make_pmat (base_generator_type * generator);
};
