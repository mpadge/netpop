/***************************************************************************
 *  Project:    netpop
 *  File:       net4-homog.c++
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
 *      2. net4:    Four the following kinds of networks of 4 nodes:
 *
 *      A =	*---*---*---*		B =	*---*---*
 *       							    |
 * 	        						    *
 *       					          *
 *          C =	*---*		         /|
 *		        |   |		D = *---* |
 *		        *---*		         \|
 *					                  *
 *      3. net_dend:    Dendritic network of 25 nodes
 *      4. net_latt:    Square lattice network of 25 nodes
 *
 *  Limitations:
 *
 *  Dependencies:       libboost
 *
 *  Compiler Options:   -std=c++11 -lboost_program_options
 *
 *  net4-homog is modified from net4, to simulate a structurally homogeneous
 *  network (that is, one with pars.k0sd = pars.alphasd = pars.k0sd = 0.0).
 ***************************************************************************/


#include "net4.h"


/************************************************************************
 ************************************************************************
 **                                                                    **
 **                         MAIN FUNCTION                              **
 **                                                                    **
 ************************************************************************
 ************************************************************************/

int main (int argc, char *argv[])
{
    double tempd, progress;
    std::string tempstr; // for simple timeout for non-linux systems
    Net4 net;
    std::ofstream out_file;
    clock_t time_start;
    time_t seed;
    base_generator_type generator(42u);

    try {
        boost::program_options::options_description generic("Generic options");
        generic.add_options()
            ("version,v", "print version std::string")
            ("help", "produce help message")    
            ;

        // The networks are identical, and only the population dynamics remain
        // random, so the only parameter is nRepeats
        boost::program_options::options_description config("Configuration");
        config.add_options()
            ("nTrials,n", boost::program_options::value <int>
             (&net.pars.nTrials)->default_value (1e6), "Number of trials")
            ("r,r", boost::program_options::value <double>
             (&net.pars.r)->default_value (0.1), "Growth rate, r")
            ("ksd,k", boost::program_options::value <double>
             (&net.pars.ksd)->default_value (0.1), "SD of k")
            ;

        // Not used here
        boost::program_options::options_description hidden("Hidden options");
        hidden.add_options()
            ("hidden-option", boost::program_options::value
             <std::vector<std::string> >(), "hidden option")
            ;

        boost::program_options::options_description cmdline_options;
        cmdline_options.add(generic).add(config).add(hidden);

        boost::program_options::options_description visible("Allowed options");
        visible.add(generic).add(config);

        boost::program_options::variables_map vm;
        store(boost::program_options::command_line_parser(argc, argv).
                options(cmdline_options).run(), vm);

        notify(vm);

        if (vm.count("help")) {
            std::cout << visible << std::endl;
            return 0;
        }

        if (vm.count("version")) {
            std::cout << "net4, version 1.0" << std::endl;
            return 0;
        }

    }
    catch(std::exception& e)
    {
        std::cout << e.what() << std::endl;
        return 1;
    }    
    std::cout << "nTrials for each alpha = " << net.pars.nTrials << 
        "; with ksd = " << net.pars.ksd << std::endl;

    int nnodes = net.get_nnodes (); // == 4

    std::string fname = "net4_results_homog_r";
    std::stringstream ss;
    ss.str ("");
    if (net.pars.r < 0.1)
        fname += "0";
    if (net.pars.r < 1.0)
        fname += "0";
    ss << round (net.pars.r * 100.0);
    fname += ss.str () + "_ksd";
    ss.str ("");
    if (net.pars.ksd < 0.01)
        fname += "0";
    if (net.pars.ksd < 0.1)
        fname += "0";
    ss << round (1000.0 * net.pars.ksd);
    fname += ss.str () + ".txt";

    time_start = clock ();

    // makepmat explicitly presumes k0=1.0 and doesn't use the value, but
    // iterate_population still uses it.
    net.pars.k0 = 1.0; 
    net.pars.k0sd = 0.0;
    net.pars.alphasd = 0.0;
    net.k0.resize (nnodes);
    for (int i=0; i<nnodes; i++) 
        net.k0 (i) = net.pars.k0;

    time (&seed);
    generator.seed (static_cast <unsigned int> (seed));

    std::cout.setf (std::ios::fixed, std::ios::floatfield);   
    std::cout.precision (4);

    out_file.open (fname.c_str(), std::ofstream::out);
    out_file << "alpha,\tconn0,\tmn0node0,\tmn0node1,\tmn0node2,\tmn0node3,\t" <<
        "mn0node_all,\tsd0node0,\tsd0node1,\tsd0node2,\tsd0node3,\t" <<
        "sd0node_all,\tmn0net,\tsd0net,\tcov001,\tcov002,\tcov003,\t" <<
        "cov012,\tcov013,\tcov023,\t" <<
        "conn1,\tmn1node0,\tmn1node1,\tmn1node2,\tmn1node3,\tmn1node_all,\t" <<
        "sd1node0,\tsd1node1,\tsd1node2,\tsd1node3,\tsd1node_all,\tmn1net,\t" <<
        "sd1net,\tcov101,\tcov102,\tcov103,\tcov112,\tcov113,\tcov123,\t" <<
        "conn2,\tmn2node0,\tmn2node1,\tmn2node2,\tmn2node3,\tmn2node_all,\t" <<
        "sd2node0,\tsd2node1,\tsd2node2,\tsd2node3,\tsd2node_all,\tmn2net,\t" <<
        "sd2net,\tcov201,\tcov202,\tcov203,\tcov212,\tcov213,\tcov223,\t" <<
        "conn3,\tmn3node0,\tmn3node1,\tmn3node2,\tmn3node3,\tmn3node_all,\t" <<
        "sd3node0,\tsd3node1,\tsd3node2,\tsd3node3,\tsd3node_all,\tmn3net,\t" <<
        "sd3net,\tcov301,\tcov302,\tcov303,\tcov312,\tcov313,\tcov323" << 
        std::endl;
    for (int i=1; i<=100; i++) 
    {
        net.pars.alpha0 = (double) i / 100.0;
        out_file << net.pars.alpha0;
        for (net.network_type = 0; net.network_type<4; net.network_type++)
        {
            net.make_pmat (&generator);
            net.iterate_population (&generator, nnodes);

            assert (net.results1.nmn_network > DOUBLE_MIN);
            
            net.results.conn_mn = net.results1.connectivity;
            for (int j=0; j<(nnodes + 1); j++) 
            {
                net.results.nmn_node [j] = net.results1.nmn_node [j];
                net.results.nsd_node [j] = net.results1.nsd_node [j];
            }
            net.results.nmn_net = net.results1.nmn_network;
            net.results.nsd_net = net.results1.nsd_network;

            out_file << ",\t" << net.results.conn_mn;
            for (int k=0; k<(nnodes + 1); k++) 
                out_file << ",\t" << net.results.nmn_node [k];
            for (int k=0; k<(nnodes + 1); k++) 
                out_file << ",\t" << net.results.nsd_node [k];
            out_file << ",\t" << net.results.nmn_net << ",\t" << 
                net.results.nsd_net;
            for (int k=0; k<(nnodes - 1); k++) 
                for (int m=(k+1); m<nnodes; m++)
                    out_file << ",\t" << net.results1.cov (k, m);
        } // end for net.network_type
        out_file << std::endl;

        progress = (double) i / 100.0;
        tempd = ((double) clock () - (double) time_start) / 
            (double) CLOCKS_PER_SEC;
        progLine (progress, tempd);
    } // end for i
    out_file.close();
    std::cout << std::endl << std::endl;

    return 0;
} // end main




/************************************************************************
 ************************************************************************
 **                                                                    **
 **                          MAKE_PMAT                                 **
 **                                                                    **
 ************************************************************************
 ************************************************************************/

void Net4::make_pmat (base_generator_type * generator)
{
    // And again, generator is passed yet not used.
    int nnodes = get_nnodes ();
    double tempd;

    if (network_type == 0) 
    {
        /*
         * 	0---1---2---3
         */
        tempd = 1.0 + pars.alpha0 + pars.alpha0 * pars.alpha0 +
            pars.alpha0 * pars.alpha0 * pars.alpha0;
        pmat (0, 0) = pmat (3, 3) = 1.0 / tempd;
        pmat (0, 1) = pmat (3, 2) = pars.alpha0 * pmat (0,0);
        pmat (0, 2) = pmat (3, 1) = pars.alpha0 * pmat (0,1);
        pmat (0, 3) = pmat (3, 0) = pars.alpha0 * pmat (0,2);

        tempd = 1.0 + 2.0 * pars.alpha0 + pars.alpha0 * pars.alpha0;
        pmat (1, 1) = pmat (2, 2) = 1.0 / tempd;
        pmat (1, 0) = pmat (2, 3) = pars.alpha0 * pmat (1, 1);
        pmat (1, 2) = pmat (2, 1) = pars.alpha0 * pmat (1, 1);
        pmat (1, 3) = pmat (2, 0) = pars.alpha0 * pmat (1, 0);
    } else if (network_type == 1) {
        /*
         *	    3
         *	    |
         *	    2
         *	   / \
         *	  0   1
         */
        tempd = 1.0 + pars.alpha0 + 2.0 * pars.alpha0 * pars.alpha0;
        pmat (0, 0) = pmat (1, 1) = pmat (3, 3) = 1.0 / tempd;
        pmat (0, 2) = pmat (1, 2) = pmat (3, 2) = pars.alpha0 * pmat (0, 0);
        pmat (0, 1) = pmat (1, 0) = pars.alpha0 * pmat (0, 2);
        pmat (0, 3) = pmat (3, 0) = pmat (0, 1);
        pmat (1, 3) = pmat (3, 1) = pmat (0, 1);

        tempd = 1.0 + 3.0 * pars.alpha0;
        pmat (2, 2) = 1.0 / tempd;
        pmat (2, 0) = pmat (2, 1) = pmat (2, 3) = pars.alpha0 * pmat (2, 2);
    } else if (network_type == 2) {
        /*
         * 	2---3
         * 	|   |
         * 	0---1
         */
        tempd = 1.0 + 2.0 * pars.alpha0 + pars.alpha0 * pars.alpha0;
        pmat (0, 0) = pmat (1, 1) = pmat (2, 2) = pmat (3, 3) = 1.0 / tempd;
        pmat (0, 1) = pmat (0, 2) = pars.alpha0 * pmat (0, 0);
        pmat (1, 0) = pmat (1, 3) = pmat (2, 0) = pmat (2, 3) = pmat (0, 1);
        pmat (3, 2) = pmat (3, 1) = pmat (0, 1);
        pmat (0, 3) = pmat (3, 0) = pars.alpha0 * pars.alpha0 * pmat (0, 0);
        pmat (1, 2) = pmat (2, 1) = pmat (0, 3);
    } else if (network_type == 3) {
        /*
         * 	     3
         * 	     |
         * 	     2
         * 	    / \
         * 	   0---1
         */
        tempd = 1.0 + 2.0 * pars.alpha0 + pars.alpha0 * pars.alpha0;
        pmat (0, 0) = pmat (1, 1) = 1.0 / tempd;
        pmat (0, 1) = pmat (0, 2) = pars.alpha0 * pmat (0, 0);
        pmat (1, 0) = pmat (1, 2) = pmat (0, 1);
        pmat (0, 3) = pmat (1, 3) = pars.alpha0 * pmat (0, 1);

        tempd = 1.0 + 3.0 * pars.alpha0;
        pmat (2, 2) = 1.0 / tempd;
        pmat (2, 0) = pmat (2, 1) = pmat (2, 3) = pars.alpha0 * pmat (2, 2);

        tempd = 1.0 + pars.alpha0 + 2.0 * pars.alpha0 * pars.alpha0;
        pmat (3, 3) = 1.0 / tempd;
        pmat (3, 2) = pars.alpha0 * pmat (3, 3);
        pmat (3, 1) = pmat (3, 0) = pars.alpha0 * pmat (3, 2);
    } // end if network type == 3

    // Calculate connectivity from full pmat values
    results1.connectivity = (double) nnodes;
    for (int i=0; i<nnodes; i++) 
        results1.connectivity -= pmat (i, i);
    results1.connectivity = results1.connectivity / (double) nnodes;
    // Then rescale to growth rate. Results from pscale = pars.r show very
    // little difference between the different networks. The value of 0.1 is a
    // random guess of a value that might enhance movement sufficiently to
    // reveal some stronger differences.
    double pscale = 2.0 * pars.r;
    for (int i=0; i<nnodes; i++) 
        pmat (i, i) = 1.0 - pscale * (1.0 - pmat (i, i));
    for (int i=0; i<(nnodes - 1); i++) 
        for (int j=(i + 1); j<nnodes; j++) 
        {
            pmat (i, j) = pscale * pmat (i, j);
            pmat (j, i) = pscale * pmat (j, i);
        }
}
