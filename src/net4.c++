/***************************************************************************
 *  Project:    netpop
 *  File:       net4.c++
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
 *
 * Nodal variances are calculated for each network type over the full range of
 * values of connectivity (= alpha).
 *
 *  Limitations:
 *
 *  Dependencies:       libboost
 *
 *  Compiler Options:   -std=c++11 -lboost_program_options
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
    int netcount;
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

        boost::program_options::options_description config("Configuration");
        config.add_options()
            ("nTrials,t", boost::program_options::value <int>
                (&net.pars.nTrials)->default_value (1000), "Number of trials")
            ("nRepeats,n", boost::program_options::value <int>
                (&net.pars.nRepeats)->default_value (1000), "Number of repeats")
            ("k0sd,s", boost::program_options::value <double>
                (&net.pars.k0sd)->default_value (0.1), "SD of k0")
            ("alphasd,a", boost::program_options::value <double>
                (&net.pars.alphasd)->default_value (0.1), "SD of alpha")
            ("ksd,k", boost::program_options::value <double>
                (&net.pars.ksd)->default_value (0.1), "SD of k")
            ("r,r", boost::program_options::value <double>
                (&net.pars.r)->default_value (0.1), "Growth rate, r")
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
    int nnodes = net.get_nnodes (); // == 4
    net.get_filename (nnodes);

    std::cout << "nTrials for each alpha = " << net.pars.nTrials <<
        "; with results averaged over " << net.pars.nRepeats << 
        " repeats." << std::endl;
    std::cout << "k0sd = " << net.pars.k0sd << "; alphasd = " << 
        net.pars.alphasd << "; ksd = " << net.pars.ksd << std::endl;

    time_start = clock ();

    net.pars.alpha0 = 0.5;
    net.pars.k0 = 0.5;

    time (&seed);
    generator.seed (static_cast <unsigned int> (seed));

    //cout.setf(0,ios::floatfield);            // floatfield not set
    std::cout.setf (std::ios::fixed, std::ios::floatfield);   
    std::cout.precision (4);

    out_file.open (net.filename.c_str(), std::ofstream::out);
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
        for (int j=0; j<4; j++) // the four network types
        {
            for (int k=0; k<(nnodes + 1); k++) 
            {
                net.results.nmn_node [k] = 0.0;
                net.results.nsd_node [k] = 0.0;
            }
            net.results.nmn_net = 0.0;
            net.results.nsd_net = 0.0;
            net.results.conn_mn = 0.0;
            netcount = 0;
            for (int k=0; k<net.pars.nRepeats; k++) 
            {
                net.fill_alpha (&generator);
                net.make_pmat (&generator);
                net.iterate_population (&generator, nnodes);
                if (net.results1.nmn_network > DOUBLE_MIN) 
                {
                    netcount++;
                    net.results.conn_mn += net.results1.connectivity;
                    for (int m=0; m<(nnodes + 1); m++) 
                    {
                        net.results.nmn_node [m] += 
                            net.results1.nmn_node [m];
                        net.results.nsd_node [m] += 
                            net.results1.nsd_node [m];
                    }
                    net.results.nmn_net += net.results1.nmn_network;
                    net.results.nsd_net += net.results1.nsd_network;
                }
            }
            if (netcount > 0) 
            {
                net.results.conn_mn = net.results.conn_mn / (double) netcount;
                for (int k=0; k<(nnodes + 1); k++) 
                {
                    net.results.nmn_node [k] = net.results.nmn_node [k] / 
                        (double) netcount;
                    net.results.nsd_node [k] = net.results.nsd_node [k] / 
                        (double) netcount;
                }
                net.results.nmn_net = net.results.nmn_net / (double) netcount;
                net.results.nsd_net = net.results.nsd_net / (double) netcount;
            } else {
                net.results.conn_mn = DOUBLE_MIN;
                for (int k=0; k<(nnodes + 1); k++) 
                {
                    net.results.nmn_node [k] = DOUBLE_MIN;
                    net.results.nsd_node [k] = DOUBLE_MIN;
                }
                net.results.nmn_net = DOUBLE_MIN;
                net.results.nsd_net = DOUBLE_MIN;
            }
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
        } // end for j
        out_file << std::endl;
        //std::cout << "."; std::cout.flush();

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
 **                          FILL_ALPHA                                **
 **                                                                    **
 ************************************************************************
 ************************************************************************/

void Net4::fill_alpha (base_generator_type * generator)
{
    int nnodes = get_nnodes ();
    double tempd;

    boost::normal_distribution<> norm_dist (0.0, 1.0);
    boost::variate_generator<base_generator_type&,
        boost::normal_distribution<> > rnorm((*generator), norm_dist);
    // Burn generator in
    for (int i=0; i<20; i++) 
        tempd = rnorm();

    k0.resize (nnodes);

    for (int i=0; i<nnodes; i++) 
    {
        k0 (i) = pars.k0 + pars.k0sd * rnorm ();
        while (k0 (i) < minqc || k0 (i) > (1.0 - minqc))
            k0 (i) = pars.k0 + pars.k0sd * rnorm ();
        for (int j=0; j<nnodes; j++)
            alpha (i, j) = 0.0;
    }

    // First set up connectivities, which obviously has to explicitly assume
    // that nnodes = 4!
    std::vector <std::pair <int, int> > connlist;
    if (network_type == 0) 
    {
        /*
         * 	0---1---2---3
         */
        connlist.push_back (std::pair <int, int> (0, 1));
        connlist.push_back (std::pair <int, int> (1, 0));
        connlist.push_back (std::pair <int, int> (1, 2));
        connlist.push_back (std::pair <int, int> (2, 1));
        connlist.push_back (std::pair <int, int> (2, 3));
        connlist.push_back (std::pair <int, int> (3, 2));
    } else if (network_type == 1) {
        /*
         *	    3
         *	    |
         *	    2
         *	   / \
         *	  0   1
         */
        connlist.push_back (std::pair <int, int> (0, 2));
        connlist.push_back (std::pair <int, int> (2, 0));
        connlist.push_back (std::pair <int, int> (1, 2));
        connlist.push_back (std::pair <int, int> (2, 1));
        connlist.push_back (std::pair <int, int> (2, 3));
        connlist.push_back (std::pair <int, int> (3, 2));
    } else if (network_type == 2) {
        /*
         * 	2---3
         * 	|   |
         * 	0---1
         */
        connlist.push_back (std::pair <int, int> (0, 1));
        connlist.push_back (std::pair <int, int> (1, 0));
        connlist.push_back (std::pair <int, int> (1, 3));
        connlist.push_back (std::pair <int, int> (3, 1));
        connlist.push_back (std::pair <int, int> (2, 3));
        connlist.push_back (std::pair <int, int> (3, 2));
        connlist.push_back (std::pair <int, int> (0, 2));
        connlist.push_back (std::pair <int, int> (2, 0));
    } else if (network_type == 3) {
        /*
         * 	     3
         * 	     |
         * 	     2
         * 	    / \
         * 	   0---1
         */
        connlist.push_back (std::pair <int, int> (0, 2));
        connlist.push_back (std::pair <int, int> (2, 0));
        connlist.push_back (std::pair <int, int> (1, 2));
        connlist.push_back (std::pair <int, int> (2, 1));
        connlist.push_back (std::pair <int, int> (2, 3));
        connlist.push_back (std::pair <int, int> (3, 2));
        connlist.push_back (std::pair <int, int> (0, 1));
        connlist.push_back (std::pair <int, int> (1, 0));
    }
    
    // Then fill the connectivities
    std::vector <std::pair <int, int> >::const_iterator itr;
    for (itr = connlist.begin(); itr < connlist.end(); itr++) 
    {
        tempd = pars.alpha0 + pars.alphasd * rnorm ();
        while (tempd < minqc || tempd > (1.0 - minqc))
            tempd = pars.alpha0 + pars.alphasd * rnorm ();
        alpha (itr -> first, itr -> second) = tempd;
    } // end for itr
    for (int i=0; i<nnodes; i++)
        alpha (i, i) = 1.0;

    connlist.resize (0);
}


/************************************************************************
 ************************************************************************
 **                                                                    **
 **                          MAKE_PMAT                                 **
 **                                                                    **
 ************************************************************************
 ************************************************************************/

void Net4::make_pmat (base_generator_type * generator)
{
    int nnodes = get_nnodes ();
    // Also has to explicitly assume here that nnodes = 4!
    double tempd;

    // Note that doing this with the full makepmat routine using the shortest
    // path algorithm (as in the 25-node versions) makes this really enormously
    // slower than the following manual construction.
    if (network_type == 0) 
    {
        /*
         * 	0---1---2---3
         */
        tempd = k0 [0] + alpha (0, 1) * k0 [1] + 
                alpha (0, 1) * alpha (1, 2) * k0 [2] +
                alpha (0, 1) * alpha (1, 2) * alpha (2, 3) * k0 [3];
        pmat (0, 0) = k0 [0] / tempd;
        pmat (0, 1) = alpha (0, 1) * k0 [1] / tempd;
        pmat (0, 2) = alpha (0, 1) * alpha (1, 2) * k0 [2] / tempd;
        pmat (0, 3) = alpha (0, 1) * alpha (1, 2) * alpha (2, 3) * 
                        k0 [3] / tempd;

        tempd = k0 [1] + alpha (1, 0) * k0 [0] + 
                alpha (1, 2) * k0 [2] + alpha (1, 2) * alpha (2, 3) * k0 [3];
        pmat (1, 0) = alpha (1, 0) * k0 [0] / tempd;
        pmat (1, 1) = k0 [1] / tempd;
        pmat (1, 2) = alpha (1, 2) * k0 [2] / tempd;
        pmat (1, 3) = alpha (1, 2) * alpha (2, 3) * k0 [3] / tempd;

        tempd = k0 [2] + alpha (2, 1) * k0 [1] + 
            alpha (2, 1) * alpha (1, 0) * k0 [0] + alpha (2, 3) * k0 [3];
        pmat (2, 0) = alpha (2, 1) * alpha (1, 0) * k0 [0] / tempd;
        pmat (2, 1) = alpha (2, 1) * k0 [1] / tempd;
        pmat (2, 2) = k0 [2] / tempd;
        pmat (2, 3) = alpha (2, 3) * k0 [3] / tempd;

        tempd = k0 [3] + alpha (3, 2) * k0 [2] + 
                alpha (3, 2) * alpha (2, 1) * k0 [1] +
                alpha (3, 2) * alpha (2, 1) * alpha (1, 0) * k0 [0];
        pmat (3, 0) = alpha (3, 2) * alpha (2, 1) * alpha (1, 0) * 
                        k0 [0] / tempd;
        pmat (3, 1) = alpha (3, 2) * alpha (2, 1) * k0 [1] / tempd;
        pmat (3, 2) = alpha (3, 2) * k0 [2] / tempd;
        pmat (3, 3) = k0 [3] / tempd;
    } else if (network_type == 1) {
        /*
         *	    3
         *	    |
         *	    2
         *	   / \
         *	  0   1
         */
        tempd = k0 [0] + alpha (0, 2) * k0 [2] + 
                alpha (0, 2) * alpha (2, 1) * k0 [1] +
                alpha (0, 2) * alpha (2, 3) * k0 [3];
        pmat (0, 0) = k0 [0] / tempd;
        pmat (0, 1) = alpha (0, 2) * alpha (2, 1) * k0 [1] / tempd;
        pmat (0, 2) = alpha (0, 2) * k0 [2] / tempd;
        pmat (0, 3) = alpha (0, 2) * alpha (2, 3) * k0 [3] / tempd;

        tempd = k0 [1] + alpha (1, 2) * k0 [2] + 
                alpha (1, 2) * alpha (2, 0) * k0 [0] +
                alpha (1, 2) * alpha (2, 3) * k0 [3];
        pmat (1, 0) = alpha (1, 2) * alpha (2, 0) * k0 [0] / tempd;
        pmat (1, 1) = k0 [1] / tempd;
        pmat (1, 2) = alpha (1, 2) * k0 [2] / tempd;
        pmat (1, 3) = alpha (1, 2) * alpha (2, 3) * k0 [3] / tempd;

        tempd = k0 [2] + alpha (2, 0) * k0 [0] + 
                alpha (2, 1) * k0 [1] + alpha (2, 3) * k0 [3];
        pmat (2, 0) = alpha (2, 0) * k0 [0] / tempd;
        pmat (2, 1) = alpha (2, 1) * k0 [1] / tempd;
        pmat (2, 2) = k0 [2] / tempd;
        pmat (2, 3) = alpha (2, 3) * k0 [3] / tempd;

        tempd = k0 [3] + alpha (3, 2) * k0 [2] + 
                alpha (3, 2) * alpha (2, 1) * k0 [1] +
                alpha (3, 2) * alpha (2, 0) * k0 [0];
        pmat (3, 0) = alpha (3, 2) * alpha (2, 0) * k0 [0] / tempd;
        pmat (3, 1) = alpha (3, 2) * alpha (2, 1) * k0 [1] / tempd;
        pmat (3, 2) = alpha (3, 2) * k0 [2] / tempd;
        pmat (3, 3) = k0 [3] / tempd;
    } else if (network_type == 2) {
        /*
         * 	2---3
         * 	|   |
         * 	0---1
         *
         ********* node 0 **********/
        if ((alpha (0, 2) * alpha (2, 3) * alpha (3, 1)) > alpha (0, 1)) 
        {
            tempd = k0 [0] + alpha (0, 2) * k0 [2] + 
                    alpha (0, 2) * alpha (2, 3) * k0 [3] +
                    alpha (0, 2) * alpha (2, 3) * alpha (3, 1) * k0 [1];
            pmat (0, 0) = k0 [0] / tempd;
            pmat (0, 1) = alpha (0, 2) * alpha (2, 3) * alpha (3, 1) * 
                            k0 [1] / tempd;
            pmat (0, 2) = alpha (0, 2) * k0 [2] / tempd;
            pmat (0, 3) = alpha (0, 2) * alpha (2, 3) * k0 [3] / tempd;
        } else if ((alpha (0, 1) * alpha (1, 3) * alpha (3, 2)) > 
                alpha (0, 2)) {
            tempd = k0 [0] + alpha (0, 1) * k0 [1] + 
                    alpha (0, 1) * alpha (1, 3) * k0 [3] +
                    alpha (0, 1) * alpha (1, 3) * alpha (3, 2) * k0 [2];
            pmat (0, 0) = k0 [0] / tempd;
            pmat (0, 1) = alpha (0, 1) * k0 [1] / tempd;
            pmat (0, 2) = alpha (0, 1) * alpha (1, 3) * alpha (3, 2) * 
                            k0 [2] / tempd;
            pmat (0, 3) = alpha (0, 1) * alpha (1, 3) * k0 [3] / tempd;
        } else if ((alpha (0, 2) * alpha (2, 3)) > 
                (alpha (0, 1) * alpha (1, 3))) {
            tempd = k0 [0] + alpha (0, 1) * k0 [1] + alpha (0, 2) * k0 [2] + 
                    alpha (0, 2) * alpha (2, 3) * k0 [3];
            pmat (0, 0) = k0 [0] / tempd;
            pmat (0, 1) = alpha (0, 1) * k0 [1] / tempd;
            pmat (0, 2) = alpha (0, 2) * k0 [2] / tempd;
            pmat (0, 3) = alpha (0, 2) * alpha (2, 3) * k0 [3] / tempd;
        } else {
            tempd = k0 [0] + alpha (0, 1) * k0 [1] + 
                    alpha (0, 2) * k0 [2] + 
                    alpha (0, 1) * alpha (1, 3) * k0 [3];
            pmat (0, 0) = k0 [0] / tempd;
            pmat (0, 1) = alpha (0, 1) * k0 [1] / tempd;
            pmat (0, 2) = alpha (0, 2) * k0 [2] / tempd;
            pmat (0, 3) = alpha (0, 1) * alpha (1, 3) * k0 [3] / tempd;
        }

        /*
         * 	2---3
         * 	|   |
         * 	0---1
         *
         ********* node 1 **********/
        if ((alpha (1, 3) * alpha (3, 2) * alpha (2, 0)) > alpha (1, 0)) 
        {
            tempd = k0 [1] + alpha (1, 3) * k0 [3] + 
                    alpha (1, 3) * alpha (3, 2) * k0 [2] +
                    alpha (1, 3) * alpha (3, 2) * alpha (2, 0) * k0 [0];
            pmat (1, 0) = alpha (1, 3) * alpha (3, 2) * alpha (2, 0) * 
                            k0 [0] / tempd;
            pmat (1, 1) = k0 [1] / tempd;
            pmat (1, 2) = alpha (1, 3) * alpha (3, 2) * k0 [2] / tempd;
            pmat (1, 3) = alpha (1, 3) * k0 [3] / tempd;
        } else if ((alpha (1, 0) * alpha (0, 2) * alpha (2, 3)) > 
                alpha (1, 3)) {
            tempd = k0 [1] + alpha (1, 0) * k0 [0] + 
                    alpha (1, 0) * alpha (0, 2) * k0 [2] +
                    alpha (1, 0) * alpha (0, 2) * alpha (2, 3) * k0 [3];
            pmat (1, 0) = alpha (1, 0) * k0 [0] / tempd;
            pmat (1, 1) = k0 [1] / tempd;
            pmat (1, 2) = alpha (1, 0) * alpha (0, 2) * k0 [2] / tempd;
            pmat (1, 3) = alpha (1, 0) * alpha (0, 2) * alpha (2, 3) * 
                k0 [3] / tempd;
        } else if ((alpha (1, 3) * alpha (3, 2)) > 
                (alpha (1, 0) * alpha (0, 2))) {
            tempd = k0 [1] + alpha (1, 3) * k0 [3] + alpha (1, 0) * k0 [0] + 
                    alpha (1, 3) * alpha (3, 2) * k0 [2];
            pmat (1, 0) = alpha (1, 0) * k0 [0] / tempd;
            pmat (1, 1) = k0 [1] / tempd;
            pmat (1, 2) = alpha (1, 3) * alpha (3, 2) * k0 [2] / tempd;
            pmat (1, 3) = alpha (1, 3) * k0 [3] / tempd;
        } else {
            tempd = k0 [1] + alpha (1, 3) * k0 [3] + alpha (1, 0) * k0 [0] + 
                    alpha (1, 0) * alpha (0, 2) * k0 [2];
            pmat (1, 0) = alpha (1, 0) * k0 [0] / tempd;
            pmat (1, 1) = k0 [1] / tempd;
            pmat (1, 2) = alpha (1, 0) * alpha (0, 2) * k0 [2] / tempd;
            pmat (1, 3) = alpha (1, 3) * k0 [3] / tempd;
        }

        /*
         * 	2---3
         * 	|   |
         * 	0---1
         *
         ********* node 2 **********/
        if ((alpha (2, 0) * alpha (0, 1) * alpha (1, 3)) > alpha (2, 3)) 
        {
            tempd = k0 [2] + alpha (2, 0) * k0 [0] + 
                    alpha (2, 0) * alpha (0, 1) * k0 [1] +
                    alpha (2, 0) * alpha (0, 1) * alpha (1, 3) * k0 [3];
            pmat (2, 0) = alpha (2, 0) * k0 [0] / tempd;
            pmat (2, 1) = alpha (2, 0) * alpha (0, 1) * k0 [1] / tempd;
            pmat (2, 2) = k0 [2] / tempd;
            pmat (2, 3) = alpha (2, 0) * alpha (0, 1) * alpha (1, 3) * 
                            k0 [3] / tempd;
        } else if ((alpha (2, 3) * alpha (3, 1) * alpha (1, 0)) > 
                alpha (2, 0)) {
            tempd = k0 [2] + alpha (2, 3) * k0 [3] + 
                    alpha (2, 3) * alpha (3, 1) * k0 [1] +
                    alpha (2, 3) * alpha (3, 1) * alpha (1, 0) * k0 [0];
            pmat (2, 0) = alpha (2, 3) * alpha (3, 1) * alpha (1, 0) * 
                k0 [0] / tempd;
            pmat (2, 1) = alpha (2, 3) * alpha (3, 1) * k0 [1] / tempd;
            pmat (2, 2) = k0 [2] / tempd;
            pmat (2, 3) = alpha (2, 3) * k0 [3] / tempd;
        } else if ((alpha (2, 3) * alpha (3, 1)) > 
                (alpha (2, 0) * alpha (0, 1))) {
            tempd = k0 [2] + alpha (2, 0) * k0 [0] + alpha (2, 3) * k0 [3] + 
                    alpha (2, 3) * alpha (3, 1) * k0 [1];
            pmat (2, 0) = alpha (2, 0) * k0 [0] / tempd;
            pmat (2, 1) = alpha (2, 3) * alpha (3, 1) * k0 [1] / tempd;
            pmat (2, 2) = k0 [2] / tempd;
            pmat (2, 3) = alpha (2, 3) * k0 [3] / tempd;
        } else {
            tempd = k0 [2] + alpha (2, 0) * k0 [0] + alpha (2, 3) * k0 [3] + 
                    alpha (2, 0) * alpha (0, 1) * k0 [1];
            pmat (2, 0) = alpha (2, 0) * k0 [0] / tempd;
            pmat (2, 1) = alpha (2, 0) * alpha (0, 1) * k0 [1] / tempd;
            pmat (2, 2) = k0 [2] / tempd;
            pmat (2, 3) = alpha (2, 3) * k0 [3] / tempd;
        }

        /*
         * 	2---3
         * 	|   |
         * 	0---1
         *
         ********* node 3 **********/
        if ((alpha (3, 2) * alpha (2, 0) * alpha (0, 1)) > alpha (3, 1)) 
        {
            tempd = k0 [3] + alpha (3, 2) * k0 [2] + 
                    alpha (3, 2) * alpha (2, 0) * k0 [0] +
                    alpha (3, 2) * alpha (2, 0) * alpha (0, 1) * k0 [1];
            pmat (3, 0) = alpha (3, 2) * alpha (2, 0) * k0 [0] / tempd;
            pmat (3, 1) = alpha (3, 2) * alpha (2, 0) * alpha (0, 1) * 
                            k0 [1] / tempd;
            pmat (3, 2) = alpha (3, 2) * k0 [2] / tempd;
            pmat (3, 3) = k0 [3] / tempd;
        } else if ((alpha (3, 1) * alpha (1, 0) * alpha (0, 2)) > 
                alpha (3, 2)) {
            tempd = k0 [3] + alpha (3, 1) * k0 [1] + 
                    alpha (3, 1) * alpha (1, 0) * k0 [0] +
                    alpha (3, 1) * alpha (1, 0) * alpha (0, 2) * k0 [2];
            pmat (3, 0) = alpha (3, 1) * alpha (1, 0) * k0 [0] / tempd;
            pmat (3, 1) = alpha (3, 1) * k0 [1] / tempd;
            pmat (3, 2) = alpha (3, 1) * alpha (1, 0) * alpha (0, 2) * 
                            k0 [2] / tempd;
            pmat (3, 3) = k0 [3] / tempd;
        } else if ((alpha (3, 2) * alpha (2, 0)) > 
                    (alpha (3, 1) * alpha (1, 0))) {
            tempd = k0 [3] + alpha (3, 1) * k0 [1] + alpha (3, 2) * k0 [2] + 
                    alpha (3, 2) * alpha (2, 0) * k0 [0];
            pmat (3, 0) = alpha (3, 2) * alpha (2, 0) * k0 [0] / tempd;
            pmat (3, 1) = alpha (3, 1) * k0 [1] / tempd;
            pmat (3, 2) = alpha (3, 2) * k0 [2] / tempd;
            pmat (3, 3) = k0 [3] / tempd;
        } else {
            tempd = k0 [3] + alpha (3, 1) * k0 [1] + alpha (3, 2) * k0 [2] + 
                    alpha (3, 1) * alpha (1, 0) * k0 [0];
            pmat (3, 0) = alpha (3, 1) * alpha (1, 0) * k0 [0] / tempd;
            pmat (3, 1) = alpha (3, 1) * k0 [1] / tempd;
            pmat (3, 2) = alpha (3, 2) * k0 [2] / tempd;
            pmat (3, 3) = k0 [3] / tempd;
        }

    } else if (network_type == 3) {
        /*
         * 	     3
         * 	     |
         * 	     2
         * 	    / \
         * 	   0---1
         *
         ********* node 0 **********/
        if ((alpha (0, 2) * alpha (2, 1)) > alpha (0, 1)) 
        {
            tempd = k0 [0] + alpha (0, 2) * alpha (2, 1) * k0 [1] + 
                    alpha (0, 2) * k0 [2] + 
                    alpha (0, 2) * alpha (2, 3) * k0 [3];
            pmat (0, 0) = k0 [0] / tempd;
            pmat (0, 1) = alpha (0, 2) * alpha (2, 1) * k0 [1] / tempd;
            pmat (0, 2) = alpha (0, 2) * k0 [2] / tempd;
            pmat (0, 3) = alpha (0, 2) * alpha (2, 3) * k0 [3] / tempd;
        } else if ((alpha (0, 1) * alpha (1, 2)) > alpha (0, 2)) {
            tempd = k0 [0] + alpha (0, 1) * k0 [1] + 
                    alpha (0, 1) * alpha (1, 2) * k0 [2] +
                    alpha (0, 1) * alpha (1, 2) * alpha (2, 3) * k0 [3];
            pmat (0, 0) = k0 [0] / tempd;
            pmat (0, 1) = alpha (0, 1) * k0 [1] / tempd;
            pmat (0, 2) = alpha (0, 1) * alpha (1, 2) * k0 [2] / tempd;
            pmat (0, 3) = alpha (0, 1) * alpha (1, 2) * alpha (2, 3) * 
                            k0 [3] / tempd;
        } else {
            tempd = k0 [0] + alpha (0, 1) * k0 [1] + alpha (0, 2) * k0 [2] + 
                    alpha (0, 2) * alpha (2, 3) * k0 [3];
            pmat (0, 0) = k0 [0] / tempd;
            pmat (0, 1) = alpha (0, 1) * k0 [1] / tempd;
            pmat (0, 2) = alpha (0, 2) * k0 [2] / tempd;
            pmat (0, 3) = alpha (0, 2) * alpha (2, 3) * k0 [3] / tempd;
        }

        /*
         * 	     3
         * 	     |
         * 	     2
         * 	    / \
         * 	   0---1
         *
         ********* node 1 **********/
        if ((alpha (1, 2) * alpha (2, 0)) > alpha (1, 0)) 
        {
            tempd = k0 [1] + alpha (1, 2) * alpha (2, 0) * k0 [0] + 
                    alpha (1, 2) * k0 [2] +
                    alpha (1, 2) * alpha (2, 3) * k0 [3];
            pmat (1, 0) = alpha (1, 2) * alpha (2, 0) * k0 [0] / tempd;
            pmat (1, 1) = k0 [1] / tempd;
            pmat (1, 2) = alpha (1, 2) * k0 [2] / tempd;
            pmat (1, 3) = alpha (1, 2) * alpha (2, 3) * k0 [3] / tempd;
        } else if ((alpha (1, 0) * alpha (0, 2)) > alpha (1, 2)) {
            tempd = k0 [1] + alpha (1, 0) * k0 [0] + 
                    alpha (1, 0) * alpha (0, 2) * k0 [2] +
                    alpha (1, 0) * alpha (0, 2) * alpha (2, 3) * k0 [3];
            pmat (1, 0) = alpha (1, 0) * k0 [0] / tempd;
            pmat (1, 1) = k0 [1] / tempd; 
            pmat (1, 2) = alpha (1, 0) * alpha (0, 2) * k0 [2] / tempd;
            pmat (1, 3) = alpha (1, 0) * alpha (0, 2) * alpha (2, 3) * 
                            k0 [3] / tempd;
        } else {
            tempd = k0 [1] + alpha (1, 0) * k0 [0] + alpha (1, 2) * k0 [2] + 
                    alpha (1, 2) * alpha (2, 3) * k0 [3];
            pmat (1, 0) = alpha (1, 0) * k0 [0] / tempd;
            pmat (1, 1) = k0 [1] / tempd;
            pmat (1, 2) = alpha (1, 2) * k0 [2] / tempd;
            pmat (1, 3) = alpha (1, 2) * alpha (2, 3) * k0 [3] / tempd;
        }

        /*
         * 	     3
         * 	     |
         * 	     2
         * 	    / \
         * 	   0---1
         *
         ********* nodes 2 & 3 **********/
        if ((alpha (2, 1) * alpha (1, 0)) > alpha (2, 0)) 
        {
            tempd = k0 [2] + alpha (2, 1) * k0 [1] + 
                    alpha (2, 1) * alpha (1, 0) * k0 [0] + 
                    alpha (2, 3) * k0 [3];
            pmat (2, 0) = alpha (2, 1) * alpha (1, 0) * k0 [0] / tempd;
            pmat (2, 1) = alpha (2, 1) * k0 [1] / tempd;
            pmat (2, 2) = k0 [2] / tempd;
            pmat (2, 3) = alpha (2, 3) * k0 [3] / tempd;
            tempd = k0 [3] + alpha (3, 2) * k0 [2] + 
                    alpha (3, 2) * alpha (2, 1) * k0 [1] +
                    alpha (3, 2) * alpha (2, 1) * alpha (1, 0) * k0 [0];
            pmat (3, 0) = alpha (3, 2) * alpha (2, 1) * alpha (1, 0) * 
                            k0 [0] / tempd;
            pmat (3, 1) = alpha (3, 2) * alpha (2, 1) * k0 [1] / tempd;
            pmat (3, 2) = alpha (3, 2) * k0 [2] / tempd;
            pmat (3, 3) = k0 [3] / tempd;
        } else if ((alpha (2, 0) * alpha (0, 1)) > alpha (2, 1)) {
            tempd = k0 [2] + alpha (2, 0) * k0 [0] + 
                    alpha (2, 0) * alpha (0, 1) * k0 [1] + 
                    alpha (2, 3) * k0 [3];
            pmat (2, 0) = alpha (2, 0) * k0 [0] / tempd;
            pmat (2, 1) = alpha (2, 0) * alpha (0, 1) * k0 [1] / tempd;
            pmat (2, 2) = k0 [2] / tempd;
            pmat (2, 3) = alpha (2, 3) * k0 [3] / tempd;
            tempd = k0 [3] + alpha (3, 2) * k0 [2] + 
                    alpha (3, 2) * alpha (2, 0) * k0 [0] +
                    alpha (3, 2) * alpha (2, 0) * alpha (0, 1) * k0 [1];
            pmat (3, 0) = alpha (3, 2) * alpha (2, 0) * k0 [0] / tempd;
            pmat (3, 1) = alpha (3, 2) * alpha (2, 0) * alpha (0, 1) * 
                            k0 [1] / tempd;
            pmat (3, 2) = alpha (3, 2) * k0 [2] / tempd;
            pmat (3, 3) = k0 [3] / tempd;
        } else {
            tempd = k0 [2] + alpha (2, 0) * k0 [0] + alpha (2, 1) * k0 [1] + 
                    alpha (2, 3) * k0 [3];
            pmat (2, 0) = alpha (2, 0) * k0 [0] / tempd;
            pmat (2, 1) = alpha (2, 1) * k0 [1] / tempd;
            pmat (2, 2) = k0 [2] / tempd;
            pmat (2, 3) = alpha (2, 3) * k0 [3] / tempd;
            tempd = k0 [3] + alpha (3, 2) * k0 [2] + 
                    alpha (3, 2) * alpha (2, 0) * k0 [0] + 
                    alpha (3, 2) * alpha (2, 1) * k0 [1];
            pmat (3, 0) = alpha (3, 2) * alpha (2, 0) * k0 [0] / tempd;
            pmat (3, 1) = alpha (3, 2) * alpha (2, 1) * k0 [1] / tempd;
            pmat (3, 2) = alpha (3, 2) * k0 [2] / tempd;
            pmat (3, 3) = k0 [3] / tempd;
        }
    } // end if network type == 3

    // Calculate connectivity from full pmat values
    results1.connectivity = (double) nnodes;
    for (int i=0; i<nnodes; i++) 
        results1.connectivity -= pmat (i, i);
    results1.connectivity = results1.connectivity / (double) nnodes;
    // Then rescale to growth rate 
    double pscale = pars.r;
    for (int i=0; i<nnodes; i++) 
        pmat (i, i) = 1.0 - pscale * (1.0 - pmat (i, i));
    for (int i=0; i<(nnodes - 1); i++) 
        for (int j=(i + 1); j<nnodes; j++) 
        {
            pmat (i, j) = pscale * pmat (i, j);
            pmat (j, i) = pscale * pmat (j, i);
        }
}
