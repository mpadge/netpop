/*
 * diff_sim4.cc
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
    int netcount, nRepeats;
    double tempd, progress, conn_mn, nmn_node [5], nsd_node [5], nmn_net, nsd_net;
    Parameters pars;
    Results results;
    std::string fname;
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
                (&pars.nTrials)->default_value (1000), "Number of trials")
            ("nRepeats,n", boost::program_options::value <int>
                (&nRepeats)->default_value (1000), "Number of repeats")
            ("k0sd,s", boost::program_options::value <double>
                (&pars.k0sd)->default_value (0.0), "SD of k0")
            ("alphasd,a", boost::program_options::value <double>
                (&pars.alphasd)->default_value (0.0), "SD of alpha")
            ("ksd,k", boost::program_options::value <double>
                (&pars.ksd)->default_value (0.1), "SD of k")
            ("r,r", boost::program_options::value <double>
                (&pars.r)->default_value (0.1), "Growth rate, r")
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
            std::cout << "diff-sim4, version 1.0" << std::endl;
            return 0;
        }

    }
    catch(std::exception& e)
    {
        std::cout << e.what() << std::endl;
        return 1;
    }    

    std::cout << "nTrials for each alpha = " << pars.nTrials <<
        "; with results averaged over " << nRepeats << " repeats." << std::endl;
    std::cout << "k0sd = " << pars.k0sd << "; alphasd = " << pars.alphasd << 
        "; ksd = " << pars.ksd << std::endl;

    time_start = clock ();

    //pars.r = 0.1;
    pars.alpha0 = 0.5;
    pars.k0 = 0.5;
    //pars.ksd = 0.1;

    time (&seed);
    generator.seed (static_cast <unsigned int> (seed));

    //cout.setf(0,ios::floatfield);            // floatfield not set
    std::cout.setf (std::ios::fixed, std::ios::floatfield);   
    std::cout.precision (4);

    std::stringstream ss;
    fname = "aaasim4_results_r0";
    if (pars.r < 0.1)
        fname += "0";
    ss.str (""); 
    //ss << floor (pars.r * 10.0); // presumes 0.1 <= r < 1
    ss << floor (pars.r * 100.0); 
    fname += ss.str () + "_ksd";
    if (pars.ksd < 0.01) 
        fname += "0";
    if (pars.ksd < 0.1) 
        fname += "0";
    ss.str (""); ss << floor (1000.0 * pars.ksd);
    fname += ss.str () + "_k0sd";
    if (pars.k0sd < 0.01) 
        fname += "0";
    if (pars.k0sd < 0.1) 
        fname += "0";
    ss.str (""); ss << floor (1000.0 * pars.k0sd);
    fname += ss.str() + "_alphasd";
    if (pars.alphasd < 0.01) 
        fname += "0";
    if (pars.alphasd < 0.1) 
        fname += "0";
    ss.str (""); ss << floor (1000.0 * pars.alphasd);
    fname += ss.str() + ".txt";

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
        pars.alpha0 = (double) i / 100.0;
        out_file << pars.alpha0;
        for (int j=0; j<4; j++) 
        {
            for (int k=0; k<5; k++) 
            {
                nmn_node [k] = 0.0;
                nsd_node [k] = 0.0;
            }
            nmn_net = 0.0;
            nsd_net = 0.0;
            conn_mn = 0.0;
            netcount = 0;
            for (int k=0; k<nRepeats; k++) 
            {
                results = runPop (pars, j, &generator);
                if (results.nmn_network > dnix) 
                {
                    netcount++;
                    conn_mn += results.connectivity;
                    for (int m=0; m<5; m++) 
                    {
                        nmn_node [m] += results.nmn_node [m];
                        nsd_node [m] += results.nsd_node [m];
                    }
                    nmn_net += results.nmn_network;
                    nsd_net += results.nsd_network;
                }
            }
            if (netcount > 0) 
            {
                conn_mn = conn_mn / (double) netcount;
                for (int k=0; k<5; k++) 
                {
                    nmn_node [k] = nmn_node [k] / (double) netcount;
                    nsd_node [k] = nsd_node [k] / (double) netcount;
                }
                nmn_net = nmn_net / (double) netcount;
                nsd_net = nsd_net / (double) netcount;
            } else {
                conn_mn = dnix;
                for (int k=0; k<5; k++) 
                {
                    nmn_node [k] = dnix;
                    nsd_node [k] = dnix;
                }
                nmn_net = dnix;
                nsd_net = dnix;
            }
            out_file << ",\t" << conn_mn;
            for (int k=0; k<5; k++) 
                out_file << ",\t" << nmn_node [k];
            for (int k=0; k<5; k++) 
                out_file << ",\t" << nsd_node [k];
            out_file << ",\t" << nmn_net << ",\t" << nsd_net;
            for (int k=0; k<3; k++) 
                for (int m=(k+1); m<4; m++)
                    out_file << ",\t" << results.cov [k] [m];
        } // end for j
        out_file << std::endl;
        //std::cout << "."; std::cout.flush();

        progress = 100.0 * (double) i / 100.0;
        std::cout << "\r[" << i << "]; progress = " << progress << "% after ";
        tempd = ((double) clock () - (double) time_start) / 
            (double) CLOCKS_PER_SEC;
        timeout (tempd);
        tempd = (tempd / progress) * (100.0 - progress);
        std::cout << "; remaining time = ";
        timeout (tempd);
        //std::cout << std::endl;
        std::cout.flush();
    } // end for i
    out_file.close();
    std::cout << std::endl << std::endl;

    return 0;
} // end main


/************************************************************************
 ************************************************************************
 **                                                                    **
 **                         SUBFUNCTIONS                               **
 **                                                                    **
 ************************************************************************
 ************************************************************************/

Results runPop (Parameters pars, int network_type, base_generator_type * generator)
{
    const int maxTrials = 100;
    int count, tempi;
    double tempd, n [4], nold [4], k0 [4], alpha [4] [4], kvals [4];
    bool flag, bigflag;
    Results results;

    boost::normal_distribution<> norm_dist (0.0, 1.0);
    boost::variate_generator<base_generator_type&,
        boost::normal_distribution<> > rnorm((*generator), norm_dist);
    // Burn generator in
    for (int i=0; i<20; i++) 
        n [0] = rnorm();

    for (int i=0; i<4; i++) 
    {
        k0 [i] = pars.k0 + pars.k0sd * rnorm ();
        while (k0 [i] < minqc || k0 [i] > (1.0 - minqc))
            k0 [i] = pars.k0 + pars.k0sd * rnorm ();
        for (int j=0; j<4; j++)
            alpha [i] [j] = 0.0;
    }

    // First set up connectivities
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
        alpha [itr -> first] [itr -> second] = tempd;
    } // end for itr
    for (int i=0; i<4; i++)
        alpha [i] [i] = 1.0;


    /********************************************************
     *********************   MAKEPMAT   *********************
     ********************************************************/

    // Doing this with the full makepmat routine using the
    // shortest path algorithm makes this really enormously
    // slower than the following manual construction.
    boost::numeric::ublas::matrix<double> pmat (4, 4);
    if (network_type == 0) 
    {
        /*
         * 	0---1---2---3
         */
        tempd = k0 [0] + alpha [0] [1] * k0 [1] + 
                alpha [0] [1] * alpha [1] [2] * k0 [2] +
                alpha [0] [1] * alpha [1] [2] * alpha [2] [3] * k0 [3];
        pmat (0, 0) = k0 [0] / tempd;
        pmat (0, 1) = alpha [0] [1] * k0 [1] / tempd;
        pmat (0, 2) = alpha [0] [1] * alpha [1] [2] * k0 [2] / tempd;
        pmat (0, 3) = alpha [0] [1] * alpha [1] [2] * alpha [2] [3] * 
                        k0 [3] / tempd;

        tempd = k0 [1] + alpha [1] [0] * k0 [0] + 
                alpha [1] [2] * k0 [2] + alpha [1] [2] * alpha [2] [3] * k0 [3];
        pmat (1, 0) = alpha [1] [0] * k0 [0] / tempd;
        pmat (1, 1) = k0 [1] / tempd;
        pmat (1, 2) = alpha [1] [2] * k0 [2] / tempd;
        pmat (1, 3) = alpha [1] [2] * alpha [2] [3] * k0 [3] / tempd;

        tempd = k0 [2] + alpha [2] [1] * k0 [1] + 
            alpha [2] [1] * alpha [1] [0] * k0 [0] + alpha [2] [3] * k0 [3];
        pmat (2, 0) = alpha [2] [1] * alpha [1] [0] * k0 [0] / tempd;
        pmat (2, 1) = alpha [2] [1] * k0 [1] / tempd;
        pmat (2, 2) = k0 [2] / tempd;
        pmat (2, 3) = alpha [2] [3] * k0 [3] / tempd;

        tempd = k0 [3] + alpha [3] [2] * k0 [2] + 
                alpha [3] [2] * alpha [2] [1] * k0 [1] +
                alpha [3] [2] * alpha [2] [1] * alpha [1] [0] * k0 [0];
        pmat (3, 0) = alpha [3] [2] * alpha [2] [1] * alpha [1] [0] * 
                        k0 [0] / tempd;
        pmat (3, 1) = alpha [3] [2] * alpha [2] [1] * k0 [1] / tempd;
        pmat (3, 2) = alpha [3] [2] * k0 [2] / tempd;
        pmat (3, 3) = k0 [3] / tempd;
    } else if (network_type == 1) {
        /*
         *	    3
         *	    |
         *	    2
         *	   / \
         *	  0   1
         */
        tempd = k0 [0] + alpha [0] [2] * k0 [2] + 
                alpha [0] [2] * alpha [2] [1] * k0 [1] +
                alpha [0] [2] * alpha [2] [3] * k0 [3];
        pmat (0, 0) = k0 [0] / tempd;
        pmat (0, 1) = alpha [0] [2] * alpha [2] [1] * k0 [1] / tempd;
        pmat (0, 2) = alpha [0] [2] * k0 [2] / tempd;
        pmat (0, 3) = alpha [0] [2] * alpha [2] [3] * k0 [3] / tempd;

        tempd = k0 [1] + alpha [1] [2] * k0 [2] + 
                alpha [1] [2] * alpha [2] [0] * k0 [0] +
                alpha [1] [2] * alpha [2] [3] * k0 [3];
        pmat (1, 0) = alpha [1] [2] * alpha [2] [0] * k0 [0] / tempd;
        pmat (1, 1) = k0 [1] / tempd;
        pmat (1, 2) = alpha [1] [2] * k0 [2] / tempd;
        pmat (1, 3) = alpha [1] [2] * alpha [2] [3] * k0 [3] / tempd;

        tempd = k0 [2] + alpha [2] [0] * k0 [0] + 
                alpha [2] [1] * k0 [1] + alpha [2] [3] * k0 [3];
        pmat (2, 0) = alpha [2] [0] * k0 [0] / tempd;
        pmat (2, 1) = alpha [2] [1] * k0 [1] / tempd;
        pmat (2, 2) = k0 [2] / tempd;
        pmat (2, 3) = alpha [2] [3] * k0 [3] / tempd;

        tempd = k0 [3] + alpha [3] [2] * k0 [2] + 
                alpha [3] [2] * alpha [2] [1] * k0 [1] +
                alpha [3] [2] * alpha [2] [0] * k0 [0];
        pmat (3, 0) = alpha [3] [2] * alpha [2] [0] * k0 [0] / tempd;
        pmat (3, 1) = alpha [3] [2] * alpha [2] [1] * k0 [1] / tempd;
        pmat (3, 2) = alpha [3] [2] * k0 [2] / tempd;
        pmat (3, 3) = k0 [3] / tempd;
    } else if (network_type == 2) {
        /*
         * 	2---3
         * 	|   |
         * 	0---1
         *
         ********* node 0 **********/
        if ((alpha [0] [2] * alpha [2] [3] * alpha [3] [1]) > alpha [0] [1]) 
        {
            tempd = k0 [0] + alpha [0] [2] * k0 [2] + 
                    alpha [0] [2] * alpha [2] [3] * k0 [3] +
                    alpha [0] [2] * alpha [2] [3] * alpha [3] [1] * k0 [1];
            pmat (0, 0) = k0 [0] / tempd;
            pmat (0, 1) = alpha [0] [2] * alpha [2] [3] * alpha [3] [1] * 
                            k0 [1] / tempd;
            pmat (0, 2) = alpha [0] [2] * k0 [2] / tempd;
            pmat (0, 3) = alpha [0] [2] * alpha [2] [3] * k0 [3] / tempd;
        } else if ((alpha [0] [1] * alpha [1] [3] * alpha [3] [2]) > 
                alpha [0] [2]) {
            tempd = k0 [0] + alpha [0] [1] * k0 [1] + 
                    alpha [0] [1] * alpha [1] [3] * k0 [3] +
                    alpha [0] [1] * alpha [1] [3] * alpha [3] [2] * k0 [2];
            pmat (0, 0) = k0 [0] / tempd;
            pmat (0, 1) = alpha [0] [1] * k0 [1] / tempd;
            pmat (0, 2) = alpha [0] [1] * alpha [1] [3] * alpha [3] [2] * 
                            k0 [2] / tempd;
            pmat (0, 3) = alpha [0] [1] * alpha [1] [3] * k0 [3] / tempd;
        } else if ((alpha [0] [2] * alpha [2] [3]) > 
                (alpha [0] [1] * alpha [1] [3])) {
            tempd = k0 [0] + alpha [0] [1] * k0 [1] + alpha [0] [2] * k0 [2] + 
                    alpha [0] [2] * alpha [2] [3] * k0 [3];
            pmat (0, 0) = k0 [0] / tempd;
            pmat (0, 1) = alpha [0] [1] * k0 [1] / tempd;
            pmat (0, 2) = alpha [0] [2] * k0 [2] / tempd;
            pmat (0, 3) = alpha [0] [2] * alpha [2] [3] * k0 [3] / tempd;
        } else {
            tempd = k0 [0] + alpha [0] [1] * k0 [1] + 
                    alpha [0] [2] * k0 [2] + 
                    alpha [0] [1] * alpha [1] [3] * k0 [3];
            pmat (0, 0) = k0 [0] / tempd;
            pmat (0, 1) = alpha [0] [1] * k0 [1] / tempd;
            pmat (0, 2) = alpha [0] [2] * k0 [2] / tempd;
            pmat (0, 3) = alpha [0] [1] * alpha [1] [3] * k0 [3] / tempd;
        }

        /*
         * 	2---3
         * 	|   |
         * 	0---1
         *
         ********* node 1 **********/
        if ((alpha [1] [3] * alpha [3] [2] * alpha [2] [0]) > alpha [1] [0]) 
        {
            tempd = k0 [1] + alpha [1] [3] * k0 [3] + 
                    alpha [1] [3] * alpha [3] [2] * k0 [2] +
                    alpha [1] [3] * alpha [3] [2] * alpha [2] [0] * k0 [0];
            pmat (1, 0) = alpha [1] [3] * alpha [3] [2] * alpha [2] [0] * 
                            k0 [0] / tempd;
            pmat (1, 1) = k0 [1] / tempd;
            pmat (1, 2) = alpha [1] [3] * alpha [3] [2] * k0 [2] / tempd;
            pmat (1, 3) = alpha [1] [3] * k0 [3] / tempd;
        } else if ((alpha [1] [0] * alpha [0] [2] * alpha [2] [3]) > 
                alpha [1] [3]) {
            tempd = k0 [1] + alpha [1] [0] * k0 [0] + 
                    alpha [1] [0] * alpha [0] [2] * k0 [2] +
                    alpha [1] [0] * alpha [0] [2] * alpha [2] [3] * k0 [3];
            pmat (1, 0) = alpha [1] [0] * k0 [0] / tempd;
            pmat (1, 1) = k0 [1] / tempd;
            pmat (1, 2) = alpha [1] [0] * alpha [0] [2] * k0 [2] / tempd;
            pmat (1, 3) = alpha [1] [0] * alpha [0] [2] * alpha [2] [3] * 
                k0 [3] / tempd;
        } else if ((alpha [1] [3] * alpha [3] [2]) > 
                (alpha [1] [0] * alpha [0] [2])) {
            tempd = k0 [1] + alpha [1] [3] * k0 [3] + alpha [1] [0] * k0 [0] + 
                    alpha [1] [3] * alpha [3] [2] * k0 [2];
            pmat (1, 0) = alpha [1] [0] * k0 [0] / tempd;
            pmat (1, 1) = k0 [1] / tempd;
            pmat (1, 2) = alpha [1] [3] * alpha [3] [2] * k0 [2] / tempd;
            pmat (1, 3) = alpha [1] [3] * k0 [3] / tempd;
        } else {
            tempd = k0 [1] + alpha [1] [3] * k0 [3] + alpha [1] [0] * k0 [0] + 
                    alpha [1] [0] * alpha [0] [2] * k0 [2];
            pmat (1, 0) = alpha [1] [0] * k0 [0] / tempd;
            pmat (1, 1) = k0 [1] / tempd;
            pmat (1, 2) = alpha [1] [0] * alpha [0] [2] * k0 [2] / tempd;
            pmat (1, 3) = alpha [1] [3] * k0 [3] / tempd;
        }

        /*
         * 	2---3
         * 	|   |
         * 	0---1
         *
         ********* node 2 **********/
        if ((alpha [2] [0] * alpha [0] [1] * alpha [1] [3]) > alpha [2] [3]) 
        {
            tempd = k0 [2] + alpha [2] [0] * k0 [0] + 
                    alpha [2] [0] * alpha [0] [1] * k0 [1] +
                    alpha [2] [0] * alpha [0] [1] * alpha [1] [3] * k0 [3];
            pmat (2, 0) = alpha [2] [0] * k0 [0] / tempd;
            pmat (2, 1) = alpha [2] [0] * alpha [0] [1] * k0 [1] / tempd;
            pmat (2, 2) = k0 [2] / tempd;
            pmat (2, 3) = alpha [2] [0] * alpha [0] [1] * alpha [1] [3] * 
                            k0 [3] / tempd;
        } else if ((alpha [2] [3] * alpha [3] [1] * alpha [1] [0]) > 
                alpha [2] [0]) {
            tempd = k0 [2] + alpha [2] [3] * k0 [3] + 
                    alpha [2] [3] * alpha [3] [1] * k0 [1] +
                    alpha [2] [3] * alpha [3] [1] * alpha [1] [0] * k0 [0];
            pmat (2, 0) = alpha [2] [3] * alpha [3] [1] * alpha [1] [0] * 
                k0 [0] / tempd;
            pmat (2, 1) = alpha [2] [3] * alpha [3] [1] * k0 [1] / tempd;
            pmat (2, 2) = k0 [2] / tempd;
            pmat (2, 3) = alpha [2] [3] * k0 [3] / tempd;
        } else if ((alpha [2] [3] * alpha [3] [1]) > 
                (alpha [2] [0] * alpha [0] [1])) {
            tempd = k0 [2] + alpha [2] [0] * k0 [0] + alpha [2] [3] * k0 [3] + 
                    alpha [2] [3] * alpha [3] [1] * k0 [1];
            pmat (2, 0) = alpha [2] [0] * k0 [0] / tempd;
            pmat (2, 1) = alpha [2] [3] * alpha [3] [1] * k0 [1] / tempd;
            pmat (2, 2) = k0 [2] / tempd;
            pmat (2, 3) = alpha [2] [3] * k0 [3] / tempd;
        } else {
            tempd = k0 [2] + alpha [2] [0] * k0 [0] + alpha [2] [3] * k0 [3] + 
                    alpha [2] [0] * alpha [0] [1] * k0 [1];
            pmat (2, 0) = alpha [2] [0] * k0 [0] / tempd;
            pmat (2, 1) = alpha [2] [0] * alpha [0] [1] * k0 [1] / tempd;
            pmat (2, 2) = k0 [2] / tempd;
            pmat (2, 3) = alpha [2] [3] * k0 [3] / tempd;
        }

        /*
         * 	2---3
         * 	|   |
         * 	0---1
         *
         ********* node 3 **********/
        if ((alpha [3] [2] * alpha [2] [0] * alpha [0] [1]) > alpha [3] [1]) 
        {
            tempd = k0 [3] + alpha [3] [2] * k0 [2] + 
                    alpha [3] [2] * alpha [2] [0] * k0 [0] +
                    alpha [3] [2] * alpha [2] [0] * alpha [0] [1] * k0 [1];
            pmat (3, 0) = alpha [3] [2] * alpha [2] [0] * k0 [0] / tempd;
            pmat (3, 1) = alpha [3] [2] * alpha [2] [0] * alpha [0] [1] * 
                            k0 [1] / tempd;
            pmat (3, 2) = alpha [3] [2] * k0 [2] / tempd;
            pmat (3, 3) = k0 [3] / tempd;
        } else if ((alpha [3] [1] * alpha [1] [0] * alpha [0] [2]) > 
                alpha [3] [2]) {
            tempd = k0 [3] + alpha [3] [1] * k0 [1] + 
                    alpha [3] [1] * alpha [1] [0] * k0 [0] +
                    alpha [3] [1] * alpha [1] [0] * alpha [0] [2] * k0 [2];
            pmat (3, 0) = alpha [3] [1] * alpha [1] [0] * k0 [0] / tempd;
            pmat (3, 1) = alpha [3] [1] * k0 [1] / tempd;
            pmat (3, 2) = alpha [3] [1] * alpha [1] [0] * alpha [0] [2] * 
                            k0 [2] / tempd;
            pmat (3, 3) = k0 [3] / tempd;
        } else if ((alpha [3] [2] * alpha [2] [0]) > 
                    (alpha [3] [1] * alpha [1] [0])) {
            tempd = k0 [3] + alpha [3] [1] * k0 [1] + alpha [3] [2] * k0 [2] + 
                    alpha [3] [2] * alpha [2] [0] * k0 [0];
            pmat (3, 0) = alpha [3] [2] * alpha [2] [0] * k0 [0] / tempd;
            pmat (3, 1) = alpha [3] [1] * k0 [1] / tempd;
            pmat (3, 2) = alpha [3] [2] * k0 [2] / tempd;
            pmat (3, 3) = k0 [3] / tempd;
        } else {
            tempd = k0 [3] + alpha [3] [1] * k0 [1] + alpha [3] [2] * k0 [2] + 
                    alpha [3] [1] * alpha [1] [0] * k0 [0];
            pmat (3, 0) = alpha [3] [1] * alpha [1] [0] * k0 [0] / tempd;
            pmat (3, 1) = alpha [3] [1] * k0 [1] / tempd;
            pmat (3, 2) = alpha [3] [2] * k0 [2] / tempd;
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
        if ((alpha [0] [2] * alpha [2] [1]) > alpha [0] [1]) 
        {
            tempd = k0 [0] + alpha [0] [2] * alpha [2] [1] * k0 [1] + 
                    alpha [0] [2] * k0 [2] + 
                    alpha [0] [2] * alpha [2] [3] * k0 [3];
            pmat (0, 0) = k0 [0] / tempd;
            pmat (0, 1) = alpha [0] [2] * alpha [2] [1] * k0 [1] / tempd;
            pmat (0, 2) = alpha [0] [2] * k0 [2] / tempd;
            pmat (0, 3) = alpha [0] [2] * alpha [2] [3] * k0 [3] / tempd;
        } else if ((alpha [0] [1] * alpha [1] [2]) > alpha [0] [2]) {
            tempd = k0 [0] + alpha [0] [1] * k0 [1] + 
                    alpha [0] [1] * alpha [1] [2] * k0 [2] +
                    alpha [0] [1] * alpha [1] [2] * alpha [2] [3] * k0 [3];
            pmat (0, 0) = k0 [0] / tempd;
            pmat (0, 1) = alpha [0] [1] * k0 [1] / tempd;
            pmat (0, 2) = alpha [0] [1] * alpha [1] [2] * k0 [2] / tempd;
            pmat (0, 3) = alpha [0] [1] * alpha [1] [2] * alpha [2] [3] * 
                            k0 [3] / tempd;
        } else {
            tempd = k0 [0] + alpha [0] [1] * k0 [1] + alpha [0] [2] * k0 [2] + 
                    alpha [0] [2] * alpha [2] [3] * k0 [3];
            pmat (0, 0) = k0 [0] / tempd;
            pmat (0, 1) = alpha [0] [1] * k0 [1] / tempd;
            pmat (0, 2) = alpha [0] [2] * k0 [2] / tempd;
            pmat (0, 3) = alpha [0] [2] * alpha [2] [3] * k0 [3] / tempd;
        }

        /*
         * 	     3
         * 	     |
         * 	     2
         * 	    / \
         * 	   0---1
         *
         ********* node 1 **********/
        if ((alpha [1] [2] * alpha [2] [0]) > alpha [1] [0]) 
        {
            tempd = k0 [1] + alpha [1] [2] * alpha [2] [0] * k0 [0] + 
                    alpha [1] [2] * k0 [2] +
                    alpha [1] [2] * alpha [2] [3] * k0 [3];
            pmat (1, 0) = alpha [1] [2] * alpha [2] [0] * k0 [0] / tempd;
            pmat (1, 1) = k0 [1] / tempd;
            pmat (1, 2) = alpha [1] [2] * k0 [2] / tempd;
            pmat (1, 3) = alpha [1] [2] * alpha [2] [3] * k0 [3] / tempd;
        } else if ((alpha [1] [0] * alpha [0] [2]) > alpha [1] [2]) {
            tempd = k0 [1] + alpha [1] [0] * k0 [0] + 
                    alpha [1] [0] * alpha [0] [2] * k0 [2] +
                    alpha [1] [0] * alpha [0] [2] * alpha [2] [3] * k0 [3];
            pmat (1, 0) = alpha [1] [0] * k0 [0] / tempd;
            pmat (1, 1) = k0 [1] / tempd; 
            pmat (1, 2) = alpha [1] [0] * alpha [0] [2] * k0 [2] / tempd;
            pmat (1, 3) = alpha [1] [0] * alpha [0] [2] * alpha [2] [3] * 
                            k0 [3] / tempd;
        } else {
            tempd = k0 [1] + alpha [1] [0] * k0 [0] + alpha [1] [2] * k0 [2] + 
                    alpha [1] [2] * alpha [2] [3] * k0 [3];
            pmat (1, 0) = alpha [1] [0] * k0 [0] / tempd;
            pmat (1, 1) = k0 [1] / tempd;
            pmat (1, 2) = alpha [1] [2] * k0 [2] / tempd;
            pmat (1, 3) = alpha [1] [2] * alpha [2] [3] * k0 [3] / tempd;
        }

        /*
         * 	     3
         * 	     |
         * 	     2
         * 	    / \
         * 	   0---1
         *
         ********* nodes 2 & 3 **********/
        if ((alpha [2] [1] * alpha [1] [0]) > alpha [2] [0]) 
        {
            tempd = k0 [2] + alpha [2] [1] * k0 [1] + 
                    alpha [2] [1] * alpha [1] [0] * k0 [0] + 
                    alpha [2] [3] * k0 [3];
            pmat (2, 0) = alpha [2] [1] * alpha [1] [0] * k0 [0] / tempd;
            pmat (2, 1) = alpha [2] [1] * k0 [1] / tempd;
            pmat (2, 2) = k0 [2] / tempd;
            pmat (2, 3) = alpha [2] [3] * k0 [3] / tempd;
            tempd = k0 [3] + alpha [3] [2] * k0 [2] + 
                    alpha [3] [2] * alpha [2] [1] * k0 [1] +
                    alpha [3] [2] * alpha [2] [1] * alpha [1] [0] * k0 [0];
            pmat (3, 0) = alpha [3] [2] * alpha [2] [1] * alpha [1] [0] * 
                            k0 [0] / tempd;
            pmat (3, 1) = alpha [3] [2] * alpha [2] [1] * k0 [1] / tempd;
            pmat (3, 2) = alpha [3] [2] * k0 [2] / tempd;
            pmat (3, 3) = k0 [3] / tempd;
        } else if ((alpha [2] [0] * alpha [0] [1]) > alpha [2] [1]) {
            tempd = k0 [2] + alpha [2] [0] * k0 [0] + 
                    alpha [2] [0] * alpha [0] [1] * k0 [1] + 
                    alpha [2] [3] * k0 [3];
            pmat (2, 0) = alpha [2] [0] * k0 [0] / tempd;
            pmat (2, 1) = alpha [2] [0] * alpha [0] [1] * k0 [1] / tempd;
            pmat (2, 2) = k0 [2] / tempd;
            pmat (2, 3) = alpha [2] [3] * k0 [3] / tempd;
            tempd = k0 [3] + alpha [3] [2] * k0 [2] + 
                    alpha [3] [2] * alpha [2] [0] * k0 [0] +
                    alpha [3] [2] * alpha [2] [0] * alpha [0] [1] * k0 [1];
            pmat (3, 0) = alpha [3] [2] * alpha [2] [0] * k0 [0] / tempd;
            pmat (3, 1) = alpha [3] [2] * alpha [2] [0] * alpha [0] [1] * 
                            k0 [1] / tempd;
            pmat (3, 2) = alpha [3] [2] * k0 [2] / tempd;
            pmat (3, 3) = k0 [3] / tempd;
        } else {
            tempd = k0 [2] + alpha [2] [0] * k0 [0] + alpha [2] [1] * k0 [1] + 
                    alpha [2] [3] * k0 [3];
            pmat (2, 0) = alpha [2] [0] * k0 [0] / tempd;
            pmat (2, 1) = alpha [2] [1] * k0 [1] / tempd;
            pmat (2, 2) = k0 [2] / tempd;
            pmat (2, 3) = alpha [2] [3] * k0 [3] / tempd;
            tempd = k0 [3] + alpha [3] [2] * k0 [2] + 
                    alpha [3] [2] * alpha [2] [0] * k0 [0] + 
                    alpha [3] [2] * alpha [2] [1] * k0 [1];
            pmat (3, 0) = alpha [3] [2] * alpha [2] [0] * k0 [0] / tempd;
            pmat (3, 1) = alpha [3] [2] * alpha [2] [1] * k0 [1] / tempd;
            pmat (3, 2) = alpha [3] [2] * k0 [2] / tempd;
            pmat (3, 3) = k0 [3] / tempd;
        }
    } // end if network type == 3

    // Calculate connectivity from full pmat values
    results.connectivity = 4.0;
    for (int i=0; i<4; i++) 
        results.connectivity -= pmat (i, i);
    results.connectivity = results.connectivity / 4.0;
    // Then rescale to growth rate (see makepmat in dend or latt routines for
    // details): NOTE that the inverse one uses s=0.31, but still divides pscale
    // by 5.  This is because the original version, which produced visually
    // clear results with distinct differences between networks, used pscale =
    // pars.r, and dividing r^-0.31 by 5 makes values again around 0.3-0.4.
    // Obviously such division is arbitrary, but just serves to dampen the
    // effect of movement, and so emphasise differences between networks.
    //double pscale = pow (pars.r, -0.31) / 5.0;
    double pscale = pars.r;
    for (int i=0; i<4; i++) 
        pmat (i, i) = 1.0 - pscale * (1.0 - pmat (i, i));
    for (int i=0; i<3; i++) 
    {
        for (int j=(i + 1); j<4; j++) 
        {
            pmat (i, j) = pscale * pmat (i, j);
            pmat (j, i) = pscale * pmat (j, i);
        }
    }

    /********************************************************
     ***************   POPULATION ITERATION   ***************
     ********************************************************/

    count = 0;
    flag = true;
    bigflag = false;
    while (flag) 
    {
        for (int j=0; j<4; j++) 
            nold [j] = k0 [j];
        for (int j=0; j<5; j++) 
        {
            results.nmn_node [j] = 0.0;
            results.nsd_node [j] = 0.0;
        }
        results.nmn_network = 0.0;
        results.nsd_network = 0.0;
        flag = false;
        for (int i=0; i<3; i++)
            for (int j=(i+1); j<4; j++)
                results.cov [i] [j] = 0.0;
        for (int i=0; i<(runin + pars.nTrials); i++) 
        {
            // Movement through network
            for (int j=0; j<4; j++) 
            {
                n [j] = 0.0;
                for (int k=0; k<4; k++) 
                    n [j] += pmat (k, j) * nold [k];
            }
            // Change carrying capacities
            for (int j=0; j<4; j++) 
            {
                kvals [j] = k0 [j] + pars.ksd * rnorm ();
                while (kvals [j] < minqc || kvals [j] > (1.0 - minqc))
                    kvals [j] = k0 [j] + pars.ksd * rnorm ();
            }
            // Population dynamic
            tempi = 0;
            for (int j=0; j<4; j++) 
            {
                nold [j] = n [j] + pars.r * n [j] * n [j] / kvals [j] - 
                        pars.r * n [j] * n [j] * n [j] / (kvals [j] * kvals [j]);
                if (nold [j] < 0.0) 
                {
                    //nold [j] = 0.0; // Set to zero after calculating mean & SD values
                    tempi++;
                }
            } // end for j
            if (tempi == 4) 
            {
                flag = true;
                break;
            } else if (i > runin) {
                tempd = 0.0;
                for (int j=0; j<4; j++) 
                {
                    tempd += nold [j];
                    results.nmn_node [j] += nold [j];
                    results.nsd_node [j] += nold [j] * nold [j];
                }
                results.nmn_network += tempd;
                results.nsd_network += tempd * tempd;
                for (int j=0; j<3; j++)
                    for (int k=(j+1); k<4; k++)
                        results.cov [j] [k] += nold [j] * nold [k];
            }
            // Then set all negative nodal abundances to zero
            for (int j=0; j<4; j++)
                if (nold [j] < 0.0)
                    nold [j] = 0.0;
        } // end for i over nTrials
        if (!flag) 
        {
            for (int j=0; j<4; j++) 
            {
                results.nmn_node [j] = results.nmn_node [j] / 
                                        (double) pars.nTrials;
                results.nsd_node [j] = results.nsd_node [j] / 
                                        (double) pars.nTrials -
                    results.nmn_node [j] * results.nmn_node [j];
            }
            for (int j=0; j<3; j++)
                for (int k=(j+1); k<4; k++)
                    results.cov [j] [k] = results.cov [j] [k] / 
                                        (double) pars.nTrials -
                        results.nmn_node [j] * results.nmn_node [k];
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
        for (int i=0; i<5; i++) 
        {
            results.nmn_node [i] = dnix;
            results.nsd_node [i] = dnix;
        }
        results.nmn_network = dnix;
        results.nsd_network = dnix;
    } else {
        results.nmn_node [4] = 0.0;
        results.nsd_node [4] = 0.0;
        for (int i=0; i<4; i++) 
        {
            results.nmn_node [4] += results.nmn_node [i];
            results.nsd_node [4] += results.nsd_node [i];
        }
        results.nmn_node [4] = results.nmn_node [4] / 4.0;
        results.nsd_node [4] = results.nsd_node [4] / 4.0;
        results.nmn_network = results.nmn_network  / (double) pars.nTrials;
        results.nsd_network = results.nsd_network / (double) pars.nTrials -
            results.nmn_network * results.nmn_network;
        results.nmn_network = results.nmn_network / 4.0; 
        // So it's on the same scale as nodal abundance.
    }

    return results;
}


void timeout(double tseconds)
{
    int hh = floor (tseconds / 3600.0);
    if (hh == 0) 
        std::cout << "00:";
    else if (hh < 10) 
        std::cout << "0" << hh << ":";
    else 
        std::cout << hh << ":";
    double trem = tseconds - (double) hh * 3600.0;
    int mm = floor (trem / 60.0);
    if (mm == 0) 
        std::cout << "00:";
    else if (mm < 10) 
        std::cout << "0" << mm << ":";
    else 
        std::cout << mm << ":";
    double ss = trem - (double) mm * 60.0;
    if (ss == 0.0) 
        std::cout << "00:";
    else if (ss < 10) 
        std::cout << "0" << ss;
    else 
        std::cout << ss;
} // end function timeout

/* R script to plot results
junk <- function (minalpha=0.1, nfiles=100)
{
setwd ("/data/Documents/analyses/Fish Barriers/c++/diffusion/results/")
dat0 <- read.csv ("aaasim4_results_basic_ksd02.txt", header=TRUE)
setwd ("/data/Documents/analyses/Fish Barriers/c++/diffusion/")
indx <- which (dat0$alpha >= minalpha)
dat0 <- dat0 [indx, ]
alpha <- dat0$alpha
nsd0 <- cbind (dat0$sd0net, dat0$sd1net, dat0$sd2net, dat0$sd3net)

lf <- list.files()

r2 <- array (NA, dim=c(nfiles, 4))
r2alpha <- r2
for (i in 1:nfiles) {
	if (i < 10) {
		fname <- paste ("aaasim4_results_ksd02_k0sd00", i, "_alphasd00", i, ".txt", sep="")
	}
	else if (i < 100) {
		fname <- paste ("aaasim4_results_ksd02_k0sd0", i, "_alphasd0", i, ".txt", sep="")
	}
	else {
		fname <- paste ("aaasim4_results_ksd02_k0sd", i, "_alphasd", i, ".txt", sep="")
	}
	if (fname %in% lf) {
		dat <- read.csv (fname, header=TRUE)
		dat <- dat [indx, ]
		nsd <- cbind (dat$sd0net, dat$sd1net, dat$sd2net, dat$sd3net)
		for (j in 1:4) {
			r2 [i, j] <- cor (nsd0 [,j], nsd [,j])
			r2alpha [i, j] <- cor (alpha, nsd [,j])
			}
		maxi <- i
	}
	else {
		maxi <- i - 1
		break	}
}
r2 <- r2 [1:maxi, ]
r2alpha <- r2alpha [1:maxi, ]
r2 <- sign (r2) * r2 ^ 2
r2alpha <- sign (r2alpha) * r2alpha ^ 2

yvals <- list ()
yvals [[1]] <- r2
yvals [[2]] <- r2alpha
ylims <- list ()
ylims [[1]] <- range (yvals [[1]], na.rm=TRUE)
ylims [[2]] <- range (yvals [[2]], na.rm=TRUE)
mts <- c ("Correlation with nework pattern", "Correlation with alpha")

x11 (width = 10, height = 5)
par (mfrow = c(1, 2), mar = c(2, 2, 1.5, 0.5), mgp = c(1, 0.3, 0), ps=10)
x <- (1:maxi) / 1000
cols <- c ("red", "orange", "lawngreen", "blue")
for (i in 1:2) {
	plot (x, yvals [[i]] [,1], "l", col=cols [1], ylim=ylims [[i]],
		xlab="SD", ylab="R2", main=mts [i])
	lines (c(-1, 100), c(0, 0), col="grey", lty=2)
	for (j in 1:4) {
		lines (x, yvals [[i]] [,j], col=cols [j])
		#points (x, r2 [,i], pch=20, col=cols [i])
		}
	if (i == 1) {
		xpos <- min (x) + 0.7 * diff (range (x))
		legend (xpos, ylims [[i]] [2], lwd=1, col=cols, bty="n",
			legend=c("A: Line", "B: Star", "C: Circle", "D: Star+"))
		}
	}
}
junk()
*/
