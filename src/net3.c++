/*
 * diff_sim3.cc
 *
 * Simulates a triangular network of three nodes
 */

#include "net3.h"


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
    double tempd, progress, conn_mn, nmn_node [4], nsd_node [4], nmn_net, nsd_net;
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
                (&pars.k0sd)->default_value (0.1), "SD of k0")
            ("alphasd,a", boost::program_options::value <double>
                (&pars.alphasd)->default_value (0.1), "SD of alpha")
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
            std::cout << "diff-sim3, version 1.0" << std::endl;
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
    fname = "aaasim3_results_r";
    if (pars.r < 0.1)
        fname += "0";
    if (pars.r < 1.0)
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
    out_file << "alpha,\tconn,\tmn0,\tmn1,\tmn2,\tmnall\t,sd0,\tsd1,\tsd2,\t" <<
        "sdall,\tmnnet,\tsdnet,\tcov01,\tcov02,\tcov12" << std::endl;
    for (int i=1; i<=100; i++) 
    {
        pars.alpha0 = (double) i / 100.0;
        out_file << pars.alpha0;
        for (int k=0; k<4; k++) 
        {
            nmn_node [k] = 0.0;
            nsd_node [k] = 0.0;
        }
        nmn_net = 0.0;
        nsd_net = 0.0;
        conn_mn = 0.0;
        netcount = 0;
        for (int j=0; j<nRepeats; j++) 
        {
            results = runPop (pars, &generator);
            if (results.nmn_network > DOUBLE_MIN) 
            {
                netcount++;
                conn_mn += results.connectivity;
                for (int k=0; k<4; k++) 
                {
                    nmn_node [k] += results.nmn_node [k];
                    nsd_node [k] += results.nsd_node [k];
                }
                nmn_net += results.nmn_network;
                nsd_net += results.nsd_network;
            }
        } // end for j
        if (netcount > 0) 
        {
            conn_mn = conn_mn / (double) netcount;
            for (int j=0; j<4; j++) 
            {
                nmn_node [j] = nmn_node [j] / (double) netcount;
                nsd_node [j] = nsd_node [j] / (double) netcount;
            }
            nmn_net = nmn_net / (double) netcount;
            nsd_net = nsd_net / (double) netcount;
        } else {
            conn_mn = DOUBLE_MIN;
            for (int j=0; j<4; j++) 
            {
                nmn_node [j] = DOUBLE_MIN;
                nsd_node [j] = DOUBLE_MIN;
            }
            nmn_net = DOUBLE_MIN;
            nsd_net = DOUBLE_MIN;
        }
        out_file << ",\t" << conn_mn;
        for (int j=0; j<4; j++) 
            out_file << ",\t" << nmn_node [j];
        for (int j=0; j<4; j++) 
            out_file << ",\t" << nsd_node [j];
        out_file << ",\t" << nmn_net << ",\t" << nsd_net;
        for (int j=0; j<2; j++) 
            for (int k=(j+1); k<3; k++)
                out_file << ",\t" << results.cov [j] [k];
        out_file << std::endl;
        //std::cout << "."; std::cout.flush();

        /*
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
        */
        progress = (double) i / 100.0;
        progLine (progress, 0);
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

Results runPop (Parameters pars, base_generator_type * generator)
{
    const int maxTrials = 100;
    int count, tempi;
    double tempd, n [3], nold [3], k0 [3], alpha [3] [3], kvals [3];
    bool flag, bigflag;
    Results results;

    boost::normal_distribution<> norm_dist (0.0, 1.0);
    boost::variate_generator<base_generator_type&,
        boost::normal_distribution<> > rnorm((*generator), norm_dist);
    // Burn generator in
    for (int i=0; i<20; i++) 
        n [0] = rnorm();

    for (int i=0; i<3; i++) 
    {
        k0 [i] = pars.k0 + pars.k0sd * rnorm ();
        while (k0 [i] < minqc || k0 [i] > (1.0 - minqc))
            k0 [i] = pars.k0 + pars.k0sd * rnorm ();
        for (int j=0; j<3; j++)
            alpha [i] [j] = 0.0;
    }

    // First set up connectivities
    std::vector <std::pair <int, int> > connlist;
    connlist.push_back (std::pair <int, int> (0, 1));
    connlist.push_back (std::pair <int, int> (1, 0));
    connlist.push_back (std::pair <int, int> (0, 2));
    connlist.push_back (std::pair <int, int> (2, 0));
    connlist.push_back (std::pair <int, int> (1, 2));
    connlist.push_back (std::pair <int, int> (2, 1));
    
    std::vector <std::pair <int, int> >::const_iterator itr;
    for (itr = connlist.begin(); itr < connlist.end(); itr++) 
    {
        tempd = pars.alpha0 + pars.alphasd * rnorm ();
        while (tempd < minqc || tempd > (1.0 - minqc))
            tempd = pars.alpha0 + pars.alphasd * rnorm ();
        alpha [itr -> first] [itr -> second] = tempd;
    } // end for itr
    for (int i=0; i<3; i++)
        alpha [i] [i] = 1.0;


    /********************************************************
     *********************   MAKEPMAT   *********************
     ********************************************************/

    boost::numeric::ublas::matrix<double> pmat (3, 3);
    tempd = k0 [0] + alpha [0] [1] * k0 [1] + alpha [0] [2] * k0 [2];
    pmat (0, 0) = k0 [0] / tempd;
    pmat (0, 1) = alpha [0] [1] * k0 [1] / tempd;
    pmat (0, 2) = alpha [0] [2] * k0 [2] / tempd;

    tempd = k0 [1] + alpha [1] [0] * k0 [0] + alpha [1] [2] * k0 [2];
    pmat (1, 0) = alpha [1] [0] * k0 [0] / tempd;
    pmat (1, 1) = k0 [1] / tempd;
    pmat (1, 2) = alpha [1] [2] * k0 [2] / tempd;

    tempd = k0 [2] + alpha [2] [1] * k0 [1] + alpha [2] [0] * k0 [0];
    pmat (2, 0) = alpha [2] [0] * k0 [0] / tempd;
    pmat (2, 1) = alpha [2] [1] * k0 [1] / tempd;
    pmat (2, 2) = k0 [2] / tempd;

    // Calculate connectivity from full pmat values
    results.connectivity = 3.0;
    for (int i=0; i<3; i++) 
        results.connectivity -= pmat (i, i);
    results.connectivity = results.connectivity / 3.0;
    // Then rescale to growth rate 
    double pscale = pars.r;
    for (int i=0; i<3; i++) 
        pmat (i, i) = 1.0 - pscale * (1.0 - pmat (i, i));
    for (int i=0; i<2; i++) 
    {
        for (int j=(i + 1); j<3; j++) 
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
        for (int j=0; j<3; j++) 
            nold [j] = k0 [j];
        for (int j=0; j<4; j++) 
        {
            results.nmn_node [j] = 0.0;
            results.nsd_node [j] = 0.0;
        }
        results.nmn_network = 0.0;
        results.nsd_network = 0.0;
        flag = false;
        for (int i=0; i<2; i++)
            for (int j=(i+1); j<3; j++)
                results.cov [i] [j] = 0.0;
        for (int i=0; i<(runin + pars.nTrials); i++) 
        {
            // Movement through network
            for (int j=0; j<3; j++) 
            {
                n [j] = 0.0;
                for (int k=0; k<3; k++) 
                    n [j] += pmat (k, j) * nold [k];
            }
            // Change carrying capacities
            for (int j=0; j<3; j++) 
            {
                kvals [j] = k0 [j] + pars.ksd * rnorm ();
                while (kvals [j] < minqc || kvals [j] > (1.0 - minqc))
                    kvals [j] = k0 [j] + pars.ksd * rnorm ();
            }
            // Population dynamic
            tempi = 0;
            for (int j=0; j<3; j++) 
            {
                nold [j] = n [j] + pars.r * n [j] * n [j] / kvals [j] - 
                        pars.r * n [j] * n [j] * n [j] / (kvals [j] * kvals [j]);
                if (nold [j] < 0.0) 
                {
                    //nold [j] = 0.0; // Set to zero after calculating mean & SD values
                    tempi++;
                }
            } // end for j
            if (tempi == 3) 
            {
                flag = true;
                break;
            } else if (i > runin) {
                tempd = 0.0;
                for (int j=0; j<3; j++) 
                {
                    tempd += nold [j];
                    results.nmn_node [j] += nold [j];
                    results.nsd_node [j] += nold [j] * nold [j];
                }
                results.nmn_network += tempd;
                results.nsd_network += tempd * tempd;
                for (int j=0; j<2; j++)
                    for (int k=(j+1); k<3; k++)
                        results.cov [j] [k] += nold [j] * nold [k];
            }
            // Then set all negative nodal abundances to zero
            for (int j=0; j<3; j++)
                if (nold [j] < 0.0)
                    nold [j] = 0.0;
        } // end for i over nTrials
        if (!flag) 
        {
            for (int j=0; j<3; j++) 
            {
                results.nmn_node [j] = results.nmn_node [j] 
                            / (double) pars.nTrials;
                results.nsd_node [j] = results.nsd_node [j] / 
                            (double) pars.nTrials -
                            results.nmn_node [j] * results.nmn_node [j];
            }
            for (int j=0; j<2; j++)
                for (int k=(j+1); k<3; k++)
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
        for (int i=0; i<4; i++) 
        {
            results.nmn_node [i] = DOUBLE_MIN;
            results.nsd_node [i] = DOUBLE_MIN;
        }
        results.nmn_network = DOUBLE_MIN;
        results.nsd_network = DOUBLE_MIN;
    } else {
        results.nmn_node [3] = 0.0;
        results.nsd_node [3] = 0.0;
        for (int i=0; i<3; i++) 
        {
            results.nmn_node [3] += results.nmn_node [i];
            results.nsd_node [3] += results.nsd_node [i];
        }
        results.nmn_node [3] = results.nmn_node [3] / 3.0;
        results.nsd_node [3] = results.nsd_node [3] / 3.0;
        results.nmn_network = results.nmn_network  / (double) pars.nTrials;
        results.nsd_network = results.nsd_network / (double) pars.nTrials -
            results.nmn_network * results.nmn_network;
        results.nmn_network = results.nmn_network / 3.0; 
        // So it's on the same scale as nodal abundance.
    }

    return results;
}
