/***************************************************************************
 *  Project:    netpop
 *  File:       net3.c++
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
    int netcount;
    double tempd, progress;
    std::string tempstr; // for simple timeout for non-linux systems
    Net3 net;
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
            std::cout << "net3, version 1.0" << std::endl;
            return 0;
        }

    }
    catch(std::exception& e)
    {
        std::cout << e.what() << std::endl;
        return 1;
    }    
    int nnodes = net.get_nnodes (); // == 3
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
    out_file << "alpha,\tconn,\tmn0,\tmn1,\tmn2,\tmnall\t,sd0,\tsd1,\tsd2,\t" <<
        "sdall,\tmnnet,\tsdnet,\tcov01,\tcov02,\tcov12" << std::endl;
    for (int i=1; i<=100; i++) 
    {
        net.pars.alpha0 = (double) i / 100.0;
        out_file << net.pars.alpha0;
        for (int k=0; k<(nnodes + 1); k++) 
        {
            net.results.nmn_node [k] = 0.0;
            net.results.nsd_node [k] = 0.0;
        }
        net.results.nmn_net = 0.0;
        net.results.nsd_net = 0.0;
        net.results.conn_mn = 0.0;
        netcount = 0;
        for (int j=0; j<net.pars.nRepeats; j++) 
        {
            net.fill_alpha (&generator);
            net.make_pmat (&generator);
            net.iterate_population (&generator, nnodes);
        
            if (net.results1.nmn_network > DOUBLE_MIN) 
            {
                netcount++;
                net.results.conn_mn += net.results1.connectivity;
                for (int k=0; k<(nnodes + 1); k++) 
                {
                    net.results.nmn_node [k] += net.results1.nmn_node [k];
                    net.results.nsd_node [k] += net.results1.nsd_node [k];
                }
                net.results.nmn_net += net.results1.nmn_network;
                net.results.nsd_net += net.results1.nsd_network;
            }
        } // end for j
        if (netcount > 0) 
        {
            net.results.conn_mn = net.results.conn_mn / (double) netcount;
            for (int j=0; j<(nnodes + 1); j++) 
            {
                net.results.nmn_node [j] = net.results.nmn_node [j] / (double) netcount;
                net.results.nsd_node [j] = net.results.nsd_node [j] / (double) netcount;
            }
            net.results.nmn_net = net.results.nmn_net / (double) netcount;
            net.results.nsd_net = net.results.nsd_net / (double) netcount;
        } else {
            net.results.conn_mn = DOUBLE_MIN;
            for (int j=0; j<(nnodes + 1); j++) 
            {
                net.results.nmn_node [j] = DOUBLE_MIN;
                net.results.nsd_node [j] = DOUBLE_MIN;
            }
            net.results.nmn_net = DOUBLE_MIN;
            net.results.nsd_net = DOUBLE_MIN;
        }
        out_file << ",\t" << net.results.conn_mn;
        for (int j=0; j<(nnodes + 1); j++) 
            out_file << ",\t" << net.results.nmn_node [j];
        for (int j=0; j<(nnodes + 1); j++) 
            out_file << ",\t" << net.results.nsd_node [j];
        out_file << ",\t" << net.results.nmn_net << ",\t" << net.results.nsd_net;
        for (int j=0; j<(nnodes - 1); j++) 
            for (int k=(j+1); k<nnodes; k++)
                out_file << ",\t" << net.results1.cov (j, k);
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
 **                          FILL_ALPHA                                **
 **                                                                    **
 ************************************************************************
 ************************************************************************/

void Net3::fill_alpha (base_generator_type * generator)
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
    // that nnodes = 3!
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

void Net3::make_pmat (base_generator_type * generator)
{
    int nnodes = get_nnodes ();
    // Also has to explicitly assume here that nnodes = 3!
    double tempd;

    tempd = k0 [0] + alpha (0, 1) * k0 [1] + alpha (0, 2) * k0 [2];
    pmat (0, 0) = k0 [0] / tempd;
    pmat (0, 1) = alpha (0, 1) * k0 [1] / tempd;
    pmat (0, 2) = alpha (0, 2) * k0 [2] / tempd;

    tempd = k0 [1] + alpha (1, 0) * k0 [0] + alpha (1, 2) * k0 [2];
    pmat (1, 0) = alpha (1, 0) * k0 [0] / tempd;
    pmat (1, 1) = k0 [1] / tempd;
    pmat (1, 2) = alpha (1, 2) * k0 [2] / tempd;

    tempd = k0 [2] + alpha (2, 1) * k0 [1] + alpha (2, 0) * k0 [0];
    pmat (2, 0) = alpha (2, 0) * k0 [0] / tempd;
    pmat (2, 1) = alpha (2, 1) * k0 [1] / tempd;
    pmat (2, 2) = k0 [2] / tempd;

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
