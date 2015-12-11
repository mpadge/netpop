/*
 * randnet_dendritic.cc
 *
 */

#include "randnet_dendritic.h"



/************************************************************************
 ************************************************************************
 **                                                                    **
 **                         MAIN FUNCTION                              **
 **                                                                    **
 ************************************************************************
 ************************************************************************/

int main (int argc, char *argv[])
{
    int tempi;
	double dtemp;
	bool flag;
	std::string fname;
	Parameters params;
	NetResults results;
	time_t seed;
	clock_t ctimer [2];
	std::ofstream out_file;
	std::ifstream in_file; 
    // in_file is for reading previous data and filling the varcheck table
	std::string linetxt;
	base_generator_type generator (42u);

	time (&seed);
	generator.seed (static_cast <unsigned int> (seed));

	ctimer [0] = clock();
	//	Output formatting stuff ...
	//cout.setf(0,ios::floatfield); 
	std::cout.setf (std::ios::fixed, std::ios::floatfield);   
	std::cout.precision (4);

    try {
        boost::program_options::options_description generic("Generic options");
        generic.add_options()
            ("version,v", "print version std::string")
            ("help", "produce help message")    
            ;

        boost::program_options::options_description config("Configuration");
        config.add_options()
            ("r,r", boost::program_options::value <double>
                (&params.r)->default_value (0.1), "Growth rate, r")
            ("sigma,s", boost::program_options::value <double>
                (&params.sigma)->default_value (0.1), "sigma")
            ("alphasd,a", boost::program_options::value <double>
                (&params.alphasd)->default_value (0.0), "SD of alpha")
            ("directed,d", boost::program_options::value <int>
                (&tempi)->default_value (0), "directed (0/1=F/T)")
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

	if (tempi == 0) 
        params.directed = false;
	else 
        params.directed = true;

    params.nnodes = 25;
    // File names for these versions have an extra "_pinv" inserted just below,
    // while the direct one is multiplied by 2 to generate pscale values of
    // the same magnitude for r=[0.1,0.2] - that is, pscale =[0.2,0.4]
    params.pscale = 0.04 / params.r;
	params.pscale = 2 * params.r;
	params.N0 = 0.0;
	params.k0 = 0.5;
	params.maxSims = 100;
    /* Maximal number of attempts at generating unique combination of structural
     * variables before a stream is accepted regardless.
	 * 
     * The older version of this worked with nTrials, which was set to 500 /
     * nradii = 500 / 7 = 71, with each network then run over shuffles = 50
     * different simulations. This new version includes simulation of networks
     * across the range of alpha0 from 0 to 1, which requires 100 values. This
     * version therefore replaces the nshuffles across the same network with
     * generation of the 101 sets of values according to alpha0.
	 * 
     * The total previous number per run was 50 * 500 = 25,000. The equivalent
     * number here thus requires nTrials = 250 / 7 = 36, compared with the
     * former 500 / 7. While an alternative is to generate entire new networks
     * for each of the alpha0 values, this slows the program enormously, so that
     * the former 50 minutes or so runtime becomes 2.5 hours!
	 */
	params.nTrials = ceil (250.0 / (double) nradii);

	std::stringstream ss;
	ss.str (""); 
    ss << round (100.0 * params.r);
	if (!params.directed) 
        fname = "results_dendritic_r";
	else 
        fname = "results_dendritic_directed_r";
	if (params.r < 0.1) 
        fname += "0";
	if (params.r < 1.0) 
        fname += "0";
	//fname += ss.str () + "_pinv_sigma";
	fname += ss.str () + "_sigma";
	ss.str (""); 
    ss << round (100.0 * params.sigma);
	if (params.sigma < 0.1) 
        fname += "0";
	if (params.sigma < 1.0) 
        fname += "0";
	fname += ss.str () + "_alphasd";
	if (params.alphasd < 0.1) 
        fname += "0";
	if (params.alphasd < 1.0) 
        fname += "0";
	ss.str (""); 
    ss << round (100.0 * params.alphasd);
	fname += ss.str () + ".txt";

	barr3 varcheck (boost::extents [13] [25] [12]); // [radius] [N.I.] [N+]
	for (int i=0; i<13; i++)
		for (int j=0; j<25; j++)
			for (int k=0; k<12; k++)
				varcheck [i] [j] [k] = false;

	out_file.open (fname.c_str(), std::ofstream::out);
	out_file << "alpha,\tconn,\tdiffq,\tradius,\tedgelength" <<
        ",\tn1,\tn2,\tn3,\tnmn,\tnsd,\t" << std::endl;
	flag = false;
	for (int i=0; i<params.nTrials; i++)
    {
		for (int j=0; j<nradii; j++)
        {
			params.pbranch = plist [j];
			results = do1trial (params, &generator, &varcheck);
			for (int k=0; k<101; k++)
            {
				out_file << (double) k / 100.0 << ",\t" << 
                    results.connectivity [k] << ",\t" << 
                    results.diff_q [k] << ",\t" << 
                    results.radius << ",\t" <<
                    results.edgelength << ",\t";
				for (int m=0; m<3; m++) 
                    out_file << results.edgedist [m] << ",\t";
				out_file << results.N_mn [k] << ",\t" << 
                    results.N_sd [k] << ",\t" << std::endl;
			}

			dtemp = 100.0 * ((double) i * (double) nradii + (double) j) /
				((double) params.nTrials * (double) nradii);
			std::cout << "\rdend: [" << i << ", " << j << "] / " <<
				nradii * params.nTrials << " = " << dtemp << "%";
			if (flag) 
            {
                // Because calling the timer when everything is zero produces junk
				ctimer [1] = clock();
				dtemp = ((double) ctimer [1] - (double) ctimer [0]) / 
                    (double) CLOCKS_PER_SEC;
				std::cout << "; [Elapsed, Each, Remaining] time = [";
				timeout (dtemp);
				dtemp = dtemp / ((double) i * (double) nradii + (double) j);
				std::cout << ", "; timeout (dtemp);
				dtemp = dtemp * ((double) params.nTrials * (double) nradii - 
                        (double) i * (double) nradii + (double) j);
				std::cout << ", "; timeout(dtemp); 
				std::cout << "]";
			}
			std::cout.flush ();
			flag = true;
		} // end for j
	} // end for i
	out_file.close();
	std::cout << std::endl;
	dtemp = ((double) clock () - (double) ctimer [0]) / (double) CLOCKS_PER_SEC;
	std::cout << "Calculation time = ";
	timeout (dtemp);
	std::cout << std::endl << std::endl;

	return 0;
} // end main


/************************************************************************
 ************************************************************************
 **                                                                    **
 **                        DO1TRIAL FUNCTION                        *****
 **                                                                    **
 ************************************************************************
 ************************************************************************/

NetResults do1trial (Parameters params, base_generator_type * generator, 
        barr3 * varcheck) {
	/*
     * This implementation creates one network and then shuffles the fixed
     * values of q & c around that network "nshuffles" times, generating
     * abundances and diversities each time. Final values for that network are
     * then averaged over these shuffled versions.
	 *
     * After the shuffling, the "0" and "5" networks are created. The former
     * cuts off all connectivities, which doesn't do anything in the present
     * case because q-values are always the same. The "5" network sets all q & c
     * values to 0.5.
	 * 
     * Extinction probabilites are calculated here from mean & SD abundances,
     * using erf functions (see en.wikipedia.org/wiki/Normal_distribution), as
	 * 	p_ext = 0.5 * (1.0 + erf(-mu / (sigma * sqrt(2.0))))
	 */
	int count, tempi;
	double tempd;
	NetStructure *NetPt = NULL;
	NetResults results;
	NetResults_OneNet results_oneNet;

	boost::numeric::ublas::matrix<double> pmat (params.nnodes, params.nnodes);
	boost::numeric::ublas::matrix<double> Nsums (2, params.nnodes + 1);

	NetPt = new NetStructure [params.nnodes];
	randnet (NetPt, params, generator);
	results_oneNet = netstats (NetPt, params);
	while (results_oneNet.radius > 12 || results_oneNet.edgelength > 24 || 
            results_oneNet.edgedist[2] > 11)
    {
		randnet (NetPt, params, generator);
		results_oneNet = netstats (NetPt, params);
	}
	if (!(*varcheck) [results_oneNet.radius] [results_oneNet.edgelength] 
            [results_oneNet.edgedist[2]])
    {
		(*varcheck) [results_oneNet.radius] [results_oneNet.edgelength] 
            [results_oneNet.edgedist[2]] = true;
    } else {
		count = 0;
		while ((*varcheck) [results_oneNet.radius] [results_oneNet.edgelength] 
                [results_oneNet.edgedist[2]])
        {
			randnet (NetPt, params, generator);
			results_oneNet = netstats (NetPt, params);
			while (results_oneNet.radius > 12 || 
                    results_oneNet.edgelength > 24 || 
                    results_oneNet.edgedist[2] > 11)
            {
				results_oneNet.radius = INT_MIN;
				results_oneNet.edgelength = INT_MIN;
				results_oneNet.edgedist[2] = INT_MIN;
				randnet (NetPt, params, generator);
				results_oneNet = netstats (NetPt, params);
			}
			count++;
			if (!(*varcheck) [results_oneNet.radius] [results_oneNet.edgelength] 
                    [results_oneNet.edgedist[2]])
            {
				(*varcheck) [results_oneNet.radius] [results_oneNet.edgelength] 
                    [results_oneNet.edgedist[2]] = true;
				break;
            }
			if (count > params.maxSims) 
                break;
		} // end while
	} // end else iterate

	results_oneNet = netstats (NetPt, params);
	for (int i=0; i<3; i++)
		results.edgedist [i] = results_oneNet.edgedist [i];
	results.radius = results_oneNet.radius;
	results.edgelength = results_oneNet.edgelength;

	// Then loop over the 101 values of alpha0
	for (int i=0; i<101; i++)
    {
		params.alpha0 = (double) i / 100.0;
		fillqc (NetPt, params, generator, true); // true makes rand-unif values
		makepmat (NetPt, &pmat);
		// Calculate the mean difference in adjacent q-values.
		tempd = 0.0;
		count = 0;
		for (int j=0; j<params.nnodes; j++)
        {
			tempi = NetPt [j].nextds;
			if (tempi > INT_MIN) 
            {
				tempd += fabs (NetPt [j].q - NetPt [tempi].q);
				count++;
			}
		}
		results.diff_q [i] = tempd / (double) count;
		// And then connectivity
		tempd = 0.0;
		for (int j=0; j<params.nnodes; j++)
        {
			NetPt [j].N = 1.0;
			tempd += pmat (j, j); 
            // Not scaled to pscale, as it actually is for movement
		}
		results.connectivity [i] = ((double) params.nnodes - tempd) / 
            (double) params.nnodes;

		Nsums = populateNet (&pmat, NetPt, params, generator);
		// Nsums params.nnodes has mean & SD total network abundance.
		results.N_mn [i] = Nsums (0, params.nnodes);
		results.N_sd [i] = Nsums (1, params.nnodes);
	} // end for i

	delete [] NetPt; NetPt = NULL;

	return results;
} // end function do1trial




/************************************************************************
 ************************************************************************
 **                                                                    **
 **                        RANDNET FUNCTION                            **
 **                                                                    **
 ************************************************************************
 ************************************************************************/

void randnet (NetStructure* NetPt, Parameters params, 
        base_generator_type * generator)
{
    int count, tempi, nodecount, xi, nodeHt;
    int branches [2];
    double randnum;
    bool xshift[3];

    boost::uniform_real<> uni_dist(0,1);
    boost::variate_generator<base_generator_type&,
        boost::uniform_real<> > runif((*generator), uni_dist);

    // Set up the null StreamPt structure:
    for (int i=0; i < params.nnodes; i++)
    {
        NetPt[i].nextus1 = INT_MIN;
        NetPt[i].nextus2 = INT_MIN;
        NetPt[i].nextds = INT_MIN;
        NetPt[i].xpos = INT_MIN;
        NetPt[i].ypos = INT_MIN;
    }
    // Then the first nodes
    NetPt[0].nextus1 = 1;
    NetPt[1].nextds = 0;
    NetPt[0].xpos = 0;
    NetPt[1].xpos = 0;
    NetPt[0].ypos = 0;
    NetPt[1].ypos = 1;

    nodeHt = 1; // Height of new nodes to be added

    // Main loop to add the remaining (nnodes - 2) nodes:
    nodecount = 2;
    while (nodecount < params.nnodes)
    {
        // Loop over all extant nodes to find top nodes (.ypos==nodeHt) ...
        for (int i=0; i < nodecount; i++) 
        { 
            if (NetPt[i].ypos == nodeHt)
            {
                xi = NetPt[i].xpos;
                /* Then make the list of potential xpos that may be filled in
                 * the y-space 1 step above TopNodeHt. Default at start is
                 * all=TRUE */
                for (int j=0; j<3; j++)
                    xshift[j] = true;
                count = 3; // Number of possible points to connect with
                /* New nodes at NodeHt+1 are added throughout, so this needs to
                   scan through the whole list. It also needs a bit to avoid
                   paths crossing over. Since each loop here adds only one new
                   path, it can only cross if one already exists for it to cross
                   over. If already existent paths simply go to -1 or +1, no new
                   path from (x,y)=(0,-1) can cross these. Thus, the only
                   existent paths that can be crossed are those that go to the
                   x=0 position. This bit of code locates the x-origins of any
                   such paths, and removes this value from the xshift list. It's
                   all a bit complicated, but that's just the way it goes ...
                   */
                for (int j=0; j<nodecount; j++)
                {
                    if (abs (NetPt[j].xpos - xi) <= 1 && 
                            NetPt[j].ypos == (nodeHt + 1))
                    {
                        xshift[NetPt[j].xpos - xi + 1] = false;
                        count = count - 1;
                        // ... and the bit to avoid crossing over ...
                        if (NetPt[j].xpos == xi)
                        {
                            tempi = NetPt[j].nextds;
                            tempi = NetPt[tempi].xpos - xi + 1;
                            xshift [tempi] = false;
                            count--;
                        }
                    }
                } // end for j
                // Then add nodes if places are available:
                if (count > 0)
                {
                    if (count == 1) 
                    { // Go to single possible place
                        tempi = INT_MIN;
                        for (int j=0; j<3; j++)
                            if (xshift [j] == true) 
                                tempi = j - 1;
                        NetPt[nodecount].ypos = nodeHt + 1;
                        NetPt[nodecount].xpos = NetPt[i].xpos + tempi;
                        NetPt[i].nextus1 = nodecount;
                        NetPt[nodecount].nextds = i;
                        nodecount++;
                    }
                    else 
                    { // Branch is possible
                        if (xshift[0] == true && xshift[1] == true)
                        {
                            branches[0] = -1; branches[1] = 0; 
                        } else if (xshift[1] == true && xshift[2] == true) {
                            branches[0] = 0; 
                            branches[1] = 1; 
                        } else {
                            branches[0] = -1; 
                            branches[1] = 1;
                        }
                        NetPt[nodecount].ypos = nodeHt + 1;
                        NetPt[nodecount].xpos = NetPt[i].xpos + branches[0];
                        NetPt[i].nextus1 = nodecount;
                        NetPt[nodecount].nextds = i;
                        nodecount++;
                        if (nodecount == params.nnodes) 
                            break;
                        if (runif () < params.pbranch)
                        {
                            NetPt[nodecount].ypos = nodeHt + 1;
                            NetPt[nodecount].xpos = NetPt[i].xpos + branches[1];
                            NetPt[i].nextus2 = nodecount;
                            NetPt[nodecount].nextds = i;
                            nodecount++;
                        } // end if 2nd branch
                    } // end else branch is possible
                    if (nodecount == params.nnodes) 
                        break;
                } // end if count > 0
                if (nodecount == params.nnodes) 
                    break;
            }	} // end for i & if .ypos==nodeHt
        nodeHt++;
    } // end while nodecount < nnodes loop
}; // end function randnet

/************************************************************************
 ************************************************************************
 **                                                                    **
 **                       FILLQC FUNCTION                              **
 **                                                                    **
 ************************************************************************
 ************************************************************************/

void fillqc (NetStructure * NetPt, Parameters params, 
        base_generator_type * generator, bool unif)
{
	/*
     * pseq can make exponential distributions centred around 0.5 (q=true), or
     * decreasing from 1 (q=false), with the latter potentially used for
     * connectivities so that most of them end up with relatively high values.
     * This, however, seems to produce much lower correlations with structural
     * network properties, so connectivities are now reverted to the q-dist.
	 *
     * To change this back, just ditch the even-to-odd conversion bit, and set 
     * > clist = pseq(count, false);
	 */
	int count;

    if (unif)
    {
        // Fill q & c with runif values
        boost::uniform_real<> uni_dist(0,1);
        boost::variate_generator<base_generator_type&,
            boost::uniform_real<> > runif((*generator), uni_dist);
        boost::normal_distribution<> norm_dist(0,1);
        boost::variate_generator<base_generator_type&,
            boost::normal_distribution<> > rnorm((*generator), norm_dist);

        if (params.alphasd == 0.0)
        {
            for (int i=0; i<params.nnodes; i++)
            {
                NetPt [i].q = runif ();
                if (NetPt[i].nextds != INT_MIN) 
                    NetPt[i].cdown = runif ();
                else 
                    NetPt[i].cdown = DOUBLE_MIN;
                if (NetPt[i].nextus1 != INT_MIN) 
                    NetPt[i].cup1 = runif ();
                else 
                    NetPt[i].cup1 = DOUBLE_MIN;
                if (NetPt[i].nextus2 != INT_MIN) 
                    NetPt[i].cup2 = runif ();
                else 
                    NetPt[i].cup2 = DOUBLE_MIN;
            } // end for i
        } else {
            for (int i=0; i<params.nnodes; i++)
            {
                NetPt [i].q = 0.5 + params.alphasd * rnorm ();
                while (NetPt [i].q <= minqc || NetPt [i].q > (1.0 - minqc))
                {
                    NetPt [i].q = 0.5 + params.alphasd * rnorm ();
                }

                if (NetPt[i].nextds != INT_MIN)
                {
                    NetPt [i].cdown = params.alpha0 + params.alphasd * rnorm ();
                    while (NetPt [i].cdown < minqc || NetPt [i].cdown > (1.0 - minqc))
                    {
                        NetPt [i].cdown = params.alpha0 + params.alphasd * rnorm ();
                    }
                } else { 
                    NetPt[i].cdown = DOUBLE_MIN;
                }
                if (NetPt[i].nextus1 != INT_MIN)
                {
                    NetPt [i].cup1 = params.alpha0 + params.alphasd * rnorm ();
                    while (NetPt [i].cup1 < minqc || NetPt [i].cup1 > (1.0 - minqc))
                    {
                        NetPt [i].cup1 = params.alpha0 + params.alphasd * rnorm ();
                    }
                } else { 
                    NetPt[i].cup1 = DOUBLE_MIN;
                }
                if (NetPt[i].nextus2 != INT_MIN)
                {
                    NetPt [i].cup2 = params.alpha0 + params.alphasd * rnorm ();
                    while (NetPt [i].cup2 < minqc || NetPt [i].cup2 > (1.0 - minqc))
                    {
                        NetPt [i].cup2 = params.alpha0 + params.alphasd * rnorm ();
                    }
                } else { 
                    NetPt[i].cup2 = DOUBLE_MIN;
                }
            } // end for i
        } // end else alphasd != 0.0
    } else { // not unif, so generate a normal distribution
        std::vector <int> indx = randseq (params.nnodes, generator);
        boost::numeric::ublas::vector<double> q (params.nnodes);
        pseq (&q, true);
        count = 0;
        for (int i=0; i<params.nnodes; i++) {
            NetPt[indx[i]].q = q (i);
            if (NetPt[i].nextds != INT_MIN) { count++;	}
            if (NetPt[i].nextus1 != INT_MIN) { count++;	}
            if (NetPt[i].nextus2 != INT_MIN) { count++;	}
        } // end for i
        // pseq for (q=true) has to have an odd-numbered length.
        if (floor((double) count / 2.0) == ((double) count / 2.0)) {
            count++;
        }
        std::vector <int> indx2 = randseq (count, generator);
        boost::numeric::ublas::vector<double> clist (count);
        pseq (&clist, true);
        boost::numeric::ublas::vector<double> clist2 (count);
        for (int i=0; i<count; i++) {
            clist2 (i) = clist (indx2 [i]);
        }
        count = 0;
        for (int i=0; i<params.nnodes; i++) {
            if (NetPt[i].nextds != INT_MIN) {
                NetPt[i].cdown = clist2 (count);
                count++;	}
            else { NetPt[i].cdown = DOUBLE_MIN;	}
            if (NetPt[i].nextus1 != INT_MIN) {
                NetPt[i].cup1 = clist2 (count);
                count++;	}
            else { NetPt[i].cup1 = DOUBLE_MIN;	}
            if (NetPt[i].nextus2 != INT_MIN) {
                NetPt[i].cup2 = clist2 (count);
                count++;	}
            else { NetPt[i].cup2 = DOUBLE_MIN;	}
        } // end for i
    }

	// Then simulate directed network by making .cdown always greater than .cup
	if (params.directed)
    {
		int nextds;
		double tempd;
		for (int i=0; i<params.nnodes; i++)
        {
			nextds = NetPt [i].nextds;
			if (nextds != INT_MIN)
            {
				if (NetPt [nextds].nextus1 == i && 
                        NetPt [i].cdown < NetPt [nextds].cup1)
                {
					tempd = NetPt [i].cdown;
					NetPt [i].cdown = NetPt [nextds].cup1;
					NetPt [i].cup1 = tempd;
				}
				else if (NetPt [nextds].nextus2 == i && 
                        NetPt [i].cdown < NetPt [nextds].cup2)
                {
					tempd = NetPt [i].cdown;
					NetPt [i].cdown = NetPt [nextds].cup2;
					NetPt [i].cup2 = tempd;
				}
			}
		} // end for i
	} // end if directed
} // end fillqc



/************************************************************************
 ************************************************************************
 **                                                                    **
 **                       NETSTATS FUNCTION                            **
 **                                                                    **
 ************************************************************************
 ************************************************************************/

NetResults_OneNet netstats (NetStructure *NetPt, Parameters params)
{
	int tempi[2];
	double tempd;
	NetResults_OneNet results_oneNet;

	// Calculate radius
	boost::numeric::ublas::matrix <double> dmat (params.nnodes, params.nnodes);
	dmat = makedmat (NetPt, params);
	tempi[0] = 9999;
	for (int k=0; k<params.nnodes; k++)
    {
		tempi[1] = 0;
		for (int m=0; m<params.nnodes; m++)
			if (dmat(k, m) > tempi[1]) 
                tempi[1] = dmat(k, m);
		if (tempi[1] < tempi[0]) 
            tempi[0] = tempi[1];
	}
	results_oneNet.radius = tempi[0];

	// Then fill edge distribution and count edgelength
	results_oneNet.edgedist [0] = 0;
	results_oneNet.edgedist [1] = 0;
	results_oneNet.edgedist [2] = 0;
	for (int i=0; i<params.nnodes; i++)
    {
		if (NetPt [i].nextds == INT_MIN)
        {
			if (NetPt [i].nextus2 > INT_MIN) 
                results_oneNet.edgedist [1]++;
			else 
                results_oneNet.edgedist [0]++;
        } else {
			if (NetPt [i].nextus2 > INT_MIN) 
                results_oneNet.edgedist [2]++;
			else if (NetPt [i].nextus1 > INT_MIN) 
                results_oneNet.edgedist [1]++;
			else 
                results_oneNet.edgedist [0]++;
		}
	} // end for i

    // When networks don't branch at all, edgelength counts can go in both
    // directions, so the "done" index prevents this happening.
	boost::numeric::ublas::vector <bool> done (params.nnodes);
	for (int i=0; i<params.nnodes; i++) 
        done (i) = false;

	results_oneNet.edgelength = 0;
	for (int i=0; i<params.nnodes; i++)
    {
		if (!done (i))
        {
			if (NetPt[i].nextds == INT_MIN && NetPt[i].nextus2 == INT_MIN)
            {
				done (i) = true;
				done (NetPt[i].nextus1) = true;
				results_oneNet.edgelength++;
            } else if (NetPt[i].nextus1 == INT_MIN) {
				done (i) = true;
				results_oneNet.edgelength++;
				tempi [0] = NetPt[i].nextds;
				if (tempi [0] > INT_MIN)
                {
					while (tempi [0] > INT_MIN && 
                            NetPt [tempi [0]].nextus2 == INT_MIN && 
                            !done (tempi [0]))
                    {
						done (tempi [0]) = true;
						results_oneNet.edgelength++;
						tempi [0] = NetPt [tempi [0]].nextds;
					} // end while
				} // end if (tempi[0] > INT_MIN)
                // This next line only arises for straight lines, and was added
                // *AFTER* the current results were generated, so nnodes=24 and
                // 25 are mixed together under nnodes=24. This requires
                // separation in analyses, through setting all networks with
                // results.edgedist[2] == 0 to results.edgelength = 25 instead
                // of 24.
				if (tempi [0] == INT_MIN) 
                    results_oneNet.edgelength++;
			} // end else if
		} // end if !done
	} // end for i

	return results_oneNet;
} // end function netstats



/************************************************************************
 ************************************************************************
 **                                                                    **
 **                     POPULATENET FUNCTION                           **
 **                                                                    **
 ************************************************************************
 ************************************************************************/

dmat populateNet (dmat * pmat, NetStructure *NetPt, Parameters params, 
        base_generator_type * generator)
{
    /* nIters is the length of the initial run-in, and also of the subsequent
     * data collection steps. Values are counted for individual nodes and for
     * total population.
	 * 
     * NOTE that mean & SD abundances are calculated across the entire range of
     * values including negative ones. The movement stage, however, can not
     * include negative values, and thus all abundance calculations are done
     * directly after implementation of population dynamics, but before
     * movement.
	 */
	const int nIters = 1000;
	double qval, tempd [2];

	boost::normal_distribution<> norm_dist(0,1);
	boost::variate_generator<base_generator_type&,
		boost::normal_distribution<> > rnorm((*generator), norm_dist);

	boost::numeric::ublas::matrix<double> oldN (2, params.nnodes);
	for (int i=0; i<params.nnodes; i++) 
        oldN (0, i) = NetPt[i].N;
    // Nsums holds values for each node, plus network means for abundance &
    // diversity.
    boost::numeric::ublas::matrix<double> Nsums (2, params.nnodes + 1);
	for (int i=0; i<Nsums.size1(); i++)
    {
		for (int j=0; j<Nsums.size2(); j++)
        {
			Nsums (i, j) = 0.0;
        }
    }

	const int runin = 100;
	for (int n=0; n<runin; n++)
    {
		/* Movement:
         * The proportion that stay at [i][i] is modified from
         * NetPt[i].N*pmat[i][i], because pmat dynamics are scaled by pscale.
         * Thus instead of sum(b)+a=1, where a=pmat[i][i], one seeks the value
         * of a' such that sum(pscale*b)+a'=1.  This is then
         * a'=1-sum(pscale*b)=1-pscale*sum(b)=1-pscale*(1-a). */
		for (int i=0; i<params.nnodes; i++)
        {
			oldN (1, i) = oldN (0, i) * 
                (1.0 - params.pscale * (1.0 - (*pmat) (i, i)));
			for (int j=0; j<params.nnodes; j++) 
            {
                if (j != i) 
                    oldN (1, i) += oldN (0, j) * params.pscale * (*pmat) (j, i);
            }
		}
		for (int i=0; i<params.nnodes; i++)
        {
			qval = NetPt[i].q + params.sigma * rnorm();
			while (qval < minqc || qval > (1.0 - minqc))
				qval = NetPt [i].q + params.sigma * rnorm();
			oldN (0, i) = oldN (1, i) + oldN (1, i) * params.r *
				(1.0 - oldN (1, i) / qval) * (oldN (1, i) / qval - 
                        params.N0 / qval);
			if (oldN (0, i) < 0.0) 
                oldN (0, i) = 0.0;
		} // end for i
	}

	for (int n=0; n<nIters; n++)
    {
		for (int i=0; i<params.nnodes; i++)
        {
			oldN (1, i) = oldN (0, i) * 
                (1.0 - params.pscale * (1.0 - (*pmat) (i, i)));
			for (int j=0; j<params.nnodes; j++)
            {
                if (j != i)
                    oldN (1, i) += oldN (0, j) * params.pscale * (*pmat) (j, i);
            }
		}
		for (int i=0; i<params.nnodes; i++)
        {
			qval = NetPt [i].q + params.sigma * rnorm();
			while (qval < minqc || qval > (1.0 - minqc))
				qval = NetPt [i].q + params.sigma * rnorm();
			oldN (0, i) = oldN (1, i) + oldN (1, i) * params.r *
				(1.0 - oldN (1, i) / qval) * (oldN (1, i) / qval - 
                        params.N0 / qval);
			//if (oldN (0, i) < 0.0) { oldN (0, i) = 0.0;	}
		} // end for i

		// Abundance calculations
		tempd [0] = 0.0;
		for (int i=0; i<params.nnodes; i++)
        {
			Nsums (0, i) += oldN (0, i);
			Nsums (1, i) += oldN (0, i) * oldN (0, i);
			tempd [0] += oldN (0, i);
        }
		Nsums (0, params.nnodes) += tempd[0];
		Nsums (1, params.nnodes) += tempd[0] * tempd[0];
		// Then only set all negative abundances to zero after calculations
		for (int i=0; i<params.nnodes; i++)
			if (oldN (0, i) < 0.0) 
                oldN (0, i) = 0.0;
	} // end for n

	// Convert Nsums to means & SDs
	for (int i=0; i<=params.nnodes; i++)
    {
		Nsums (0, i) = Nsums (0, i) / (double) nIters;
		Nsums (1, i) = Nsums (1, i) / (double) nIters - 
            Nsums (0, i) * Nsums (0, i);
		Nsums (1, i) = sqrt (Nsums (1, i));
	}

	return Nsums;
} // end function popnet


/************************************************************************
 ************************************************************************
 **                                                                    **
 **                   POPULATENET_EQ FUNCTION                          **
 **                                                                    **
 ************************************************************************
 ************************************************************************/

dmat populateNet_eq (dmat * pmat, Parameters params)
{
    /* This one just calculates "equilibrium" values as the eigenvector of pmat,
     * rather than implementing population dynamics. As such, variances in
     * abundance and diversity can not be calculated here, only mean values.
     * Nevertheless, the same (2,nnodes+2) structure is returned, but with the
     * final values (nnodes+1) not filled.*/
	const int maxiter = 10000;
	const double tol = 1.0e-4;
	int count;
	double qval, tempd[2], tolsum[2];

	boost::numeric::ublas::matrix<double> oldN (2, params.nnodes);
	for (int i=0; i<params.nnodes; i++) 
        oldN (0, i) = 1.0;

	// Nsums holds values for each node, plus network means for abundance & diversity.
	boost::numeric::ublas::matrix<double> Nsums (2, params.nnodes + 1);
	for (int i=0; i<Nsums.size1(); i++)
		for (int j=0; j<Nsums.size2(); j++)
			Nsums (i, j) = 0.0;

	tolsum[0] = 99999.9; tolsum[1] = 0.0;
	count = 0;
	while (fabs(tolsum[1] - tolsum[0]) > tol)
    {
		tolsum[0] = tolsum[1];
		tolsum[1] = 0.0;
		for (int i=0; i<params.nnodes; i++)
        {
			oldN (1, i) = 0.0;
			for (int j=0; j<params.nnodes; j++)
            {
				if (j == i) 
					oldN (1, i) += oldN (0, i) * (1.0 - params.pscale *
						(1.0 - (*pmat) (i, i)));
				else 
					oldN (1, i) += oldN (0, j) * params.pscale * (*pmat) (j, i);
			} // end for j
			tolsum[1] += oldN (1, i);
		} // end for i
		for (int i=0; i<params.nnodes; i++) 
            oldN (0, i) = oldN (1, i);
		count++;
		if (count > maxiter)
        {
			std::cout<<"ERROR: popnet_eq failed to converge!"<<std::endl;
			break;
        }
	} // end while

    /* Abundance calculation - produces mean and SD across all nodes for the
     * individual network, so aggregated values are just mean values of these
     * two (with mean = nnodes at al times).*/
	for (int i=0; i<2; i++) 
        Nsums (i, params.nnodes) = 0.0;
	for (int i=0; i<params.nnodes; i++)
    {
		Nsums (0, i) = oldN (1, i);
		Nsums (0, params.nnodes) += oldN (1, i);
		Nsums (1, params.nnodes) += oldN (1, i) * oldN (1, i);
	}
	Nsums (0, params.nnodes) = Nsums (0, params.nnodes) / (double) params.nnodes;
	Nsums (1, params.nnodes) = Nsums (1, params.nnodes) / (double) params.nnodes -
		Nsums (0, params.nnodes) * Nsums (0, params.nnodes);

	return Nsums;
} // end function popnet_eq


/************************************************************************
 ************************************************************************
 **                                                                    **
 **                       MAKEPMAT FUNCTION                            **
 **                                                                    **
 ************************************************************************
 ************************************************************************/

void makepmat(const NetStructure* NetPt, dmat * pmat)
{
/*
 * The original traced all paths explicitly (in the old main.cc), but
 * this new versions uses Dikstra's SP from each node to just find
 * the shortest path to all others. The full pmat is then calculated
 * on the basis of all shortest paths.
 */
	int here, there, lenQ, plen, nnodes = (*pmat).size1();
	double dmin;

	boost::numeric::ublas::vector<double> dist(nnodes);
	boost::numeric::ublas::vector<int> Q(nnodes);
	boost::numeric::ublas::vector<int> prev(nnodes);
	boost::numeric::ublas::vector<bool> done(nnodes);

	for (int i=0; i<nnodes; i++)
		for (int j=0; j<nnodes; j++)
			(*pmat) (i, j) = 0.0;

	for (int i=0; i<nnodes; i++)
    {
		for (int j=0; j<nnodes; j++)
        {
			dist (j) = 1e99;
			Q (j) = j;
			prev (j) = INT_MIN;
			done (j) = false;
        }
		here = i;
		dist (here) = 0.0;
		lenQ = 19;

		while (lenQ > 0)
        {
			dmin = 1e99; here = INT_MIN;
			for (int j=0; j<nnodes; j++) {
				if (!done[j] && dist[j] < dmin)
                {
					dmin = dist[j];
					here = j;
                }
            }
			if (here == INT_MIN) 
                std::cout << "ERROR: here == nix!" << std::endl;
			Q (here) = INT_MIN;
			done (here) = true;
			// Neighbour relaxation, i.e. update of shortest distances
			// for this dendritic version, this is done explicitly for the 3 possibilites.
			if (NetPt[here].nextds > INT_MIN)
            {
				if (dist (here) == 0.0) 
                    dmin = 1.0 / NetPt[here].cdown;
				else 
                    dmin = dist (here) / NetPt[here].cdown;
				there = NetPt[here].nextds;
				if (dist(there) > dmin)
                {
					dist (there) = dmin;
					prev (there) = here;
                }
            }
			if (NetPt[here].nextus1 > INT_MIN)
            {
				if (dist (here) == 0.0) 
                    dmin = 1.0 / NetPt[here].cup1;
				else 
                    dmin = dist (here) / NetPt[here].cup1;
				there = NetPt[here].nextus1;
				if (dist (there) > dmin)
                {
					dist (there) = dmin;
					prev (there) = here;
                }
            }
			if (NetPt[here].nextus2 > INT_MIN)
            {
				if (dist (here) == 0.0) 
                    dmin = 1.0 / NetPt[here].cup2;
				else 
                    dmin = dist (here) / NetPt[here].cup2;
				there = NetPt[here].nextus2;
				if (dist (there) > dmin)
                {
					dist (there) = dmin;
					prev (there) = here;
                }
            }
			lenQ = 0;
			for (int j=0; j<nnodes; j++) 
                if (Q(j) > INT_MIN) 
                    lenQ++;
		} // end while (lenq > 0)

		// Then trace back the prev nodes to contstruct pmat
		for (int j=0; j<nnodes; j++)
        {
			here = j;
			plen = 0;
			while (prev[here] != INT_MIN)
            {
				plen++;
				here = prev (here);
            }
			if (plen == 0) 
                (*pmat) (i, j) = NetPt[j].q;
			else 
                (*pmat) (i, j) = NetPt[j].q / dist[j];
		} // end for j
	} // end for i - the main loop
	
	// Then normalise the pmat
	for (int i=0; i<nnodes; i++)
    {
		dmin = 0.0;
		for (int j=0; j<nnodes; j++) 
            dmin += (*pmat) (i, j);
		for (int j=0; j<nnodes; j++) 
            (*pmat) (i, j) = (*pmat) (i, j) / dmin;
	}
} // end function makepmat


/************************************************************************
 ************************************************************************
 **                                                                    **
 **                       MAKEDMAT FUNCTION                            **
 **                                                                    **
 ************************************************************************
 ************************************************************************/

imat makedmat(const NetStructure* NetPt, Parameters params)
{
/*
 * Slightly modified version of makepmat that extracts the shortest
 * inter-nodal distances between all pairs of points
 */
	int here, there, lenQ, plen, dmin;

	boost::numeric::ublas::vector<int> dist(params.nnodes);
	boost::numeric::ublas::vector<int> Q(params.nnodes);
	boost::numeric::ublas::vector<int> prev(params.nnodes);
	boost::numeric::ublas::vector<bool> done(params.nnodes);

	boost::numeric::ublas::matrix<int> dmat(params.nnodes, params.nnodes);
	for (int i=0; i<params.nnodes; i++)
		for (int j=0; j<params.nnodes; j++)
			dmat(i, j) = 0;

	for (int i=0; i<params.nnodes; i++)
    {
		for (int j=0; j<params.nnodes; j++)
        {
			dist(j) = 9999;
			Q(j) = j;
			prev(j) = INT_MIN;
			done(j) = false;
        }
		here = i;
		dist(here) = 0;
		lenQ = 19;

		while (lenQ > 0)
        {
			dmin = 9999; here = INT_MIN;
			for (int j=0; j<params.nnodes; j++) 
            { 
                if (!done(j) && dist(j) < dmin)
                {
                    dmin = dist(j); 
                    here = j;
                }
            }
			Q(here) = INT_MIN;
			done(here) = true;
			// Neighbour relaxation, i.e. update of shortest distances
			// for this dendritic version, this is done explicitly for the 3 possibilites.
			dmin = dist(here) + 1;
			if (NetPt[here].nextds > INT_MIN)
            {
				there = NetPt[here].nextds;
				if (dist(there) > dmin)
                {
					dist(there) = dmin;
					prev(there) = here;
                }
            }
			if (NetPt[here].nextus1 > INT_MIN)
            {
				there = NetPt[here].nextus1;
				if (dist(there) > dmin)
                {
					dist(there) = dmin;
					prev(there) = here;
                }
            }
			if (NetPt[here].nextus2 > INT_MIN)
            {
				there = NetPt[here].nextus2;
				if (dist(there) > dmin)
                {
					dist(there) = dmin;
					prev(there) = here;
                }
            }
			lenQ = 0;
			for (int j=0; j<params.nnodes; j++) 
                if (Q(j) > INT_MIN) 
                    lenQ++;
		} // end while (lenq > 0)

		// Then trace back the prev nodes to contstruct pmat
		for (int j=0; j<params.nnodes; j++)
        {
			if (prev(j) == INT_MIN)
            {
				dmat(j, j) = 0;
            } else {
				here = j;
				plen = 0;
				while (prev(here) != INT_MIN)
                {
					plen++;
					here = prev(here);
                }
				if (plen == 0) 
                    dmat(i, j) = 0;
				else 
                    dmat(i, j) = plen;
            }
		} // end for j
	} // end for i - the main loop

	// Then convert all dmat values to the shortest of the paired distances
	for (int i=0; i<params.nnodes; i++) 
    { 
        for (int j=0; j<params.nnodes; j++)
        {
            if (dmat(i, j) < dmat(j, i)) 
                dmat(j, i) = dmat(i, j);
            else if (dmat(i, j) > dmat(j, i)) 
                dmat(i, j) = dmat(j, i);
        } // end for j
    } // end for i

	return dmat;
} // end function makedmat

