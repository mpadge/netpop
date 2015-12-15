/*
 * randnet_lattice.cc
 * 
 * Normally distributed connectivities and habitat qualities,
 * all using fixed mean values of 0.5, with SDs both controlled
 * by the same parameter of alphasd, EXCEPT if alphasd == 0, in
 * which case, runif () values are used.
 */

#include "randnet_lattice.h"


/************************************************************************
 ************************************************************************
 **                                                                    **
 **                         MAIN FUNCTION                              **
 **                                                                    **
 ************************************************************************
 ************************************************************************/

int main (int argc, char *argv[])
{
	double dtemp, progress;
	std::string fname;
	Parameters params;
	NetResults results;
	time_t seed;
	clock_t time_start;
	std::ofstream out_file;
	std::ifstream in_file; // For reading previous data and filling the varcheck table
	std::string linetxt;
	base_generator_type generator (42u);

	time (&seed);
	generator.seed (static_cast <unsigned int> (seed));

	time_start = clock();
	//	Output formatting stuff ...
	//std::cout.setf(0,ios::floatfield);            // floatfield not set
	std::cout.setf (std::ios::fixed, std::ios::floatfield);   // floatfield set to fixed
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
	//std::cout << "[r, sigma] = [" << params.r << ", " << params.sigma << "]" << std::endl;

	params.nnodes [0] = 25;
	params.nnodes [1] = 5; // y-dimension
	params.nnodes [2] = 5; // x-dimension
    // NOTE: See comments in randnet_dendritic for the following expression for
    // pscale:
    //params.pscale = 0.04 / params.r;
    params.pscale = 1.0 * params.r;
	params.N0 = 0.0;
	params.k0 = 0.5;
	params.maxSims = 100;
    /* Maximal number of attempts at generating unique combination of structural
     * variables before a stream is accepted regardless.
	 * 
     * See _dendritic version for expanation of the changes from former version
     * with nshuffles around the same network, to this version with iteration of
     * 101 values of 0 <= alpha0 <= 1.  This _lattice file had nTrials = 500
     * with nshuffles = 50 for each, so just changes to nTrials = 250 with the
     * 101 values of alpha0.
	 */
	params.nTrials = 250; 
    /* Important NOTE: The configurational space is limited, to around 6 radii
     * in [3,8] times 15 numbers of junction nodes in [6,20] times 9 values of
     * N.I. in [0,8] = 810 points, so nTrials can't approach this value! */

	std::stringstream ss;
	ss.str (""); 
    ss << round (100.0 * params.r);
	fname = "results_lattice_r";
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

	barr3 varcheck (boost::extents [9] [9] [23]); // [radius] [N.I.] [N+]
	for (int i=0; i<9; i++)
		for (int j=0; j<9; j++)
			for (int k=0; k<23; k++)
				varcheck [i] [j] [k] = false;

	out_file.open (fname.c_str(), std::ofstream::out);
	out_file << "alpha,\tconn,\tdiffq,\tradius,\tnloops,\tedgelength,\t" <<
		"n1,\tn2,\tn3,\tn4,\tn5,\tn6,\tn7,\tn8,\tnmn,\tnsd" << std::endl;
	for (int i=0; i<params.nTrials; i++)
    {
		results = do1trial (params, &generator, &varcheck);
		for (int j=0; j<101; j++)
        {
			out_file << (double) j / 100.0 << ",\t" << results.connectivity [j] << ",\t" << 
				results.diff_q [j] << ",\t" << results.radius << ",\t" <<
				results.nloops << ",\t" << results.edgelength << ",\t";
			for (int k=0; k<8; k++)
				out_file << results.edgedist [k] << ",\t";
			out_file << results.N_mn [j] << ",\t" << results.N_sd [j] << std::endl;
		} // end for j

        dtemp = ((double) clock () - (double) time_start) / 
            (double) CLOCKS_PER_SEC;
        progress = (double) i / (double) params.nTrials;
        progLine (progress, dtemp);
	} // end for i
	out_file.close();

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
	 */
	int count, tempi, n38;
	double tempd, diffq, n_mn, n_sd, n_mn0, n_sd0, connectivity;
	NetStructure *NetPt = NULL;
	NetResults results;
	NetResults_OneNet results_oneNet;

	boost::numeric::ublas::matrix<double> pmat (params.nnodes[0], 
            params.nnodes[0]);
	boost::numeric::ublas::matrix<double> Nsums (2, params.nnodes[0] + 1);

	NetPt = new NetStructure[params.nnodes[0]];
	randnet (NetPt, params, generator);
	results_oneNet = netstats (NetPt, params);
	n38 = 0;
	for (int i=2; i<9; i++) 
        n38 += results_oneNet.edgedist [i];
	while (results_oneNet.radius > 8 || results_oneNet.edgelength > 8 || n38 > 22)
    {
		randnet (NetPt, params, generator);
		results_oneNet = netstats (NetPt, params);
		n38 = 0;
		for (int i=2; i<9; i++) 
            n38 += results_oneNet.edgedist [i];
	}
	if (!(*varcheck) [results_oneNet.radius] [results_oneNet.edgelength] [n38])
    {
		(*varcheck) [results_oneNet.radius] [results_oneNet.edgelength] [n38] = true;
    } else {
		count = 0;
		while ((*varcheck) [results_oneNet.radius] [results_oneNet.edgelength] 
                [n38])
        {
			randnet (NetPt, params, generator);
			results_oneNet = netstats (NetPt, params);
			n38 = 0;
			for (int i=2; i<9; i++) 
                n38 += results_oneNet.edgedist [i];
			while (results_oneNet.radius > 8 || results_oneNet.edgelength > 8 || 
                    n38 > 22)
            {
				randnet (NetPt, params, generator);
				results_oneNet = netstats (NetPt, params);
				n38 = 0;
				for (int i=2; i<9; i++) 
                    n38 += results_oneNet.edgedist [i];
			}
			count++;
			if (!(*varcheck) [results_oneNet.radius] [results_oneNet.edgelength] 
                    [n38])
            {
				(*varcheck) [results_oneNet.radius] [results_oneNet.edgelength] 
                    [n38] = true;
				break;
            }
			if (count > params.maxSims) 
                break;
		} // end while
	} // end else iterate

	results_oneNet = netstats (NetPt, params);
	for (int i=0; i<8; i++)
		results.edgedist [i] = results_oneNet.edgedist [i];
	results.radius = results_oneNet.radius;
	results.edgelength = results_oneNet.edgelength;

	for (int i=0; i<101; i++)
    {
		params.alpha0 = (double) i / 100.0;
		fillqc (NetPt, params, generator, true); // true makes rand-unif values
		makepmat (NetPt, &pmat);
        // Calculate mean difference in adjacent q-values. For lattice networks,
        // this counts each neighbour pair twice, but the result is identical to
        // if they were only counted once.
		tempd = 0.0;
		count = 0;
		for (int j=0; j<params.nnodes [0]; j++)
        {
			for (int k=0; k<8; k++)
            {
				tempi = NetPt [j].nextnode [k];
				if (tempi > INT_MIN)
                {
					tempd += fabs (NetPt [j].q - NetPt [tempi].q);
					count++;
                }
			} // end for k
		} // end for j
		results.diff_q [i] = tempd / (double) count;
		// And then connectivity:
		tempd = 0.0;
		for (int j=0; j<params.nnodes [0]; j++)
        {
			NetPt [j].N = 1.0;
			tempd += pmat (j, j); 
            // Not scaled to pscale, as it actually is for movement
		}
		results.connectivity [i] = ((double) params.nnodes [0] - tempd) / 
            (double) params.nnodes [0];

		Nsums = populateNet (&pmat, NetPt, params, generator);
		results.N_mn [i] = Nsums (0, params.nnodes [0]);
		results.N_sd [i] = Nsums (1, params.nnodes [0]);
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
	/*
     * Each node has .nextnode, .xpos, and .ypos. xpos & ypos count
     * incrementally from bottom left, across the horizontal, then up to top
     * right.  Directions [0..7] go clockwise from 0 = SW. The xshift and yshift
     * are then the values needed to move *away* from these directions.
	 *
     * Networks are built with a single node index mapped onto the 2-D grid.
     * While using two indices would make many aspects simpler, it's still much
     * easier with one index not just to trace connections, but to store the
     * lists of connected nodes. The main drawback is the extensive fiddling at
     * the outset to make the nextnode[nnodes][8] matrix.
	 */
	const int xshift [8] = {1,1,1,0,-1,-1,-1,0};
	const int yshift [8] = {1,0,-1,-1,-1,0,1,1};
	const int nodeback[8] = {4,5,6,7,0,1,2,3}; 
    // nodeback are reverse direction indices of [0..7]
	int hmi, tempi[2], tempj[2], fromi, toi;
	double tempd;
	bool connected, found[2];

	boost::uniform_real<> uni_dist(0,1);
	boost::variate_generator<base_generator_type&,
		boost::uniform_real<> > runif((*generator), uni_dist);

	for (int i=0; i<params.nnodes[0]; i++)
		for (int j=0; j<8; j++)
			NetPt[i].conn[j] = DOUBLE_MIN;
	/*tempd = 0.0;
	for (int i=0; i<params.nnodes[0]; i++)
    {
		NetPt[i].q = minqc + (1.0 - 2.0 * minqc) * runif();
		tempd += NetPt[i].q;
		for (int j=0; j<8; j++)
			NetPt[i].conn[j] = DOUBLE_MIN;
    }
	// Standardise qualities to mean value of 0.5
	for (int i=0; i<params.nnodes[0]; i++)
		NetPt[i].q = NetPt[i].q * 0.5 * (double) params.nnodes[0] / tempd;
	*/

	// Then set x & y positions on the square grid.
	// nnodes[1] is x; nnodes[2] is y
	tempi[0] = 0;
	for (int i=0; i<params.nnodes[2]; i++)
    {
		for (int j=0; j<params.nnodes[1]; j++)
        {
			NetPt[tempi[0]].xpos = j;
			NetPt[tempi[0]].ypos = i;
			tempi[0]++;
        }
    }

	// Then set up the nextnode matrix that indexes all connections
	// Start by making them all nix, then only set those that actually connect.
	for (int i=0; i<params.nnodes[0]; i++) 
        for (int j=0; j<8; j++)
            NetPt[i].nextnode[j] = INT_MIN;
	// Connections are filled individually in numeric order starting at [0]=SW.
	for (int i=1; i<params.nnodes[2]; i++) 
    { 
        for (int j=1; j<params.nnodes[1]; j++) 
        {
            tempi[0] = i * params.nnodes[1] + j; // from index
            NetPt[tempi[0]].nextnode[0] = tempi[0] - params.nnodes[1] - 1;
        }
    }
	for (int i=0; i<params.nnodes[2]; i++) 
    { 
        for (int j=1; j<params.nnodes[1]; j++) 
        { // [1] = W
            tempi[0] = i * params.nnodes[1] + j;
            NetPt[tempi[0]].nextnode[1] = tempi[0] - 1;
        }
    }
	for (int i=0; i<(params.nnodes[2]-1); i++) 
    { 
        for (int j=1; j<params.nnodes[1]; j++) 
        { // [2] = NW
            tempi[0] = i * params.nnodes[1] + j; 
            NetPt[tempi[0]].nextnode[2] = tempi[0] + params.nnodes[1] - 1;
        }
    }
	for (int i=0; i<(params.nnodes[2]-1); i++) 
    { 
        for (int j=0; j<params.nnodes[1]; j++) 
        { // [3] = N
            tempi[0] = i * params.nnodes[1] + j;
            NetPt[tempi[0]].nextnode[3] = tempi[0] + params.nnodes[1];
        }
    }
  	for (int i=0; i<(params.nnodes[2]-1); i++) 
    { 
        for (int j=0; j<(params.nnodes[1]-1); j++) { 
            // [4] = NE
            tempi[0] = i * params.nnodes[1] + j;
            NetPt[tempi[0]].nextnode[4] = tempi[0] + params.nnodes[1] + 1;
        }
    }
	for (int i=0; i<params.nnodes[2]; i++) { for 
        (
         int j=0; j<(params.nnodes[1]-1); j++) 
        { // [5] = E
            tempi[0] = i * params.nnodes[1] + j; 
            NetPt[tempi[0]].nextnode[5] = tempi[0] + 1;
        }
    }
	for (int i=1; i<params.nnodes[2]; i++) { for 
        (
         int j=0; j<(params.nnodes[1]-1); j++) { 
            // [6] = SE
            tempi[0] = i * params.nnodes[1] + j; 
            NetPt[tempi[0]].nextnode[6] = tempi[0] - params.nnodes[1] + 1;
        }
    }
	for (int i=1; i<params.nnodes[2]; i++) { for 
        (
         int j=0; j<params.nnodes[1]; j++) 
        { // [7] = S
            tempi[0] = i * params.nnodes[1] + j; 
            NetPt[tempi[0]].nextnode[7] = tempi[0] - params.nnodes[1];
        }
    }

	// Then the hold Matrix which stores each disconnected segment until they become connected.
	boost::numeric::ublas::matrix<int> holdMatrix(params.nnodes[0], 
            params.nnodes[0]);
	for (int i=0; i<params.nnodes[0]; i++)
		for (int j=0; j<params.nnodes[0]; j++)
			holdMatrix(i, j) = INT_MIN;
	hmi = 0;

	connected = false;
	while (!connected)
    {
		fromi = (int) (runif() * params.nnodes[0]);
		tempi[0] = 0; // Counts the number of available positions
		for (int i=0; i<8; i++)
			if (NetPt[fromi].nextnode[i] > INT_MIN && NetPt[fromi].conn[i] == DOUBLE_MIN) 
                tempi[0]++;
		while (tempi[0] == 0) 
        { // Just in case a full position was picked.
			fromi = (int) (runif() * params.nnodes[0]);
			tempi[0] = 0; // Counts the number of available positions
			for (int i=0; i<8; i++) 
				if (NetPt[fromi].nextnode[i] > INT_MIN && 
                        NetPt[fromi].conn[i] == DOUBLE_MIN) 
                    tempi[0]++;
        }
		// tempi[0] is now the number of accessible nodes, so that ... */
		tempi[0] = (int) (runif() * tempi[0]);
		// will directly reference one of these, and ...
		tempj[0] = 0; 
		for (int i=0; i<8; i++)
        {
			if (tempj[0] == tempi[0] && NetPt[fromi].nextnode[i] > INT_MIN && 
                    NetPt[fromi].conn[i] == DOUBLE_MIN)
            {
				tempi[1] = i;
				break;
            }
			if (NetPt[fromi].nextnode[i] > INT_MIN && NetPt[fromi].conn[i] == DOUBLE_MIN) 
                tempj[0]++;
		} // end for i
		// tempi[1] is now the new [1..8] direction in terms of xshift, etc.
		toi = NetPt[fromi].nextnode[tempi[1]];
		NetPt[fromi].conn[tempi[1]] = minqc + (1.0 - 2.0 * minqc) * runif();
		NetPt[toi].conn[nodeback[tempi[1]]] = minqc + (1.0 - 2.0 * minqc) * runif();

		// Switch nextnodes off to stop crossing over
		if (tempi[1] == 0)
        {
			NetPt[NetPt[fromi].nextnode[7]].nextnode[2] = INT_MIN;
			NetPt[NetPt[fromi].nextnode[1]].nextnode[6] = INT_MIN;
        } else if (tempi[1] == 2) {
			NetPt[NetPt[fromi].nextnode[1]].nextnode[4] = INT_MIN;
			NetPt[NetPt[fromi].nextnode[3]].nextnode[0] = INT_MIN;
        } else if (tempi[1] == 4) {
			NetPt[NetPt[fromi].nextnode[3]].nextnode[6] = INT_MIN;
			NetPt[NetPt[fromi].nextnode[5]].nextnode[2] = INT_MIN;
        } else if (tempi[1] == 6) {
			NetPt[NetPt[fromi].nextnode[5]].nextnode[0] = INT_MIN;
			NetPt[NetPt[fromi].nextnode[7]].nextnode[4] = INT_MIN;
        }

		// Then do the holdMat stuff to check is the network is connected
		if (hmi == 0)
        {
			if (holdMatrix(0, 0) == INT_MIN)
            {
				holdMatrix(0, 0) = fromi;
				holdMatrix(0, 1) = toi;
				hmi++;
            } else {
				tempi[0] = 0; // tempi index the holdMat[0] vector
				found[0] = false; found[1] = false;
				while (holdMatrix(0, tempi[0]) > INT_MIN && 
                        tempi[0] < params.nnodes[0])
                {
					if (holdMatrix(0, tempi[0]) == fromi)
						found[0] = true;
					if (holdMatrix(0, tempi[0]) == toi)
						found[1] = true;
					tempi[0]++;
                }
				if (found[0] && !found[1])
                {
					holdMatrix(0, tempi[0]) = toi;
                } else if (!found[1] && found[1]) {
					holdMatrix(0, tempi[0]) = fromi;
                } else if (!found[0] && !found[1]) {
					holdMatrix(1, 0) = fromi;
					holdMatrix(1, 1) = toi;
					hmi++;
                } else if (found[0] && found[1] && tempi[0] == params.nnodes[0]) {
					connected = true;
                }
			} // end else holdMatrix[0][0] != nix
        } else { // hmi != 0
			tempi[0] = INT_MIN;
			for (int i=0; i<hmi; i++)
            {
				tempj[1] = INT_MIN; // Dummy counter at this point
				found[0] = false; found[1] = false;
				for (int j=0; j<params.nnodes[0]; j++)
                {
					if (holdMatrix(i, j) > INT_MIN)
                    {
						if (holdMatrix(i, j) == fromi) 
                            found[0] = true;
						else if (holdMatrix(i, j) == toi) 
                            found[1] = true;
                    }
                } // end for j
				tempj[0] = 0;
				while (holdMatrix(i, tempj[0]) > INT_MIN && 
                        tempj[0] < params.nnodes[0])
					tempj[0]++;
				if (found[0] && found[1]) 
                { 
                    break;
                } else if (found[0] && !found[1]) {
					holdMatrix(i, tempj[0]) = toi;
					tempi[0] = i; 
                    break;
                } else if (!found[0] && found[1]) {
					holdMatrix(i, tempj[0]) = fromi;
					tempi[0] = i; 
                    break;
                }
			} // end for i
			if (tempi[0] == INT_MIN && tempj[1] == INT_MIN && !found[0] && !found[1])
            {
				holdMatrix(hmi, 0) = fromi;
				holdMatrix(hmi, 1) = toi;
				hmi++;
            } else if ((found[0] && !found[1]) || (!found[0] && found[1]))
            {
				/* 
                 * Value was added to an existing holdMat vector, so now see if
                 * that modified vector shares entires with previous columns,
                 * and if so, pack em on the end and delete the previously
                 * isolated one.
				 *
                 * Because this is check for every single addition, the whole
                 * column apart from the newly added value should be unique, and
                 * is thus added in its entirety.
				 *
                 * Indices are used here such that the new value was stored
                 * above at tempi[0]tempj[0], and, when this value is found in
                 * another vector, then all values in tempi[0] from 0 to
                 * (tempj[0]-1) are added onto that other vector. The latter is
                 * then indexed by tempi[1]tempj[1].
				 */
				tempi[1] = INT_MIN;
				for (int i=0; i<hmi; i++)
                {
					if (i != tempi[0])
                    {
						for (int j=0; j<params.nnodes[0]; j++)
                        {
							if (holdMatrix(i, j) > INT_MIN && 
								holdMatrix(i, j) == holdMatrix(tempi[0], tempj[0]))
                            {
									tempi[1] = i; break;
                            }
                        } // end for j
						if (tempi[1] > INT_MIN) 
                            break;	
                    }
                } // end for i
				if (tempi[1] > INT_MIN)
                {
					// Pack the new tempi[0] vector onto this tempi[1] one.
					tempj[1] = 0;
					while (holdMatrix(tempi[1], tempj[1]) > INT_MIN && 
                            tempj[1] < params.nnodes[0])
						tempj[1]++;
					for (int i=0; i<(tempj[0]); i++)
						holdMatrix(tempi[1], tempj[1] + i) = 
                            holdMatrix(tempi[0], i);
                    // Then delete the tempi[0] vector, and if it's not the hmi
                    // vector, then move [mhi] to tempi[0]
					if (tempi[0] < (hmi - 1))
                    {
						for (int i=0; i<params.nnodes[0]; i++)
							holdMatrix(tempi[0], i) = holdMatrix(hmi - 1, i);
						tempi[0] = hmi - 1;
                    }
					for (int i=0; i<params.nnodes[0]; i++)
						holdMatrix(tempi[0], i) = INT_MIN;
					hmi--;
                    // And then finally check new holdMat vector to see if it
                    // has all nodes ...
                    tempi[0] = 0;
					while (holdMatrix(tempi[1], tempi[0]) > INT_MIN && 
                            tempi[0] < params.nnodes[0])
						tempi[0]++;
					if (tempi[0] == params.nnodes[0]) 
                        connected = true;
				} // end if tempi[1] > nix
			} // end else tempi!=nix && tempj!=nix - that is, end shuffling holdMat
		} // end else analyse holdMat
		// Then check is all are connected
		if (hmi == 1)
        {
			tempi[0] = 0;
			for (int i=0; i<params.nnodes[0]; i++)
				if (holdMatrix(0, i) > INT_MIN) 
                    tempi[0]++;
			if (tempi[0] >= params.nnodes[0]) 
                connected = true;
        } // end if hmi == 1
	} // end while connected == false

	// And finally nextnode to nix for all non-connected nodes:
	for (int i=0; i<params.nnodes[0]; i++) 
        for (int j=0; j<8; j++)
            if (NetPt[i].conn[j] == DOUBLE_MIN) 
                NetPt[i].nextnode[j] = INT_MIN;
}; // end function randnet


/************************************************************************
 ************************************************************************
 **                                                                    **
 **                       FILLQC FUNCTION                              **
 **                                                                    **
 ************************************************************************
 ************************************************************************/

void fillqc(NetStructure* NetPt, Parameters params, 
        base_generator_type * generator, bool unif)
{
	/*
     * pseq can make exponential distributions centred around 0.5 (q=true), or
     * decreasing from 1 (q=false), with the latter potentially used for
     * connectivities so that most of them end up with relatively high values.
     * This, however, seems to produce much lower correlations with structural
     * network properties, so connectivities are now reverted to the q-dist.
	 *
     * To change this back, just ditch the even-to-odd conversion bit, and set >
     * clist = pseq(count, false);
	 * 
     * NOTE very importantly - this entire routine has been disabled here, so
     * that q & c values are simply drawn from a uniform distribution in [0, 1].
	 */
	int count;

    if (unif)
    {
        boost::uniform_real<> uni_dist(0,1);
        boost::variate_generator<base_generator_type&,
            boost::uniform_real<> > runif((*generator), uni_dist);
        boost::normal_distribution<> norm_dist(0,1);
        boost::variate_generator<base_generator_type&,
            boost::normal_distribution<> > rnorm((*generator), norm_dist);

        if (params.alphasd == 0.0)
        {
            for (int i=0; i<params.nnodes[0]; i++)
            {
                NetPt [i].q = runif ();
                for (int j=0; j<8; j++)
                    if (NetPt[i].nextnode[j] > INT_MIN)
                        NetPt[i].conn[j] = runif ();
            } // end for i
        } else { // if alphasd != 0.0
            for (int i=0; i<params.nnodes[0]; i++)
            {
                NetPt [i].q = params.k0 + params.alphasd * rnorm ();
                while (NetPt [i].q < minqc || NetPt [i].q > (1.0 - minqc))
                    NetPt [i].q = params.k0 + params.alphasd * rnorm ();
                for (int j=0; j<8; j++)
                {
                    if (NetPt[i].nextnode[j] > INT_MIN)
                    {
                        NetPt [i].conn [j] = params.alpha0 + params.alphasd * rnorm ();
                        while (NetPt [i].conn [j] < minqc || 
                                NetPt [i].conn [j] > (1.0 - minqc))
                            NetPt [i].conn [j] = params.alpha0 + 
                                params.alphasd * rnorm ();
                    }
                } // end for j
            } // end for i
        } // end else alpha != 0.0
    } else { // not unif, so generate a normal distribution
        std::vector <int> indx = randseq (params.nnodes[0], generator);
        boost::numeric::ublas::vector<double> q (params.nnodes[0]);
        pseq (&q, true);
        count = 0;
        for (int i=0; i<params.nnodes[0]; i++)
        {
            NetPt[indx[i]].q = q(i);
            for (int j=0; j<8; j++)
                if (NetPt[i].nextnode[j] > INT_MIN) 
                    count++;
        }
        // pseq for (q=true) has to have an odd-numbered length.
        if (floor((double) count / 2.0) == ((double) count / 2.0))
            count++;
        std::vector <int> indx2 = randseq (count, generator);
        boost::numeric::ublas::vector<double> clist (count);
        pseq (&clist, true);
        boost::numeric::ublas::vector<double> clist2 (count);
        for (int i=0; i<count; i++)
            clist2 (i) = clist (indx2 [i]);
        count = 0;
        for (int i=0; i<params.nnodes[0]; i++)
        {
            for (int j=0; j<8; j++)
            {
                if (NetPt[i].nextnode[j] > INT_MIN)
                {
                    NetPt[i].conn[j] = clist2(count);
                    count++;
                }
            } // end for j
        } // end for i
    } // end else !unif

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
	int tempi[4];
	double tempd;
	NetResults_OneNet results_oneNet;
	
	// Calculate radius
	boost::numeric::ublas::matrix<int> dmat(params.nnodes[0], params.nnodes[0]);
	dmat = makedmat(NetPt, params);

	tempi [0] = 9999;
	for (int i=0; i<params.nnodes [0]; i++)
    {
		tempi [1] = 0;
		for (int j=0; j<params.nnodes [0]; j++)
			if (dmat (i, j) > tempi [1]) 
                tempi [1] = dmat (i, j);
		if (tempi [1] < tempi [0]) 
            tempi [0] = tempi [1];
	}
	results_oneNet.radius = tempi [0];

	/* Then fill edge distribution and calculate nodal isolation. For the latter,
	 * 	tempi[1] is always where trace came from;
	 * 	tempi[2] is the node being traced; and
	 * 	tempi[3] holds the next node (if there turns out to be only 2 edges).*/
	for (int i=0; i<8; i++) 
        results_oneNet.edgedist [i] = 0;
	results_oneNet.edgelength = 0;
	for (int i=0; i<params.nnodes [0]; i++)
    {
		tempi [0] = 0;
		for (int j=0; j<8; j++)
        {
			if (NetPt [i].nextnode [j] > INT_MIN)
            {
				tempi [0]++;
				tempi [1] = i; // Used for tracing edgelengths below
				tempi [2] = NetPt [i].nextnode [j];
			}
		} // end for j
		results_oneNet.edgedist [tempi [0] - 1]++;

		if (tempi [0] == 1) 
        { // trace edges away from terminal nodes
			results_oneNet.edgelength++;
			tempi [0] = 0;
			for (int j=0; j<8; j++) 
            {
				if (NetPt [tempi [2]].nextnode [j] > INT_MIN) 
                {
					tempi [0]++;
					if (NetPt [tempi [2]].nextnode [j] != tempi [1]) 
						tempi [3] = NetPt [tempi [2]].nextnode [j];
				}
			} // end for j
			while (tempi [0] == 2) 
            {
				results_oneNet.edgelength++;
				tempi [1] = tempi [2];
				tempi [2] = tempi [3];
				tempi [0] = 0;
				for (int j=0; j<8; j++) 
                {
					if (NetPt [tempi [2]].nextnode [j] > INT_MIN) 
                    {
						tempi [0]++;
						if (NetPt [tempi [2]].nextnode [j] != tempi [1])
							tempi [3] = NetPt [tempi [2]].nextnode [j];
					}
				} // end for j
			} // end while
		}
	}
	results_oneNet.nloops = countloops (NetPt, params);

	return results_oneNet;
} // end function netstats



/************************************************************************
 ************************************************************************
 **                                                                    **
 **                      COUNTLOOPS FUNCTION                           **
 **                                                                    **
 ************************************************************************
 ************************************************************************/

int countloops(const NetStructure *NetPt, Parameters params)
{
/*
 * Counts number of loops of sizes 3 and 4, presuming any loops of larger
 * sizes will actually encompass smaller ones. This works because 4-loops
 * cannot contain any 3-loops as sub-graphs, but 5 loops etc can.
 */
	int loopcount, next1, next2, next3;
	bool same;

	boost::numeric::ublas::matrix<int> looplist(params.nnodes[0] * 10, 4);
	for (int i=0; i<(params.nnodes[0]*10); i++)
		for (int j=0; j<4; j++)
			looplist(i, j) = INT_MIN;
	loopcount = 0;
	boost::numeric::ublas::vector<int> oneloop(4);

	// Count 3-loops first
	for (int i=0; i<params.nnodes[0]; i++) 
    {
		for (int dirj=0; dirj<8; dirj++) 
        {
			next1 = NetPt[i].nextnode[dirj];
			if (next1 > INT_MIN && next1 != i) 
            {
				for (int dirk=0; dirk<8; dirk++) 
                {
					next2 = NetPt[next1].nextnode[dirk];
					if (next2 > INT_MIN && !(next2 == i || next2 == next1)) 
                    {
						for (int j=0; j<8; j++) 
                        {
							if (NetPt[next2].nextnode[j] == i) 
                            {
								// Add loop to looplist
								oneloop(0) = i;
								oneloop(1) = next1;
								oneloop(2) = next2;
								oneloop(3) = 9999;
								oneloop = sort(oneloop,3);
								if (loopcount > 0) 
                                {
									for (int k=0; k<loopcount; k++) 
                                    {
										same = true;
										for (int m=0; m<3; m++) 
                                        {
											if (looplist(k, m) != oneloop(m)) 
                                            {
												same = false;	
                                            }
										} // end for m
										if (same) 
                                            break;
									} // end for k
								} // end if loopcount
								if (!same || loopcount == 0) 
                                {
									for (int k=0; k<3; k++) 
                                    {
										looplist(loopcount, k) = oneloop(k);
									}
									loopcount++;
								}
							} // end if nextnode[j] == i
						} // end for j over 8 directions
					} // end if next2 > INT_MIN
				} // end for dirk over 8 directions
			} // end if .nextnode[dirj] > INT_MIN
		} // end for dirj over 8 directions
	} // end for i over nnodes

	// Then count 4-loops
	for (int i=0; i<params.nnodes[0]; i++) 
    {
		for (int dirj=0; dirj<8; dirj++) 
        {
			next1 = NetPt[i].nextnode[dirj];
			if (next1 > INT_MIN && next1 != i) 
            {
				for (int dirk=0; dirk<8; dirk++) 
                {
					next2 = NetPt[next1].nextnode[dirk];
					if (next2 > INT_MIN && !(next2 == i || next2 == next1)) 
                    {
						for (int dirm=0; dirm<8; dirm++) 
                        {
							next3 = NetPt[next2].nextnode[dirm];
							if (next3 > INT_MIN && 
                                    !(next3 == i || next3 == next2 || 
                                        next3 == next1)) 
                            {
								for (int j=0; j<8; j++) 
                                {
									if (NetPt[next3].nextnode[j] == i) 
                                    {
										// Add loop to looplist
										oneloop(0) = i;
										oneloop(1) = next1;
										oneloop(2) = next2;
										oneloop(3) = next3;
										oneloop = sort(oneloop,4);
										for (int k=0; k<loopcount; k++) 
                                        {
											same = true;
											for (int m=0; m<4; m++) 
                                            {
												if (looplist(k, m) != oneloop(m)) 
                                                {
													same = false;	
                                                }
											} // end for m
											if (same) 
                                                break;
										} // end for k
										if (!same) 
                                        {
											for (int k=0; k<4; k++) 
												looplist(loopcount, k) = oneloop(k);
											loopcount++;
										} // end if !same
									} // end if nextnode[j] == i
								} // end for j over 8 directions
							} // end if next3 > INT_MIN
						} // end for dirm over 8 directions
					} // end if next2 > INT_MIN
				} // end for dirk over 8 directions
			} // end if .nextnode[dirj] > INT_MIN
		} // end for dirj over 8 directions
	} // end for i over nnodes

	return --loopcount;
} // end function countloops


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
     * total population.*/
	const int nIters = 1000;
	double qval, tempd[2];

	boost::normal_distribution<> norm_dist(0,1);
	boost::variate_generator<base_generator_type&,
		boost::normal_distribution<> > rnorm((*generator), norm_dist);

	boost::numeric::ublas::matrix<double> oldN (2, params.nnodes [0]);
	for (int i=0; i<params.nnodes [0]; i++) 
        oldN (0, i) = NetPt[i].N;
	// Nsums holds values for each node, plus network means for abundance & diversity.
	boost::numeric::ublas::matrix<double> Nsums (2, params.nnodes [0] + 1);
	for (int i=0; i<Nsums.size1(); i++)
		for (int j=0; j<Nsums.size2(); j++)
			Nsums (i, j) = 0.0;

	const int runin = 100;
	for (int n=0; n<runin; n++) 
    {
		/* Movement:
         * The proportion that stay at [i][i] is modified from
         * NetPt[i].N*pmat[i][i], because pmat dynamics are scaled by pscale.
         * Thus instead of sum(b)+a=1, where a=pmat[i][i], one seeks the value
         * of a' such that sum(pscale*b)+a'=1.  This is then
         * a'=1-sum(pscale*b)=1-pscale*sum(b)=1-pscale*(1-a). */
		for (int i=0; i<params.nnodes [0]; i++) 
        {
			oldN (1, i) = oldN (0, i) * 
                (1.0 - params.pscale * (1.0 - (*pmat) (i, i)));
			for (int j=0; j<params.nnodes [0]; j++) 
                if (j != i)
                    oldN (1, i) += oldN (0, j) * params.pscale * (*pmat) (j, i);
		}
		for (int i=0; i<params.nnodes [0]; i++) 
        {
			qval = NetPt[i].q + params.sigma * rnorm();
			while (qval < minqc || qval > (1.0 - minqc))
				qval = NetPt[i].q + params.sigma * rnorm();
			oldN (0, i) = oldN (1, i) + oldN (1, i) * params.r *
				(1.0 - oldN (1, i) / qval) * 
                (oldN (1, i) / qval - params.N0 / qval);
			if (oldN (0, i) < 0.0) 
                oldN (0, i) = 0.0;
		} // end for i
	} // end for n

	for (int n=0; n<nIters; n++) 
    {
		for (int i=0; i<params.nnodes [0]; i++) 
        {
			oldN (1, i) = oldN (0, i) * (1.0 - params.pscale * (1.0 - (*pmat) (i, i)));
			for (int j=0; j<params.nnodes[0]; j++) 
                if (j != i)
                    oldN (1, i) += oldN (0, j) * params.pscale * (*pmat) (j, i);
		} // end for i
		for (int i=0; i<params.nnodes [0]; i++) 
        {
			qval = NetPt[i].q + params.sigma * rnorm();
			while (qval < minqc || qval > (1.0 - minqc))
				qval = NetPt[i].q + params.sigma * rnorm();
			oldN (0, i) = oldN (1, i) + oldN (1, i) * params.r *
				(1.0 - oldN (1, i) / qval) * (oldN (1, i) / qval - params.N0 / qval);
			//if (oldN (0, i) < 0.0) { oldN (0, i) = 0.0;	}
		} // end for i

		// Abundance calculations
		tempd [0] = 0.0;
		for (int i=0; i<params.nnodes [0]; i++) 
        {
			Nsums (0, i) += oldN (0, i);
			Nsums (1, i) += oldN (0, i) * oldN (1, i);
			tempd [0] += oldN (0, i);
        }
		Nsums (0, params.nnodes [0]) += tempd [0];
		Nsums (1, params.nnodes [0]) += tempd [0] * tempd [0];
		// Then only set all negative abundances to zero after calculations
		for (int i=0; i<params.nnodes [0]; i++)
			if (oldN (0, i) < 0.0) 
                oldN (0, i) = 0.0;
	} // end for n

	// Convert Nsums to means & SDs
	for (int i=0; i<=params.nnodes[0]; i++) 
    {
		Nsums (0, i) = Nsums (0, i) / (double) nIters;
		Nsums (1, i) = Nsums (1, i) / (double) nIters - Nsums (0, i) * Nsums (0, i);
		Nsums (1, i) = sqrt (Nsums (1, i));
	}

	return Nsums;
}; // end function popnet


/************************************************************************
 ************************************************************************
 **                                                                    **
 **                   POPULATENET_EQ FUNCTION                          **
 **                                                                    **
 ************************************************************************
 ************************************************************************/

dmat populateNet_eq (dmat * pmat, Parameters params)
{
	const int maxiter = 10000;
	const double tol = 1.0e-4;
	int count;
	double qval, tempd[2], tolsum[2];

	boost::numeric::ublas::matrix<double> oldN (2, params.nnodes[0]);
	for (int i=0; i<params.nnodes[0]; i++) 
        oldN (0, i) = 1.0;

	// Nsums holds values for each node, plus network means for abundance & diversity.
	boost::numeric::ublas::matrix<double> Nsums (2, params.nnodes[0] + 1);
	for (int i=0; i<Nsums.size1(); i++) 
		for (int j=0; j<Nsums.size2(); j++)
			Nsums (i, j) = 0.0;

	tolsum[0] = 99999.9; tolsum[1] = 0.0;
	count = 0;
	while (fabs(tolsum[1] - tolsum[0]) > tol) 
    {
		tolsum[0] = tolsum[1];
		tolsum[1] = 0.0;
		for (int i=0; i<params.nnodes[0]; i++) 
        {
			oldN (1, i) = 0.0;
			for (int j=0; j<params.nnodes[0]; j++) 
            {
				if (j == i)
					oldN (1, i) += oldN (0, i) * (1.0 - params.pscale *
						(1.0 - (*pmat) (i, i)));
				else
					oldN (1, i) += oldN (0, j) * params.pscale * (*pmat) (j, i);
			} // end for j
			tolsum[1] += oldN (1, i);
		} // end for i
		for (int i=0; i<params.nnodes[0]; i++) 
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
        Nsums (i, params.nnodes[0]) = 0.0;
	for (int i=0; i<params.nnodes[0]; i++) 
    {
		Nsums (0, i) = oldN (1, i);
		Nsums (0, params.nnodes[0]) += oldN (1, i);
		Nsums (1, params.nnodes[0]) += oldN (1, i) * oldN (1, i);
	}
	Nsums (0, params.nnodes[0]) = Nsums (0, params.nnodes[0]) / (double) params.nnodes[0];
	Nsums (1, params.nnodes[0]) = Nsums (1, params.nnodes[0]) / (double) params.nnodes[0] -
		Nsums (0, params.nnodes[0]) * Nsums (0, params.nnodes[0]);

	return Nsums;
}; // end function popnet_eq



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
			dist(j) = 1e99;
			Q(j) = j;
			prev(j) = INT_MIN;
			done(j) = false;
        } // end for j
		here = i;
		dist(here) = 0.0;
		lenQ = 19;

		while (lenQ > 0) 
        {
			dmin = 1e99; here = INT_MIN;
			for (int j=0; j<nnodes; j++) 
            { 
                if (!done(j) && dist(j) < dmin) 
                {
				dmin = dist(j); here = j;	
                }	
            }
			if (here == INT_MIN)
				std::cout << "ERROR: here == nix!" << std::endl;
			Q(here) = INT_MIN;
			done(here) = true;
			// Neighbour relaxation, i.e. update of shortest distances
			for (int j=0; j<8; j++) 
            {
				if (NetPt[here].nextnode[j] > INT_MIN) 
                {
					if (dist(here) == 0.0) 
                        dmin = 1.0 / NetPt[here].conn[j];
					else 
                        dmin = dist(here) / NetPt[here].conn[j];
					there = NetPt[here].nextnode[j];
					if (dmin < dist(there)) 
                    {
						dist(there) = dmin;
						prev(there) = here;	
                    }	
                }	
            } // end for j
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
			while (prev(here) != INT_MIN) 
            {
				plen++;
				here = prev(here);	
            }
			if (plen == 0) 
                (*pmat) (i, j) = NetPt[j].q;
			else 
                (*pmat) (i, j) = NetPt[j].q / dist(j);
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
	} // end for i
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

	boost::numeric::ublas::vector<int> dist(params.nnodes[0]);
	boost::numeric::ublas::vector<int> Q(params.nnodes[0]);
	boost::numeric::ublas::vector<int> prev(params.nnodes[0]);
	boost::numeric::ublas::vector<bool> done(params.nnodes[0]);

	boost::numeric::ublas::matrix<int> dmat(params.nnodes[0], params.nnodes[0]);
	for (int i=0; i<params.nnodes[0]; i++) 
		for (int j=0; j<params.nnodes[0]; j++)
			dmat(i, j) = 0;

	for (int i=0; i<params.nnodes[0]; i++) 
    {
		for (int j=0; j<params.nnodes[0]; j++) 
        {
			dist(j) = 9999;
			Q(j) = j;
			prev(j) = INT_MIN;
			done(j) = false;
        }
		here = i;
		dist[here] = 0.0;
		lenQ = 19;

		while (lenQ > 0) 
        {
			dmin = 9999; here = INT_MIN;
			for (int j=0; j<params.nnodes[0]; j++) 
            { 
                if (!done(j) && dist(j) < dmin) 
                {
                    dmin = dist(j); 
                    here = j;	
                }	
            } // end for j
			if (here == INT_MIN)
				std::cout << "ERROR: here == nix!" << std::endl;
			Q(here) = INT_MIN;
			done(here) = true;
			// Neighbour relaxation, i.e. update of shortest distances
			dmin = dist(here) + 1;
			for (int j=0; j<8; j++) 
            {
				if (NetPt[here].nextnode[j] > INT_MIN) 
                {
					there = NetPt[here].nextnode[j];
					if (dmin < dist(there)) 
                    {
						dist(there) = dmin;
						prev(there) = here;	
                    }	
                }	
            } // end for j
			lenQ = 0;
			for (int j=0; j<params.nnodes[0]; j++) 
                if (Q(j) > INT_MIN) 
                    lenQ++;
		} // end while (lenq > 0)

		// Then trace back the prev nodes to contstruct pmat
		for (int j=0; j<params.nnodes[0]; j++) 
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
	for (int i=0; i<params.nnodes[0]; i++) 
    { 
        for (int j=0; j<params.nnodes[0]; j++) 
        {
            if (dmat(i, j) < dmat(j, i)) 
                dmat(j, i) = dmat(i, j);
            else if (dmat(i, j) > dmat(j, i)) 
                dmat(i, j) = dmat(j, i);
        } // end for j
    } // end for i

	return dmat;
} // end function makedmat


