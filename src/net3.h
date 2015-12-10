/*
 * diff_sim3.h
 */


#include <iostream>
#include <math.h>
#include <iomanip> // has setprecision, etc.
#include <fstream> // file in & out
#include <cstdlib> // has rand and srand
#include <time.h> 

#include <boost/random/linear_congruential.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include "boost/multi_array.hpp"

#include <boost/program_options.hpp>

// This is a typedef for a random number generator.
// Try boost::mt19937 or boost::ecuyer1988 instead of boost::minstd_rand
typedef boost::minstd_rand base_generator_type;
typedef boost::numeric::ublas::vector<double> dvec;
typedef boost::numeric::ublas::matrix<double> dmat;
typedef boost::numeric::ublas::vector<int> ivec;
typedef boost::numeric::ublas::matrix<int> imat;
typedef std::vector <double> stdvec;

struct Parameters 
{
    int nTrials;
    double k0, k0sd, ksd, alpha0, alphasd, r;
};

struct Results 
{
    double connectivity, nmn_node [4], nsd_node [4], nmn_network, nsd_network, 
           cov [3] [3];
};

void timeout (double tseconds);
Results runPop (Parameters pars, base_generator_type * generator);

const int inix = -9999, runin = 100;
const double minqc = 0.001, dnix = -9999.0;
// minqc is important to prevent dividing by sometimes *really* small numbers
