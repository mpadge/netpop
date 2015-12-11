/*
 * utils.cc
 *
 */

#include "utils_dendlatt.h"


/************************************************************************
 ************************************************************************
 **                                                                    **
 **                       RANDSEQ FUNCTION                             **
 **                                                                    **
 ************************************************************************
 ************************************************************************/

std::vector <int> randseq (int n, base_generator_type * generator)
{
    /* NOTE that this could be done with pointers, but the pair vector has to be
     * copied anyway, and one more copy of what is always a relative short
     * vector is unlikely to really make any difference.
	 */
	boost::uniform_real<> uni_dist(0,1);
	boost::variate_generator<base_generator_type&,
		boost::uniform_real<> > runif((*generator), uni_dist);

	std::vector <std::pair <double, int> > junkvec;
	for (int i=0; i<n; i++)
		junkvec.push_back (std::pair <double, int> (runif(), i));
	std::sort (junkvec.begin(), junkvec.end());
	std::vector <int> vecout;
	std::vector <std::pair <double, int> >::const_iterator itr;
	for (itr = junkvec.begin(); itr < junkvec.end(); itr++)
    {
		vecout.push_back (itr -> second);
	}
	return(vecout);
} // end function randseq


/************************************************************************
 ************************************************************************
 **                                                                    **
 **                         PSEQ FUNCTION                              **
 **                                                                    **
 ************************************************************************
 ************************************************************************/

void pseq (dvec * pvals, bool q)
{
    /* Generates a sequence of evenly-spaced doubles along an inverted
     * exponential distribution. That is, there are more of them near the middle
     * (0.5) than at the ends. The const k = 0.3 controls how much the
     * distribution "bends." Increasing k just a little (->0.5) makes it almost
     * a straight line, or a uniform distribution, but decreasing it below 0.2
     * makes it bend so much that almost all values are in the middle (0.5).
     * k=0.3 seems just right.
	 *
	 * n is the length, which must be an odd number!
	 *
     * Connectivity distribution (q = false) differs from quality (q = true). It
     * is basically only half a distribution, so that most values are ~1, while
     * there are increasingly fewer lower values. The distribution in this case
     * uses exp(-x/k), rather than exp(-(x-0.5)/k^2).
	 *
	 * The effect of k can be examined with this R function:
		getp <- function(k = 0.3, n = 25)
		{
			x <- (1:n) / (n + 1)
			y <- exp(-(x - 0.5)^2 / k^2)
			n2 <- floor(n / 2)
			p <- rep(NA,n)
			p[n2 + 1] <- 0.5
			for (i in 1:n2) {
				p[n2 + 1 - i] = p[n2 +2 - i] - y[i]
				p[n2 + 1 + i] = p[n2 + i] + y[i]
			}
			p <- 0.05 + 0.9 * (p - min(p)) / (max(p) - min(p))
			plot(x,p,"l",col="red")
			points(x,p,pch=20,col="red")
			lines(c(0.5,0.5),c(0,1),col="lawngreen")
			lines(c(0,1),c(0.5,0.5),col="lawngreen")
        }
        x11()
        par(mfrow=c(2,2))
        getp(0.01); title(main="0.01")
        getp(0.1); title(main="0.1")
        getp(0.2); title(main="0.2")
        getp(0.3); title(main="0.3")
	 */
	const double k = 0.2;

	int nhalf, n = (*pvals).size();
	double tempd[2];

	boost::numeric::ublas::vector <double> refdist (n);

	if (q)
    {
		for (int i=0; i<n; i++)
        {
			refdist (i) = ((double) i + 1.0) / ((double) n + 1.0) - 0.5;
			refdist (i) = exp(-(refdist (i) * refdist (i)) / (k * k));
		}

		nhalf = floor (n / 2);
		(*pvals) (nhalf) = 0.5;
		for (int i=0; i<nhalf; i++)
        {
			(*pvals) (nhalf - i - 1) = (*pvals) (nhalf - i) - refdist (i);
			(*pvals) (nhalf + i + 1) = (*pvals) (nhalf + i) + refdist (i);
		}
    } else {
		for (int i=0; i<n; i++)
        {
			refdist (i) = ((double) i + 1.0) / ((double) n + 1.0);
			refdist (i) = exp(-(refdist (i) * refdist (i)) / k);
		}
		(*pvals) (0) = 1.0;
		for (int i=1; i<n; i++)
			(*pvals) (i) = (*pvals) (i - 1) - refdist (n - i);
	} // end else connectivity distribution

	// Normalise pvals to [0.05,0.95]
	tempd [0] = (*pvals) (0);
	tempd [1] = (*pvals) (n - 1);
	for (int i=0; i<n; i++) 
		(*pvals) (i) = minqc + (1.0 - 2.0 * minqc) * 
			(((*pvals) (i) - tempd [0]) / (tempd [1] - tempd [0]));
} // end function pseq


/************************************************************************
 ************************************************************************
 **                                                                    **
 **                         SORT FUNCTION                              **
 **                                                                    **
 ************************************************************************
 ************************************************************************/

// Just a lazy bubble sort here, coz they're only 3 or 4 long
ivec sort(ivec sortvec, int veclen)
{
	int tempi;

	for(int i=0; i<veclen; i++) 
    {
		for (int j=(i + 1); j<veclen; j++)	
        {
			if (sortvec(j) < sortvec(i)) 
            {
				tempi = sortvec(i);
				sortvec(i) = sortvec(j);
				sortvec(j) = tempi;
			}
		} // end for j
	} // end for i
	return sortvec;
} // end function sort

/************************************************************************
 ************************************************************************
 **                                                                    **
 **                          TIMEOUT FUNCTION                          **
 **                                                                    **
 ************************************************************************
 ************************************************************************/

void timeout(double tseconds)
{
	int hh = floor(tseconds / 3600.0);
	if (hh == 0) 
        std::cout << "00:";
	else if (hh < 10) 
        std::cout << "0" << hh << ":";
	else 
        std::cout << hh << ":";
	double trem = tseconds - (double) hh * 3600.0;
	int mm = floor(trem / 60.0);
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


/************************************************************************
 ************************************************************************
 **                                                                    **
 **                          DUMPSEED FUNCTION                         **
 **                                                                    **
 ************************************************************************
 ************************************************************************/

void dumpseed()
{
	// from: http://sci.tuomastonteri.fi/programming/cplus/wrapper-to-boost-random
	uint64_t tseed;
	std::ifstream urandom;
	urandom.open ("/dev/urandom");
	urandom.read (reinterpret_cast<char*> (&tseed), sizeof (tseed));
	urandom.close ();
	std::cout << "seed = " << tseed << std::endl;
}
