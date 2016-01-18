#include <Rcpp.h>
#include <fstream>

using namespace Rcpp;

// [[Rcpp::export]]
NumericVector iteratePop (NumericMatrix pmat, 
             int nt=1e5, double ksd=0.1, double rr=0.1)
{
    const int nnodes = 4;
    bool flag;
    int count = 0;
    const double minqc = 0.001;
    double tempd, nold [nnodes], n [nnodes], k [nnodes], nmn = 0.0, nsd = 0.0;

    for (int i=0; i<nnodes; i++)
        nold [i] = 1.0;

    //std::ofstream out_file ("aaajunk.txt", std::ios::out);

    for (int i=0; i<nt; i++)
    {
        // Movement through network; transforms nold -> n
        for (int j=0; j<nnodes; j++) 
        {
            n [j] = 0.0;
            for (int k=0; k<nnodes; k++) 
                n [j] += pmat (k, j) * nold [k];
        }
        // Change carrying capacities
        for (int j=0; j<nnodes; j++)
        {
            k [j] = 1.0 + ksd * R::rnorm (0, 1);
            while (k [j] < minqc | k [j] > (2.0 - minqc))
                k [j] = 1.0 + ksd * R::rnorm (0, 1);
        }
        // Population dynamic; transforms n -> nold
        tempd = 0.0;
        flag = false;
        for (int j=0; j<nnodes; j++)
        {
            nold [j] = n [j] + rr * n [j] * n [j] / k [j] -
                rr * n [j] * n [j] * n [j] / (k [j] * k [j]);
            tempd += nold [j];
            if (nold [j] < 0.0)
                flag = true;
        }
        if (!flag)
        {
            count++;
            nmn += tempd;
            nsd += tempd * tempd;
        }
        // Then set all negative nodal abundances to zero
        for (int j=0; j<nnodes; j++)
            if (nold [j] < 0.0)
                nold [j] = 0.0;
    } // end for i
    nmn = nmn / (double) count;
    nsd = nsd / (double) count - nmn * nmn;
    //out_file.close ();

    NumericVector out (2);
    out [0] = nmn;
    out [1] = nsd;

    return out;
}
