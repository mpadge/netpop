# The boost random number generator doesn't seem as good as the R one, as can be
# seen by simply comparing rnorm histograms---R is perfectly smooth (for
# sufficiently large N>1e6), while boost is still notably jittery. These test
# files implement the four-node network using Rcpp which allows the C++ routine
# to call the R rnorm routine.

require (Rcpp)
#pmat <- make.pmat ()
#sourceCpp ("net4.cpp")
#iteratePop (pmat)

do1trial <- function (nt=1e5, ksd=0.1, rr=0.1)
{
    pb <- txtProgressBar (max=1, style = 3) # shows start and end positions
    nmn <- nsd <- array (NA, dim=c(101, 4))
    alpha <- 0:100 / 100
    for (a in 1:length (alpha))
    {
        for (net in 1:4)
        {
            pmat <- make.pmat (net.type = net - 1, alpha = alpha [a])
            results <- iterate.pop (pmat, nt=nt, ksd=ksd, rr=rr)
            nmn [a, net] <- results [1]
            nsd [a, net] <- results [2]
        }
        setTxtProgressBar(pb, a / 101)
    }
    close (pb)

    ylims <- range (nsd)
    plot (NULL, NULL, xlim=c(0,1), ylim=ylims, xlab="", ylab="",
          xaxt="n", yaxt="n")
    cols <- rainbow (4)
    for (i in 1:4)
        lines (alpha, nsd [,i], col=cols [i])
    return (list (nmn, nsd))
}

iterate.pop <- function (pmat, nt=1e4, ksd=0.1, rr=0.1)
{
    minqc <- 0.001

    nold <- rep (1, 4)
    nmn.net <- nsd.net <- 0

    for (i in 1:nt)
    {
        # Movement through network; transforms nold -> n
        n <- colSums (array (nold, dim=c(4, 4)) * pmat)
        # Change carrying capacities
        k <- 1 + ksd * rnorm (4)
        while (any (k < minqc) | any (k > (2 - minqc)))
            k <- 1 + ksd * rnorm (4)
        # Population dynamic; transforms n -> nold
        nold <- n + rr * n ^ 2 / k - rr * n ^ 3 / (k ^ 2)
        #if (any (nold < 0))
        #    cat ("Warning: n[", which (nold < 0), "] < 0\n")
        tempd <- sum (nold)
        nmn.net <- nmn.net + tempd
        nsd.net <- nsd.net + tempd ^ 2
        # Then set all negative nodal abundances to zero
        nold [which (nold < 0)] <- 0
    } # end for i
    nmn.net <- nmn.net / nt
    nsd.net <- nsd.net / nt - nmn.net ^ 2

    return (c (nmn.net, nsd.net))
}

make.pmat <- function (net.type = 0, alpha0 = 0.5)
{
    pmat <- array (NA, dim=c(4, 4))

    if (net.type == 0) 
    {
        #
        #   1---2---3---4
        #
        tempd <- 1 + alpha0 + alpha0 ^ 2 + alpha0 ^ 3
        pmat [1, 1] <- pmat [4, 4] <- 1 / tempd
        pmat [1, 2] <- pmat [4, 3] <- alpha0 * pmat [1,1]
        pmat [1, 3] <- pmat [4, 2] <- alpha0 * pmat [1,2]
        pmat [1, 4] <- pmat [4, 1] <- alpha0 * pmat [1,3]

        tempd <- 1 + 2 * alpha0 + alpha0 ^ 2
        pmat [2, 2] <- pmat [3, 3] <- 1 / tempd;
        pmat [2, 1] <- pmat [3, 4] <- alpha0 * pmat [2, 2]
        pmat [2, 3] <- pmat [3, 2] <- alpha0 * pmat [2, 2]
        pmat [2, 4] <- pmat [3, 1] <- alpha0 * pmat [2, 1]
    } else if (net.type == 1) {
        #
        #	    4
        #	    |
        #	    3
        #	   / \
        #	  1   2
        #
        tempd <- 1 + alpha0 + 2 * alpha0 ^2
        pmat [1, 1] <- pmat [2, 2] <- pmat [4, 4] <- 1.0 / tempd
        pmat [1, 3] <- pmat [2, 3] <- pmat [4, 3] <- alpha0 * pmat [1, 1]
        pmat [1, 2] <- pmat [2, 1] <- alpha0 * pmat [1, 3]
        pmat [1, 4] <- pmat [4, 1] <- pmat [1, 2]
        pmat [2, 4] <- pmat [4, 2] <- pmat [1, 2]

        tempd <- 1 + 3 * alpha0
        pmat [3, 3] <- 1 / tempd
        pmat [3, 1] <- pmat [3, 2] <- pmat [3, 4] <- alpha0 * pmat [3, 3]
    } else if (net.type == 2) {
        #
        # 	3---4
        # 	|   |
        # 	1---2
        #
        tempd <- 1 + 2 * alpha0 + alpha0 ^ 2
        pmat [1, 1] <- pmat [2, 2] <- pmat [3, 3] <- pmat [4, 4] <- 1 / tempd
        pmat [1, 2] <- pmat [1, 3] <- alpha0 * pmat [1, 1]
        pmat [2, 1] <- pmat [2, 4] <- pmat [3, 1] <- pmat [3, 4] <- pmat [1, 2]
        pmat [4, 3] <- pmat [4, 2] <- pmat [1, 2]
        pmat [1, 4] <- pmat [4, 1] <- alpha0 ^ 2 * pmat [1, 1]
        pmat [2, 3] <- pmat [3, 2] <- pmat [1, 4]
    } else if (net.type == 3) {
        #
        # 	     4
        # 	     |
        # 	     3
        # 	    / \
        # 	   1---2
        #
        tempd <- 1 + 2 * alpha0 + alpha0 ^ 2 
        pmat [1, 1] <- pmat [2, 2] <- 1.0 / tempd
        pmat [1, 2] <- pmat [1, 3] <- alpha0 * pmat [1, 1]
        pmat [2, 1] <- pmat [2, 3] <- pmat [1, 2]
        pmat [1, 4] <- pmat [2, 4] <- alpha0 * pmat [1, 2]

        tempd <- 1 + 3 * alpha0
        pmat [3, 3] <- 1 / tempd
        pmat [3, 1] <- pmat [3, 2] <- pmat [3, 4] <- alpha0 * pmat [3, 3]

        tempd <- 1 + alpha0 + 2 * alpha0 ^ 2
        pmat [4, 4] <- 1 / tempd
        pmat [4, 3] <- alpha0 * pmat [4, 4]
        pmat [4, 2] <- pmat [4, 1] <- alpha0 * pmat [4, 3]
    } # end if network type == 3

    return (pmat)
}
