# Compares different R routines for multi-dimensional numerical integration. The
# R2Cuba packages has a number of routines, which are compared with the single
# adaptive integration routine from "cubature." This latter is revealed to be
# superior to any others, as quantified in terms of computation times, numbers of
# iterations, absolute errors, and chi-squared probabilities of results being
# accurate. The latter are unbiased for the R2Cuba routines, unlike absolute errors,
# as explained in the manual.
#
# The only real comparison is between "cuhre" and "adaptIntegrate," with "cuhre"
# producing lower errors (~0.93), yet taking slightly longer to compute (~1.44), and
# requiring slightly more iterations (~1.45). Because speed is likely more important
# than extreme precision, and because the gain in accuracy is not as large as the
# loss in speed, adaptIntegrate seems to be the better routine to use. This
# script simply serves to justify that choice.

getvals <- function (tol=1e-4)
{
    require (cubature)
    require (R2Cuba)
    # The function is the third components of the mean density.
    fn <- function (arg)
    {
        ni <- arg [1]
        nj <- arg [2]
        k <- arg [3]
        alpha <- arg [4]
        (ni + alpha * nj) * (ni * k) ^ 2 / ((ni + alpha * nj) ^ 2 * (ni + k) ^ 2)
    }
    dn <- runif (1)
    dk <- runif (1)
    dalpha <- runif (1)
    alpha <- runif (1)
    while (alpha < dalpha / 2 | (alpha + dalpha / 2) > 1)
    {
        dalpha <- runif (1)
        alpha <- runif (1)
    }
    lims.lo <- c (1-dn/2, 1-dn/2, -dk/2, alpha-dalpha/2)
    lims.hi <- c (1+dn/2, 1+dn/2, dk/2, alpha+dalpha/2)
    vals <- niters <- err <- times <- NULL
    probs <- NA # chi-squared probs, not valid for adaptIntegrate

    args1 <- list (f=fn, lowerLimit=lims.lo, upperLimit=lims.hi, tol=tol)
    args2 <- list (ndim=4, ncomp=1, integrand=fn, lower=lims.lo, upper=lims.hi,
                   rel.tol=tol, flags=list(verbose=0))
    val.names <- c ("integral", rep ("value", 4))
    niter.names <- c ("functionEvaluations", rep ("neval", 4))
    err.names <- c ("error", rep ("abs.error", 4))
    f <- c ("adaptIntegrate", "cuhre", "divonne", "suave", "vegas")
    args <- list (args1, args2, args2, args2, args2)
    for (i in 1:5) {
        pt <- proc.time ()
        y <- do.call (f [i], args [[i]])
        times <- c (times, (proc.time () - pt) [3])
        vals <- c (vals, y [[val.names [i]]])
        niters <- c (niters, y [[niter.names [i]]])
        err <- c (err, y [[err.names [i]]])
        if (i > 1) probs <- c (probs, y$prob)
    }
    return (list (times=times, vals=vals, niters=niters, err=err, probs=probs))
}

getallvals <- function (nrpts=20)
{
    pb <- txtProgressBar (max=1, style = 3) # shows start and end positions
    results <- getvals ()
    for (i in 1:(nrpts - 1)) {
        temp <- getvals ()
        for (j in 1:length (results)) {
            results [[j]] <- rbind (results [[j]], temp [[j]])
        }
        setTxtProgressBar(pb, (i+1) / nrpts)
    }
    close (pb)
    return (lapply (results, colMeans))
}
results <- getallvals ()

plotvals <- function (results)
{
    x11 ()
    par (mfrow=c(2,3))
    x <- c ("aI", "c", "d", "s", "v")
    for (i in 1:length (results)) {
        barplot (results [[i]], width=0.9, space=0.5, bty="l", names.arg=x)
        title (main=names(results) [i])
    }
    plot.new ()
    err.ratio <- results$err [2] / results$err [1]
    text (0.5, 0.8, labels="Error ratio")
    text (0.5, 0.7, labels="for c/aI")
    text (0.5, 0.6, labels=paste ("=", formatC (err.ratio, format="f", digits=2)))
    time.ratio <- results$times [1] / results$times [2]
    text (0.5, 0.4, labels="Time ratio")
    text (0.5, 0.3, labels="for aI/c")
    text (0.5, 0.2, labels=paste ("=", formatC (time.ratio, format="f", digits=2)))
    return (c (err.ratio, time.ratio))
}
ratios <- plotvals (results)
