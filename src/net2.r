# NOTE that arguments passed to adaptIntegrate have to be exactly those ones
# used. Passing additional arguments that are not used in a particular function
# fails to generate the correct values! 
#
# Note that only O(2) versions include the movement-to-growth scaling parameter
# of s. Note further that "switch-test.r" contains various speed tests for
# different ways of handling different cases of s. This reveals that there is no
# speed advantage other than through constructing entirely separate routines for
# the case of s=0, and there are thus duplicates of calcmean.o2 and calcmean2.o2
# for this case, along with get1val.

calcmean.o2 <- function (r=0.1, rs=1, dn=0, dk=0, dalpha=0, alpha=0, tol=1e-5)
{
    # rs = r^s
    require (cubature)
    # Uses adaptIntegrate. Comparison with speed and accuracy of other potential
    # routines is detailed in "numeric_integration.r".

    fn1 <- function (arg)
    {
        ni <- arg [1]
        k <- arg [2]
        ni * k / (ni + k)
    }
    fn2 <- function (arg, r=0.1, rs=0)
    {
        ni <- arg [1]
        nj <- arg [2]
        k <- arg [3]
        alpha <- arg [4]
        pterm <- (ni + (2*alpha*(1-rs) + alpha^2*(1-rs)^2*nj/ni + 
                        rs^2*alpha^2) * nj) / (ni + alpha*nj)^2
        pterm * (ni * k / (ni + k)) ^ 2
    }
    lims.lo <- c (1-dn/2, 1-dn/2, -dk/2, alpha-dalpha/2)
    lims.hi <- c (1+dn/2, 1+dn/2, dk/2, alpha+dalpha/2)
    norms <- c (dn * dk, dn ^ 2 * dk * dalpha)
    # These indices are essential to get adaptIntegrate to give correct results!
    indxs <- list (c (1, 3), 1:4)
    extras <- list (NULL, list(r=r, rs=rs))
    results <- list (y=list(), maxerror=0, niters=0, calctime=0)
    f <- c ("fn1", "fn2")
    for (i in 1:length (f))
    {
        # Note that proc.time is the appropriate way to time this, because
        # the recommended alternative of using system.time just calls this function
        # exactly as done here anyway.
        pt0 <- proc.time ()
        if (is.null (extras [[i]]))
            y <- do.call (adaptIntegrate, list (f=f [i], 
                                        lowerLimit=lims.lo [indxs [[i]] ],
                                        upperLimit=lims.hi [indxs [[i]] ],
                                        tol=tol))
        else
            y <- do.call (adaptIntegrate, list (f=f [i], 
                                        lowerLimit=lims.lo [indxs [[i]] ],
                                        upperLimit=lims.hi [indxs [[i]] ], 
                                        r=extras[[i]]$r, rs=extras[[i]]$rs,
                                        tol=tol))
        results$calctime <- results$calctime + as.numeric ((proc.time () - pt0) [3])
        if (y$error > results$maxerror) results$maxerror <- y$error
        results$niters <- results$niters + y$functionEvaluations
        if (norms [i] > 0) results$y [[i]] <- y$integral / norms [i]
        else results$y [[i]] <- y$integral # will be zero regardless
    }
    results$value <- 2 + r*(1-r) * results$y [[1]] - r^3 * results$y [[2]]
    return (results)
} # end calcmean

calcmean.o2.s0 <- function (r=0.1, dn=0, dk=0, dalpha=0, alpha=0, tol=1e-5)
{
    require (cubature)
    # Uses adaptIntegrate. Comparison with speed and accuracy of other potential
    # routines is detailed in "numeric_integration.r".

    fn1 <- function (arg)
    {
        ni <- arg [1]
        k <- arg [2]
        ni * k / (ni + k)
    }
    fn2 <- function (arg)
    {
        ni <- arg [1]
        nj <- arg [2]
        k <- arg [3]
        alpha <- arg [4]
        pterm <- (ni + alpha^2 * nj) / (ni + alpha*nj) ^ 2
        pterm * (ni * k / (ni + k)) ^ 2
    }
    lims.lo <- c (1-dn/2, 1-dn/2, -dk/2, alpha-dalpha/2)
    lims.hi <- c (1+dn/2, 1+dn/2, dk/2, alpha+dalpha/2)
    norms <- c (dn * dk, dn ^ 2 * dk * dalpha)
    # These indices are essential to get adaptIntegrate to give correct results!
    indxs <- list (c (1, 3), 1:4)
    results <- list (y=list(), maxerror=0, niters=0, calctime=0)
    f <- c ("fn1", "fn2")
    for (i in 1:length (f))
    {
        # Note that proc.time is the appropriate way to time this, because
        # the recommended alternative of using system.time just calls this function
        # exactly as done here anyway.
        pt0 <- proc.time ()
        y <- do.call (adaptIntegrate, list (f=f [i], 
                                    lowerLimit=lims.lo [indxs [[i]] ],
                                    upperLimit=lims.hi [indxs [[i]] ],
                                    tol=tol))
        results$calctime <- results$calctime + as.numeric ((proc.time () - pt0) [3])
        if (y$error > results$maxerror) results$maxerror <- y$error
        results$niters <- results$niters + y$functionEvaluations
        if (norms [i] > 0) results$y [[i]] <- y$integral / norms [i]
        else results$y [[i]] <- y$integral # will be zero regardless
    }
    results$value <- 2 + r*(1-r) * results$y [[1]] - r^3 * results$y [[2]]
    return (results)
} # end calcmean.s0

calcmean.o3 <- function (r=0.1, dn=0, dk=0, dalpha=0, alpha=0, tol=1e-5)
{
    # This is currently WRONG because it doesn't have the necessary indxs
    require (cubature)

    fn1 <- function (arg)
    {
        ni <- arg [1]
        k <- arg [3]
        ni ^ 2 * k / (ni + k) ^ 2
    }
    fn2 <- function (arg)
    {
        ni <- arg [1]
        nj <- arg [2]
        k <- arg [3]
        alpha <- arg [4]
        ni^4 * k^2 * (ni + alpha^2 * nj) / ((ni + k)^4 * (ni + alpha*nj)^2)
    }
    fn3 <- function (arg)
    {
        ni <- arg [1]
        nj <- arg [2]
        k <- arg [3]
        alpha <- arg [4]
        ni^6 * k^3 * (ni + alpha^3*nj) / ((ni + k)^6 * (ni + alpha*nj)^3)
    }
    lims.lo <- c (1-dn/2, 1-dn/2, -dk/2, alpha-dalpha/2)
    lims.hi <- c (1+dn/2, 1+dn/2, dk/2, alpha+dalpha/2)
    norms <- c (dn * dk, dn ^ 2 * dk * dalpha, dn ^ 2 * dk * dalpha)
    results <- list (y=list(), maxerror=0, niters=0, calctime=0)
    f <- c ("fn1", "fn2", "fn3")
    for (i in 1:3)
    {
        pt0 <- proc.time ()
        y <- do.call (adaptIntegrate, list (f=f [i], lowerLimit=lims.lo,
                                            upperLimit=lims.hi, tol=tol))
        results$calctime <- results$calctime + as.numeric ((proc.time () - pt0) [3])
        if (y$error > results$maxerror) results$maxerror <- y$error
        results$niters <- results$niters + y$functionEvaluations
        if (norms [i] > 0) results$y [[i]] <- y$integral / norms [i]
        else results$y [[i]] <- y$integral 
    }
    results$value <- 2 + r*(1-r) * results$y [[1]] - 2*r^3 * results$y [[2]] -
        r^4 * results$y [[3]]
    return (results)
} # end calcmean.o3

calcmean2.o2 <- function (r=0.1, rs=1, dn=0, dk=0, dalpha=0, alpha=0, tol=1e-5)
{
    require (cubature)

    fn1 <- function (arg)
    {
        ni <- arg [1]
        nj <- arg [2]
        (ni + nj) ^ 2
    }
    fn2 <- function (arg)
    {
        ni <- arg [1]
        k <- arg [2]
        (ni + 1) * ni * k / (ni + k)
    }
    fn3 <- function (arg, r=0.1, rs=1)
    {
        ni <- arg [1]
        nj <- arg [2]
        k <- arg [3]
        alpha <- arg [4]
        pterm <- (ni + (2*alpha*(1-rs) + alpha^2*(1-rs)^2*nj/ni + 
                        rs^2*alpha^2) * nj) / (ni + alpha*nj)^2
        (ni + nj) * pterm * (ni * k / (ni + k)) ^ 2
    }
    fn4 <- function (arg)
    {
        ni <- arg [1]
        k <- arg [2]
        (ni * k / (ni + k)) ^ 2
    }
    fn5 <- function (arg, r=0.1, rs=1)
    {
        ni <- arg [1]
        nj <- arg [2]
        k <- arg [3]
        alpha <- arg [4]
        pterm <- (ni + (2*alpha*(1-rs) + alpha^2*(1-rs)^2*nj/ni + 
                        rs^2*alpha^2) * nj) / (ni + alpha*nj)^2
        pterm * (ni * k / (ni + k)) ^ 3
    }
    fn6 <- function (arg, r=0.1, rs=1)
    {
        ni <- arg [1]
        nj <- arg [2]
        k <- arg [3]
        alpha <- arg [4]
        pterm <- (ni + (2*alpha*(1-rs) + alpha^2*(1-rs)^2*nj/ni + 
                        rs^2*alpha^2) * nj) / (ni + alpha*nj)^2
        pterm ^ 2 * (ni * k / (ni + k)) ^ 4
    }
    lims.lo <- c (1-dn/2, 1-dn/2, -dk/2, alpha-dalpha/2)
    lims.hi <- c (1+dn/2, 1+dn/2, dk/2, alpha+dalpha/2)
    n1 <- dn ^ 2
    n2 <- dn * dk
    n3 <- dn ^ 2 * dk * dalpha
    norms <- c (n1, n2, n3, n2, n3, n3)
    indxs <- list (1:2, c(1,3), 1:4, c(1,3), 1:4, 1:4)
    extras <- list (NULL, NULL, list(r=r, rs=rs), NULL, 
                    list(r=r, rs=rs), list(r=r, rs=rs))
    results <- list (y=list(), maxerror=0, niters=0, calctime=0)
    f <- c ("fn1", "fn2", "fn3", "fn4", "fn5", "fn6")
    for (i in 1:length (f))
    {
        pt0 <- proc.time ()
        if (is.null (extras [[i]]))
            y <- do.call (adaptIntegrate, list (f=f [i], 
                                        lowerLimit=lims.lo [indxs [[i]] ],
                                        upperLimit=lims.hi [indxs [[i]] ],
                                        tol=tol))
        else
            y <- do.call (adaptIntegrate, list (f=f [i], 
                                        lowerLimit=lims.lo [indxs [[i]] ],
                                        upperLimit=lims.hi [indxs [[i]] ], 
                                        r=extras[[i]]$r, rs=extras[[i]]$rs,
                                        tol=tol))
        results$calctime <- results$calctime + as.numeric ((proc.time () - pt0) [3])
        if (y$error > results$maxerror) results$maxerror <- y$error
        results$niters <- results$niters + y$functionEvaluations
        if (norms [i] > 0) results$y [[i]] <- y$integral / norms [i]
        else results$y [[i]] <- y$integral
    }
    results$value <- results$y[[1]] + 2*r*(1-r) * results$y[[2]] -
            2*r^3 * results$y[[3]] + r^2*(1-r)^2 * results$y[[4]] -
            2*r^4*(1-r) * results$y[[5]] + r^6 * results$y[[6]]
    return (results)
} # end calcmean2.o2

calcmean2.o2.s0 <- function (r=0.1, dn=0, dk=0, dalpha=0, alpha=0, tol=1e-5)
{
    require (cubature)

    fn1 <- function (arg)
    {
        ni <- arg [1]
        nj <- arg [2]
        (ni + nj) ^ 2
    }
    fn2 <- function (arg)
    {
        ni <- arg [1]
        k <- arg [2]
        (ni + 1) * ni * k / (ni + k)
    }
    fn3 <- function (arg)
    {
        ni <- arg [1]
        nj <- arg [2]
        k <- arg [3]
        alpha <- arg [4]
        pterm <- (ni + alpha^2 * nj) / (ni + alpha*nj) ^ 2
        (ni + nj) * pterm * (ni * k / (ni + k)) ^ 2
    }
    fn4 <- function (arg)
    {
        ni <- arg [1]
        k <- arg [2]
        (ni * k / (ni + k)) ^ 2
    }
    fn5 <- function (arg)
    {
        ni <- arg [1]
        nj <- arg [2]
        k <- arg [3]
        alpha <- arg [4]
        pterm <- (ni + alpha^2 * nj) / (ni + alpha*nj) ^ 2
        pterm * (ni * k / (ni + k)) ^ 3
    }
    fn6 <- function (arg)
    {
        ni <- arg [1]
        nj <- arg [2]
        k <- arg [3]
        alpha <- arg [4]
        pterm <- (ni + alpha^2 * nj) / (ni + alpha*nj) ^ 2
        pterm ^ 2 * (ni * k / (ni + k)) ^ 4
    }
    lims.lo <- c (1-dn/2, 1-dn/2, -dk/2, alpha-dalpha/2)
    lims.hi <- c (1+dn/2, 1+dn/2, dk/2, alpha+dalpha/2)
    n1 <- dn ^ 2
    n2 <- dn * dk
    n3 <- dn ^ 2 * dk * dalpha
    norms <- c (n1, n2, n3, n2, n3, n3)
    indxs <- list (1:2, c(1,3), 1:4, c(1,3), 1:4, 1:4)
    results <- list (y=list(), maxerror=0, niters=0, calctime=0)
    f <- c ("fn1", "fn2", "fn3", "fn4", "fn5", "fn6")
    for (i in 1:length (f))
    {
        pt0 <- proc.time ()
        y <- do.call (adaptIntegrate, list (f=f [i], 
                                    lowerLimit=lims.lo [indxs [[i]] ],
                                    upperLimit=lims.hi [indxs [[i]] ],
                                    tol=tol))
        results$calctime <- results$calctime + as.numeric ((proc.time () - pt0) [3])
        if (y$error > results$maxerror) results$maxerror <- y$error
        results$niters <- results$niters + y$functionEvaluations
        if (norms [i] > 0) results$y [[i]] <- y$integral / norms [i]
        else results$y [[i]] <- y$integral
    }
    results$value <- results$y[[1]] + 2*r*(1-r) * results$y[[2]] -
            2*r^3 * results$y[[3]] + r^2*(1-r)^2 * results$y[[4]] -
            2*r^4*(1-r) * results$y[[5]] + r^6 * results$y[[6]]
    return (results)
} # end calcmean2.o2.s0

calcmean2.o3 <- function (r=0.1, dn=0, dk=0, dalpha=0, alpha=0, tol=1e-5)
{
    require (cubature)

    fn1 <- function (arg)
    {
        ni <- arg [1]
        k <- arg [3]
        (ni + 1) * ni ^ 2 * k / (ni + k) ^ 2
    }
    fn2 <- function (arg)
    {
        ni <- arg [1]
        nj <- arg [2]
        k <- arg [3]
        alpha <- arg [4]
        (ni + nj) * ni^4 * k^2 * (ni + alpha ^ 2 * nj) / 
            ((ni + k)^4 * (ni + alpha * nj) ^ 2)
    }
    fn3 <- function (arg)
    {
        ni <- arg [1]
        nj <- arg [2]
        k <- arg [3]
        alpha <- arg [4]
        (ni + nj) * ni^6 * k^3 * (ni + alpha^3*nj) / ((ni+k)^6 * (ni+alpha*nj)^3)
    }
    fn4 <- function (arg)
    {
        ni <- arg [1]
        k <- arg [3]
        ni^4 * k^2 / (ni + k) ^ 4
    }
    fn5 <- function (arg)
    {
        ni <- arg [1]
        nj <- arg [2]
        k <- arg [3]
        alpha <- arg [4]
        ni^6 * k^3 * (ni+alpha^2*nj) / ((ni+k)^6 * (ni+alpha*nj)^2)
    }
    fn6 <- function (arg)
    {
        ni <- arg [1]
        nj <- arg [2]
        k <- arg [3]
        alpha <- arg [4]
        ni^8 * k^4 * (ni + alpha^3*nj) / ((ni+k)^8 * (ni+alpha*nj)^3)
    }
    fn7 <- function (arg)
    {
        ni <- arg [1]
        nj <- arg [2]
        k <- arg [3]
        alpha <- arg [4]
        ni^8 * k^4 * (ni+alpha^2*nj)^2 / ((ni+k)^8 * (ni + alpha*nj)^4)
    }
    fn8 <- function (arg)
    {
        ni <- arg [1]
        nj <- arg [2]
        k <- arg [3]
        alpha <- arg [4]
        ni^(10) * k^5 * (ni+alpha^2*nj) * (ni+alpha^3*nj) /
            ((ni + k) ^ (10) * (ni + alpha * nj) ^ 5)
    }
    fn9 <- function (arg)
    {
        ni <- arg [1]
        nj <- arg [2]
        k <- arg [3]
        alpha <- arg [4]
        ni^(12) * k^6 * (ni + alpha^3*nj)^2 / ((ni+k)^(12) * (ni+alpha*nj)^6)
    }
    lims.lo <- c (1-dn/2, 1-dn/2, -dk/2, alpha-dalpha/2)
    lims.hi <- c (1+dn/2, 1+dn/2, dk/2, alpha+dalpha/2)
    n1 <- dn * dk
    n2 <- dn ^ 2 * dk * dalpha
    norms <- c (n1, n2, n2, n1, n2, n2, n2, n2, n2)
    results <- list (y=list(), maxerror=0, niters=0, calctime=0)
    f <- c ("fn1", "fn2", "fn3", "fn4", "fn5", "fn6", "fn7", "fn8", "fn9")
    for (i in 1:length (f)) {
        pt0 <- proc.time ()
        y <- do.call (adaptIntegrate, list (f=f [i], lowerLimit=lims.lo,
                                            upperLimit=lims.hi, tol=tol))
        results$calctime <- results$calctime + as.numeric ((proc.time () - pt0) [3])
        if (y$error > results$maxerror) results$maxerror <- y$error
        results$niters <- results$niters + y$functionEvaluations
        if (norms [i] > 0) results$y [[i]] <- y$integral / norms [i]
        else results$y [[i]] <- y$integral 
    }
    results$value <- 4 - 2*r*(1-r)*results$y [[1]] - 4*r^3*results$y [[2]] -
        2*r^4*results$y [[3]] + r^2*(1-r)^2*results$y [[4]] - 
        4*r^4*(1-r)*results$y [[5]] - 2*r^5*(1-r)*results$y [[6]] + 
        4*r^6*results$y [[7]] + 4*r^7*results$y [[8]] +
        r^8*results$y [[9]]
    return (results)
}


get1val <- function (r=0.1, rs=1, dn=0, dk=0, dalpha=0, alpha=0, tol=1e-5, ord=2)
{
    if (ord == 2) {
        d1 <- calcmean.o2 (r, rs, dn, dk, dalpha, alpha, tol=tol)
        d2 <- calcmean2.o2 (r, rs, dn, dk, dalpha, alpha, tol=tol)
    } else {
        d1 <- calcmean.o3 (r, dn, dk, dalpha, alpha, tol=tol)
        d2 <- calcmean2.o3 (r, dn, dk, dalpha, alpha, tol=tol)
    }
    sigma <- d2$value - d1$value ^ 2
    niters <- d1$niters + d2$niters
    calctime <- d1$calctime + d2$calctime
    error <- max (c (d1$maxerror, d2$maxerror))
    results <- list (pars=list (dn=dn, dk=dk, dalpha=dalpha, alpha=alpha,
                                rs=rs),
                     error=error, niters=niters, calctime=calctime, sigma=sigma)
    return (results)
} # end calcmean2.o3

get1val.s0 <- function (r=0.1, dn=0, dk=0, dalpha=0, alpha=0, tol=1e-5, ord=2)
{
    if (ord == 2) {
        d1 <- calcmean.o2.s0 (r, dn, dk, dalpha, alpha, tol=tol)
        d2 <- calcmean2.o2.s0 (r, dn, dk, dalpha, alpha, tol=tol)
    } else {
        d1 <- calcmean.o3 (r, dn, dk, dalpha, alpha, s=0, tol=tol)
        d2 <- calcmean2.o3 (r, dn, dk, dalpha, alpha, s=0, tol=tol)
    }
    sigma <- d2$value - d1$value ^ 2
    niters <- d1$niters + d2$niters
    calctime <- d1$calctime + d2$calctime
    error <- max (c (d1$maxerror, d2$maxerror))
    results <- list (pars=list (dn=dn, dk=dk, dalpha=dalpha, alpha=alpha, s=0),
                     error=error, niters=niters, calctime=calctime, sigma=sigma)
    return (results)
} # end calcmean2.o3

get1slope <- function (r=0.1, rs=1, dn=0, dk=0, dalpha=0, tol=1e-5, ord=2, 
                       regr.len=50, plot=FALSE) {
    alo <- dalpha / 2
    ahi <- 1 - dalpha / 2
    alpha <- seq (alo, ahi, length.out=regr.len)
    slope <- slope.rel <- calctime <- sigma <- NA

    pt0 <- proc.time ()
    if (rs == 1)
        sigma <- sapply (alpha, function (a) { get1val.s0 (r=r, dn=dn, dk=dk, 
                        dalpha=dalpha, alpha=a, tol=tol, ord=ord)$sigma })
    else
        sigma <- sapply (alpha, function (a) { get1val (r=r, rs=rs, dn=dn, dk=dk, 
                        dalpha=dalpha, alpha=a, tol=tol, ord=ord)$sigma })
    #calctime <- format (.POSIXct ((proc.time () - pt0) [3], tz="GMT"), "%H:%M:%S")
    if (plot) {
        plot (alpha, sigma, "l", col="blue")
        lines (alpha, lm(sigma ~ alpha)$fitted.values, col="red")
    }
    calctime <- as.numeric ((proc.time () - pt0) [3])
    slope <- as.numeric (lm (sigma ~ alpha)$coefficients [2])
    if (!all (is.na (sigma))) {
        if (mean (sigma, na.rm=T) > 0) {
            slope.rel <- slope / mean (sigma, na.rm=T)
        } else slope.rel <- 0
    }
    return (list (slope=slope, slope.rel=slope.rel, calctime=calctime))
}

getslopes <- function (r=0.1, s=0, tol=1e-3, ord=2)
{
    dn0 <- (0:99) / 100
    dk0 <- (0:100) / 100
    kmax <- 1 / (1 + r)
    indx <- which (dk0 < kmax)
    dk0 <- dk0 [indx]
    regr.len <- 50

    rs <- r ^ s

    # The following use mapply to do neat evaluation, but there is no way of
    # including a progress indicator, so explicit loops are used. (There is
    # "pbapply" but this only works for apply, sapply, and lapply.)
    #dk <- rep (dk0, length (dn0))
    #dn <- unlist (lapply (dn0, function (x) rep (x, length (dk0))))
    #dat.in <- data.frame (cbind (dk, dn))
    #names (dat.in) <- c ("dk", "dn")
    #dat.out <- mapply (function (dk, dn) 
    #                   get1slope (r=r, dn=dn, dk=dk, dalpha=dn, tol=tol), 
    #                   dat.in$dk, dat.in$dalpha)
    #
    # Note that the maximal value of dk for a cubic logistic is
    # dk = 1 + r/2 - sqrt (r^2 + 4*r)

    arr <- array (NA, dim=c(length (dn0), length (dk0)))
    results <- list (r=r, s=s, dn=dn0, dk=dk0, slope=arr, slope.rel=arr,
                     regr.len=regr.len)
    fname <- paste ("slopes_o", ord, "_r0", r * 10, "_s", s, ".rdat", sep="")
    ptm <- proc.time ()
    pb <- txtProgressBar (max=1, style = 3) # shows start and end positions, and %
    for (i in 1:length (dn0)) {
        slopes <- lapply (dk0, function (k) {
                          get1slope (r=r, rs=rs, dn=dn0 [i], dk=k,
                                     dalpha=dn0 [i], tol=tol, ord=ord,
                                     regr.len=regr.len)
                     })
        results$slope [i, ] <- sapply (slopes, function (s) s$slope)
        results$slope.rel [i, ] <- sapply (slopes, function (s) s$slope.rel)
        setTxtProgressBar(pb, i / length (dn0))
        save (results, file=fname)
    }
    close(pb)
    results$calctime <- format (.POSIXct ((proc.time () - ptm) [3], tz="GMT"), "%H:%M:%S")
    cat ("Done; calculation time = ", results$calctime, "\n", sep="")
    save (results, file=fname)
}

getslopes.1d <- function (r=0.1, s=0, dn=0.5, tol=1e-3, ord=2)
{
    dk0 <- (0:100) / 100
    kmax <- 1 / (1 + r)
    indx <- which (dk0 < kmax)
    dk0 <- dk0 [indx]
    regr.len <- 50

    rs <- r ^ s

    vec <- rep (NA, length (dk0))
    results <- list (r=r, s=s, dn=dn, dk=dk0, slope=vec, slope.rel=vec,
                     regr.len=regr.len, tol=tol)
    fname <- paste ("slopes_o", ord, "_r0", r * 10, "_s", s, "_1D.rdat", sep="")
    ptm <- proc.time ()
    pb <- txtProgressBar (max=1, style = 3) # shows start and end positions, and %
    for (i in 1:length (dk0)) {
        slopes <- get1slope (r=r, rs=rs, dn=dn, dk=dk0 [i], 
                             dalpha=dn, tol=tol, ord=ord, regr.len=regr.len)
        results$slope [i] <- slopes$slope
        results$slope.rel [i] <- slopes$slope.rel
        prog <- i / length (dk0)
        setTxtProgressBar(pb, prog)
        save (results, file=fname)
    }
    close(pb)
    results$calctime <- format (.POSIXct ((proc.time () - ptm) [3], tz="GMT"), "%H:%M:%S")
    cat ("Done; calculation time = ", results$calctime, "\n", sep="")
    save (results, file=fname)
}

plotslopes <- function (r=0.1, rel.slopes=FALSE)
{
    fname <- paste ("./results/slopes_o3_r0", r * 10, ".rdat", sep="")
    load (fname)
    if (rel.slopes) { slopes <- results$slope.rel    }
    else { slopes <- results$slope   }
    
    ym = 0.5 * (1 - results$dn) * (2 + r - sqrt((2 + r) ^ 2 - 4))
    ymi <- sapply (ym, function (i) which.min (abs (i - results$dk)))
    ds <- dim (slopes)
    for (i in 1:ds [1]) slopes [i, ymi [i]:ds [2]] <- NA

    x11 ()
    layout (matrix(1:2,1,2),widths=c(0.9,0.1))
    plot (NULL, NULL, xlim=range(results$dn), ylim=range(results$dk),
        xlab="Predictable variation (n)", ylab="Unpredictable variation (k)")
    title (main=paste ("r = ", r, ", max slope = ", 
                       max (slopes, na.rm=TRUE), sep=""))
    ncols <- 30
    levs <- pretty (slopes, n=ncols)
    cols <- terrain.colors (ncols)
    .filled.contour(as.double (results$dn), as.double (results$dk), as.matrix (slopes),
            as.double (levs), cols)
    contour (results$dn, results$dk, slopes, levels=c(-0.1,0,0.1), 
             col=c("blue", "red", "orange"), add=TRUE)
    # The maximal value of dk for a cubic logistic is
    # dk = 0.5 * n * (2 + r - sqrt((2 + r) ^ 2 - 4))
    ym = 0.5 * (1 - results$dn) * (2 + r - sqrt((2 + r) ^ 2 - 4))
    yp = 0.5 * (1 + results$dn) * (2 + r - sqrt((2 + r) ^ 2 - 4))
    lines (results$dn, ym, col="red", lty=2)
    lines (results$dn, yp, col="red")


    # legend
    xps <- c (0, 0.3, 0.2)
    zlims <- range (slopes, na.rm=TRUE)
    par (mar=c(2,0,2,0),mgp=c(0,0,0), ps=8)
    plot (NULL, NULL, xlim=c(0,1), ylim=zlims,
          xlab="", ylab="", xaxt="n", yaxt="n", frame=FALSE)
    tc <- rev (cols)
    coly <- abs(diff(zlims)) / ncols # Height of each colour band             
    for (i in 1:ncols) {                                                      
        tcol <- ncols - i + 1                                             
        bottom <- zlims[1] + (i - 1) * coly                               
        top <- bottom + coly                                              
        polygon(c(xps[1],xps[1],xps[2],xps[2]),c(bottom,top,top,bottom),  
                    col=tc[tcol],border=NA)
    }
    lines (c (xps[1], xps[2]), c (zlims[2], zlims[2]), col="black")                  
    lines (c (xps[1], xps[1]), c (zlims[1], zlims[2]), col="black")                  
    lines (c (xps[1], xps[2]), c (zlims[1], zlims[1]), col="black")                  
    lines (c (xps[2], xps[2]), c (zlims[1], zlims[2]), col="black") 
    levs <- pretty (zlims, 5)
    for (i in 1:length (levs)) {
        lines (c (xps[1], xps [2]), rep (levs [i], 2))
        text (xps [3], levs [i], levs [i], pos=4)
    }
}


# *****************************************************************
# This function takes about 6 minutes to execute, and plots the effect of different
# adaptIntegrate tolerances on calculation times and errors. Decreasing tolerances
# cause roughly linear (or maybe quadratic) increases in calculation times, while
# errors decrease roughly logarithmically (can be seen clearly in the final of the 3
# plots).
#
# A tolerance of 1e-5 has errors around 4-5 times as much as with tol=1e-6, and is
# around 10 times faster. Thus tol=1e5 could be used for development, with final
# runs done with tol=1e-6.

gettimes <- function (n=100) {
    tols <- 10 ^ -(3:6)
    tols <- 10 ^ (-(30:60) / 10)
    tmn <- tsd <- emn <- esd <- rep (NA, length (tols))
    pt0 <- proc.time ()
    for (i in 1:length (tols)) {
        junk <- sapply (1:n, function (j) {
                        val <- get1val (tol=tols [i])
                        c (val$calctime, val$error) })
        cat ("\r", i, "/", length (tols), sep="")
        tmn [i] <- mean (junk [1,])
        tsd [i] <- sd (junk [1,])
        emn [i] <- mean (junk [2,])
        esd [i] <- sd (junk [2,])
    }
    calctime <- format (.POSIXct ((proc.time () - pt0) [3], tz="GMT"), "%H:%M:%S")
    cat (". Calculation time = ", calctime, "\n", sep="")

    tlo <- tmn - tsd
    thi <- tmn + tsd
    ylims <- range (c (tlo, thi))
    x11 (width=12)
    par (mfrow=c(1,3))
    plot (tols, tmn, "l", col="red", lwd=2, ylim=ylims, log="x",
          main="Calculation Time")
    points (tols, tmn, pch=19, col="red")
    lines (tols, tlo, col="red", lwd=2, lty=2)
    lines (tols, thi, col="red", lwd=2, lty=2)

    elo <- emn - esd
    ehi <- emn + esd
    ylims <- range (c (elo, ehi))
    plot (tols, emn, "l", col="blue", lwd=2, ylim=ylims, log="x", main="Error")
    points (tols, emn, pch=19, col="blue")
    lines (tols, elo, col="blue", lwd=2, lty=2)
    lines (tols, ehi, col="blue", lwd=2, lty=2)

    plot (tmn, emn, pch=19, col="lawngreen", log="y",
          xlab="Calculation Time", ylab="Error")
    i5 <- which (tols == 1e-5)
    i6 <- which (tols == 1e-6)
    tratio <- tmn [i6] / tmn [i5]
    eratio <- emn [i6] / emn [i5]
    cat ("tol [1e-5 -> 1e-6] increases calculation time by ", tratio,
         " and decreases error by ", eratio, "\n", sep="")
}
