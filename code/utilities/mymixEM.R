mymixEM <- function (x, lambda = NULL, mu = NULL, sigma = NULL, k = 2, mean.constr = NULL, 
    sd.constr = NULL, epsilon = 1e-08, maxit = 1000, maxrestarts = 20, 
    verb = FALSE, fast = FALSE, ECM = FALSE, arbmean = TRUE, 
    arbvar = TRUE, mean.limits = NULL) 
{
    require(mixtools)
    getMap <- function(xrange, murange) {
        xmin <- min(xrange)
        xmax <- max(xrange)
        mumin <- min(murange)
        mumax <- max(murange)
        f <- function(x) {
            ## writeLines(paste(
            ##     "xmin:",xmin,"\n",
            ##     "xmax:",xmax,"\n",
            ##     "mumin:",mumin,"\n",
            ##     "mumax:",mumax))
            if (xmin < x & x < xmax) {
                (mumax-mumin)/(xmax-xmin)*(x-xmin)+mumin
            } else if (x >= xmax) {
                mumax
            } else {                        #else if (x <= xmin)
                mumin
            }
        }
    }
    warn <- options(warn = -1)
    x <- as.vector(x)
    tmp <- normalmix.init(x = x, lambda = lambda, mu = mu, s = sigma, 
        k = k, arbmean = arbmean, arbvar = arbvar)
    lambda <- tmp$lambda
    mu <- tmp$mu
    sigma <- tmp$s
    k <- tmp$k
    arbvar <- tmp$arbvar
    arbmean <- tmp$arbmean
    ## Verify that the mean.limits structure matches with mu and mean.constr
    if (!is.null(mean.constr) && !is.null(mean.limits)) {
        stop("mean.constr and mean.limits may not both be specified")
    }
    if (!is.null(mean.limits)) {
        xmin <- min(x)
        xmax <- max(x)
        xrange <- c(xmin, xmax)
        mean.limits <- matrix(as.vector(mean.limits), ncol=2, byrow=T)
        ## Replace any NA entries with xmin in the first column and xmax in the second
        mean.limits[,1] <- ifelse(is.na(mean.limits[,1]),
                                  xmin,
                                  mean.limits[,1])
        mean.limits[,2] <- ifelse(is.na(mean.limits[,2]),
                                  xmax,
                                  mean.limits[,2])
        ## Verify that the lower bound is less than or equal to the upper bound
        if (!all(apply(mean.limits, 1, function(row) {row[1] <= row[2]}))) {
            stop("Not all lower bounds are less than the respective upper bounds")
        }
        ## Verify that all mean limit values are within the range of the data
        if (!all(xmin <= mean.limits & mean.limits <= xmax)) {
            stop("Not all specified mean limits are within the data range")
        }
        maps <- apply(mean.limits, 1, function(murange) { #returns a list
            getMap(xrange, murange)
        })
        print(mean.limits)
        ## Generate the mapping to redirect log liklihood calculations by the algorithm
        loglikmap <- function(mu) {
            sapply(seq(mu), function(i) {
                maps[[i]](mu[1])
            })
        }
    } else {
        loglikmap <- function(mu) {mu}
    }
    pseudo.log.lik <- function(n, k, x, mu, sigma, lambda) {
        .C(mixtools:::C_normpost, as.integer(n), as.integer(k), 
           as.double(x), as.double(loglikmap(mu)), as.double(sigma), 
           as.double(lambda), res2 = double(n * k), double(3 * k),
           post = double(n * k), loglik = double(1), 
           PACKAGE = "mixtools")                               
    }
    if (fast == TRUE && k == 2 && arbmean == TRUE) {
        a <- normalmixEM2comp(x, lambda = lambda[1], mu = mu, 
            sigsqrd = sigma^2, eps = epsilon, maxit = maxit, 
            verb = verb)
    }
    else {
        z <- parse.constraints(mean.constr, k = k, allsame = !arbmean)
        meancat <- z$category
        meanalpha <- z$alpha
        z <- parse.constraints(sd.constr, k = k, allsame = !arbvar)
        sdcat <- z$category
        sdalpha <- z$alpha
        ECM <- ECM || any(meancat != 1:k) || any(sdcat != 1)
        n <- length(x)
        notdone <- TRUE
        restarts <- 0
        while (notdone) {
            notdone <- FALSE
            tmp <- normalmix.init(x = x, lambda = lambda, mu = mu, 
                s = sigma, k = k, arbmean = arbmean, arbvar = arbvar)
            lambda <- tmp$lambda
            mu <- tmp$mu
            k <- tmp$k
            sigma <- tmp$s
            var <- sigma^2
            diff <- epsilon + 1
            iter <- 0
            postprobs <- matrix(nrow = n, ncol = k)
            mu <- rep(mu, k)[1:k]
            sigma <- rep(sigma, k)[1:k]
            z <- pseudo.log.lik(n, k, x, mu, sigma, lambda)
            postprobs <- matrix(z$post, nrow = n)
            res <- matrix(z$res2, nrow = n)
            ll <- obsloglik <- z$loglik
            while (diff > epsilon && iter < maxit) {
                lambda <- colMeans(postprobs)
                mu[meancat == 0] <- meanalpha[meancat == 0]
                if (max(meancat) > 0) {
                  for (i in 1:max(meancat)) {
                    w <- which(meancat == i)
                    if (length(w) == 1) {
                      mu[w] <- sum(postprobs[, w] * x)/(n * lambda[w])
                    }
                    else {
                      tmp <- t(postprobs[, w]) * (meanalpha[w]/sigma[w]^2)
                      mu[w] <- meanalpha[w] * sum(t(tmp) * x)/sum(tmp * 
                        meanalpha[w])
                    }
                  }
                }
                if (ECM) {
                  z <- pseudo.log.lik(n, k, x, mu, sigma, lambda)
                  postprobs <- matrix(z$post, nrow = n)
                  res <- matrix(z$res2, nrow = n)
                  lambda <- colMeans(postprobs)
                }
                sigma[sdcat == 0] <- sdalpha[sdcat == 0]
                if (max(sdcat) > 0) {
                  for (i in 1:max(sdcat)) {
                    w <- which(sdcat == i)
                    if (length(w) == 1) {
                      sigma[w] <- sqrt(sum(postprobs[, w] * res[, 
                        w])/(n * lambda[w]))
                    }
                    else {
                      tmp <- t(postprobs[, w])/sdalpha[w]
                      sigma[w] <- sdalpha[w] * sqrt(sum(t(tmp) * 
                        res[, w])/(n * sum(lambda[w])))
                    }
                  }
                  if (any(sigma < 1e-08)) {
                    notdone <- TRUE
                    cat("One of the variances is going to zero; ", 
                      "trying new starting values.\n")
                    restarts <- restarts + 1
                    lambda <- mu <- sigma <- NULL
                    if (restarts > maxrestarts) {
                      stop("Too many tries!")
                    }
                    break
                  }
                }
                z <- pseudo.log.lik(n, k, x, mu, sigma, lambda)
                postprobs <- matrix(z$post, nrow = n)
                res <- matrix(z$res2, nrow = n)
                newobsloglik <- z$loglik
                diff <- newobsloglik - obsloglik
                obsloglik <- newobsloglik
                ll <- c(ll, obsloglik)
                iter <- iter + 1
                if (verb) {
                  cat("iteration =", iter, " log-lik diff =", 
                    diff, " log-lik =", obsloglik, "\n")
                  print(rbind(lambda, mu, sigma))
                }
            }
        }
        if (iter == maxit) {
            cat("WARNING! NOT CONVERGENT!", "\n")
        }
        cat("number of iterations=", iter, "\n")
        if (arbmean == FALSE) {
            scale.order = order(sigma)
            sigma.min = min(sigma)
            postprobs = postprobs[, scale.order]
            colnames(postprobs) <- c(paste("comp", ".", 1:k, 
                sep = ""))
            a = list(x = x, lambda = lambda[scale.order], mu = mu, 
                sigma = sigma.min, scale = sigma[scale.order]/sigma.min, 
                loglik = obsloglik, posterior = postprobs, all.loglik = ll, 
                restarts = restarts, ft = "normalmixEM")
        }
        else {
            colnames(postprobs) <- c(paste("comp", ".", 1:k, 
                sep = ""))
            a = list(x = x, lambda = lambda, mu = mu, sigma = sigma, 
                loglik = obsloglik, posterior = postprobs, all.loglik = ll, 
                restarts = restarts, ft = "normalmixEM")
        }
    }
    class(a) = "mixEM"
    options(warn)
    ## Map to modified values to the actual values
    a$mu <-  loglikmap(a$mu)
    ## ordering <- order(a$mu)
    ## a$mu <- a$mu[ordering]
    ## a$lambda <- a$lambda[ordering]
    ## a$sigma <- a$sigma[ordering]
    a
}
