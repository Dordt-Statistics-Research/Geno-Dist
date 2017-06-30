library(parallel)
library(mixtools)

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


mymixEM <- function (x, lambda = NULL, mu = NULL, sigma = NULL, k = 2, mean.constr = NULL, 
    sd.constr = NULL, epsilon = 1e-08, maxit = 1000, maxrestarts = 20, 
    verb = FALSE, fast = FALSE, ECM = FALSE, arbmean = TRUE, 
    arbvar = TRUE, mean.limits = NULL) 
{
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

plot_data <- function(model, ...) {
    gaussian <- function(x, lambda, mu, sigma) {
        mixnorm <- function(x){
            sum(lambda/sqrt(2*pi*sigma^2)*exp(-(x-mu)^2 / (2*sigma^2)))
        }
        sapply(x, mixnorm)
    }
    
    plot.mymixEM(model, whichplots=2, ...)
    curve( gaussian(x, model$lambda, model$mu, model$sigma),
          from=min(model$x),
          to=max(model$x),
          add=TRUE,
          col="purple")
}

plot.mymixEM <- function (x, whichplots = 1, loglik = 1 %in% whichplots, density = 2 %in% 
    whichplots, xlab1 = "Iteration", ylab1 = "Log-Likelihood", 
    main1 = "Observed Data Log-Likelihood", col1 = 1, lwd1 = 2, 
    xlab2 = NULL, ylab2 = NULL, main2 = NULL, col2 = NULL, lwd2 = 2, 
    alpha = 0.05, marginal = FALSE, ...) 
{
    def.par <- par(ask = (loglik + density > 1), "mar")
    mix.object <- x
    if (!inherits(mix.object, "mixEM")) 
        stop("Use only with \"mixEM\" objects!")
    if (loglik) {
        plot(mix.object$all.loglik, xlab = xlab1, ylab = ylab1, 
            main = main1, type = "l", lwd = lwd1, col = col1, 
            ...)
    }
    if (density) {
        if (mix.object$ft == "logisregmixEM") {
            if (ncol(mix.object$x) != 2) {
                stop("The predictors must have 2 columns!")
            }
            if (sum((mix.object$y == 1) + (mix.object$y == 0)) != 
                length(mix.object$y)) {
                stop("The response must be binary!")
            }
            k = ncol(mix.object$beta)
            x = mix.object$x[, 2]
            if (is.null(main2)) {
                main2 <- "Most Probable Component Membership"
            }
            if (is.null(xlab2)) {
                xlab2 <- "Predictor"
            }
            if (is.null(ylab2)) {
                ylab2 <- "Response"
            }
            if (is.null(col2)) {
                col2 <- 2:(k + 1)
            }
            plot(x, mix.object$y, main = main2, xlab = xlab2, 
                ylab = ylab2, col = col2[apply(mix.object$posterior, 
                  1, which.max)], ...)
            a = cbind(x, mix.object$y)
            a = a[order(a[, 1]), ]
            for (i in 1:k) {
                lines(a[, 1], plogis(mix.object$beta[1, i] + 
                  mix.object$beta[2, i] * a[, 1]), col = col2[i])
            }
        }
        if (mix.object$ft == "normalmixEM") {
            k <- ncol(mix.object$posterior)
            x <- sort(mix.object$x)
            a <- hist(x, plot = FALSE)
            maxy <- max(max(a$density), 0.3989 * mix.object$lambda/mix.object$sigma)
            if (is.null(main2)) {
                main2 <- "Density Curves"
            }
            if (is.null(xlab2)) {
                xlab2 <- "Data"
            }
            if (is.null(col2)) {
                col2 <- 2:(k + 1)
            }
            hist(x, prob = TRUE, main = main2, xlab = xlab2, 
                ylim = c(0, maxy), ...)
            if (length(mix.object$mu) == 1) {
                arbvar <- TRUE
                mix.object$sigma <- mix.object$scale * mix.object$sigma
                arbmean <- FALSE
            }
            if (length(mix.object$mu) == k && length(mix.object$sigma) == 
                1) {
                arbmean <- TRUE
                arbvar <- FALSE
            }
            if (length(mix.object$sigma) == k && length(mix.object$mu) == 
                k) {
                arbmean <- TRUE
                arbvar <- TRUE
            }
            gaussian <- function(x, lambda, mu, sigma) {
                mixnorm <- function(x){
                    sum(lambda/sqrt(2*pi*sigma^2)*exp(-(x-mu)^2 / (2*sigma^2)))
                }
                sapply(x, mixnorm)
            }
            for (i in 1:k) {
                j <- order(mix.object$mu)[i]
                lines(x, mix.object$lambda[j] * dnorm(x, mean = mix.object$mu[j * 
                  arbmean + (1 - arbmean)], sd = mix.object$sigma[j * 
                  arbvar + (1 - arbvar)]), col = col2[i], lwd = lwd2)
            }
            curve( gaussian(x, mix.object$lambda, mix.object$mu, mix.object$sigma),
                  from=min(x),
                  to=max(x),
                  add=TRUE,
                  col=1)
        }
        if (mix.object$ft == "repnormmixEM") {
            x <- as.vector(as.matrix(x))
            k <- ncol(mix.object$posterior)
            x <- sort(mix.object$x)
            a <- hist(x, plot = FALSE)
            maxy <- max(max(a$density), 0.3989 * mix.object$lambda/mix.object$sigma)
            if (is.null(main2)) {
                main2 <- "Density Curves"
            }
            if (is.null(xlab2)) {
                xlab2 <- "Data"
            }
            if (is.null(col2)) {
                col2 <- 2:(k + 1)
            }
            hist(x, prob = TRUE, main = main2, xlab = xlab2, 
                ylim = c(0, maxy), ...)
            if (length(mix.object$mu) == 1) {
                arbvar <- TRUE
                mix.object$sigma = mix.object$scale * mix.object$sigma
                arbmean <- FALSE
            }
            if (length(mix.object$mu) == k && length(mix.object$sigma) == 
                1) {
                arbmean <- TRUE
                arbvar <- FALSE
            }
            if (length(mix.object$sigma) == k && length(mix.object$mu) == 
                k) {
                arbmean <- TRUE
                arbvar <- TRUE
            }
            for (i in 1:k) {
                lines(x, mix.object$lambda[i] * dnorm(x, mean = mix.object$mu[i * 
                  arbmean + (1 - arbmean)], sd = mix.object$sigma[i * 
                  arbvar + (1 - arbvar)]), col = col2[i], lwd = lwd2)
            }
        }
        if (mix.object$ft == "regmixEM.mixed") {
            x.1 = mix.object$x
            n = sum(sapply(x.1, nrow))
            x.1.sum = sum(sapply(1:length(x.1), function(i) length(x.1[[i]][, 
                1])))
            if (x.1.sum == n) {
                x = lapply(1:length(x.1), function(i) matrix(x.1[[i]][, 
                  -1], ncol = 1))
            }
            else {
                x = x.1
            }
            post.beta(x = x, y = mix.object$y, p.beta = mix.object$posterior.beta, 
                p.z = mix.object$posterior.z)
        }
        if (mix.object$ft == "mvnormalmixEM") {
            x = mix.object$x
            if (ncol(x) != 2) {
                stop("The data must have 2 columns!")
            }
            post = apply(mix.object$posterior, 1, which.max)
            k <- ncol(mix.object$posterior)
            if (is.list(mix.object$sigma)) {
                sigma = mix.object$sigma
            }
            else {
                sigma = lapply(1:k, function(i) mix.object$sigma)
            }
            if (is.list(mix.object$mu)) {
                mu = mix.object$mu
            }
            else {
                mu = lapply(1:k, function(i) mix.object$mu)
            }
            if (is.null(xlab2)) {
                xlab2 <- "X.1"
            }
            if (is.null(ylab2)) {
                ylab2 <- "X.2"
            }
            if (is.null(col2)) {
                col2 <- 2:(k + 1)
            }
            if (marginal == FALSE) {
                if (is.null(main2)) {
                  main2 <- "Density Curves"
                }
                plot(x, col = col2[post], xlab = xlab2, ylab = ylab2, 
                  main = main2, ...)
                lapply(1:k, function(i) points(mu[[i]][1], mu[[i]][2], 
                  pch = 19))
                for (i in 1:k) {
                  for (j in 1:length(alpha)) {
                    ellipse(mu = mu[[i]], sigma = sigma[[i]], 
                      alpha = alpha[j], col = col2[i])
                  }
                }
            }
            else {
                if (is.null(main2)) {
                  main2 <- ""
                }
                x <- mix.object$x[, 1]
                y <- mix.object$x[, 2]
                xhist <- hist(x, plot = FALSE)
                yhist <- hist(y, plot = FALSE)
                top <- max(c(xhist$counts, yhist$counts))
                xrange <- range(x)
                yrange <- range(y)
                nf <- layout(matrix(c(2, 0, 1, 3), 2, 2, byrow = TRUE), 
                  c(4, 1), c(1, 4), TRUE)
                layout.show(nf)
                par(mar = c(3, 3, 1, 1))
                plot(mix.object$x[, 1], mix.object$x[, 2], col = col2[post], 
                  xlab = xlab2, ylab = ylab2, main = main2, ...)
                lapply(1:k, function(i) points(mu[[i]][1], mu[[i]][2], 
                  pch = 19))
                for (i in 1:k) {
                  for (j in 1:length(alpha)) {
                    ellipse(mu = mu[[i]], sigma = sigma[[i]], 
                      alpha = alpha[j], col = (i + 1))
                  }
                }
                par(mar = c(0, 3, 1, 1))
                barplot(xhist$counts, axes = FALSE, ylim = c(0, 
                  top), space = 0, ...)
                par(mar = c(3, 0, 1, 1))
                barplot(yhist$counts, axes = FALSE, xlim = c(0, 
                  top), space = 0, horiz = TRUE, ...)
            }
        }
        if (mix.object$ft == "regmixEM") {
            if (ncol(mix.object$x) != 2) {
                stop("The predictors must have 2 columns!")
            }
            post <- apply(mix.object$posterior, 1, which.max)
            k <- ncol(mix.object$posterior)
            x <- mix.object$x[, 2]
            y <- mix.object$y
            n <- length(y)
            if (is.null(main2)) {
                main2 <- "Most Probable Component Membership"
            }
            if (is.null(xlab2)) {
                xlab2 <- "Predictor"
            }
            if (is.null(ylab2)) {
                ylab2 <- "Response"
            }
            if (is.null(col2)) {
                col2 <- 2:(k + 1)
            }
            plot(x, y, main = main2, xlab = xlab2, ylab = ylab2, 
                type = "n", ...)
            a = cbind(mix.object$x[, 2], mix.object$y, post)
            for (i in 1:k) {
                xy = subset(cbind(a, mix.object$posterior[, i]), 
                  a[, 3] == i)[, -3]
                xy = matrix(xy, ncol = 3)
                points(xy[, 1], xy[, 2], col = col2[i])
                if (is.matrix(mix.object$beta) == FALSE) {
                  abline(coef = mix.object$beta)
                  beta = matrix(mix.object$beta, ncol = k, nrow = 2)
                }
                else {
                  abline(coef = mix.object$beta[, i], col = col2[i])
                  beta = mix.object$beta
                }
                out = lm(y ~ x, weights = mix.object$posterior[, 
                  i])
                fit = beta[1, i] + beta[2, i] * x
                out.aov = anova(out)
                MSE = out.aov$Mean[2]
                xy.f = cbind(x, y, fit)
                xy.sort = xy.f[order(xy.f[, 1]), ]
                x.new = seq(from = min(x), to = max(x), length.out = 100)
                y.new = beta[1, i] + beta[2, i] * x.new
                s.h <- sqrt(MSE * (1/n + (x.new - mean(xy.sort[, 
                  1]))^2/var(xy.sort[, 1])/(n - 1)))
                for (j in 1:length(alpha)) {
                  W = sqrt(qf(1 - alpha[j], 2, n - 2))
                  upper = y.new + W * s.h
                  lower = y.new - W * s.h
                  lines(x.new, upper, col = (i + 1))
                  lines(x.new, lower, col = (i + 1))
                }
            }
        }
        if (mix.object$ft == "expRMM_EM") {
            plotexpRMM(mix.object, ...)
        }
        if (mix.object$ft == "weibullRMM_SEM") {
            plotweibullRMM(mix.object, ...)
        }
    }
    par(def.par)
}



## ## Run tests
## data1 <- c(rnorm(100,0,1),rnorm(100,1,1))
## myres1 <- mymixEM(data1, mean.limits=c(-1,1,0,2), sd.constr=c("a","a"))
## myres2 <- mymixEM(data1, mean.limits=c(-10,0,0,10), sd.constr=c("a","a"))
## mixres1 <- normalmixEM(data1, sd.constr=c("a","a"))
## pdf("dev6.pdf")
## par(mfrow=c(1,3))
## plot_data(myres1)
## plot_data(myres2)
## plot_data(mixres1)
## dev.off()


load("../../data/exploratory/minus_data_by_variety.RData")
names(minus_data_by_variety) <- lapply(minus_data_by_variety, function(x){names(x)<-as.character(x$plot_id[1])})
data <- minus_data_by_variety[["17-LDH-STN-SAG-1057"]]$height_diff_cm
mixes <- normalmixEM(data, mean.constr=c(0,NA))
mine <- mymixEM(data, mean.limits=c(-20,NA,NA,-20))
par(mfrow=c(1,2))
plot_data(mixes, main2="Mean constrained to 0")
plot_data(mine, main2="Limit one mean above -20 and the other below -20")
