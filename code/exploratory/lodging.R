dirh <- "/home/jsnvndrw/MyWD/genodist"
setwd(dirh)

## data <- read.csv("data/provided/17ldh_0315minus0302_DEM.csv")

## sub <- data[data$plot_id==data$plot_id[1],]
## hist(sub[[2]])
## plot(density(sub[[2]]))

library("mclust")
library(parallel)
library(mixtools)

## model <- Mclust(sub[[2]])
## model$parameters

## model <- Mclust(sub[[2]], G=2)
## plot.Mclust(model, what="density")
## lines(density(sub[[2]]), col="red")


plot.real_v_mclust<- function(x,G) {
    model <- Mclust(x,G=G)
    plot.Mclust(model, what="density")
    lines(density(sub[[2]]), col="red")
}

extract_data <- function(data) {
    stopifnot(require(parallel))
    ids <- sort(unique(data$plot_id))
    data_sets <- mclapply(ids, function(id) {data[data$plot_id==id,]}, mc.cores=20)
}

## minus_data_by_variety <- extract_data(data)

## minus_models_2_comp <- mclapply(minus_data_by_variety, function(sub) {Mclust(sub[[2]], 2)})
## minus_models_3_comp <- mclapply(minus_data_by_variety, function(sub) {Mclust(sub[[2]], 3)})


load("data/exploratory/minus_data_by_variety.RData")
## load("data/exploratory/minus_models_2_comp.RData")

## mixtools_model <- normalmixEM(minus_data_by_variety[[1]]$height_diff_cm)
## mixtools_model_fix <- normalmixEM(minus_data_by_variety[[1]]$height_diff_cm, mean.constr=c(0,NA))
## mixtools_model_fix_2 <- normalmixEM(minus_data_by_variety[[1]]$height_diff_cm, mean.constr=c(NA,0))

## mixtools_model_fix <- normalmixEM(minus_data_by_variety[[1]]$height_diff_cm, mu=c(0,1), mean.constr=c(0,NA))
## mixtools_model_fix <- normalmixEM(minus_data_by_variety[[1]]$height_diff_cm, mu=c(1,0), mean.constr=c(NA,0))

## mixtools_model_fix <- normalmixEM(minus_data_by_variety[[1]]$height_diff_cm, mu=c(0,1), mean.constr=c(NA,NA), verb=T)


## pdf("dev.pdf",width=7,height=5)
## plot(mixtools_model_fix <- normalmixEM(minus_data_by_variety[[1]]$height_diff_cm, mean.constr=c(NA,0)), whichplot=2)
## plot(mixtools_model_fix <- normalmixEM(minus_data_by_variety[[1]]$height_diff_cm, mean.constr=c(0,NA)), whichplot=2)
## plot(mixtools_model_fix <- normalmixEM(minus_data_by_variety[[1]]$height_diff_cm), whichplot=2)
## dev.off()

names(minus_data_by_variety) <- lapply(minus_data_by_variety, function(x){names(x)<-as.character(x$plot_id[1])})

gaussian <- function(x, lambda, mu, sigma) {
    mixnorm <- function(x){
        sum(lambda/sqrt(2*pi*sigma^2)*exp(-(x-mu)^2 / (2*sigma^2)))
    }
    #temp <- Map(norm, x, lambda, mu, sigma)
    sapply(x, mixnorm)
}

BIC <- function(n, k, l) {
    k*log(n) - 2*l
}

plot.NULLmixEM <- function(model, ...) {
    plot(x=
}

plot_data <- function(data, mean.constr = NULL, sd.constr = NULL, ...){
    ## invisible(models <- mclapply(1:20, function(x) {normalmixEM(data, mean.constr = mean.constr, sd.constr = sd.constr)}, mc.cores=20))
    invisible(models <- mclapply(1:20, function(x) {
        tryCatch(
        {
            normalmixEM(data, mean.constr = mean.constr, sd.constr = sd.constr)
        },
        error=function(cond) {
            NULL
            message("JASON: ERROR")
        })
    }))
    print("here")
    ## browser()
    ## Remove models that aren't of class mixEM. Examples include elements with class try-error
    ## Then only 
    models <- models[sapply(models, inherits, "mixEM")]
    if (length(models) == 0) {
        model <- NULL
    } else {
        model <- models[[which.max(sapply(models, function(model) {model$loglik}))]]
    }
    xmin <- min(data)
    xmax <- max(data)
    if (is.null(model)) {
        hist(data)
    } else {
        plot(model, whichplot=2, xlab2 = "Height Difference (cm)", ylab2 = "Density", ...)
        curve( gaussian(x, model$lambda, model$mu, model$sigma),
              from=xmin,
              to=xmax,
              add=TRUE,
              col="purple")
    }
}
 
sample_plots <- c(1009,1010,1029,1044,1046,1057,5071,5072,5073,6093,8005,8006,8007,11072,11073,11074)
## sample_plots <- rep(8006,16)

doit <- function() {
    pdf("data/exploratory/samples10.pdf",width=14,height=8.5)
    ##pdf("graphs/prelimresults.pdf", width=11, height=8.5)
    par(mfrow=c(1,4))
    for (sample in sample_plots) {         
        name <- paste0("17-LDH-STN-SAG-", sample)
        print(name)
        invisible(data <- minus_data_by_variety[[name]]$height_diff_cm)
        print(length(data))   
        plot_data(data, main2 = name, sub = "Unconstrained")
        plot_data(data, mean.constr=c(NA,0), main2 = name, sub = "Mean Constrained")
        plot_data(data, sd.constr=c("a","a"), main2 = name, sub = "SD Constrained")
        plot_data(data, mean.const=c(NA,0), sd.constr=c("a","a"), main2 = name, sub = "Both Constrained")
    }
    dev.off()
}


do_it <- function() {
    out <- tryCatch(
        {
            # Just to highlight: if you want to use more than one 
            # R expression in the "try" part then you'll have to 
            # use curly brackets.
            # 'tryCatch()' will return the last evaluated expression 
            # in case the "try" part was completed successfully

            message("This is the 'try' part")

            normalmixEM(minus_data_by_variety[["17-LDH-STN-SAG-8006"]]$height_diff_cm, mean.constr=c(NA,0))

            # The return value of `readLines()` is the actual value 
            # that will be returned in case there is no condition 
            # (e.g. warning or error). 
            # You don't need to state the return value via `return()` as code 
            # in the "try" part is not wrapped insided a function (unlike that
            # for the condition handlers for warnings and error below)
        },
        error=function(cond) {
            message("THIS IS MY ERROR MESSAGE")
            message("Here's the original error message:")
            message(cond)
            # Choose a return value in case of error
            return("JASON THIS IS AN ERROR")
        },
        warning=function(cond) {
            message("THIS IS MY WARNGIN MESSAGE")
            message("Here's the original warning message:")
            message(cond)
            # Choose a return value in case of warning
            return("JASON THIS IS A WARNING")
        },
        finally={
        # NOTE:
        # Here goes everything that should be executed at the end,
        # regardless of success or error.
        # If you want more than one expression to be executed, then you 
        # need to wrap them in curly brackets ({...}); otherwise you could
        # just have written 'finally=<expression>' 
            message("Some ending message")
            message("Some other message at the end")
        }
    )    
    return(out)
}
