load("../../data/exploratory/minus_data_by_variety.RData")
names(minus_data_by_variety) <- lapply(minus_data_by_variety, function(x){names(x)<-as.character(x$plot_id[1])})

sample_plots <- c(1009,1010,1029,1044,1046,1057,5071,5072,5073,6093,8005,8006,8007,11072,11073,11074)
sample_plots <- c(1057)




projdir <- "/home/jsnvndrw/MyWD/genodist"
projfile <- function(filename) {
    projdir <- gsub("/*$", "", projdir) #remove trailing slash
    paste0(projdir, "/", filename)                     #return full file path
}

source.util <- function(filename) {
    source(paste0(projfile("code/utilities/"), filename))
}

source.util("mymixEM.R")
source.util("plot.mixEM.R")


plot_data <- function(){
    mtext(paste("BIC =",round(BIC(length(data), free_params, model$loglik))))
}

best_model <- function(x, n, fun, ...) {
    if (n < 1)
        stop("Choose n >= 1")
    models <- lapply(1:n, function(i) {
        fun(x, ...)
    })
    logliks <- mclapply(models, function(model) {
        model$loglik
    }, mc.cores=20)
    if (all(!sapply(models, inherits, "mixEM"))) {
        browser()
    }
    models[[which.max(logliks)]]
}

doit <- function() {
    pdf(projfile("data/exploratory/comp-0-lim-(3).pdf"),width=14,height=8.5)
    par(mfrow=c(1,2))
    for (sample in sample_plots) {         
        name <- paste0("17-LDH-STN-SAG-", sample)
        print(name)
        invisible(data <- minus_data_by_variety[[name]]$height_diff_cm)
        print(length(data))

        mine <- best_model(data, 20, mymixEM, mean.limits=c(-15,NA,NA,-15))
        theirs <- best_model(data, 20, normalmixEM, mean.constr=c(NA,0))
        
        plot(
            theirs,
            main2=name, sub="Mean Constrained to 0", whichplot=2)
        mtext(paste("Liklihood =", theirs$loglik))
        plot(
            mine,
            main2=name, sub="Mean Limited Around -15", whichplot=2)
        mtext(paste("Liklihood =", mine$loglik))
    }
    dev.off()
}
