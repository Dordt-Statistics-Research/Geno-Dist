dirh <- "/home/jsnvndrw/MyWD/genodist"
setwd(dirh)

data <- read.csv("data/provided/17ldh_0315minus0302_DEM.csv")

sub <- data[data$plot_id==data$plot_id[1],]
hist(sub[[2]])
plot(density(sub[[2]]))

library("mclust")
library(parallel)
library(mixtools)

model <- Mclust(sub[[2]])
model$parameters

model <- Mclust(sub[[2]], G=2)
plot.Mclust(model, what="density")
lines(density(sub[[2]]), col="red")


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

minus_data_by_variety <- extract_data(data)

minus_models_2_comp <- mclapply(minus_data_by_variety, function(sub) {Mclust(sub[[2]], 2)})
minus_models_3_comp <- mclapply(minus_data_by_variety, function(sub) {Mclust(sub[[2]], 3)})


load("data/exploratory/minus_data_by_variety.RData")
load("data/exploratory/minus_models_2_comp.RData")

mixtools_model <- normalmixEM(minus_data_by_variety[[1]]$height_diff_cm)
mixtools_model_fix <- normalmixEM(minus_data_by_variety[[1]]$height_diff_cm, mean.constr=c(0,NA))
mixtools_model_fix_2 <- normalmixEM(minus_data_by_variety[[1]]$height_diff_cm, mean.constr=c(NA,0))

mixtools_model_fix <- normalmixEM(minus_data_by_variety[[1]]$height_diff_cm, mu=c(0,1), mean.constr=c(0,NA))
mixtools_model_fix <- normalmixEM(minus_data_by_variety[[1]]$height_diff_cm, mu=c(1,0), mean.constr=c(NA,0))

mixtools_model_fix <- normalmixEM(minus_data_by_variety[[1]]$height_diff_cm, mu=c(0,1), mean.constr=c(NA,NA), verb=T)


pdf("dev.pdf",width=7,height=5)
plot(mixtools_model_fix <- normalmixEM(minus_data_by_variety[[1]]$height_diff_cm, mean.constr=c(NA,0)), whichplot=2)
plot(mixtools_model_fix <- normalmixEM(minus_data_by_variety[[1]]$height_diff_cm, mean.constr=c(0,NA)), whichplot=2)
plot(mixtools_model_fix <- normalmixEM(minus_data_by_variety[[1]]$height_diff_cm), whichplot=2)
dev.off()

names(minus_data_by_variety) <- lapply(minus_data_by_variety, function(x){names(x)<-as.character(x$plot_id[1])})

gaussian <- function(x, lambda, mu, sigma) {
    mixnorm <- function(x){
        sum(lambda/sqrt(2*pi*sigma^2)*exp(-(x-mu)^2 / (2*sigma^2)))
    }
    #temp <- Map(norm, x, lambda, mu, sigma)
    sapply(x, mixnorm)
}

plot_data <- function(data, ...){
    invisible(model <- normalmixEM(data, ...))
    xmin <- min(data)
    xmax <- max(data)
    plot(model, whichplot=2)
    curve( gaussian(x, model$lambda, model$mu, model$sigma),
          from=xmin,
          to=xmax,
          add=TRUE,
          col="purple", ...)
}

sample_plots <- c(1,1009,1010,1029,1044,1046,1057,5071,5072,5073,6093,8005,8006,8007,11072,11073,11074)


pdf("data/exploratory/samples7.pdf",width=14,height=8.5)
par(mfrow=c(1,4))
for (sample in sample_plots) {
##for (i in 1:20) {
    name <- paste0("17-LDH-STN-SAG-", sample)
    ##name <- paste0("17-LDH-STN-SAG-", sample_plots[9])
    print(name)
    invisible(data <- minus_data_by_variety[[name]]$height_diff_cm)
    print(length(data))

    plot_data(data)
    plot_data(data, mean.constr=c(0,NA), sd.constr=c("a","a"))
    plot_data(data, mean.constr=c(NA,0), sd.constr=c("a","a"))
    plot_data(data, sd.constr=c("a","a"))
}
dev.off()
