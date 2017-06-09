dirh <- "/home/jsnvndrw/MyWD/genodist"

data <- read.csv("../../data/provided/17ldh_0315minus0302_DEM.csv")

sub <- data[data$plot_id==data$plot_id[1],]
hist(sub[[2]])
plot(density(sub[[2]]))

library("mclust")
library(parallel)

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

