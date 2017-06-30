library(parallel)
library(mixtools)

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

load("../../data/exploratory/minus_data_by_variety.RData")
names(minus_data_by_variety) <- lapply(minus_data_by_variety, function(x){names(x)<-as.character(x$plot_id[1])})
data <- minus_data_by_variety[["17-LDH-STN-SAG-1057"]]$height_diff_cm
mixes <- normalmixEM(data, mean.constr=c(0,NA))
mine <- mymixEM(data, mean.limits=c(-20,NA,NA,-20))
par(mfrow=c(1,2))
plot(mixes, main2="Mean constrained to 0", whichplot=2)
plot(mine, main2="Limit one mean above -20 and the other below -20", whichplot=2)
