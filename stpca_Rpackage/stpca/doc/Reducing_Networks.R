## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)


## ----results='hide', message=FALSE, warning=FALSE------------------------
library(stpca) 


## ------------------------------------------------------------------------
data(demoNet)
names(demoNet)

## ---- echo=FALSE, message=FALSE, warning=FALSE---------------------------
library(SSN)

## ---- echo=FALSE---------------------------------------------------------
plot(demoNet, xlab="x-coordinate", ylab="y-coordinate",
     main="Observed and prediction locations")
plot(demoNet, PredPointsID="preds", add=TRUE, color.palette="blue")

## ---- echo=FALSE---------------------------------------------------------
plot(demoNet, "Sim_Values", xlab="x-coordinate", ylab="y-coordinate", 
     main="Response variable in simulated network")

river <- as.SpatialLines(demoNet)
coords <- getSSNdata.frame(demoNet)[,c("NEAR_X", "NEAR_Y", "strata")]
plot(summary(river)$bbox[,1], summary(river)$bbox[,2], type="n",
     xlab="x-coordinate", ylab="y-coordinate", main="Strata in simulated network",
     xlim=c(-5, 2), ylim=c(0, 14))
lines(river, col="darkgray", lwd=1.5)
points(coords[,1], coords[,2], col=coords[,3], pch=19)


## ------------------------------------------------------------------------
# get the data
x <- importSSN(system.file("demoSSN/demoNet.ssn", package = "stpca"),
               predpts = "preds")

# create sampled networks
test.samp <- sampNet(x, nsim = 3, prop.sites = c(0.5),
                     samp.scheme = c("random"),
                     sub.filepath = paste(tempdir(),"/subset1.ssn", sep = ""),
                     predpts = "preds", response.col = "Sim_Values",
                     addfunccol = "addfunccol", siteID = "locID")


## ------------------------------------------------------------------------
names(test.samp)

## ---- echo=FALSE---------------------------------------------------------
names(test.samp$random.50)

## ---- echo=FALSE---------------------------------------------------------
names(test.samp$all)

## ---- eval=FALSE---------------------------------------------------------
#  # get prediction points coordinates
#  library(SSN)
#  
#  # get data as a data.frame
#  x <- getSSNdata.frame(demoNet, "Obs")
#  
#  # get prediction locations
#  preds.coords <- getSSNdata.frame(demoNet, "preds")[,c("NEAR_X", "NEAR_Y")]
#  
#  # create sampled networks
#  test.samp <- sampNet(x, nsim = 3, prop.sites = c(0.5),
#                       samp.scheme = c("random"),
#                       data.col = "Sim_Values", coords.col = c(9, 10),
#                       preds.coords = preds.coords, siteID = "locID")
#  

## ---- eval=FALSE---------------------------------------------------------
#  set.seed(618999)
#  data(demoNet)
#  
#  # Observed values for the simulated network.
#  x <- getSSNdata.frame(demoNet, "Obs")
#  
#  # Prediction locations for the simulated network
#  preds.coords <- getSSNdata.frame(demoNet, "preds")[,c("NEAR_X", "NEAR_Y")]
#  
#  # Create the sampled networks
#  sampDemo <- sampNet(x, nsim = 5, prop.sites = c(0.8, 0.7),
#                      samp.scheme = c("random", "stratified", "weighted"),
#                      data.col = "Sim_Values", coords.col = c(9, 10),
#                      preds.coords = preds.coords, siteID = "locID",
#                      strat.col = "strata", weight.col = "weight")
#  

## ------------------------------------------------------------------------
data(sampDemo)
plot_sampNet(samp.results = sampDemo)

