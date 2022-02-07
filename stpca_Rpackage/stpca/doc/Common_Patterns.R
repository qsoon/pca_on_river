## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)


## ----results='hide', message=FALSE, warning=FALSE------------------------
library(stpca) 


## ---- echo=FALSE---------------------------------------------------------
plot(1:10,1:10, type="n", xaxt="n", yaxt="n", bty="n", xlab="", ylab="")
lines(x=c(3, 5.5, 8), y=c(9, 6, 9))
arrows(x0=c(5.5), y0=c(6), x1=c(5.5), y1=c(1))
text(x=c(3.5, 7.5, 6), y=c(7.5, 7.5, 3), 
     labels=c("a1", "a2", "a3"))

## ------------------------------------------------------------------------
# load the demo dataframe containing spatiotemporal observations and 
# inspect the first few rows and columns
data(demoY)
head(demoY)[,1:7]


## ------------------------------------------------------------------------
data(demoNet)

plot(demoNet, "Sim_Values", 
     xlab="x-coordinate", ylab="y-coordinate",
     cex=1.5)

## ----results='hide', message=FALSE, warning=FALSE------------------------
library(SSN)

# get the filepath for the SSN object, including the .ssn extension
demoNet.path <- system.file("demoSSN/demoNet.ssn", package="stpca")

# import the SSN object.  Please refer to the importSSN() help file 
# in the SSN package for more details
demoNet <- importSSN(demoNet.path)



## ---- echo=FALSE---------------------------------------------------------
plot(1:25, demoY[,1], type="n", ylim=c(min(demoY), max(demoY)),
     xlab="Date", ylab="concentration (units)", main="Time series for 30 monitoring sites")
for(i in 1:ncol(demoY)) {
  lines(1:25, demoY[,i], col=i)
}

## ----results='hide', messages=FALSE--------------------------------------
spatial.wt <- createWeightS(ssndata = demoNet.path, 
                            afvcol = "addfunccol")

## ----messages=FALSE------------------------------------------------------
# display the first few rows and columns
head(spatial.wt)[,1:7]

## ------------------------------------------------------------------------

temporal.wt <- createWeightT(n = 25, rho = 0.75)
head(temporal.wt)[,1:7]

## ------------------------------------------------------------------------
data(demoYmiss) 
head(demoYmiss[,1:7])

## ----warning=FALSE, message=FALSE, results='hide'------------------------
Smode.impute <- completeData(x = demoYmiss, pca.mode = "Smode",
                             estim.plot = FALSE, 
                             mi.plot = FALSE, 
                             mean.plot = FALSE)

## ----warning=FALSE, message=FALSE----------------------------------------
names(Smode.impute)
names(Smode.impute$res.MIPCA)

## ----warning=FALSE, message=FALSE, results='hide'------------------------
Smode.impute <- completeData(x = demoYmiss, pca.mode = "Smode", 
                             ncp.min = 2, ncp.max = 4,
                             estim.plot = TRUE, 
                             mi.plot = TRUE,
                             mean.plot = TRUE,
                             response.var = "determinand (units)",
                             nbsim = 10, nboot = 10)

## ------------------------------------------------------------------------

Smode.pca.uw <- stpca(x = demoY,
                      pca.mode = "Smode", 
                      pca.wt = c("unweighted"))


## ------------------------------------------------------------------------
Smode.pca.uw$unweighted$var.prop[1]

## ------------------------------------------------------------------------
Smode.pca.uw$unweighted$cumvar[1:3]

## ------------------------------------------------------------------------
mean(colnames(spatial.wt) == colnames(demoY))


## ---- message=FALSE------------------------------------------------------

Smode.pca.all <- stpca(x = demoY,
                       pca.mode = "Smode", 
                       spatial.wt = spatial.wt, 
                       temporal.wt = temporal.wt,
                       pca.wt = c("unweighted", "spatial", 
                                  "temporal", "spatiotemporal"))

## ----results='hide'------------------------------------------------------
# first impute missing values
Tmode.impute <- completeData(x = demoYmiss, pca.mode = "Tmode",
                             mean.plot = FALSE, estim.plot = FALSE,
                             mi.plot = FALSE, 
                             nbsim = 10, nboot = 10)

## ------------------------------------------------------------------------
# transpose the completed data so that columns are monitoring sites and
# rows are time points
Tmode.completedata <- t(Tmode.impute$res.MIPCA$res.imputePCA)

# check column names are in same format/order as "spatial.wt"
mean(colnames(Tmode.completedata) == colnames(spatial.wt))
colnames(Tmode.completedata)[1:5]
colnames(spatial.wt)[1:5]

# need to remove "x" from column names
colnames(Tmode.completedata) <- 1:30

## ------------------------------------------------------------------------

Tmode.pca.uw <- stpca(x = Tmode.completedata,
                      pca.mode = "Tmode", 
                      pca.wt = c("unweighted"),
                      scree = FALSE)


## ------------------------------------------------------------------------

Tmode.pca.all <- stpca(x = Tmode.completedata,
                       pca.mode = "Tmode", 
                       spatial.wt = spatial.wt, 
                       temporal.wt = temporal.wt,
                       pca.wt = c("unweighted", "spatial", 
                                  "temporal", "spatiotemporal"),
                       scree = FALSE)


## ----message=FALSE, warning=FALSE, results='hide'------------------------
# load maptools library, required for reading in a .shp file
library(maptools)
# set the filepath for the .shp file
river.path <- paste(demoNet.path, "edges.shp", sep="/")
# read in the shapefile for the river network
river.dat <- readShapeLines(river.path)

# the coordinate information for the monitoring sites is stored here in the SSN 
# object for the demo dataset.  First, import the SSN object
library(SSN)
demo <- importSSN(demoNet.path)

# pull out the "Obs"" dataframe from the SSN object to get the coordinates for 
# the monitoring sites
coords <- getSSNdata.frame(demo, "Obs")[, c("NEAR_X", "NEAR_Y")]

## ------------------------------------------------------------------------
# produce the map plot
plot_stpca(x = Smode.pca.uw, plots = "map",
           river = river.dat, coords = coords)

## ------------------------------------------------------------------------
plot_stpca(x = Smode.pca.uw, plots = "biplot", 
           biplot.factor = 100)


## ---- tidy=TRUE----------------------------------------------------------

plot_stpca(x = Smode.pca.uw, plots  = "glyph", 
           river = river.dat, coords = coords, r1 = 10)

## ------------------------------------------------------------------------
plot_stpca(x = Smode.pca.uw, plots = "ts")

## ------------------------------------------------------------------------
plot_stpca(x = Smode.pca.uw, plots = "meanPM")

## ---- tidy=TRUE----------------------------------------------------------
plot_stpca(x = Tmode.pca.all, plots = "map",
          river = river.dat, coords = coords)

## ------------------------------------------------------------------------
plot_stpca(Tmode.pca.all, plots="biplot", biplot.factor=100)

## ------------------------------------------------------------------------
plot_stpca(x = Tmode.pca.all, plots = "glyph", 
           river = river.dat, coords = coords)

## ---- tidy=TRUE----------------------------------------------------------
plot_stpca(x = Tmode.pca.uw, plots = "ts")
plot_stpca(x = Tmode.pca.uw, plots = "meanPM")

