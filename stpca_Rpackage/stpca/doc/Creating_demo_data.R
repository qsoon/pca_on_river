## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)


## ---- echo=FALSE, warning=FALSE, message=FALSE---------------------------
library(stpca)
library(SSN)
data(demoNet)
plot(demoNet, "Sim_Values", 
     xlab="x-coordinate", ylab="y-coordinate",
     cex=1.5)


## ---- echo=FALSE---------------------------------------------------------
data(demoY)
plot(1:25, demoY[,1], type="l", ylim=c(min(demoY), max(demoY)),
     xlab="Time", ylab="Sim_Values", main="simulated time series")
for(i in 1:30) {
  lines(1:25, demoY[,i], col=i)
}


