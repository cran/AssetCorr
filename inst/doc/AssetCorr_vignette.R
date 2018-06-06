## ----setup, include=FALSE,warning=TRUE-----------------------------------
library(knitr)

## ------------------------------------------------------------------------
library(AssetCorr)
set.seed(111)

#number of obligors: 1000
#intra asset correlation: 0.3
#length of the default time series: 20
#probability of default: 0.01

D1=defaultTimeseries(1000,0.3,20,0.01)
N1=rep(1000,20)

## ------------------------------------------------------------------------
intraAMM(D1,N1)

## ------------------------------------------------------------------------
Output<-intraAMM(D1,N1, B=1000, CI_Boot = 0.95, plot=TRUE)

Output$Original
Output$Bootstrap
Output$CI_Boot

## ------------------------------------------------------------------------
Output<-intraAMM(D1,N1, DB=c(100,100))

Output$Original
Output$Bootstrap
Output$Double_Bootstrap


## ------------------------------------------------------------------------
intraAMM(D1,N1, JC=TRUE)


## ------------------------------------------------------------------------
Output<-intraALL(D1,N1, B=500, plot=TRUE,Adjust = 0.0001, Estimator = c("AMM","FMM","CMM","JDP1","JDP2","AMLE") )
Output


## ------------------------------------------------------------------------
library(mvtnorm)
set.seed(2)
#number of obligors: 1000
#intra asset correlation 1: 0.3
#intra asset correlation 2: 0.1
#inter correlation of the systematic factors: 0.5
#length of the default time series: 20
#probability of default: 0.01

Psi=rmvnorm(20,sigma=matrix(c(1,0.5,0.5,1),2))
PDcond1=pnorm((qnorm(0.01)-sqrt(0.3)*Psi[,1])/sqrt(1-0.3))
PDcond2=pnorm((qnorm(0.01)-sqrt(0.1)*Psi[,2])/sqrt(1-0.1))

D1=rbinom(20,1000,PDcond1)
D2=rbinom(20,1000,PDcond2)
N1=N2=rep(1000,20)

## ------------------------------------------------------------------------
rho1=intraAMM(D1,N1)$Original
rho2=intraAMM(D2,N2)$Original

interCov(D1,N1,D2,N2,rho1,rho2)

## ------------------------------------------------------------------------
#Single bootstrap Correction
rho1=intraAMM(D1,N1)$Original
rho2=intraAMM(D2,N2)$Original

Output<- interCov(D1,N1,D2,N2,rho1,rho2, B=1000, CI_Boot = 0.95, plot=TRUE)

Output$Original
Output$Bootstrap
Output$CI_Boot

#Double bootstrap correction

Output<- interCov(D1,N1,D2,N2,rho1,rho2, DB=c(100,100))

Output$Original
Output$Bootstrap
Output$Double_Bootstrap

#Furthermore, a Jackknife correction would be possible

