pkgname <- "ROptEst"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
options(pager = "console")
library('ROptEst')

assign(".oldSearch", search(), pos = 'CheckExEnv')
cleanEx()
nameEx("0ROptEst-package")
### * 0ROptEst-package

flush(stderr()); flush(stdout())

### Name: ROptEst-package
### Title: Optimally robust estimation
### Aliases: ROptEst-package ROptEst
### Keywords: package

### ** Examples

library(ROptEst)

## Example: Rutherford-Geiger (1910); cf. Feller~(1968), Section VI.7 (a)
x <- c(rep(0, 57), rep(1, 203), rep(2, 383), rep(3, 525), rep(4, 532), 
       rep(5, 408), rep(6, 273), rep(7, 139), rep(8, 45), rep(9, 27), 
       rep(10, 10), rep(11, 4), rep(12, 0), rep(13, 1), rep(14, 1))

## ML-estimate from package distrMod
MLest <- MLEstimator(x, PoisFamily())
MLest
## confidence interval based on CLT
confint(MLest)

## compute optimally (w.r.t to MSE) robust estimator (unknown contamination)
robest <- roptest(x, PoisFamily(), eps.upper = 0.1, steps = 3)
estimate(robest)
## check influence curve
checkIC(pIC(robest))
## plot influence curve
plot(pIC(robest))
## confidence interval based on LAN - neglecting bias
confint(robest)
## confidence interval based on LAN - including bias
confint(robest, method = symmetricBias())



cleanEx()
nameEx("asAnscombe-class")
### * asAnscombe-class

flush(stderr()); flush(stdout())

### Name: asAnscombe-class
### Title: Asymptotic Anscombe risk
### Aliases: asAnscombe-class eff eff,asAnscombe-method
###   show,asAnscombe-method
### Keywords: classes

### ** Examples

new("asAnscombe")



cleanEx()
nameEx("asAnscombe")
### * asAnscombe

flush(stderr()); flush(stdout())

### Name: asAnscombe
### Title: Generating function for asAnscombe-class
### Aliases: asAnscombe
### Keywords: robust

### ** Examples

asAnscombe()

## The function is currently defined as
function(eff = .95, biastype = symmetricBias(), normtype = NormType()){ 
    new("asAnscombe", eff = eff, biastype = biastype, normtype = normtype) }



cleanEx()
nameEx("asL1-class")
### * asL1-class

flush(stderr()); flush(stdout())

### Name: asL1-class
### Title: Asymptotic mean absolute error
### Aliases: asL1-class
### Keywords: classes

### ** Examples

new("asMSE")



cleanEx()
nameEx("asL1")
### * asL1

flush(stderr()); flush(stdout())

### Name: asL1
### Title: Generating function for asMSE-class
### Aliases: asL1
### Keywords: robust

### ** Examples

asL1()

## The function is currently defined as
function(biastype = symmetricBias(), normtype = NormType()){ 
         new("asL1", biastype = biastype, normtype = normtype) }



cleanEx()
nameEx("asL4-class")
### * asL4-class

flush(stderr()); flush(stdout())

### Name: asL4-class
### Title: Asymptotic mean power 4 error
### Aliases: asL4-class
### Keywords: classes

### ** Examples

new("asMSE")



cleanEx()
nameEx("asL4")
### * asL4

flush(stderr()); flush(stdout())

### Name: asL4
### Title: Generating function for asL4-class
### Aliases: asL4
### Keywords: robust

### ** Examples

asL4()

## The function is currently defined as
function(biastype = symmetricBias(), normtype = NormType()){ 
         new("asL4", biastype = biastype, normtype = normtype) }



cleanEx()
nameEx("cniperCont")
### * cniperCont

flush(stderr()); flush(stdout())

### Name: cniperCont
### Title: Generic Functions for Computation and Plot of Cniper
###   Contamination and Cniper Points.
### Aliases: cniperCont cniperCont-methods
###   cniperCont,IC,IC,L2ParamFamily,ContNeighborhood,asMSE-method
###   cniperPoint cniperPoint-methods
###   cniperPoint,L2ParamFamily,ContNeighborhood,asMSE-method
###   cniperPointPlot cniperPointPlot-methods
###   cniperPointPlot,L2ParamFamily,ContNeighborhood,asMSE-method
### Keywords: robust

### ** Examples

## cniper contamination
P <- PoisFamily(lambda = 4)
RobP1 <- InfRobModel(center = P, neighbor = ContNeighborhood(radius = 0.1))
IC1 <- optIC(model=RobP1, risk=asMSE())
RobP2 <- InfRobModel(center = P, neighbor = ContNeighborhood(radius = 1))
IC2 <- optIC(model=RobP2, risk=asMSE())
cniperCont(IC1 = IC1, IC2 = IC2, L2Fam = P, 
           neighbor = ContNeighborhood(radius = 0.5), 
           risk = asMSE(),
           lower = 0, upper = 8, n = 101)

## cniper point plot
cniperPointPlot(P, neighbor = ContNeighborhood(radius = 0.5), 
                risk = asMSE(), lower = 0, upper = 10)

## cniper point
cniperPoint(P, neighbor = ContNeighborhood(radius = 0.5), 
            risk = asMSE(), lower = 0, upper = 4)
cniperPoint(P, neighbor = ContNeighborhood(radius = 0.5), 
            risk = asMSE(), lower = 4, upper = 8)



cleanEx()
nameEx("getL1normL2deriv")
### * getL1normL2deriv

flush(stderr()); flush(stdout())

### Name: getL1normL2deriv
### Title: Calculation of L1 norm of L2derivative
### Aliases: getL1normL2deriv getL1normL2deriv-methods
###   getL1normL2deriv,UnivariateDistribution-method
###   getL1normL2deriv,RealRandVariable-method
### Keywords: robust

### ** Examples

##



cleanEx()
nameEx("getL2normL2deriv")
### * getL2normL2deriv

flush(stderr()); flush(stdout())

### Name: getL2normL2deriv
### Title: Calculation of L2 norm of L2derivative
### Aliases: getL2normL2deriv
### Keywords: robust

### ** Examples

##



cleanEx()
nameEx("getMaxIneff")
### * getMaxIneff

flush(stderr()); flush(stdout())

### Name: getMaxIneff
### Title: getMaxIneff - computation of the maximal inefficiency of an IC
### Aliases: getMaxIneff
### Keywords: robust

### ** Examples

N0 <- NormLocationFamily(mean=2, sd=3)
## L_2 family + infinitesimal neighborhood
neighbor <- ContNeighborhood(radius = 0.5)
N0.Rob1 <- InfRobModel(center = N0, neighbor = neighbor)
## OBRE solution (ARE 95%)
N0.ICA <- optIC(model = N0.Rob1, risk = asAnscombe(.95))
## OMSE solution radius 0.5
N0.ICM <- optIC(model=N0.Rob1, risk=asMSE())
## RMX solution 
N0.ICR <- radiusMinimaxIC(L2Fam=N0, neighbor=neighbor,risk=asMSE())

getMaxIneff(N0.ICA,neighbor)
getMaxIneff(N0.ICM,neighbor)
getMaxIneff(N0.ICR,neighbor)

N0ls <- NormLocationScaleFamily()
ICsc <- makeIC(list(sin,cos),N0ls)
getMaxIneff(ICsc,neighbor)




cleanEx()
nameEx("getReq")
### * getReq

flush(stderr()); flush(stdout())

### Name: getReq
### Title: getReq - computation of the radius interval where IC1 is better
###   than IC2
### Aliases: getReq
### Keywords: robust

### ** Examples

N0 <- NormLocationFamily(mean=2, sd=3)
## L_2 family + infinitesimal neighborhood
neighbor <- ContNeighborhood(radius = 0.5)
N0.Rob1 <- InfRobModel(center = N0, neighbor = neighbor)
## OBRE solution (ARE 95%)
N0.ICA <- optIC(model = N0.Rob1, risk = asAnscombe(.95))
## MSE solution
N0.ICM <- optIC(model=N0.Rob1, risk=asMSE())
## RMX solution
N0.ICR <- radiusMinimaxIC(L2Fam=N0, neighbor=neighbor,risk=asMSE())

getReq(asMSE(),neighbor,N0.ICA,N0.ICM,n=1)
getReq(asMSE(),neighbor,N0.ICA,N0.ICM,n=30)
getReq(asL1(),neighbor,N0.ICA,N0.ICM,n=30)
getReq(asL4(),neighbor,N0.ICA,N0.ICM,n=30)
getReq(asMSE(),neighbor,N0.ICA,N0.ICR,n=30)
getReq(asL1(),neighbor,N0.ICA,N0.ICR,n=30)
getReq(asL4(),neighbor,N0.ICA,N0.ICR,n=30)
getReq(asMSE(),neighbor,N0.ICM,N0.ICR,n=30)

### when to use MAD and when Qn 
##  for Qn, see C. Croux, P. Rousseeuw (1993). Alternatives to the Median 
##      Absolute Deviation, JASA 88(424):1273-1283
L2M <- NormScaleFamily()
IC.mad <- makeIC(function(x)sign(abs(x)-qnorm(.75)),L2M)
d.qn <- (2^.5*qnorm(5/8))^-1
IC.qn <- makeIC(function(x) d.qn*(1/4 - pnorm(x+1/d.qn) + pnorm(x-1/d.qn)), L2M)
getReq(asMSE(), neighbor, IC.mad, IC.qn)
# => MAD is better once r > 0.5144 (i.e. for more than 2 outliers for n = 30)



cleanEx()
nameEx("leastFavorableRadius")
### * leastFavorableRadius

flush(stderr()); flush(stdout())

### Name: leastFavorableRadius
### Title: Generic Function for the Computation of Least Favorable Radii
### Aliases: leastFavorableRadius leastFavorableRadius-methods
###   leastFavorableRadius,L2ParamFamily,UncondNeighborhood,asGRisk-method
### Keywords: robust

### ** Examples

N <- NormLocationFamily(mean=0, sd=1) 
leastFavorableRadius(L2Fam=N, neighbor=ContNeighborhood(),
                     risk=asMSE(), rho=0.5)



cleanEx()
nameEx("lowerCaseRadius")
### * lowerCaseRadius

flush(stderr()); flush(stdout())

### Name: lowerCaseRadius
### Title: Computation of the lower case radius
### Aliases: lowerCaseRadius lowerCaseRadius-methods
###   lowerCaseRadius,L2ParamFamily,ContNeighborhood,asMSE,ANY-method
###   lowerCaseRadius,L2ParamFamily,TotalVarNeighborhood,asMSE,ANY-method
###   lowerCaseRadius,L2ParamFamily,ContNeighborhood,asMSE,onesidedBias-method
###   lowerCaseRadius,UnivariateDistribution,ContNeighborhood,asMSE,onesidedBias-method
###   lowerCaseRadius,L2ParamFamily,ContNeighborhood,asMSE,asymmetricBias-method
### Keywords: robust

### ** Examples

lowerCaseRadius(BinomFamily(size = 10), ContNeighborhood(), asMSE())
lowerCaseRadius(BinomFamily(size = 10), TotalVarNeighborhood(), asMSE())



cleanEx()
nameEx("optIC")
### * optIC

flush(stderr()); flush(stdout())

### Name: optIC
### Title: Generic function for the computation of optimally robust ICs
### Aliases: optIC optIC-methods optIC,InfRobModel,asRisk-method
###   optIC,InfRobModel,asUnOvShoot-method
###   optIC,FixRobModel,fiUnOvShoot-method
### Keywords: robust

### ** Examples

B <- BinomFamily(size = 25, prob = 0.25) 

## classical optimal IC
IC0 <- optIC(model = B, risk = asCov())
plot(IC0) # plot IC
checkIC(IC0, B)



cleanEx()
nameEx("optRisk")
### * optRisk

flush(stderr()); flush(stdout())

### Name: optRisk
### Title: Generic function for the computation of the minimal risk
### Aliases: optRisk optRisk-methods optRisk,L2ParamFamily,asCov-method
###   optRisk,InfRobModel,asRisk-method
###   optRisk,FixRobModel,fiUnOvShoot-method
### Keywords: robust

### ** Examples

optRisk(model = NormLocationScaleFamily(), risk = asCov())



cleanEx()
nameEx("radiusMinimaxIC")
### * radiusMinimaxIC

flush(stderr()); flush(stdout())

### Name: radiusMinimaxIC
### Title: Generic function for the computation of the radius minimax IC
### Aliases: radiusMinimaxIC radiusMinimaxIC-methods
###   radiusMinimaxIC,L2ParamFamily,UncondNeighborhood,asGRisk-method
### Keywords: robust

### ** Examples

N <- NormLocationFamily(mean=0, sd=1) 
radIC <- radiusMinimaxIC(L2Fam=N, neighbor=ContNeighborhood(), 
                         risk=asMSE(), loRad=0.1, upRad=0.5)
checkIC(radIC)



cleanEx()
nameEx("roptest")
### * roptest

flush(stderr()); flush(stdout())

### Name: roptest
### Title: Optimally robust estimation
### Aliases: roptest
### Keywords: robust

### ** Examples

#############################
## 1. Binomial data
#############################
## generate a sample of contaminated data
ind <- rbinom(100, size=1, prob=0.05) 
x <- rbinom(100, size=25, prob=(1-ind)*0.25 + ind*0.9)

## ML-estimate
MLest <- MLEstimator(x, BinomFamily(size = 25))
estimate(MLest)
confint(MLest)

## compute optimally robust estimator (known contamination)
robest1 <- roptest(x, BinomFamily(size = 25), eps = 0.05, steps = 3)
estimate(robest1)
confint(robest1, method = symmetricBias())
## neglecting bias
confint(robest1)
plot(pIC(robest1))
qqplot(x, robest1, cex.pch=1.5, exp.cex2.pch = -.25,
       exp.fadcol.pch = .55, jit.fac=.9)

## compute optimally robust estimator (unknown contamination)
robest2 <- roptest(x, BinomFamily(size = 25), eps.lower = 0, eps.upper = 0.2, steps = 3)
estimate(robest2)
confint(robest2, method = symmetricBias())
plot(pIC(robest2))

## total variation neighborhoods (known deviation)
robest3 <- roptest(x, BinomFamily(size = 25), eps = 0.025, 
                   neighbor = TotalVarNeighborhood(), steps = 3)
estimate(robest3)
confint(robest3, method = symmetricBias())
plot(pIC(robest3))

## total variation neighborhoods (unknown deviation)
robest4 <- roptest(x, BinomFamily(size = 25), eps.lower = 0, eps.upper = 0.1, 
                   neighbor = TotalVarNeighborhood(), steps = 3)
estimate(robest4)
confint(robest4, method = symmetricBias())
plot(pIC(robest4))


#############################
## 2. Poisson data
#############################
## Example: Rutherford-Geiger (1910); cf. Feller~(1968), Section VI.7 (a)
x <- c(rep(0, 57), rep(1, 203), rep(2, 383), rep(3, 525), rep(4, 532), 
       rep(5, 408), rep(6, 273), rep(7, 139), rep(8, 45), rep(9, 27), 
       rep(10, 10), rep(11, 4), rep(12, 0), rep(13, 1), rep(14, 1))

## ML-estimate
MLest <- MLEstimator(x, PoisFamily())
estimate(MLest)
confint(MLest)

## compute optimally robust estimator (unknown contamination)
robest <- roptest(x, PoisFamily(), eps.upper = 0.1, steps = 3)
estimate(robest)
confint(robest, symmetricBias())
plot(pIC(robest))
qqplot(x, robest, cex.pch=1.5, exp.cex2.pch = -.25,
       exp.fadcol.pch = .55, jit.fac=.9)
 
## total variation neighborhoods (unknown deviation)
robest1 <- roptest(x, PoisFamily(), eps.upper = 0.05, 
                  neighbor = TotalVarNeighborhood(), steps = 3)
estimate(robest1)
confint(robest1, symmetricBias())
plot(pIC(robest1))


#############################
## 3. Normal (Gaussian) location and scale
#############################
## 24 determinations of copper in wholemeal flour
library(MASS)
data(chem)
plot(chem, main = "copper in wholemeal flour", pch = 20)

## ML-estimate
MLest <- MLEstimator(chem, NormLocationScaleFamily())
estimate(MLest)
confint(MLest)

## compute optimally robust estimator (known contamination)
## takes some time -> you can use package RobLox for normal 
## location and scale which is optimized for speed
robest <- roptest(chem, NormLocationScaleFamily(), eps = 0.05, steps = 3)
estimate(robest)
confint(robest, symmetricBias())
plot(pIC(robest))
## plot of relative and absolute information; cf. Kohl (2005)
infoPlot(pIC(robest))

qqplot(chem, robest, cex.pch=1.5, exp.cex2.pch = -.25,
       exp.fadcol.pch = .55, withLab = TRUE, which.Order=1:4,
       exp.cex2.lbl = .12,exp.fadcol.lbl = .45,
       nosym.pCI = TRUE, adj.lbl=c(1.7,.2),
       exact.pCI = FALSE, log ="xy")

## finite-sample correction
if(require(RobLox)){
    n <- length(chem)
    r <- 0.05*sqrt(n)
    r.fi <- finiteSampleCorrection(n = n, r = r)
    fsCor <- r.fi/r
    robest <- roptest(chem, NormLocationScaleFamily(), eps = 0.05, 
                      fsCor = fsCor, steps = 3)
    estimate(robest)
}

## compute optimally robust estimator (unknown contamination)
## takes some time -> use package RobLox!
robest1 <- roptest(chem, NormLocationScaleFamily(), eps.lower = 0.05, 
                   eps.upper = 0.1, steps = 3)
estimate(robest1)
confint(robest1, symmetricBias())
plot(pIC(robest1))
## plot of relative and absolute information; cf. Kohl (2005)
infoPlot(pIC(robest1))



### * <FOOTER>
###
cat("Time elapsed: ", proc.time() - get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
