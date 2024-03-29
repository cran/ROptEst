###############################################################################
## Optimally robust IC for infinitesimal robust model and asymptotic risks
###############################################################################
setMethod("optIC", signature(model = "InfRobModel", risk = "asRisk"),
    function(model, risk, z.start = NULL, A.start = NULL, upper = 1e4,
             lower = 1e-4, OptOrIter = "iterate",
             maxiter = 50, tol = .Machine$double.eps^0.4,
             warn = TRUE, noLow = FALSE, verbose = NULL, ...,
             .withEvalAsVar = TRUE, withMakeIC = FALSE,
             returnNAifProblem = FALSE, modifyICwarn = NULL){
        if(missing(verbose)|| is.null(verbose))
           verbose <- getRobAStBaseOption("all.verbose")
        L2derivDim <- numberOfMaps(model@center@L2deriv)
        if(missing(warn)|| is.null(warn)) warn <- TRUE
        #L2Fam <- model@center
        #model@center <- moveL2Fam2RefParam(L2Fam)
        if(L2derivDim == 1){
            res <- getInfRobIC(L2deriv = model@center@L2derivDistr[[1]],
                        neighbor = model@neighbor, risk = risk, 
                        symm = model@center@L2derivDistrSymm[[1]],
                        Finfo = model@center@FisherInfo, trafo = trafo(model@center@param), 
                        upper = upper, lower = lower, maxiter = maxiter, tol = tol, warn = warn,
                        noLow = noLow, verbose = verbose, ...)
            linfo <- length(res$info)
            res$info <- if(linfo>1) c(rep("optIC",linfo),res$info) else c("optIC", res$info)
            res <- c(res, modifyIC = getModifyIC(L2FamIC = model@center,
                                                 neighbor = model@neighbor,
                                                 risk = risk, warn = warn, verbose = verbose))
            if(returnNAifProblem) if(!is.null(res$problem)) if(res$problem) return(NA)
            IC.o <- generateIC(model@neighbor, model@center, res)
        }else{
            if(is(model@center@distribution, "UnivariateDistribution")){
                if((length(model@center@L2deriv) == 1) & is(model@center@L2deriv[[1]], "RealRandVariable")){
                    L2deriv <- model@center@L2deriv[[1]]
                    L2derivSymm <- model@center@L2derivSymm
                    L2derivDistrSymm <- model@center@L2derivDistrSymm
                }else{
                    L2deriv <- diag(dimension(model@center@L2deriv)) %*% model@center@L2deriv
                    L2deriv <- RealRandVariable(Map = L2deriv@Map, Domain = L2deriv@Domain)
                    nrvalues <- numberOfMaps(L2deriv)
                    if(numberOfMaps(model@center@L2deriv) != nrvalues){
                        L1 <- vector("list", nrvalues)
                        L2 <- vector("list", nrvalues)
                        for(i in 1:nrvalues){
                            L1[[i]] <- NonSymmetric()
                            L2[[i]] <- NoSymmetry()
                        }
                        L2derivSymm <- new("FunSymmList", L1)
                        L2derivDistrSymm <- new("DistrSymmList", L2)
                    }
                }
                res <- getInfRobIC(L2deriv = L2deriv, neighbor = model@neighbor,
                            risk = risk,  Distr = model@center@distribution, 
                            DistrSymm = model@center@distrSymm, L2derivSymm = L2derivSymm,
                            L2derivDistrSymm = L2derivDistrSymm, Finfo = model@center@FisherInfo, 
                            trafo = trafo(model@center@param), z.start = z.start, A.start = A.start, 
                            upper = upper, lower = lower, OptOrIter = OptOrIter,
                            maxiter = maxiter, tol = tol, warn = warn,
                            verbose = verbose, ...,.withEvalAsVar = .withEvalAsVar)
                if(returnNAifProblem) if(!is.null(res$problem)) if(res$problem) return(NA)
                linfo <- length(res$info)
                res$info <- if(linfo>1) c(rep("optIC",linfo),res$info) else c("optIC", res$info)
                res <- c(res, modifyIC = getModifyIC(L2FamIC = model@center,
                                                 neighbor = model@neighbor,
                                                 risk = risk, warn = warn, verbose = verbose))
                IC.o <- generateIC(model@neighbor, model@center, res)
            }else{
                stop("not yet implemented")
            }
        }
        #IC.o <- moveICBackFromRefParam(IC.o,L2Fam)
        if(withMakeIC) IC.o <- makeIC(IC.o, model, ...)
        return(IC.o)
    })

###############################################################################
## Optimally robust IC for robust model with fixed neighborhood 
## and asymptotic under-/overshoot risk
###############################################################################
setMethod("optIC", signature(model = "InfRobModel", risk = "asUnOvShoot"),
    function(model, risk, upper = 1e4, lower = 1e-4, maxiter = 50,
             tol = .Machine$double.eps^0.4, withMakeIC = FALSE,
             warn = TRUE, verbose = NULL, modifyICwarn = NULL, ...){
        L2derivDistr <- model@center@L2derivDistr[[1]]
        if(missing(warn)|| is.null(warn)) warn <- TRUE
        if((length(model@center@L2derivDistr) == 1) & is(L2derivDistr, "UnivariateDistribution")){
            if(identical(all.equal(model@neighbor@radius, 0), TRUE)){
               return(optIC(model@center, risk = asCov(), withMakeIC = withMakeIC))
            }else{
               res <- getInfRobIC(L2deriv = L2derivDistr, 
                        neighbor = model@neighbor, risk = risk, 
                        symm = model@center@L2derivDistrSymm[[1]],
                        Finfo = model@center@FisherInfo, trafo = trafo(model@center@param), 
                        upper = upper, maxiter = maxiter, tol = tol, warn = warn)
               linfo <- length(res$info)
               if(is(model@neighbor, "ContNeighborhood")){
                  res$info <- c(rep("optIC",linfo+1),res$info,
                      "Optimal IC for 'InfRobModel' with 'ContNeighborhood'!!!")
               }else{
                  res$info <- if(linfo>1) c(rep("optIC",linfo),
                                            res$info) else c("optIC", res$info)
               }
               res <- c(res, modifyIC = getModifyIC(L2FamIC = model@center,
                                                 neighbor = model@neighbor,
                                                 risk = risk, warn = warn, verbose = verbose))
               IC.o <- generateIC(TotalVarNeighborhood(radius = model@neighbor@radius), model@center, res)
               if(withMakeIC) IC.o <- makeIC(IC.o, model, ...)
               return(IC.o)
           }    
        }else{
            stop("restricted to 1-dimensional parameteric models")
        }
    })

###############################################################################
## Optimally robust IC for robust model with fixed neighborhood 
## and finite-sample under-/overshoot risk
###############################################################################
setMethod("optIC", signature(model = "FixRobModel", risk = "fiUnOvShoot"),
    function(model, risk, sampleSize, upper = 1e4, lower = 1e-4, maxiter = 50,
             tol = .Machine$double.eps^0.4, withMakeIC = FALSE, warn = TRUE, Algo = "A",
             cont = "left", verbose = NULL, modifyICwarn = NULL, ...){
        if(missing(verbose)|| is.null(verbose))
           verbose <- getRobAStBaseOption("all.verbose")
        if(missing(warn)|| is.null(warn)) warn <- TRUE
        if(!identical(all.equal(sampleSize, trunc(sampleSize)), TRUE))
            stop("'sampleSize' has to be an integer > 0")
        if(is(model@center@distribution, "UnivariateDistribution")){
            res <- getFixRobIC(Distr = model@center@distribution, 
                        neighbor = model@neighbor, risk = risk, 
                        sampleSize = sampleSize, upper = upper, maxiter = maxiter, 
                        tol = tol, warn = warn, Algo = Algo, cont = cont)
            linfo <- length(res$info)
            if(is(model@neighbor, "ContNeighborhood")){
               res$info <- c(rep("optIC",linfo+1),res$info,
                   "Optimal IC for 'InfRobModel' with 'ContNeighborhood'!!!")
            }else{
               res$info <- if(linfo>1) c(rep("optIC",linfo),
                                         res$info) else c("optIC", res$info)
            }
            res <- c(res, modifyIC = getModifyIC(L2FamIC = model@center,
                                                 neighbor = model@neighbor,
                                                 risk = risk, warn = warn, verbose = verbose))
            IC.o <- generateIC(TotalVarNeighborhood(radius = model@neighbor@radius), model@center, res)
            if(withMakeIC) IC.o <- makeIC(IC.o, model, ...)
            return(IC.o)
        }else{
            stop("restricted to 1-dimensional parametric models")
        }
    })
