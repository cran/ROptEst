\name{roptest}
\alias{roptest}
\alias{roptest.old}
\title{ Optimally robust estimation }
\description{
  Function to compute optimally robust estimates for L2-differentiable
  parametric families via k-step construction.
}
\usage{
roptest(x, L2Fam, eps, eps.lower, eps.upper, fsCor = 1, initial.est, 
        neighbor = ContNeighborhood(), risk = asMSE(), steps = 1L, 
        distance = CvMDist, startPar = NULL, verbose = NULL,
        OptOrIter = "iterate",
        useLast = getRobAStBaseOption("kStepUseLast"),
        withUpdateInKer = getRobAStBaseOption("withUpdateInKer"),
        IC.UpdateInKer = getRobAStBaseOption("IC.UpdateInKer"),
        withICList = getRobAStBaseOption("withICList"),
        withPICList = getRobAStBaseOption("withPICList"),
        na.rm = TRUE, initial.est.ArgList, ...,
        withLogScale = TRUE, ..withCheck = FALSE, withTimings = FALSE,
        withMDE = NULL, withEvalAsVar = NULL, withMakeIC = FALSE,
        modifyICwarn = NULL, E.argList = NULL, diagnostic = FALSE)
roptest.old(x, L2Fam, eps, eps.lower, eps.upper, fsCor = 1, initial.est,
        neighbor = ContNeighborhood(), risk = asMSE(), steps = 1L,
        distance = CvMDist, startPar = NULL, verbose = NULL,
        OptOrIter = "iterate",
        useLast = getRobAStBaseOption("kStepUseLast"),
        withUpdateInKer = getRobAStBaseOption("withUpdateInKer"),
        IC.UpdateInKer = getRobAStBaseOption("IC.UpdateInKer"),
        withICList = getRobAStBaseOption("withICList"),
        withPICList = getRobAStBaseOption("withPICList"),
        na.rm = TRUE, initial.est.ArgList, ...,
        withLogScale = TRUE)
}
\arguments{
  \item{x}{ sample }
  \item{L2Fam}{ object of class \code{"L2ParamFamily"} }
  \item{eps}{ positive real (0 < \code{eps} <= 0.5): amount of gross errors. 
          See details below. }
  \item{eps.lower}{ positive real (0 <= \code{eps.lower} <= \code{eps.upper}): 
        lower bound for the amount of gross errors. See details below. }
  \item{eps.upper}{ positive real (\code{eps.lower} <= \code{eps.upper} <= 0.5): 
        upper bound for the amount of gross errors. See details below. }
  \item{fsCor}{ positive real: factor used to correct the neighborhood radius;
        see details. }
  \item{initial.est}{ initial estimate for unknown parameter. If missing, a
        minimum distance estimator is computed. }
  \item{neighbor}{ object of class \code{"UncondNeighborhood"} }
  \item{risk}{ object of class \code{"RiskType"} }
  \item{steps}{ positive integer: number of steps used for k-steps construction }
  \item{distance}{ distance function used in \code{MDEstimator}, which in turn
        is used as (default) starting estimator. }
  \item{startPar}{ initial information used by \code{optimize} resp. \code{optim};
    i.e; if (total) parameter is of length 1, \code{startPar} is 
    a search interval, else it is an initial parameter value; if \code{NULL}
    slot \code{startPar} of \code{ParamFamily} is used to produce it;
    in the multivariate case, \code{startPar} may also be of class \code{Estimate},
    in which case slot \code{untransformed.estimate} is used.}
  \item{verbose}{ logical: if \code{TRUE}, some messages are printed }
  \item{useLast}{ which parameter estimate (initial estimate or
    k-step estimate) shall be used to fill the slots \code{pIC},
    \code{asvar} and \code{asbias} of the return value. }
  \item{OptOrIter}{character; which method to be used for determining Lagrange
  multipliers \code{A} and \code{a}: if (partially) matched to \code{"optimize"},
  \code{getLagrangeMultByOptim} is used; otherwise: by default, or if matched to
  \code{"iterate"} or to \code{"doubleiterate"},
  \code{getLagrangeMultByIter} is used. More specifically,
  when using \code{getLagrangeMultByIter}, and if argument \code{risk} is of
  class \code{"asGRisk"}, by default and if matched to \code{"iterate"}
  we use only one (inner) iteration, if matched to \code{"doubleiterate"}
  we use up to \code{Maxiter} (inner) iterations.}
  \item{withUpdateInKer}{if there is a non-trivial trafo in the model with matrix \eqn{D}, shall
  the parameter be updated on \eqn{{\rm ker}(D)}{ker(D)}?}
  \item{IC.UpdateInKer}{if there is a non-trivial trafo in the model with matrix \eqn{D},
  the IC to be used for this; if \code{NULL} the result of \code{getboundedIC(L2Fam,D)} is taken;
  this IC will then be projected onto \eqn{{\rm ker}(D)}{ker(D)}.}
  \item{withPICList}{logical: shall slot \code{pICList} of return value
  be filled?}
  \item{withICList}{logical: shall slot \code{ICList} of return value
  be filled?}
  \item{na.rm}{logical: if  \code{TRUE}, the estimator is evaluated at \code{complete.cases(x)}.}
  \item{initial.est.ArgList}{a list of arguments to be given to argument \code{start} if the latter
   is a function; this list by default already starts with two unnamed items,
   the sample \code{x}, and the model \code{L2Fam}.}
  \item{\dots}{ further arguments }
  \item{withLogScale}{logical; shall a scale component (if existing and found
   with name \code{scalename}) be computed on log-scale and backtransformed
   afterwards? This avoids crossing 0.}
   \item{..withCheck}{logical: if \code{TRUE}, debugging info is issued.}
  \item{withTimings}{logical: if \code{TRUE}, separate (and aggregate)
   timings for the three steps evaluating the starting value, finding
   the starting influence curve, and evaluating the k-step estimator is
   issued.}
  \item{withMDE}{ logical or \code{NULL}: Shall a minimum distance estimator be used as
       starting estimator---in addition to the function given in slot
       \code{startPar} of the L2 family? If \code{NULL} (default), the content
       of slot \code{.withMDE} in the L2 family is used instead to take
       this decision.}
  \item{withEvalAsVar}{logical or \code{NULL}: if \code{TRUE} (default), tells R
  to evaluate the asymptotic variance or if \code{FALSE} just to produces a call
  to do so. If \code{withEvalAsVar} is \code{NULL} (default),  the content
       of slot \code{.withEvalAsVar} in the L2 family is used instead to take
       this decision.}
  \item{withMakeIC}{logical; if \code{TRUE} the [p]IC is passed through
        \code{makeIC} before return.}
  \item{modifyICwarn}{logical: should a (warning) information be added if
     \code{modifyIC} is applied and hence some optimality information could
     no longer be valid? Defaults to \code{NULL} in which case this value
     is taken from \code{RobAStBaseOptions}.}
  \item{E.argList}{\code{NULL} (default) or a list of arguments to be passed
                   to calls to \code{E} from (a) \code{MDEstimator}
                   (here this additional argument is only used if
                   \code{initial.est} is missing), (b) \code{getStartIC},
                   and (c) \code{kStepEstimator}. Potential clashes with
        arguments of the same name in \code{\dots} are resolved by inserting
        the items of argument list \code{E.argList} as named items, so
        in case of collisions the item of \code{E.argList} overwrites the
        existing one from \code{\dots}.}
  \item{diagnostic}{ logical; if \code{TRUE},
    diagnostic information on the performed integrations is gathered and
    shipped out as attributes \code{kStepDiagnostic} (for the kStepEstimator-step)
    and \code{diagnostic} for the remaining steps of the return value
    of \code{roptest}. }
}
\details{
  Computes the optimally robust estimator for a given L2 differentiable
  parametric family. The computation uses a k-step construction with an 
  appropriate initial estimate; cf. also \code{\link[RobAStBase]{kStepEstimator}}. 
  Valid candidates are e.g. Kolmogorov(-Smirnov) or von Mises minimum 
  distance estimators (default); cf. Rieder (1994) and Kohl (2005).

  Before package version 0.9, this computation was done with the code of
  function \code{roptest.old} (with the same formals). From package version
  0.9 on, this function uses the modularized function \code{\link{robest}}
  internally.

  If the amount of gross errors (contamination) is known, it can be 
  specified by \code{eps}. The radius of the corresponding infinitesimal 
  contamination neighborhood is obtained by multiplying \code{eps} 
  by the square root of the sample size. 

  If the amount of gross errors (contamination) is unknown, try to find a 
  rough estimate for the amount of gross errors, such that it lies 
  between \code{eps.lower} and \code{eps.upper}.

  In case \code{eps.lower} is specified and \code{eps.upper} is missing, 
  \code{eps.upper} is set to 0.5. In case \code{eps.upper} is specified and
  \code{eps.lower} is missing, \code{eps.lower} is set to 0.

  If neither \code{eps} nor \code{eps.lower} and/or \code{eps.upper} is 
  specified, \code{eps.lower} and \code{eps.upper} are set to 0 and 0.5, 
  respectively.

  If \code{eps} is missing, the radius-minimax estimator in sense of 
  Rieder et al. (2001, 2008), respectively Section 2.2 of Kohl (2005) is returned.

  Finite-sample and higher order results suggest that the asymptotically
  optimal procedure is to liberal. Using \code{fsCor} the radius can be
  modified - as a rule enlarged - to obtain a more conservative estimate.
  In case of normal location and scale there is function 
  \code{\link[RobLox]{finiteSampleCorrection}} which returns a finite-sample 
  corrected (enlarged) radius based on the results of large Monte-Carlo
  studies.

  The logic in argument \code{initial.est} is as follows: It can be
  a numeric vector of the length of the unknow parameter or a function or
  it can be missing. If it is missing, one consults argument \code{startPar}
  for a search interval (if a one dimensional unknown parameter) or a starting
  value for the search (if the dimension of the unknown parameter is larger
  than one). If \code{startPar} is missing, too, it takes the value from
  the corresponding slot of argument \code{L2Fam}. Then, if argument \code{withMDE}
  is \code{TRUE} a Minimum-Distance estimator is computed as initial value
  \code{initial.est} with distance as specified in argument \code{distance}
  and possibly further arguments as passed through \code{\dots}.

  In the next step, the value of \code{initial.est} (either if not missing
  from beginning or as computed through the MDE) is then passed on to
  \code{kStepEstimator.start} which  then takes out the essential
  information for the sequel, i.e., a numeric vector of the estimate.

  At this initial value the optimal influence curve is computed through
  interface \code{getStartIC}, which in turn, depending on the risk calls
  \code{optIC}, \code{radiusMinimaxIC}, or computes the IC
  from precomputed grid values in case of \code{risk} being of class
  \code{interpolRisk}. With the obtained optimal IC, \code{kStepEstimator}
  is called.

  The default value of argument \code{useLast} is set by the
  global option \code{kStepUseLast} which by default is set to 
  \code{FALSE}. In case of general models \code{useLast} 
  remains unchanged during the computations. However, if 
  slot \code{CallL2Fam} of \code{IC} generates an object of 
  class \code{"L2GroupParamFamily"} the value of \code{useLast} 
  is changed to \code{TRUE}.
  Explicitly setting \code{useLast} to \code{TRUE} should
  be done with care as in this situation the influence curve
  is re-computed using the value of the one-step estimate
  which may take quite a long time depending on the model.

  If \code{useLast} is set to \code{TRUE} the computation of \code{asvar}, 
  \code{asbias} and \code{IC} is based on the k-step estimate.

  Timings for the steps run through in \code{roptest} are available
  in attributes \code{timings}, and for the step of the
  \code{kStepEstimator} in \code{kStepTimings}.

  One may also use the arguments \code{startCtrl}, \code{startICCtrl}, and
  \code{kStepCtrl} of function \code{\link{robest}}. This allows for individual
  settings of \code{E.argList}, \code{withEvalAsVar}, and
  \code{withMakeIC} for the different steps. If any of the three arguments
  \code{startCtrl}, \code{startICCtrl}, and \code{kStepCtrl} is used, the
  respective attributes set in the correspondig argument are used and, if
  colliding with arguments directly passed to \code{roptest}, the directly
  passed ones are ignored.

  Diagnostics on the involved integrations are available if argument
  \code{diagnostic} is \code{TRUE}. Then there are attributes \code{diagnostic}
  and \code{kStepDiagnostic} attached to the return value, which may be inspected
  and assessed through \code{\link[distrEx:distrExIntegrate]{showDiagnostic}} and
  \code{\link[distrEx:distrExIntegrate]{getDiagnostic}}.
}
\value{Object of class \code{"kStepEstimate"}.  In addition, it has
       an attribute \code{"timings"} where computation time is stored.}
\references{
  Kohl, M. (2005) \emph{Numerical Contributions to the Asymptotic Theory of Robustness}. 
  Bayreuth: Dissertation. \url{https://epub.uni-bayreuth.de/id/eprint/839/2/DissMKohl.pdf}.

  Kohl, M. and Ruckdeschel, P. (2010): R package distrMod: 
  Object-Oriented Implementation of Probability Models. 
  J. Statist. Softw.  \bold{35}(10), 1--27. \doi{10.18637/jss.v035.i10}.
  
  Kohl, M. and Ruckdeschel, P., and Rieder, H. (2010):
  Infinitesimally Robust Estimation in General Smoothly Parametrized Models.
  \emph{Stat. Methods Appl.}, \bold{19}, 333--354. 
  \doi{10.1007/s10260-010-0133-0}.

  Rieder, H. (1994) \emph{Robust Asymptotic Statistics}. New York: Springer.
  \doi{10.1007/978-1-4684-0624-5}.

  Rieder, H., Kohl, M. and Ruckdeschel, P. (2008) The Costs of not Knowing
  the Radius. Statistical Methods and Applications \bold{17}(1) 13-40.
  \doi{10.1007/s10260-007-0047-7}.
  
  Rieder, H., Kohl, M. and Ruckdeschel, P. (2001) The Costs of not Knowing
  the Radius. Appeared as discussion paper Nr. 81. 
  SFB 373 (Quantification and Simulation of Economic Processes),
  Humboldt University, Berlin; also available under
  \doi{10.18452/3638}
}
\author{Matthias Kohl \email{Matthias.Kohl@stamats.de},\cr
        Peter Ruckdeschel \email{peter.ruckdeschel@uni-oldenburg.de}}
%\note{}
\seealso{ \code{\link[RobLox]{roblox}}, 
          \code{\link[distrMod]{L2ParamFamily-class}}
          \code{\link[RobAStBase]{UncondNeighborhood-class}},
          \code{\link[distrMod]{RiskType-class}} }
\examples{
## Don't run to reduce check time on CRAN
\dontrun{
#############################
## 1. Binomial data
#############################
## generate a sample of contaminated data
set.seed(123)
ind <- rbinom(100, size=1, prob=0.05)
x <- rbinom(100, size=25, prob=(1-ind)*0.25 + ind*0.9)

## ML-estimate
MLest <- MLEstimator(x, BinomFamily(size = 25))
estimate(MLest)
confint(MLest)

## compute optimally robust estimator (known contamination)
robest1 <- roptest(x, BinomFamily(size = 25), eps = 0.05, steps = 3)
robest1.0 <- roptest.old(x, BinomFamily(size = 25), eps = 0.05, steps = 3)
identical(robest1,robest1.0)
estimate(robest1)
confint(robest1, method = symmetricBias())
## neglecting bias
confint(robest1)
plot(pIC(robest1))
tmp <- qqplot(x, robest1, cex.pch=1.5, exp.cex2.pch = -.25,
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
tmp <- qqplot(x, robest, cex.pch=1.5, exp.cex2.pch = -.25,
              exp.fadcol.pch = .55, jit.fac=.9)
 
## total variation neighborhoods (unknown deviation)
robest1 <- roptest(x, PoisFamily(), eps.upper = 0.05, 
                  neighbor = TotalVarNeighborhood(), steps = 3)
estimate(robest1)
confint(robest1, symmetricBias())
plot(pIC(robest1))
}

#############################
## 3. Normal (Gaussian) location and scale
#############################

## this example of a two dimensional parameter
## to be estimated will need more time than 
## 5 seconds to run 
## you can find it in 
## system.file("scripts", "examples_taking_longer.R", 
##              package="ROptEst")

}
\concept{k-step construction}
\concept{optimally robust estimation}
\keyword{robust}
