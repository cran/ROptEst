\name{cniperCont}
\alias{cniperCont}
\alias{cniperPoint}
\alias{cniperPointPlot}
\title{ Functions for Computation and Plot of Cniper Contamination
        and Cniper Points. }
\description{
  These functions and their methods can be used to determine cniper
  contamination as well as cniper points. That is, under which (Dirac) 
  contamination is the risk of one procedure larger than the risk of some
  other procedure.
}
\usage{
cniperCont(IC1, IC2, data = NULL, ...,
           neighbor, risk, lower=getdistrOption("DistrResolution"),
           upper=1-getdistrOption("DistrResolution"), n = 101,
           with.automatic.grid = TRUE, scaleX = FALSE, scaleX.fct,
           scaleX.inv, scaleY = FALSE, scaleY.fct = pnorm, scaleY.inv=qnorm,
           scaleN = 9, x.ticks = NULL, y.ticks = NULL, cex.pts = 1,
           cex.pts.fun = NULL, col.pts = par("col"), pch.pts = 19,
           cex.npts = 0.6, cex.npts.fun = NULL, col.npts = "red", pch.npts = 20,
           jit.fac = 1, jit.tol = .Machine$double.eps, with.lab = FALSE,
           lab.pts = NULL, lab.font = NULL, alpha.trsp = NA, which.lbs = NULL,
           which.Order  = NULL, which.nonlbs = NULL, attr.pre = FALSE,
           return.Order = FALSE, withSubst = TRUE)

cniperPoint(L2Fam, neighbor, risk, lower, upper)

cniperPointPlot(L2Fam, data=NULL, ..., neighbor, risk= asMSE(),
                        lower=getdistrOption("DistrResolution"),
                        upper=1-getdistrOption("DistrResolution"), n = 101,
                        withMaxRisk = TRUE, with.automatic.grid = TRUE,
                           scaleX = FALSE, scaleX.fct, scaleX.inv,
                           scaleY = FALSE, scaleY.fct = pnorm, scaleY.inv=qnorm,
                           scaleN = 9, x.ticks = NULL, y.ticks = NULL,
                           cex.pts = 1, cex.pts.fun = NULL, col.pts = par("col"),
                           pch.pts = 19,
                           cex.npts = 1, cex.npts.fun = NULL, col.npts = par("col"),
                           pch.npts = 19,
                           jit.fac = 1, jit.tol = .Machine$double.eps,
                           with.lab = FALSE,
                           lab.pts = NULL, lab.font = NULL, alpha.trsp = NA,
                           which.lbs = NULL, which.nonlbs = NULL,
                           which.Order  = NULL, attr.pre = FALSE, return.Order = FALSE,
                           withSubst = TRUE, withMakeIC = FALSE)
}
\arguments{
  \item{IC1}{ object of class \code{IC} }
  \item{IC2}{ object of class \code{IC} }
  \item{L2Fam}{ object of class \code{L2ParamFamily} }
  \item{neighbor}{ object of class \code{Neighborhood} }
  \item{risk}{ object of class \code{RiskType} }
  \item{\dots}{ additional parameters (in particular to be passed on to \code{plot}). }
  \item{data}{data to be plotted in}
  \item{lower, upper}{ the lower and upper end points of the 
          contamination interval (in prob-scale). }
  \item{n}{ number of points between \code{lower} and \code{upper}}
  \item{withMaxRisk}{logical; if \code{TRUE}, for risk comparison
     uses the maximal risk of the classically optimal IC \eqn{\psi}{psi} in all
     situations with contamination in Dirac points 'no larger' than
     the respective evaluation point and the optimally-robust
     IC \eqn{\eta}{eta} at its least favorable contamination situation
     ('over all real Dirac contamination points'). This is the default and
     was the  behavior prior to package version 0.9).
     If \code{FALSE} it uses exactly the situation
     with Dirac contamination in the evaluation point for both ICs
     \eqn{\psi}{psi} and \eqn{\eta}{eta} which amounts to calling \code{cniperCont}
     with \code{IC1=psi}, \code{IC2=eta}.}
  \item{with.automatic.grid}{logical; should a grid be plotted alongside
      with the ticks of the axes, automatically? If \code{TRUE} a respective
      call to \code{grid} in argument \code{panel.first} is ignored. }
  \item{scaleX}{logical; shall X-axis be rescaled (by default according to the cdf of
          the underlying distribution)?}
  \item{scaleY}{logical; shall Y-axis be rescaled (by default according to a probit scale)?}
  \item{scaleX.fct}{an isotone, vectorized function mapping the domain of the IC(s)
            to [0,1]; if \code{scaleX} is \code{TRUE} and \code{scaleX.fct} is
            missing, the cdf of the underlying observation distribution.}
  \item{scaleX.inv}{the inverse function to \code{scale.fct}, i.e., an isotone,
            vectorized function mapping [0,1] to the domain of the IC(s)
            such that for any \code{x} in the domain,
            \code{scaleX.inv(scaleX.fct(x))==x}; if \code{scaleX} is \code{TRUE}
            and \code{scaleX.inv} is
            missing, the quantile function of the underlying observation distribution.}
  \item{scaleY.fct}{an isotone, vectorized function mapping for each coordinate the
            range of the respective coordinate of the IC(s)
            to [0,1]; defaulting to the cdf of \eqn{{\cal N}(0,1)}{N(0,1)}.}
 \item{scaleY.inv}{an isotone, vectorized function mapping for each coordinate
            the range [0,1] into the range of the respective coordinate of the IC(s);
            defaulting to the quantile function of  \eqn{{\cal N}(0,1)}{N(0,1)}.}
  \item{scaleN}{integer; defaults to 9; on rescaled axes, number of x
                and y ticks if drawn automatically;}
  \item{x.ticks}{numeric; defaults to NULL; (then ticks are chosen automatically);
                 if non-NULL, user-given x-ticks (on original scale);}
  \item{y.ticks}{numeric; defaults to NULL; (then ticks are chosen automatically);
                 if non-NULL, user-given y-ticks (on original scale);}
  \item{cex.pts}{size of the points of the second argument plotted (vectorized);}
  \item{cex.pts.fun}{rescaling function for the size of the points to be plotted;
        either \code{NULL} (default), then \code{log(1+abs(x))} is used for
        the rescaling, or a function which is then used for the
        rescaling.}
  \item{col.pts}{color of the points of the second argument plotted (vectorized);}
  \item{pch.pts}{symbol of the points of the second argument plotted (vectorized);}
  \item{col.npts}{color of the non-labelled points of the \code{data} argument
     plotted (vectorized);}
  \item{pch.npts}{symbol of the non-labelled points of the \code{data} argument
         plotted (vectorized);}
  \item{cex.npts}{size of the non-labelled points of the \code{data} argument
                 plotted (vectorized);}
  \item{cex.npts.fun}{rescaling function for the size of the non-labelled points
  to be plotted; either \code{NULL} (default), then \code{log(1+abs(x))}
  is used for each of the rescalings, or a function which is then used
  for each of the rescalings.}
  \item{with.lab}{logical; shall labels be plotted to the observations?}
  \item{lab.pts}{character or NULL; labels to be plotted to the observations; if \code{NULL}
                 observation indices;}
  \item{lab.font}{font to be used for labels}
  \item{alpha.trsp}{alpha transparency to be added ex post to colors
        \code{col.pch} and \code{col.lbl}; if one-dim and NA all colors are
        left unchanged. Otherwise, with usual recycling rules \code{alpha.trsp}
        gets shorted/prolongated to length the data-symbols to be plotted.
        Coordinates of this vector \code{alpha.trsp} with NA are left unchanged,
        while for the remaining ones, the alpha channel in rgb space is set
        to the respective coordinate value of \code{alpha.trsp}. The non-NA
        entries must be integers in [0,255] (0 invisible, 255 opaque).}
  \item{jit.fac}{jittering factor used in case of a \code{DiscreteDistribution}
                    for plotting points of the second argument in a jittered fashion.}
  \item{jit.tol}{jittering tolerance used in case of a \code{DiscreteDistribution}
                    for plotting points of the second argument in a jittered fashion.}
  \item{which.lbs}{either an integer vector with the indices of the observations
          to be plotted into graph or \code{NULL} --- then no observation is excluded}
  \item{which.nonlbs}{indices of the observations which should be plotted but
        not labelled; either an integer vector with the indices of the observations
          to be plotted into graph or \code{NULL} --- then all non-labelled
          observations are plotted.}
  \item{which.Order}{we order the observations (descending) according to the norm given by
           \code{normtype(object)}; then \code{which.Order}
           either is an integer vector with the indices of the \emph{ordered}
           observations (remaining after a possible reduction by argument \code{which.lbs})
           to be plotted into graph or \code{NULL} --- then no (further) observation
           is excluded.}
  \item{attr.pre}{logical; do graphical attributes for plotted data refer
                  to indices prior (\code{TRUE}) or posterior to selection
                  via arguments \code{which.lbs}, \code{which.Order}, \code{which.nonlbs}
                  (\code{FALSE})? }
  \item{return.Order}{logical; if \code{TRUE}, an order vector
    is returned; more specifically, the order of the (remaining) observations
    given by their original index is returned (remaining means: after a possible
    reduction by argument \code{which.lbs}, and ordering is according to the norm given by
           \code{normtype(object)});
   otherwise we return \code{invisible()} as usual.}
  \item{withSubst}{logical; if \code{TRUE} (default) pattern substitution for
      titles and lables is used; otherwise no substitution is used. }
  \item{withMakeIC}{logical; if \code{TRUE} the [p]IC is passed through
        \code{makeIC} before return.}
}
\details{
  In case of \code{cniperCont} the difference between the risks of two ICs 
  is plotted.

  The function \code{cniperPoint} can be used to determine cniper
  points. That is, points such that the optimally robust estimator
  has smaller minimax risk than the classical optimal estimator under 
  contamination with Dirac measures at the cniper points. 

  As such points might be difficult to find, we provide the
  function \code{cniperPointPlot} which can be used to obtain a plot
  of the risk difference; in this function the usual arguments for
  \code{plot} can be used. For arguments \code{col}, \code{lwd},
  vectors can be used; then the first coordinate is taken for the
  curve, the second one for the balancing line. For argument \code{lty},
  a list can be used; its first component is then taken for the
  curve, the second one for the balancing line.

  If argument \code{withSubst} is \code{TRUE}, in all title 
  and axis lable arguments of \code{cniperCont} and \code{cniperPointPlot}, 
  the following patterns are substituted:
  \describe{
  \item{\code{"\%C"}}{class of argument \code{L2Fam} (for  \code{cniperPointPlot})}
  \item{\code{"\%A"}}{deparsed argument  \code{L2Fam} (for  \code{cniperPointPlot})}
  \item{\code{"\%C1"}}{class of argument \code{IC1} (for  \code{cniperCont})}
  \item{\code{"\%A1"}}{deparsed argument  \code{IC1} (for  \code{cniperCont})}
  \item{\code{"\%C2"}}{class of argument \code{IC2} (for  \code{cniperCont})}
  \item{\code{"\%A2"}}{deparsed argument  \code{IC2} (for  \code{cniperCont})}
  \item{\code{"\%D"}}{time/date-string when the plot was generated}
  }

  For more details about cniper contamination and cniper points we refer 
  to Section~3.5 of Kohl et al. (2008) as well as Ruckdeschel (2004) and 
  the Introduction of Kohl (2005).
}
\value{The cniper point is returned by \code{cniperPoint}.
In case of \code{cniperPointPlot}, we return
an S3 object of class \code{c("plotInfo","DiagnInfo")}, i.e., a list
  containing the information needed to produce the
  respective plot, which at a later stage could be used by different
  graphic engines (like, e.g. \code{ggplot}) to produce the plot
  in a different framework. A more detailed description will follow in
  a subsequent version.

}
\references{
  M. Kohl (2005). Numerical Contributions to the Asymptotic Theory of Robustness.
  Dissertation. University of Bayreuth. \url{https://epub.uni-bayreuth.de/id/eprint/839/2/DissMKohl.pdf}.

  M. Kohl, P. Ruckdeschel, and H. Rieder (2010). Infinitesimally Robust Estimation in 
  General Smoothly Parametrized Models. Statistical Methods and Applications \emph{19}(3): 333-354.
  \doi{10.1007/s10260-010-0133-0}.

  P. Ruckdeschel (2004). Higher Order Asymptotics for the MSE of M-Estimators
  on Shrinking Neighborhoods. Unpublished Manuscript.
}
\author{Matthias Kohl \email{Matthias.Kohl@stamats.de}}
%\note{}
%\seealso{ ~~objects to See Also as \code{\link{help}}, ~~~ }
\examples{
## cniper contamination
P <- PoisFamily(lambda = 4)
RobP1 <- InfRobModel(center = P, neighbor = ContNeighborhood(radius = 0.1))
IC1 <- optIC(model=RobP1, risk=asMSE())
RobP2 <- InfRobModel(center = P, neighbor = ContNeighborhood(radius = 1))
IC2 <- optIC(model=RobP2, risk=asMSE())
cniperCont(IC1 = IC1, IC2 = IC2,
           neighbor = ContNeighborhood(radius = 0.5), 
           risk = asMSE(),
           lower = 0, upper = 8, n = 101)

## cniper point plot
cniperPointPlot(P, neighbor = ContNeighborhood(radius = 0.5), 
                risk = asMSE(), lower = 0, upper = 10)

## Don't run to reduce check time on CRAN
\donttest{
## cniper point
cniperPoint(P, neighbor = ContNeighborhood(radius = 0.5), 
            risk = asMSE(), lower = 0, upper = 4)
cniperPoint(P, neighbor = ContNeighborhood(radius = 0.5), 
            risk = asMSE(), lower = 4, upper = 8)
}
}
\concept{cniper contamination}
\concept{cniper point}
\keyword{robust}
