\name{GarchDistributions}

\alias{GarchDistributions}

\alias{dsnorm}
\alias{psnorm}
\alias{qsnorm}
\alias{rsnorm}

\alias{dstd}
\alias{pstd}
\alias{qstd}
\alias{rstd}

\alias{dsstd}
\alias{psstd}
\alias{qsstd}
\alias{rsstd}

\alias{dged}
\alias{pged}
\alias{qged}
\alias{rged}

\alias{dsged}
\alias{psged}
\alias{qsged}
\alias{rsged}


\title{GARCH Distributions}


\description{
    
    A collection and description of functions 
    to compute density, distribution function, 
    quantile function and to generate random 
    variates for the skew normal, the skew 
    Student-t, and skew generalized error distribution. 
    \cr
    
    The functions are:
    
    \tabular{ll}{
	\code{[dpqr]norm} \tab Normal distribution from R's base package, \cr
	\code{[dpqr]snorm} \tab Skew Normal distribution, \cr
	\code{[dpqr]std} \tab Symmetric Student-t Distribution, \cr
	\code{[dpqr]sstd} \tab Skew Student-t Distribution, \cr
	\code{[dpqr]ged} \tab Symmetric GED distribution, \cr
	\code{[dpqr]sged} \tab Skew GED distribution. }
}


\usage{
dsnorm(x, mean = 0, sd = 1, xi = 1.5)
psnorm(q, mean = 0, sd = 1, xi = 1.5)
qsnorm(p, mean = 0, sd = 1, xi = 1.5)
rsnorm(n, mean = 0, sd = 1, xi = 1.5)

dstd(x, mean = 0, sd = 1, nu = 5)
pstd(q, mean = 0, sd = 1, nu = 5)
qstd(p, mean = 0, sd = 1, nu = 5)
rstd(n, mean = 0, sd = 1, nu = 5)

dsstd(x, mean = 0, sd = 1, nu = 5, xi = 1.5)
psstd(q, mean = 0, sd = 1, nu = 5, xi = 1.5)
qsstd(p, mean = 0, sd = 1, nu = 5, xi = 1.5)
rsstd(n, mean = 0, sd = 1, nu = 5, xi = 1.5)

dged(x, mean = 0, sd = 1, nu = 2)
pged(q, mean = 0, sd = 1, nu = 2)
qged(p, mean = 0, sd = 1, nu = 2)
rged(n, mean = 0, sd = 1, nu = 2)

dsged(x, mean = 0, sd = 1, nu = 2, xi = 1.5)
psged(q, mean = 0, sd = 1, nu = 2, xi = 1.5)
qsged(p, mean = 0, sd = 1, nu = 2, xi = 1.5)
rsged(n, mean = 0, sd = 1, nu = 2, xi = 1.5)
}


\arguments{

    \item{mean, sd, nu, xi}{
        location parameter \code{mean},
        scale parameter \code{sd},
        shape parameter \code{nu}, 
        skewness parameter \code{xi}.
        }
    \item{n}{
        number of observations.
        } 
    
    \item{p}{
        a numeric vector of probabilities.
        }
    \item{x, q}{
        a numeric vector of quantiles.
        }
    
}


\value{

    All values are numeric vectors:
    \code{d*} returns the density,
    \code{p*} returns the distribution function,
    \code{q*} returns the quantile function, and
    \code{r*} generates random deviates.
    
}

\details{
  
    \bold{Symmetric Normal Distibution:}
    \cr\cr
    The functions for the normal distribution are part of R's
    base package. The functions for the symmetric Student-t 
    distribution are rescaled in such a way that they have unit 
    variance in contrast to the Student-t family \code{dt}, \code{pt}, 
    \code{qt} and \code{rt} which are part of R's base package.
    The generalized error distribution functions are defined as 
    described by Nelson (1991).
    \cr
    
    \bold{Skew Normal Distribution:}
    \cr\cr
    The skew normal distribution functions are defined as described
    by Fernandez and Steel (2000).
    
}


\examples{
## SOURCE("fBasics.A0-SPlusCompatibility")
## SOURCE("fSeries.C1-HeavisideFunction")
## SOURCE("fSeries.C2-GarchDistributions")

## snorm -
   xmpSeries("\nStart: Skew Normal Distribuion:  > ")
   par(mfrow = c(2, 2), cex = 0.75)
   set.seed(1953)
   r = rsnorm(n = 1000, mean = 1, sd = 0.5, xi = 1.5)
   plot(r, type = "l", main = "snorm: xi = 1.5")
   # Plot empirical density and compare with true density:
   hist(r, n = 25, probability = TRUE, border = "white", col = "steelblue4")
   x = seq(-4, 6, 0.1)
   lines(x, dsnorm(x = x, mean = 1, sd = 0.5, xi = 1.5))
   # Plot df and compare with true df:
   plot(sort(r), (1:1000/1000), main = "Probability", col = "steelblue4")
   lines(x, psnorm(x, mean = 1, sd = 0.5, xi = 1.5))
   # Compute quantiles:
   qsnorm(psnorm(q = -4:6, mean = 1, sd = 0.5, xi = 1.5), 
     mean = 1, sd = 0.5, xi = 1.5) 
     
## sstd -
   xmpSeries("\nStart: Skew Student-t Distribution:  > ")
   par(mfrow = c(2, 2), cex = 0.75)
   set.seed(1953)
   r = rsstd(n = 1000, nu = 4, xi = 1.5)
   # Print Variance:
   var(r)
   plot(r, type = "l", main = "sstd: xi = 1.5")
   # Plot empirical density and compare with true density:
   hist(r, n = 30, xlim = c(-5, 5), probability = TRUE, 
     border = "white", col = "steelblue4")
   x = seq(-5, 5, 0.1)
   lines(x, dsnorm(x = x, xi = 1.5))
   # Plot df and compare with true df:
   plot(sort(r), (1:1000/1000), main = "Probability", col = "steelblue4")
   lines(x, psstd(x, xi = 1.5))
   # Compute quantiles:
   qsstd(psstd(q = -5:5, xi = 1.5), xi = 1.5) 
   
## sged -
   xmpSeries("\nStart: Skew Generalized Error Distribuion:  > ")
   par(mfrow = c(2, 2), cex = 0.75)
   set.seed(1953)
   r = rsged(n = 1000, mean = 1, sd = 0.5, xi = 1.5)
   plot(r, type = "l", main = "sged: xi = 1.5")
   # Plot empirical density and compare with true density:
   hist(r, n = 25, probability = TRUE, border = "white", col = "steelblue4")
   x = seq(-1, 5, 0.1)
   lines(x, dsged(x = x, mean = 1, sd = 0.5, xi = 1.5))
   # Plot df and compare with true df:
   plot(sort(r), (1:1000/1000), main = "Probability", col = "steelblue4")
   lines(x, psged(x, mean = 1, sd = 0.5, xi = 1.5))
   # Compute quantiles:
   qsged(psged(q = -1:5, mean = 1, sd = 0.5, xi = 1.5), 
     mean = 1, sd = 0.5, xi = 1.5) 
}


\author{

    Diethelm Wuertz for the Rmetrics \R-port.
    
}


\references{

Nelson D.B. (1991);
    \emph{Conditional Heteroscedasticity in Asset Returns: A New Approach},
    Econometrica, 59, 347--370.
    
Fernandez C., Steel M.F.J. (2000); 
    \emph{On Bayesian Modelling of Fat Tails and Skewness},
    Preprint, 31 pages. 
    
}


\seealso{
  
    \code{HeavisideFunction},
    \code{GarchDistributionFits}.
    
}


\keyword{distribution}

