\name{ChaoticTimeSeries}

\alias{ChaoticTimeSeries}

\alias{henonSim}
\alias{ikedaSim}
\alias{logisticSim}
\alias{lorentzSim}
\alias{roesslerSim}

\title{Chaotic Time Series Modelling}


\description{

	A collection and description of functions to 
	investigate the chaotic behavior of time series 
	processes. Included are functions to simulate 
	different types of chaotic time series maps.
	\cr
	
	Chaotic Time Series Maps:
		
	\tabular{ll}{
	\code{henonSim} \tab Simulates data from theHenon Map, \cr
	\code{ikedaSim} \tab simulates data from the Ikeda Map, \cr
	\code{logisticSim} \tab simulates data from the Logistic Map, \cr
	\code{lorentzSim} \tab simulates data from the Lorentz Map, \cr
	\code{roesslerSim} \tab simulates data from the Roessler Map. }  
	
	Sorry, currently are implemented only functions to simulate
	chaotic time maps.
	
}


\usage{
henonSim(n = 1000, n.skip = 100, parms = c(a = 1.4, b = 0.3), 
	start = runif(2), doplot = FALSE)
ikedaSim(n = 1000, n.skip = 100, parms = c(a = 0.4, b = 6.0, c = 0.9), 
	start = runif(2), doplot = FALSE)
logisticSim(n = 1000, n.skip = 100, parms = c(r = 4), start = runif(1), 
	doplot = FALSE)
lorentzSim(times = seq(0, 40, by = 0.01), parms = c(sigma = 16, r = 45.92, 
	b = 4), start = c(-14, -13, 47), doplot = TRUE, \dots)
roesslerSim(times = seq(0, 100, by = 0.01), parms = c(a = 0.2, b = 0.2, c = 8.0),
	start = c(-1.894, -9.920, 0.0250), doplot = TRUE, \dots)
}


\arguments{
  
  	\item{n, n.skip}{
        [henonSim][ikedaSim][logisticSim] - \cr
        the number of chaotic time series points to be generated and the 
        number of initial values to be skipped from the series. 
        }
  	\item{parms}{
		the parameter vector characterizing the chaotic map.
		}
  	\item{start}{
		the vector of start values to initiate the chaotic map.
		}
  	\item{doplot}{
		a logical value. Should a plot be displayed? By default FALSE.
		}
	\item{times}{
        [lorentzSim][roesslerSim] - \cr
        the sequence of time series points at which to generate the map. 
    	}
	\item{\dots}{
		arguments to be passed.
		}

}


\value{

	All functions return invisible a vector of time series data. 
		
}


\references{

Brock, W.A., Dechert W.D., Sheinkman J.A. (1987); 
	\emph{A Test of Independence Based on the Correlation 
	Dimension}, 
	SSRI no. 8702, Department of Economics, University of 
	Wisconsin, Madison.
	
}


\author{

	Diethelm Wuertz for the Rmetrics \R-port.
	
}


\seealso{

	\code{RandomInnovations}.
	
}


\examples{
## SOURCE("fBasics.A0-SPlusCompatibility")
## SOURCE("fBasics.B4-TestsClass")
## SOURCE("fSeries.A4-TseriesTests")

## bdsTest -
   xmpSeries("\nNext: Brock-Dechert-Sheinkman Test for iid Series >")
   # iid Time Series:
   par(mfrow = c(3, 1))
   x = rnorm(100)
   plot(x, type = "l", main = "iid Time Series")
   bdsTest(x, m = 3)
   # Non Identically Distributed Time Series:
   x = c(rnorm(50), runif(50))
   plot(x, type = "l", main = "Non-iid Time Series")
   bdsTest(x, m = 3)  
   # Non Independent Innovations from Quadratic Map:
   x = rep(0.2, 100)
   for (i in 2:100) x[i] = 4*(1-x[i-1])*x[i-1]
   plot(x, type = "l", main = "Quadratic Map")
   bdsTest(x, m = 3)
   
## tnnTest -
   xmpSeries("\nNext: Teraesvirta NN test for Neglected Nonlinearity >")
   # Time Series Non-linear in "mean" regression 
   par(mfrow = c(2, 1))
   n = 1000
   x = runif(1000, -1, 1)  
   tnnTest(x)
   # Generate time series which is nonlinear in "mean"
   x[1] = 0.0
   for (i in (2:n)) {
     x[i] = 0.4*x[i-1] + tanh(x[i-1]) + rnorm (1, sd = 0.5) }
   plot(x, main = "Teraesvirta Test", type = "l")
   tnnTest(x)
   
## wnnTest -
   xmpSeries("\nNext: White NN test for Neglected Nonlinearity >")
   # Time Series Non-Linear in "mean" Regression
   par(mfrow = c(2, 1))
   n = 1000
   x = runif(1000, -1, 1)
   wnnTest(x)
   # Generate time series which is nonlinear in "mean"
   x[1] = 0.0
   for (i in (2:n)) {
     x[i] = 0.4*x[i-1] + tanh(x[i-1]) + rnorm (1, sd = 0.5) }
   plot(x, main = "White Test", type = "l")
   wnnTest(x)                     
}


\keyword{models}

