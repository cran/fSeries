\name{LongMemoryModelling}

\alias{LongMemoryModelling}

\alias{fgnSim}

\title{Long Memory Behaviour of Time Series}


\description{

	A collection and description of functions to investigate the 
	long memory behavior of an univariate time series process.
	\cr
	
	The functions and methods are:
	
	\tabular{ll}{
	\code{fgnSim} \tab Simulates Fractional Gaussian Noise. }
	
	Sorry, not yet implemented are the functions to fit the Hurst
	exponent.
	
}


\usage{
fgnSim(n = 1000, H = 0.7, method = c("beran", "durbin", "paxson"), 
mean = 0, std = 1)
}


\arguments{
	
	\item{H}{
		the Hurst exponent, a numeric value between 0.5 and 1,
		by default 0.7.
		}
	\item{mean, std}{
		mean and standard deviation of the random innovations.
		By default a zero mean and a unit standard deviation is
		assumed.
		}
	\item{method}{
		the method how to generate the time series sequence, one
		of the following character strings: \code{"beran"}, 
		\code{"durbin"}, or \code{"paxson"}.
		}
	\item{n}{
		number of data points to be simulated, a numeric value, 
		by default 1000.
		}
		
}


\value{

	returns a numeric vector of length \code{n}, the fractional 
	Gaussian noise series.
	
}


\authors{

	From Vern Paxson, code as listed in the Appendix of his paper 1995,\cr
	From Beran, code as listed in the Appendix of his Book, \cr
	Diethelm Wuertz for the Rmetrics \R-port.
	
}


\examples{
## fgnSim -
   par(mfrow = c(3, 1), cex = 0.75)  
   # Beran's Method:
   plot(fgnSim(n = 200, H = 0.75), type = "l",  
   	 ylim = c(-3, 3), xlab = "time", ylab = "x(t)", main = "Beran")
   # Durbin's Method:
   plot(fgnSim(n = 200, H = 0.75, method = "durbin"), type = "l",
   	 ylim = c(-3, 3), xlab = "time", ylab = "x(t)", main = "Durbin")
   # Paxson's Method:
   plot(fgnSim(n = 200, H = 0.75, method = "paxson"), type = "l",
     ylim = c(-3, 3), xlab = "time", ylab = "x(t)", main = "Paxson")
}


\references{

Paxson V. (1995); 
	\emph{Fast Approximation of Self-Similar Network Traffic},
	Berkeley. 
	
}


\keyword{models}

