\name{randomInnovations}

\alias{randomInnovations}

\alias{set.lcgseed}
\alias{get.lcgseed}
\alias{runif.lcg}
\alias{rnorm.lcg}
\alias{rt.lcg}

\title{Generator for random Innovations}


\description{

	A collection of functions to generate portable random innovations.
	The function runs under R and SPlus and generates the same
	sequence of random numbers. Supported are uniform, normal and 
	Student-t distributed random numbers.
		
	The functions are:
	
	\tabular{rll}{
		1  	\tab 	\code{set.lcgseed}	 	
			\tab	Set initial random seed, \cr
		2  	\tab 	\code{get.lcgseed}		
			\tab	Get the current valus of the random seed, \cr
		3  	\tab 	\code{runif.lcg}		
			\tab	Uniform linear congruational generator, \cr
		4  	\tab 	\code{rnorm.lcg}		
			\tab	Normal linear congruational generator, \cr
		5  	\tab 	\code{rt.lcg}			
			\tab	Student-t linear congruential generator. }
	
}


\usage{
set.lcgseed(seed = 4711)
get.lcgseed()
runif.lcg(n, min = 0, max = 1)
rnorm.lcg(n, mean = 0, sd = 1)
rt.lcg(n, df)	
}


\arguments{

  	\item{df}{
    	number of degrees of freedom, a positive integer, maybe 
    	non-integer.
    	}
    \item{mean, sd}{
  		means and standard deviation of the normal distributed 
  		innovations.
  		}
  	\item{min, max}{
    	lower and upper limits of the uniform distributed 
    	innovations.
    	}
    \item{seed}{
  		an integer value, the random number seed.
  		}
  	\item{n}{
  		an integer, the number of random innovations to be 
  		generated.
  		}
  			
}


\details{
	
	A simple portable random number generator for use in R and 
	SPlus. We recommend to use this generator only for comparisons 
	of calculations in R and Splus.

	The generator is a linear congruential generator with parameters
	\code{LCG(a=13445, c=0, m=2^31-1, X=0)}. It is a simple random 
	number generator which passes the bitwise randomness test.

}


\value{

  	A vector of generated random innovations. The value of the
  	current seed is stored in the variable \code{lcg.seed}.
  	
}


\references{

Altman, N.S. (1988);
	\emph{Bitwise Behavior of Random Number Generators},
 	 SIAM J. Sci. Stat. Comput., 9(5), September, 941--949.

}


\author{

	Diethelm Wuertz for this \R-port.
	
}


\seealso{

	\code{\link{runif}},
	\code{\link{rnorm}},
	\code{\link{rt}}.
	
}


\examples{
## set.lcgseed -
   xmpSeries("\nStart: Set Initial Seed >")
   set.lcgseed(seed = 65890)
 
## runif.lcg - rnorm.lcg - rt.lcg -
   xmpSeries("\nNext: Create Random Numbers >")
   cbind(runif.lcg(10), rnorm.lcg(10), rt.lcg(10, df = 4))

## get.lcgseed -
   xmpSeries("\nNext: What is the current value of the seed? >")
   get.lcgseed()  
   
## Note, to overwrite rnorm, use
   # rnorm = rnorm.lcg
   # Going back to rnorm
   # rm(rnorm)
}


\keyword{programming}

