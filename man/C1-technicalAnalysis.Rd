\name{technicalAnalysis}

\alias{technicalAnalysis}

\alias{emaTA}
\alias{biasTA}
\alias{medpriceTA}
\alias{typicalpriceTA}
\alias{wcloseTA}

\alias{rocTA}
\alias{oscTA}
\alias{momTA}
\alias{macdTA}
\alias{cdsTA}

\alias{cdoTA}
\alias{vohlTA}
\alias{vorTA}
\alias{fpkTA}
\alias{fpdTA}

\alias{spdTA}
\alias{apdTA}
\alias{wprTA}
\alias{rsiTA}


\title{Tools for the Technical Analysis}


\description{

	A collection of functions for the technical analysis of stock markets. 
	This collection provides a set of the most common technical indicators 
	together with functions for open--high--close--low charts, and three 
	functions for benchmark measurements.
	\cr
	
	The functions are:
	
	\tabular{rll}{
	 	1	\tab 	\code{emaTA}			
	 		\tab 	Exponential Moving Average, \cr
	 	2	\tab	\code{biasTA}			
	 		\tab 	Bias Indicator, \cr
	 	3	\tab	\code{medpriceTA}		
	 		\tab 	Medium Price Indicator, \cr
	 	4	\tab	\code{typicalpriceTA} 	
	 		\tab 	Typical Price Indicator, \cr
	 	5	\tab	\code{wcloseTA}		
	 		\tab 	Weighted Close Indicator, \cr
	 	6	\tab	\code{rocTA}		
	 		\tab 	Rate of Change, \cr
	 	7	\tab	\code{oscTA}			
	 		\tab 	Oscillator Indicator, \cr
	 	8	\tab	\code{momTA}			
	 		\tab 	Momentum Indicator, \cr
	 	9	\tab	\code{macdTA}			
	 		\tab 	MACD Indicator, \cr
	   10	\tab	\code{cdsTA}		
	   		\tab 	MACD Signal Line, \cr
	   11	\tab	\code{cdoTA}			
	   		\tab 	MACD Oscillator, \cr
	   12	\tab	\code{vohlTA}			
	   		\tab 	High/Low Volatility, \cr
	   13	\tab	\code{vorTA}			
	   		\tab 	Volatility Ratio, \cr
	   14	\tab	\code{fpkTA}		
	   		\tab 	Fast Percent K, \cr
	   15	\tab	\code{fpdTA}			
	   		\tab 	Fast Percent D, \cr
	   16	\tab	\code{spdTA}		
	   		\tab 	Slow Percent D, \cr
	   17	\tab	\code{apdTA}		
	   		\tab 	Averaged Percent D, \cr
	   18	\tab	\code{wprTA}		
	   		\tab 	William's Percent R, \cr
	   19	\tab	\code{rsiTA}		
	   		\tab 	Relative Strength Index. }
	
}


\usage{
emaTA(x, lambda, startup = 0)
biasTA(x, lag)
rocTA(x, lag)
oscTA(x, lag1, lag2)
momTA(x, lag)
macdTA(x, lag1, lag2)
cdsTA(x, lag1, lag2, lag3)
cdoTA(x, lag1, lag2, lag3)
vohlTA(high, low)
vorTA(high, low)
fpkTA(close, high, low, lag)
fpdTA(close, high, low, lag1, lag2)
spdTA(close, high, low, lag1, lag2, lag3)
apdTA(close, high, low, lag1, lag2, lag3, lag4)
wprTA(close, high, low, lag)
rsiTA(close, lag)

medpriceTA(high, low)
typicalpriceTA(high, low, close)
wcloseTA(high, low, close)
}


\arguments{

  	\item{lag, lag1, lag2, lag3, lag4}{
		integer values, time lags.
		}
  	\item{lambda}{
		a numeric value, the decay length of the exponential
		moving average.
		}
	 	\item{startup}{
		an integer value, the startup position of the exponential
		moving average, by default 0.
		}
	\item{x, high, low, close}{
		a numeric vector of prices, either opening, closing, or
		high and low values.
		For \code{ohlcPlot} a multivariate time series object of 
		class \code{mts}.
		}

}


\value{

	\code{*TA}
	\cr	
	The technical Indicators return the following numeric vectors:

	\code{emaTA}  returns the Exponential Moving Average, EMA \cr
	\code{biasTA} returns the EMA-Bias, \cr
	\code{rocTA}  returns the Rate of Change Indicator, \cr
	\code{oscTA}  returns the EMA Oscillator Indicator, \cr
	\code{momTA}  returns the Momentum Oscillator, \cr
	\code{macdTA} returns the MACD Oscillator, \cr
	\code{cdsTA}  returns the MACD Signal Line, \cr
	\code{cdo}    returns the MACD Oscillator, \cr
	\code{vohlTA} returns the High/Low Volatility Oscillator, \cr
	\code{vorTA}  returns Volatility Ratio Oscillator, \cr
	\code{fpkTA}  returns the Fast Percent-K Stochastics Indicator, \cr
	\code{fpdTA}  returns the Fast Percent-D Stochastics Indicator, \cr
	\code{spdTA}  returns the Slow Percent-D Stochastics Indicator, \cr
	\code{apdTA}  returns the Averaged Percent-D Stochastics Indicator, \cr
	\code{wprTA}  returns the Williams Percent-R Stochastics Indicator, \cr
	\code{rsiTA}  returns the Relative Strength Index Stochastics Indicator, \cr
	\code{medpriceTA} returns the Medium Price, \cr
	\code{typicalpriceTA} returns the Typical Price, \cr
	\code{wcloseTA} returns the Weighted Closing Price. \cr

}


\author{

	Diethelm Wuertz for this \R-port.
}


\examples{
## Currently no examples ...
}
   

\keyword{math}

