\name{fSeriesData}

\alias{fSeriesData}

\alias{nyseres}
\alias{recession}
\alias{spc1970}
\alias{spcindis}

\alias{dem2gbp}
\alias{cac40}


\title{fSeries Data Sets} 


\description{

    A collection of Data Sets used in the examples of the \code{fSeries}
	library.
	\cr
  
  	\code{recession}
  	\cr
  	Data sets used in the regression analysis.
  	Recession data from the US,
  	3 Month Tbills data from US FED,
  	10 Year Tbonds data from US FED,
  	Stock-Watson experimental recession index.
  	\cr	
  	
  	\code{spc1970} and 
  	\code{spcindis}
  	\cr
  	These data are indicators for daily SP500 Stock Data together
  	with an set of trading indicators for technical analysis.
  	\cr
  	
  	\code{dem2gbp}
  	\cr
  	The file "dem2gbp" contains daily observations of the 
	Deutschmark / British Pound foreign exchange log returns. 
	This data set has been promoted as an informal benchmark 
	for GARCH time-series software validation. See McCullough and 
	Renfro [1991], and Brooks, Burke, and Persand (2001) for details.
	The nominal returns are expressed in percent, as published in 
	Bollerslev and Ghysels (2001). The data set is available from 
	the \emph{Journal of Business and Economic Statistics}, (JBES), 
	\emph{ftp://www.amstat.org}. A text file has one column of 
	data listing the percentual log-returns of the DEM/GBP exchange 
	rates. The sample period is from January 3, 1984, to December 
	31, 1991, for a total of 1975 daily observations of FX exchange 
	rates.
	
	\code{cac40}
	\cr
	Daily CAC40 index returns with realized volatility. The data 
	cover the period January 1995 until December 1999, and have 1249 
	observations. The first column of the file \code{CAC40} lists 
	the averaged return of the Index, and the second column 
	\code{CAC40VOL} lists the realized volatility. The data are
	used a s a benchmark for GARCH modeling by Laurent and Peters
	(2002).
 
}


\format{
	
    All files are in CSV Excel spreadsheet format. 
    
}


\references{

Brooks C., Burke S.P, Persand G. (2001);
	\emph{Benchmarks and the Accuracy of GARCH Model Estimation},
	International Journal of Forecasting 17, 45--56.

McCullough B.D., Renfro C.G. (1998);
	\emph{Benchmarks and Software Standards: A Case Study of GARCH 
	Procedures},
	Journal of Economic and Social Measurement 25, 59--71. 
	
Laurent S., Peters J.P. (2002), 
	\emph{G@RCH 2.2: An Ox Package for Estimating and Forecasting  
	Various ARCH Models}, 
	Journal of Economic Surveys 16, 447--485. 

}


\examples{
## NYSE Residuals:
   data(nyseres)
   ts.plot(nyseres, xlab = "Index", ylab = "log-Returns",
   	 main = "NYSE: log-Returns")
} 
 
 
\keyword{datasets}

