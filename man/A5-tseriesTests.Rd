\name{tseriesTests}

\alias{tseriesTests}

\alias{tsTest}
\alias{adfTest}
\alias{bdsTest}
\alias{jbTest}
\alias{kpssTest}
\alias{poTest}
\alias{ppTest}
\alias{tnnTest}
\alias{wnnTest}

\title{Time Series Tests}


\description{

	A collection of functions for testing various aspects of time series, 
	including unit roots, independence, normality of the distribution, 
	trend stationary, cointegration and neglected nonlinearities.
	\cr
	
	The family of time series tests includes the following 
	hypothesis tests:
		
	\tabular{rll}{
	 	1	\tab 	\code{adfTest} 	
	 		\tab  	Augmented Dickey--Fuller test for unit roots, \cr
		2	\tab 	\code{bdsTest}	
			\tab  	Brock--Dechert--Scheinkman test for iid series, \cr
		3	\tab 	\code{jbTest} 	
			\tab  	Jarque--Bera Test for normality, \cr
		4	\tab 	\code{kpssTest} 	
			\tab  	KPSS test for level or trend stationary, \cr
		5	\tab 	\code{poTest} 	
			\tab  	Phillips--Ouliaris test for cointegration, \cr
		6	\tab 	\code{ppTest} 	
			\tab  	Phillips--Perron test for unit roots, \cr	
		7	\tab 	\code{tnnTest} 	
			\tab  	Teraesvirta NN test for neglected nonlinearity, \cr
		8	\tab 	\code{wnnTest} 	
			\tab  	White NN test for neglected nonlinearity.}
	
}


\note{

	There is nothing really new in this package. 
	The benefit you will get from this help page, that we bring together
	a broad collection of time series tests from several \R-packages 
	which may be of interest for economists and financial engineers.
	
	On the other hand the user can still use the underlying function
	calls from the imported \R-packages.
	
	The output of the various hypothesis tests is an object of class
	\code{htest}. The associated \code{print} method gives an unique 
	report about the test results.
	
}


\usage{
adfTest(x, alternative = c("stationary", "explosive"), k)
bdsTest(x, m = 3, eps = seq(0.5*sd(x), 2*sd(x), length = 4), trace = FALSE)
jbTest(x)
kpssTest(x, null = c("Level", "Trend"), lshort = TRUE)
poTest(x, demean = TRUE, lshort = TRUE) 
ppTest(x, alternative = c("stationary", "explosive"), 
	type = c("Z(alpha)", "Z(t_alpha)"), lshort = TRUE)
tnnTest(x, lag = 1, type = c("Chisq", "F"), scale = TRUE, \dots) 
wnnTest(x, lag = 1, qstar = 2, q = 10, range = 4, 
	type = c("Chisq", "F"), scale = TRUE, \dots) 
}


\arguments{
  
  	\item{alternative, null}{
		In the case of \code{alternative}, the argument indicates an 
		alternative hypothesis. The choice must be one of the elements 
		specified in the input vector, by default the first entry is 
		selected.
		In the case of \code{null}, the argument indicates for the
		\code{kpssTest} the null hypothesis and must be one of 
		\code{"Level"} (default) or \code{"Trend"}. You can specify 
		just the initial letter.
		}
  	\item{demean}{
		a logical indicating whether an intercept is included in
		the cointegration regression or not.
		}
  	\item{eps}{
		a numeric vector of epsilon values for close points. The
		BDS test is computed for each element of \code{eps}. It should 
		be set in terms of the standard deviation of \code{x}.
		}
  	\item{k}{
		the lag order to calculate the test statistic.
		}
  	\item{lag}{
		an integer which specifies the model order in terms of lags.
		}
  	\item{lshort}{
		a logical indicating whether the short or long version
		of the truncation lag parameter is used.
		}
  	\item{m}{
		an integer indicating that the BDS test statistic is computed
		for embedding dimensions \code{2}, \dots, \code{m}.
		}
  	\item{q}{
		an integer representing the number of phantom hidden units
		used to compute the test statistic.
		}
  	\item{qstar}{
		the test is conducted using \code{qstar} principal
		components of the phantom hidden units. The first principal
		component is omitted since in most cases it appears to be collinear
		with the input vector of lagged variables. This strategy preserves
		power while still conserving degrees of freedom.
		}
  	\item{range}{
		the input to hidden unit weights are initialized uniformly
		over \code{[-range/2, range/2]}.
		}
  	\item{scale}{
		a logical indicating whether the data should be scaled 
		before computing the test statistic. The default arguments to
		\code{\link{scale}} are used.
		}
  	\item{trace}{
		a logical indicating whether some informational output is
		traced.
		}
  	\item{type}{
    	indicates which variant of a test will be computed. The choice 
    	must be one of the elements specified in the input vector, by 
    	default the first entry is selected.
    	} 
    \item{x}{
		a numeric vector or time series.
		}	
    \item{\dots}{
		further arguments to be passed from or to methods.
		}	

}


\details{

	\bold{Augmented Dickey--Fuller Test for Unit Roots:}
	\cr\cr
	The general regression equation which incorporates a constant and a
	linear trend is used and the t-statistic for a first order
	autoregressive coefficient equals one is computed. The number of lags
	used in the regression is \code{k}. The default value of
	\code{trunc((length(x)-1)^(1/3))} corresponds to the suggested upper
	bound on the rate at which the number of lags, \code{k}, should be
	made to grow with the sample size for the general \code{ARMA(p,q)}
	setup. Note that for \code{k} equals zero the standard Dickey-Fuller
	test is computed. The p-values are interpolated from Table 4.2, p. 103
	of Banerjee et al. (1993). If the computed statistic is outside the
	table of critical values, then a warning message is generated. 
	Note, that missing values are not allowed. 
	\cr
	\code{[tseries:adf.test]}
	\cr
	
	\bold{Brock--Dechert--Sheinkman Test:}
	\cr\cr
	The \code{bdsTest} test examines the \emph{spatial dependence} of the 
	observed series. To do this, the series is embedded in \code{m}-space 
	and the dependence of \code{x} is examined by counting \emph{near} 
	points. Points for which the distance is less than \code{eps} are 
	called near. The BDS test statistic is asymptotically standard Normal.
  	Note, that missing values are not allowed.
	There is a special print method for objects of class \code{"bdsTest"}
  	which by default uses 4 digits to format real numbers. 
  	\cr
  	\code{[tseries:bds.test]}
  	\cr
	
	\bold{Jarque--Bera Test:}
	\cr\cr
	The \code{jbTest} evaluates the hypothesis that the series has a 
	normal distribution with unspecified mean and variance, against the 
	alternative that the series does not have a normal distribution. The 
	test is based on the sample skewness and kurtosis. For a true normal 
	distribution, the sample skewness should be near 0 and the sample 
	kurtosis should be near 3. The Jarque-Bera test determines whether 
	the sample skewness and kurtosis are unusually different than their 
	expected values, as measured by a chi-square statistic. The Jarque-Bera 
	test is an asymptotic test, and should not be used with small samples. 
	You may want to use the \code{lillieTest} for small samples. 
	\cr
	\code{[jb.test]}
	\cr
	
	\bold{Kwiatkowski--Phillips--Schmidt--Shin Test:}
	\cr\cr
	To estimate \code{sigma^2} the Newey-West estimator is used.
  	If \code{lshort} is \code{TRUE}, then the truncation lag parameter is
  	set to \code{trunc(3*sqrt(n)/13)}, otherwise
  	\code{trunc(10*sqrt(n)/14)} is used.  The p-values are interpolated
  	from Table 1 of Kwiatkowski et al. (1992). If the computed statistic
  	is outside the table of critical values, then a warning message is
  	generated. 
  	\cr
  	\code{[tseries:kpss.test]}
  	\cr
	
	\bold{Phillips--Ouliaris Test:}
	\cr\cr
	The \code{poTest} Z(alpha) statistic for a unit root in the
  	residuals of the cointegration regression is computed, see also
  	\code{ppTest}. The unit root is estimated from a regression of
  	the first variable (column) of \code{x} on the remaining variables of
  	\code{x} without a constant and a linear trend.  To estimate
  	\code{sigma^2} the Newey-West estimator is used.  If \code{lshort} is
  	\code{TRUE}, then the truncation lag parameter is set to
  	\code{trunc(n/100)}, otherwise \code{trunc(n/30)} is used.  The
  	p-values are interpolated from Table Ia and Ib, p. 189 of Phillips and
  	Ouliaris (1990).  If the computed statistic is outside the table of
  	critical values, then a warning message is generated. The dimension 
  	of \code{x} is restricted to six variables. 
  	\cr
  	\code{[tseries:po.test]}
  	\cr  
  	
	\bold{Phillips--Perron Test:}
	\cr\cr	
	The general regression equation which incorporates a constant and a
  	linear trend is used and the \code{Z(alpha)} or \code{Z(t_alpha)}
  	statistic for a first order autoregressive coefficient equals one are
  	computed.  To estimate \code{sigma^2} the Newey-West estimator is
  	used.  If \code{lshort} is \code{TRUE}, then the truncation lag
  	parameter is set to \code{trunc(4*(n/100)^0.25)}, otherwise
  	\code{trunc(12*(n/100)^0.25)} is used.  The p-values are interpolated
  	from Table 4.1 and 4.2, p. 103 of Banerjee et al. (1993).  If the
  	computed statistic is outside the table of critical values, then a
  	warning message is generated. 
  	\cr
  	\code{[tseries:pp.test]}
  	\cr
  
  	\bold{Teraesvirta Neural Network Test:}
  	\cr\cr	
  	The null is the hypotheses of linearity in \code{mean}. This test 
  	uses a Taylor series expansion of the activation function to 
  	arrive at a suitable test statistic. If \code{type} equals \code{"F"}, 
  	then the F-statistic instead of the Chi-Squared statistic
	is used in analogy to the classical linear regression. 
  	Missing values are not allowed. 
  	\cr 
  	\code{[tseries:teraesvirta.test]}
  	\cr
	
	\bold{White Neural Network Test:}
	\cr\cr
	The null is the hypotheses of linearity in ``mean''. This
	type of test is consistent against arbitrary nonlinearity
	in mean. If \code{type} equals \code{"F"}, then the F-statistic
	instead of the Chi-Squared statistic is used in analogy to the
	classical linear regression. 
	\cr
	\code{[tseries:white.test]}
	
}


\value{

	A list with class \code{"htest"} containing the following 
	components:

	\item{statistic}{
		the value of the test statistic.
		}
	\item{parameter}{
		the lag order.
		}
	\item{p.value}{
		the p-value of the test.
		}
	\item{method}{
		a character string indicating what type of test was
		performed.
		} 
	\item{data.name}{
		a character string giving the name of the data.
		}
	\item{alternative}{
		a character string describing the alternative
		hypothesis.
		}
		
}


\references{
Banerjee A., Dolado J.J., Galbraith J.W., Hendry D.F. (1993);
	\emph{Cointegration, Error Correction, and the Econometric 
	Analysis of Non-Stationary Data},
	Oxford University Press, Oxford.

Brock, W.A., Dechert W.D., Sheinkman J.A. (1987); 
	\emph{A Test of Independence Based on the Correlation 
	Dimension}, 
	SSRI no. 8702, Department of Economics, University of 
	Wisconsin, Madison.

Conover, W.J. (1980);
 	\emph{Practical Nonparametric Statistics}, 
 	New York, Wiley. 

Cromwell J.B., Labys W.C., Terraza M. (1994);
	\emph{Univariate Tests for Time Series Models},
	Sage, Thousand Oaks, CA, pages 32--36.
	
Kwiatkowski D., Phillips P.C.B, Schmidt P., Shin Y. (1992);
	\emph{Testing the Null Hypothesis of Stationarity against 
	the Alternative of a Unit Root},
	Journal of Econometrics 54, 159--178.
		
Lee T.H., White H., Granger C.W.J. (1993); 
	\emph{Testing for neglected nonlinearity in time series models},
	Journal of Econometrics 56, 269--290.
	
Perron P. (1988);
	\emph{Trends and Random Walks in Macroeconomic Time Series},
	Journal of Economic Dynamics and Control 12, 297--332.

Phillips P.C.B., Ouliaris S. (1990);
	\emph{Asymptotic Properties of Residual Based Tests for Cointegration},
	Econometrica 58, 165--193.
  	
Said S.E., Dickey D.A. (1984);
	\emph{Testing for Unit Roots in Autoregressive-Moving Average 
	Models of Unknown Order},
	Biometrika 71, 599--607.

Teraesvirta T., Lin C.F., Granger C.W.J. (1993);
	\emph{Power of the Neural Network Linearity Test},
	Journal of Time Series Analysis 14, 209--220.
	
}


\author{

	Adrian Trapletti for the tests from \R's tseries package,\cr
	Blake LeBaron for the \code{bds} C program,\cr 
	Diethelm Wuertz for this \R-port.
}


\seealso{

	\code{classicalTests},
	\code{\link{lmTests}}.
	
}


\examples{
## adfTest - 
   xmpSeries("\nStart: Augmented Dickey-Fuller Test for Unit Roots >")
   # A time series which contains no unit-root:
   x = rnorm(1000)  
   adfTest(x)
   # A time series which contains a unit-root:
   y = diffinv(x)
   adfTest(y)
   
## bdsTest -
   xmpSeries("\nNext: Brock-Dechert-Sheinkman Test for iid Series >")
   # iid Time Series:
   par(mfrow = c(3,1))
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
   
## jbTest -
   xmpSeries("\nNext: Jarque-Bera Test for Normality >")
   # Null:
   x = rnorm(100)
   jbTest(x)
   # Alternative:
   x = runif(100)
   jbTest(x)
   
## kpssTest -
   xmpSeries("\nNext: KPSS test for Level or Trend Stationary >")
   # Time Series is Level Stationary:
   x = rnorm(1000)
   kpssTest(x)
   # Time Series has Unit Root:
   y = cumsum(x)
   kpssTest(y)
   # Time Series is Trend Stationary:
   x = 0.3*(1:1000)+rnorm(1000)
   kpssTest(x, null = "Trend")
   
## poTest -
   xmpSeries("\nNext: Phillips-Ouliaris Test for Cointegration >")
   # Non-Cointegrated Case:
   x = ts(diffinv(matrix(rnorm(2000), 1000, 2)))
   poTest(x)
   # Cointegrated Case:
   x = diffinv(rnorm(1000))
   y = 2.0 - 3.0*x + rnorm(x, sd=5)
   z = ts(cbind(x, y))
   poTest(z)
   
## ppTest -
   xmpSeries("\nNext: Phillips-Perron Test for Unit Roots >")
   # The Time Series has no Unit Root:
   x = rnorm(1000)
   ppTest(x)
   # The Time Series has Unit Root:
   y = cumsum(x)  
   ppTest(y)
   
## tnnTest -
   xmpSeries("\nNext: Teraesvirta NN test for Neglected Nonlinearity >")
   # Time Series Non-linear in ``mean'' regression 
   par(mfrow = c(2, 1))
   n = 1000
   x = runif(1000, -1, 1)  
   tnnTest(x)
   # Generate time series which is nonlinear in ``mean''
   x[1] = 0.0
   for (i in (2:n)) {
     x[i] = 0.4*x[i-1] + tanh(x[i-1]) + rnorm (1, sd = 0.5) }
   plot(x, main = "Teraesvirta Test", type = "l")
   tnnTest(x)
   
## wnnTest -
   xmpSeries("\nNext: White NN test for Neglected Nonlinearity >")
   # Time Series Non-Linear in ``mean'' Regression
   par(mfrow = c(2, 1))
   n = 1000
   x = runif(1000, -1, 1)
   wnnTest(x)
   # Generate time series which is nonlinear in ``mean''
   x[1] = 0.0
   for (i in (2:n)) {
     x[i] = 0.4*x[i-1] + tanh(x[i-1]) + rnorm (1, sd = 0.5) }
   plot(x, main = "White Test", type = "l")
   wnnTest(x)                     
}


\keyword{htest}

