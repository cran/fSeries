
#
# Example: 
#	European Stock Market Data
#
# Author:
#	(C) 2002, Diethelm Wuertz, GPL
#


# ------------------------------------------------------------------------------


# Load Data:
	
	data(EuStockMarkets)
	colnames(EuStockMarkets)
	class(EuStockMarkets)
	Continue = readline("Press any key > ")
		
# Time Series Attributes:
	
	start(EuStockMarkets)
	end(EuStockMarkets)
	frequency(EuStockMarkets)
	Continue = readline("Press any key > ")

	
# Log Returns:
	
	x = diff(log(EuStockMarkets[,"SMI"]))

	
# Fit AR(2)
	
	fit = armaFit(x ~ ar(2))
	par(mfrow = c(2, 2))
	summary(fit, doplot = TRUE)
	Continue = readline("Press any key > ") 
	par(mfrow = c(2, 1))
	ts.plot(x)
	predict(fit)
	Continue = readline("Press any key > ") 

	
# Fit ARMA(2,1) with ARMA Model Approach
	
	fit = armaFit(x ~ arma(2,1))
	par(mfrow = c(2, 2))
	summary(fit, doplot = TRUE)
	Continue = readline("Press any key > ") 
	par(mfrow = c(2, 1))
	predict(fit, 10)
	Continue = readline("Press any key > ") 

	
# Fit ARMA(2,1) with ARIMA Model Approach
	
	fit = armaFit(x ~ arima(2, 0, 1))
	summary(fit, doplot = FALSE)
	Continue = readline("Press any key > ") 
	predict(fit, 10)
	par(ask = FALSE)
	
	