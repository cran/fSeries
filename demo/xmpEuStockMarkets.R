
#
# Example: 
#	European Stock Market Data
#
# Description:
#	The examples show the analysis of the European Stock Market Data
#	using ARMA models.
#
# Author:
#	(C) 2002, Diethelm Wuertz, GPL
#


################################################################################
# Data:
	
	# Load Data:
	data(EuStockMarkets)
	colnames(EuStockMarkets)
	class(EuStockMarkets)
	###	
	
	
	# Show Time Series Attributes:
	start(EuStockMarkets)
	end(EuStockMarkets)
	frequency(EuStockMarkets)
	###

	
	# Compute Log Returns:
	x = diff(log(EuStockMarkets[,"SMI"]))
	###
	
	
################################################################################
# Fit AR(2)
	
	# Fit:
	fit = armaFit(x ~ ar(2))
	par(mfrow = c(2, 2))
	summary(fit, doplot = TRUE)
	par(mfrow = c(2, 1))
	ts.plot(x)
	predict(fit)
	### 


################################################################################	
# Fit ARMA(2,1) 
	
	# Fit:
	fit = armaFit(x ~ arma(2,1))
	par(mfrow = c(2, 2))
	summary(fit, doplot = TRUE)
	par(mfrow = c(2, 1))
	predict(fit, 10)
	###

	
################################################################################
# Fit ARMA(2,1) 
	
	# Fit
	fit = armaFit(x ~ arima(2, 0, 1))
	summary(fit, doplot = FALSE)
	predict(fit, 10)
	par(ask = FALSE)
	###
	
	
################################################################################

