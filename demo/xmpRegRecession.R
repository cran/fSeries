
#
# Example: 
#	Model US Recession by several regression models
#
# Author:
#	(C) 2002, Diethelm Wuertz, GPL
#


# ------------------------------------------------------------------------------


# Read the Dataset from File:
	
	# Column 1: Date - CCYYMM
	# Column 2: Recession 0 | 1
	# Column 3: 3 Month Treasury Bills
	# Column 4: 10 Years Treasury Bonds
	data(recession)
		
# Build the Data Set for 6 Months Regression Predictions:	
	
	horizon = 6
	n = length(recession[, 1])
	# Response and Predictors:
	response = recession[,2][(1+horizon):n]
	predictor1 = recession[,3][1:(n-horizon)]
	predictor2 = recession[,4][1:(n-horizon)]
	StockWatson = recession[,5][(1+horizon):n]
	# Date/Time:
	ccyy = floor(recession[,1]/100)
	mm = recession[, 1] - 100 * ccyy
	time = ccyy + (mm-0.5)/12  # mid-month
	time = time[1:(n-horizon)]
	# Data Frame:
	data = data.frame(cbind(response, predictor1, predictor2))	

# Fit Regression Models:

	family = binomial(link = probit)
	glm.fit = regFit(response ~ predictor1 + predictor2, 
		family = family, data = data, method = "GLM")
		print(glm.fit)
	gam.fit = regFit(response ~ predictor1 + predictor2, 
		family = family, data = data, method = "GAM")
		print(gam.fit)
	ppr.fit = regFit(response ~ predictor1 + predictor2, 
		data = data, method = "PPR", nterms = 2)
		print(ppr.fit)
	mars.fit = regFit(response ~ predictor1 + predictor2, 
		data = data, method = "MARS")
		print(mars.fit)
	polymars.fit = regFit(response ~ predictor1 + predictor2, 
		data = data, method = "POLYMARS")
		print(polymars.fit)
	nnet.fit = regFit(response ~ predictor1 + predictor2, data = data, 
		method = "NNET", linout = TRUE, trace = FALSE, maxit = 1000, 
		size = 6)
		print(nnet.fit)
		
# Write Plot Function:
	myPlot = function(time, response, in.sample, StockWatson, title) {
		plot(time, response, type = "n", main = title, col = "steelblue3")
		grid()
		lines(time, response, type = "h", col = "steelblue3")
		# The Benchmark - Compare with StockWatson:
		lines(time, StockWatson, col = "red", lty = 3)
		# Set Negative Values to Zero:
		in.sample = (in.sample+abs(in.sample))/2
		# Set Values Larger than One to One:
		in.sample = -(1-in.sample+abs(1-in.sample))/2 + 1
		lines(time, in.sample) 
		}
		
# Compare with Fitted Values:
	par(mfrow = c(3, 2), cex = 0.6)
	myPlot(time, response, 
		glm.fit@fit$fitted.values, StockWatson, "GLM")
	myPlot(time, response, 
		gam.fit@fit$fitted.values, StockWatson, "GAM")
	myPlot(time, response, 
		ppr.fit@fit$fitted.values, StockWatson, "PPR")
	myPlot(time, response, 
		mars.fit@fit$fitted.values, StockWatson, "MARS")
	myPlot(time, response, 
		polymars.fit@fit$fitted.values, StockWatson, "POLYMARS")
	myPlot(time, response, 
		nnet.fit@fit$fitted.values, StockWatson, "NNET")

