
#
# Example:
#	Perform a Teraesvirta Neural Network Test and 
#	a White Neural Network Test for Nonlinearity
#
# Author:
#	(C) 2002, Diethelm Wuertz, GPL
#


# ------------------------------------------------------------------------------

	
# Teraesvirta Neural Network Test for Nonlinearity	
	
	# Time Series:
	par(mfrow = c(2,2), cex = 0.6)
	n = 1000
	x = runif(1000, -1, 1)  
	x = as.ts(x)
	plot(x)
	tnnTest(x)
	Continue = readline("Press any key > ")
	# Generate time series which is nonlinear in mean:
	x[1] = 0.0
	for (i in (2:n)) {
		x[i] = 0.4*x[i-1] + tanh(x[i-1]) + rnorm (1, sd = 0.5) }
	x = as.ts(x)
	plot(x)
	tnnTest(x)
	Continue = readline("Press any key > ")

# White Neural Network Test for Nonlinearity:

	# Time Series:
	n = 1000
	x = runif(1000, -1, 1)  
	x = as.ts(x)
	plot(x)
	wnnTest(x)
	Continue = readline("Press any key > ")	
	# Generate time series which is nonlinear in mean:
	x[1] = 0.0
	for (i in (2:n)) 
		x[i] = 0.4*x[i-1] + tanh(x[i-1]) + rnorm (1, sd = 0.5) 
	x = as.ts(x)
	plot(x)
	wnnTest(x)

