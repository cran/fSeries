
#
# Title: 
#	Parameter Fit of a Density Function
#   
# Description: 
#	Fits the parameters of a distribution function by the
#	maximum log-likelihood approach using the "nlm"
#	optimizer from R's base package.

	
# ------------------------------------------------------------------------------


# Simulate NORM Series:

	x = rnorm(1000)
	
# Estimate the Parameters:

	fit = dFit(x, density = dnorm, parm = c(mean = mean(x), sd = mad(x)), 
		trace = TRUE, doplot = TRUE, labels = TRUE, title = NA, 
		description = NA) 
		
# Print the Results:
		
	print(fit)
	
	
# ------------------------------------------------------------------------------


# Simulate SGED Series:

	x = rsged(n = 1000, mean = 1, sd = 1/2, nu = 1.5, xi = 1.5)
	
# Estimate the Parameters:

	fit = dFit(x, density = dsged, parm = c(mean = 0, sd = 1,
		nu = 1, xi = 0), trace = TRUE, doplot = TRUE, labels = TRUE, 
		title = NA, description = NA) 
		
# Print the Results:
		
	print(fit)
	
	
# ------------------------------------------------------------------------------

