
#
# Example: 
#	Compute the absolute moments of a density function. 
#

	
# ------------------------------------------------------------------------------

	
# Try GED:

	absM = absMoments(1:5, density = "dged", nu = 1)

# Integrate and compare:

	moments = function(x, n, nu) { abs(x)^n * dged(x = x, nu = nu) }
	intM = NULL
	for (i in 1:5)
		intM = c(intM, integrate(moments, -Inf, +Inf, nu = 1, n = i)$value)
	
# Try Standardized Laplace Distribution:

	dlaplace = function(x) { exp(-sqrt(2)*abs(x)) / sqrt(2) }
	laplaceM = absMoments(1:5, density = "dlaplace")
	absMoments.error
	
# Print
	data.frame(absM, intM, laplaceM)
	
	