
#
# Example: 
#	Write a function which displays and tests the properties
#	of a distribution function specified by name. The components
#	which should be plotted and tested are:
#	1. Generate and plot a series of random numbers
#   2. Compute mean and variance of random numbers
#   3. Plot a histogram and compare with true density
#   4. Integrate the density
#   5. Plot empirical and true probability
#   6. Test quantile function and probability function
#


# ------------------------------------------------------------------------------


distTest = 
function(dist = "norm", x = seq(-4, 4, length = 256), n = 5000, 
breaks = 30, ...)
{	# A function implemented by Diethelm Wuertz

	# Description
	#	Test function for distributions
	
	# FUNCTION:
	
	# Random Deviates:
	rdist = match.fun(paste("r", dist, sep = "")) 
	r = rdist(n, ...)
	plot(r, type = "l", main = dist, col = "steelblue4") 
	abline(h = mean(r), lty = 3, col = "grey")
	plot(r[-1], r[-n], main = "Scatterplot", col = "steelblue4")
	abline(h = mean(r), lty = 3, col = "grey")
	abline(v = mean(r), lty = 3, col = "grey")
	Mean = mean(r)
	Variance = var(r)

	# Density:
	ddist = match.fun(paste("d", dist, sep = "")) 
	d = ddist(x, ...)
	hist(r, breaks = breaks, probability = TRUE, 
		xlim = range(x), ylim = c(0, max(d)), 
		col = "steelblue4", border = "white")
	lines(x, d, type = "l") 
	Norm = integrate(ddist, -Inf, Inf, ...)$value
	
	# Probability:
	pdist = match.fun(paste("p", dist, sep = "")) 
	p = pdist(x, ...)
	plot(x, p, type = "l", main = dist)
	points(sort(r), (1:length(r))/length(r), 
		col = "steelblue4") 
		
	# Quantiles:
	qdist = match.fun(paste("q", dist, sep = "")) 	
	Quart1 = pdist(qdist(1/4, ...), ...)
	Quart3 = pdist(qdist(3/4, ...), ...)	

	# Test Results:
	test = c(Mean = Mean, Variance = Variance, Norm = Norm, 
		Quart1 = Quart1, Quart3 = Quart3)
		
	# Return Value:
	round(test, 2)
}

