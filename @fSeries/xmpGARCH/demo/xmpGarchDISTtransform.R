
#
# Title: 
#	Shift and Scale a Distribution
#   
# Description: 
#   This example shows how to scale and shift the
#	density function, the distribution function, the
#	quantile function and random deviates from the
#	view of the normal distribution. The approach 
#	can easily be extended to other distributions.	


# ------------------------------------------------------------------------------

	
# What is the scaled/shifted density?

	mean = 1; sd = 2
	Location = -4:4
	Density1 = dnorm(x = Location, mean = mean, sd = sd)	
	Density2 = dnorm(x = (Location-mean)/sd, mean = 0, sd = 1) / sd	
	cbind(Location, Density1, Density2)
	
# What is the scaled/shifted probability?

	mean = 1; sd = 2
	Quantiles = -4:4
	Probability1 = pnorm(q = Quantiles, mean = mean, sd = sd)
	Probability2 = pnorm(q = (Quantiles-mean)/sd, mean = 0, sd = 1)
	cbind(Quantiles, Probability1, Probability2)
	
# What are the scaled/shifted quantiles?
	
	mean = 1; sd = 2
	Probability = (1:9)*0.1
	Quantiles1 = qnorm(p = Probability, mean = mean, sd = sd)
	Quantiles2 = qnorm(p = Probability, mean = 0, sd = 1) * sd + mean
	cbind(Probability, Quantiles1, Quantiles2)
	
# What are the scaled/shifted random deviates?	

	mean = 1; sd = 2
	N = 10
	set.seed(4711)
	Random1 = rnorm(n = N, mean = mean, sd = sd)
	set.seed(4711)
	Random2 = rnorm(n = N, mean = 0, sd = 1) * sd + mean
	cbind(1:10, Random1, Random2)
	
