
#
# Example: 
#	Plot the skew GED Distribution in the limit of the 
#	- Laplace Distribution, mu=1, in the limit of the 
#   - Gaussian Distribution, mu=2, and in the limit of the
#   - Uniform distribution, mu->Inf.
#   

	
# ------------------------------------------------------------------------------


# Investigate the Special Cases:
    
	# nu = 1   : Laplace Distribution
    # nu = 2   : Gaussian Distribution
    # nu = Inf : Uniform Distribution

	# Settings:
	mean = 0; sd = 1; dx = 0.001
	x = seq(-5+mean, 5+mean, by = dx)


# ------------------------------------------------------------------------------
	

# Plot Skew Normal:

	# Density Function:
	par(mfrow = c(3, 2), cex = 0.5)
	d1 = dsnorm(x, mean = mean, sd = sd, xi = 2/3)
	d2 = dsnorm(x, mean = mean, sd = sd, xi = 1)
	plot(x, d1, type = "n", xlim = c(-5, 5), ylim = c(0, 0.5), 
		xlab = "z", ylab = "f(z)", main = "Skew GED: Density")
	lines(x, d1, col = "blue")
	lines(x, d2, col = "red")
	text(3.6, 0.15, "xi=2/3", col = "blue")
	text(2.8, 0.25, "xi=1", col = "red")
	
	# Distribution Function:
	p1 = psged(x, mean = mean, sd = sd, xi = 2/3)
	p2 = psged(x, mean = mean, sd = sd, xi = 1)
	plot(x, p1, type = "n", xlim = c(-5, 5),
		xlab = "z", ylab = "F(z)", main = "Skew GED: Distribution")
	lines(x, p1, col = "blue")
	lines(x, p2, col = "red")
	text(1.8, 0.20, "nu=1", col = "blue")
	text(2.5, 0.50, "nu=2", col = "red")
	

# ------------------------------------------------------------------------------
	

# Plot Skew GED:

	# Density Function:
	xi = 3/2
	d1 = dsged(x, mean = mean, sd = sd, nu = 1, xi = xi)
	d2 = dsged(x, mean = mean, sd = sd, nu = 2, xi = xi)
	d3 = dsged(x, mean = mean, sd = sd, nu = 10, xi = xi)
	plot(x, d1, type = "n", xlim = c(-5, 5), ylim = c(0, 0.8), 
		xlab = "z", ylab = "f(z)", main = "skew GED: Density")
	lines(x, d1, col = "blue")
	lines(x, d2, col = "red")
	lines(x, d3, col = "green")
	text(3.6, 0.15, "nu=10", col = "green")
	text(2.4, 0.40, "nu=2",  col = "red")
	text(1.7, 0.64, "nu=1",  col = "blue")

	# Distribution Function:
	xi = 3/2
	p1 = psged(x, mean = mean, sd = sd, nu = 1, xi = xi)
	p2 = psged(x, mean = mean, sd = sd, nu = 2, xi = xi)
	p3 = psged(x, mean = mean, sd = sd, nu = 10, xi = xi)
	plot(x, p1, type = "n", xlim = c(-5, 5),
		xlab = "z", ylab = "F(z)", main = "skew GED: Distribution")
	lines(x, p1, col = "blue")
	lines(x, p2, col = "red")
	lines(x, p3, col = "green")
	text(2.0, 0.20, "nu=1", col = "blue")
	text(2.7, 0.50, "nu=2", col = "red")
	text(3.8, 0.80, "nu=10", col = "green")
	
	
# ------------------------------------------------------------------------------


par(mfrow = c(2, 2))

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
	abline(h = mean(r), col = "grey")
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
	
distTest("norm")

distTest("t", df=4)