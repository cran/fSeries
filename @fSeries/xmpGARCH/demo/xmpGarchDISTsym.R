
#
# Example 1:
# 	Plot the symmetric GED and Student-t Distribution,
# 	use the distTest() function
#

#
# Example 2: 
#	Plot the GED Distribution in the limit of the 
#	- Laplace Distribution, mu=1, in the limit of the 
#   - Gaussian Distribution, mu=2, and in the limit of the
#   - Uniform distribution, mu->Inf.
#   

#
# Example 3: 
#	Plot the Student-t Distribution for
#	- nu = 2.5 resulting in a quite fat tailed distribution, 
#   - nu = 5.0 resulting in a less fat tailed distributuin, and 
#   - nu = 10.0 resulting in a almost Normal distribution.
#   


# ------------------------------------------------------------------------------


# Example 1:

	# Load Function:
	source("src/library/xmpGARCH/demo/fun-distTest.R")
	
	# Test:
	par(mfrow = c(2,2), cex = 0.5)
	distTest("norm")
	distTest("ged")
	distTest("std")
	
	# Change Parameters:
	distTest("norm", breaks = 20, mean = 1, sd = 0.5)
	distTest("ged", breaks = 50, nu = 1)
	
	
# ------------------------------------------------------------------------------


# Example 2: Investigate the Special Cases:
    
	# nu = 1   : Laplace Distribution
    # nu = 2   : Gaussian Distribution
    # nu = Inf : Uniform Distribution
    
    # Figure 1

	# Settings:
	mean = 0; sd = 1; dx = 0.001
	x = seq(-5+mean, 5+mean, by = dx)

	# 1. Plot Density Function:
	par(mfrow = c(2, 2), cex = 0.75)
	d1 = dged(x, mean = mean, sd = sd, nu = 1)
	d2 = dged(x, mean = mean, sd = sd, nu = 2)
	d3 = dged(x, mean = mean, sd = sd, nu = 10)
	plot(x, d1, type = "n", xlim = c(-5, 5), ylim = c(0, 0.7), 
		xlab = "z", ylab = "f(z)", main = "GED: Density")
	lines(x, d1, col = "blue")
	lines(x, d2, col = "red")
	lines(x, d3, col = "green")
	text(3.6, 0.15, "nu=10", col = "green")
	text(2.4, 0.40, "nu=2",  col = "red")
	text(1.7, 0.64, "nu=1",  col = "blue")

	# 2. Plot Distribution Function:
	p1 = pged(x, mean = mean, sd = sd, nu = 1)
	p2 = pged(x, mean = mean, sd = sd, nu = 2)
	p3 = pged(x, mean = mean, sd = sd, nu = 10)
	plot(x, p1, type = "n", xlim = c(-5, 5),
		xlab = "z", ylab = "F(z)", main = "GED: Distribution")
	lines(x, p1, col = "blue")
	lines(x, p2, col = "red")
	lines(x, p3, col = "green")
	text(2.0, 0.20, "nu=1", col = "blue")
	text(2.7, 0.50, "nu=2", col = "red")
	text(3.8, 0.80, "nu=10", col = "green")
	
	# 3. Plot Moment Ratios:
	plot(c(0, 11), c(0, 15), type = "n", 
		xlab = "n", ylab = "M(n)/M(n-1)", main = "absMoment Ratio")
	for (nu in c(3/4, 1, 2)) {
		# True Values:
		m1 = absMoments(n = 1:10, density = "dged", nu = nu)
		m2 = absMoments(n = 2:11, density = "dged", nu = nu)
		cat("\n nu: ", nu, "\n")
		print(diff(m2/m1))
		points(1:10, m2/m1)
		lines (1:10, m2/m1) 
		text(10,  2.4, "nu=2")
		text(10,  6.6, "nu=1")
		text(10, 11.6, "nu=3/4") }

	# 4. Plot M4 (Kurtosis):
	# The Kurtosis: M4/M2^2 - 3
	# Figure 4:
	nu = seq(0.1, 10, length = 100)
	M4 = absMoments(n = 4, density = "dged", nu = nu)
	# Plot:
	plot(nu, log(M4), type = "l", ylim = c(0,10), 
		xlab = "nu", ylab = "ln M(4)", main = "GED: 4th Moment")
	abline(v = 0, lty = 3)
	abline(h = log(3), lty = 3)
	points(2, log(3))
	points(1, log(6))
	text(3.1, log(4.8), "Normal")
	text(2.1, log(9.3), "Laplace")
	text(9.1, log(1.2), "Uniform")
	
		
# ------------------------------------------------------------------------------


# Example 3: Plot the Student-t Distribution for

	#	- nu = 2.5 resulting in a quite fat tailed distribution, 
	#   - nu = 5.0 resulting in a less fat tailed distributuin, and 
	#   - nu = 10.0 resulting in a almost Normal distribution. 
	
	# Figure 2

	# Settings:
	mean = 0; sd = 1; dx = 0.001
	x = seq(-10+mean, 10+mean, by = dx)

	# 1. Plot Density Function:
	par(mfrow = c(2, 2), cex = 0.75)
	d1 = dstd(x, mean = mean, sd = sd, nu = 2.5)
	d2 = dstd(x, mean = mean, sd = sd, nu = 5)
	d3 = dstd(x, mean = mean, sd = sd, nu = 10)
	plot(x, d1, type = "n", xlim = c(-5, 5), ylim = c(0, 0.9), 
		xlab = "z", ylab = "f(z)", main = "Student-t: Density")
	lines(x, d1, col = "blue")
	lines(x, d2, col = "red")
	lines(x, d3, col = "green")
	text(3.6, 0.15, "nu=10", col = "green")
	text(2.4, 0.40, "nu=5",  col = "red")
	text(1.7, 0.64, "nu=2.5",  col = "blue")
	
	# 2. Plot Distribution Function:
	p1 = pstd(x, mean = mean, sd = sd, nu = 2.5)
	p2 = pstd(x, mean = mean, sd = sd, nu = 5)
	p3 = pstd(x, mean = mean, sd = sd, nu = 10)
	plot(x, p1, type = "n", xlim = c(-5, 5),
		xlab = "z", ylab = "F(z)", main = "Student-t Distribution")
	lines(x, p1, col = "blue")
	lines(x, p2, col = "red")
	lines(x, p3, col = "green")
	text(1.2, 0.20, "nu=2.5",  col = "blue")
	text(1.7, 0.50, "nu=5",  col = "red")
	text(2.6, 0.80, "nu=10", col = "green")
	
	# 3. Plot Random Deviates:
	set.seed(0.4711)
	r = rstd(n = 10000, mean = 1, sd = 2, nu = 2.5)
	mean(r)
	sqrt(var(r))
	hist(r, n = 500, probability = TRUE, xlim = c(-5, 7), 
		main = "10000 Random Deviates")
	x = seq(-5, 7, length = 201)
	lines(x, dstd (x, mean = 1, sd = 2, nu = 2.5), lwd = 2, col = "red")
	lines(x, dnorm(x, mean = 1, sd = 2), col = "blue")

	
	# 4. Plot the Kurtosis:
	nu = seq(2.1, 20, by = 0.1)
	M4 = M2 = rep(NA, times = length(nu))
	for (i in 1:length(nu)) {
		M4[i] = absMoments(n = 4, density = "dstd", nu[i])
		M2[i] = absMoments(n = 2, density = "dstd", nu[i]) }
	Kurtosis = M4/M2^2 - 3
	plot(nu, Kurtosis, type = "l", xlim = c(0, 20), ylim = c(0, 10),
		col = "steelblue4", main = "Kurtosis")
	abline(v = 2, lty = 3)
	abline(h = 0, lty = 2)
	
	
# ------------------------------------------------------------------------------

