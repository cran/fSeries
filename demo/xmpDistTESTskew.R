
#
# Example:
#	Skew Normal, Student-t and GED Distribution
#
# Description:
#   Test if the Skew NORMAL, STUDENT-t, and GED have the desired properties:
#	1. Fullfills the density the relation f(x|xi) = f(-x|1/xi)?
#	2. Is the density normalized?
#   3. Is the first moment equal to "mean"?
#   4. Is the sqrt(variance) equal to "sd"?
#   5. Is the computed distribution the same as the integrated?
#   6. Is the quantile function the inverse of the distribution function?
#   7. Have the random deviates proper mean and variance?
#   8. Fits the true density to the histogram?
#
# Author:
#	(C) 2002, Diethelm Wuertz, GPL
#


################################################################################
# SKEW NORMAL SPECIFICATION TESTS


	# Test if the Normal Distribution has the desired properties:
	# 1. Fullfills the density the relation f(x|xi) = f(-x|1/xi)?
	# 2. Is the density normalized?
	# 3. Is the first moment equal to "mean"?
	# 4. Is the sqrt(variance) equal to "sd"?
	# 5. Is the computed distribution the same as the integrated?
	# 6. Is the quantile function the inverse of the distribution function?
	# 7. Have the random deviates proper mean and variance?
	# 8. Fits the true density to the histogram?
	###
	
	
	# Test snorm Settings:
	par(mfrow = c(1, 1))
	mean = 0.5; sd = 1.5; xi = 1.5
	x = -4:4; q = x; n = 100000
	###
		
	
	# Test snorm Density Function:
	# 1. Fullfills the density the relation f(x|xi) = f(-x|1/xi)?
	Density1 = dsnorm( x, mean = 0, sd = 1, xi = xi)
	Density2 = dsnorm(-x, mean = 0, sd = 1, xi = 1/xi)
	c(mean = mean, sd = sd, xi = xi)
	cbind(x, Density1, Density2)
	# 2. Is the density normalized?
	f = function(x, mean, sd, nu, xi) { 
	  dsnorm(x, mean = mean, sd = sd, xi = xi) }
	Norm = integrate(f, -Inf, Inf, mean = mean, sd = sd, xi = xi)$value
	c(mean = mean, sd = sd, xi = xi)
	c(norm = 1, Norm = Norm)	
	# 3. Is the first moment equal to "mean"?
	f = function(x, mean, sd, nu, xi) { 
	  x * dsnorm(x, mean = mean, sd = sd, xi = xi) }
	mu = integrate(f, -Inf, Inf, mean = mean, sd = sd, xi = xi)$value
	c(mean = mean, sd = sd, xi = xi)
	c(mean = mean, mu = mu)	
	# 4. Is the variance equal to "sd"?
	f = function(x, mean, sd, nu, xi) { 
		(x-mean)^2 * dsnorm(x, mean = mean, sd = sd, xi = xi) }
	sigma = sqrt(integrate(f, -Inf, Inf, mean = mean, sd = sd, xi = xi)$value)
	c(mean = mean, sd = sd, xi = xi)
	c(sd = sd, sigma = sigma)
	###
	
	
	# Test snorm Distribution Function:
	# 5. Is the computed distribution the same as the integrated?	
	Probability1 = psnorm(q = q, mean = mean, sd = sd, xi = xi)
	Probability2 = NULL
	for (i in 1:length(q)) Probability2 = c(Probability2, integrate(dsnorm, 
		-Inf, q[i], mean = mean, sd = sd, xi = xi)$value)
	c(mean = mean, sd = sd, xi = xi)
	cbind(q, Probability1, Probability2)
	###
	
	
	# Test snorm Quantile Function:	
	# 6. Is the quantile function the inverse of the distribution function?
	Probability = psnorm(q = q, mean = mean, sd = sd, xi = xi)
	Quantile = qsnorm(p = Probability, mean = mean, sd = sd, xi = xi)
	c(mean = mean, sd = sd, xi = xi)
	cbind(q, Probability, Quantile)
	
	
	# Test snorm Random Variates:
	# 7. Have the random deviates proper mean and variance?
	r = rsnorm(n = n, mean = mean, sd = sd, xi = xi)
	c(mean = mean, sd = sd, xi = xi)
	c(mu = mean(r), sigma =  sqrt(var(r)))	
	# 8. Fits the true density to the histogram?
	x1 = x[1] * sd + mean
	x2 = x[length(x)] * sd + mean
	hist(r, n = 100, probability = TRUE, xlim = c(x1, x2))
	z = seq(x1, x2, length = 801)
	lines(z, dsnorm(z, mean = mean, sd = sd, xi = xi), col = "red")


################################################################################
# SYMMETRIC-STUDENT-t SPECIFICATION TESTS


	# Test if the Symmetric tudent-t has the desired properties:
	# 1. Fullfills the density the relation f(x|xi) = f(-x|1/xi)?
	# 2. Is the density normalized?
	# 3. Is the first moment equal to "mean"?
	# 4. Is the sqrt(variance) equal to "sd"?
	# 5. Is the computed distribution the same as the integrated?
	# 6. Is the quantile function the inverse of the distribution function?
	# 7. Have the random deviates proper mean and variance?
	# 8. Fits the true density to the histogram?
	###
	
	
	# Test STD Settings:
	par(mfrow = c(1, 1))
	mean = 0.5; sd = 1.5; nu = 2.5
	x = -4:4; q = x; n = 100000
	###
	
		
	# Test STD density function:
	# 1. Fullfills the density the relation f(x) = f(-x)?
	Density1 = dstd( x, mean = 0, sd = 1, nu = nu)
	Density2 = dstd(-x, mean = 0, sd = 1, nu = nu)
	c(mean = mean, sd = sd, nu = nu)
	cbind(x, Density1, Density2)
	# 2. Is the density normalized?
	f = function(x, mean, sd, nu) { 
	  dstd(x, mean = mean, sd = sd, nu = nu) }
	Norm = integrate(f, -Inf, Inf, mean = mean, sd = sd, nu = nu)$value
	c(mean = mean, sd = sd, nu = nu)
	c(norm = 1, Norm = Norm)
	# 3. Is the first moment equal to "mean"?
	f = function(x, mean, sd, nu) { 
	  x * dstd(x, mean = mean, sd = sd, nu = nu) }
	mu = integrate(f, -Inf, Inf, mean = mean, sd = sd, nu = nu)$value
	c(mean = mean, sd = sd, nu = nu)
	c(mean = mean, mu = mu)
	# 4. Is the variance equal to "sd"?
	f = function(x, mean, sd, nu) { 
	  (x-mean)^2 * dstd(x, mean = mean, sd = sd, nu = nu) }
	sigma = sqrt(integrate(f, -Inf, Inf, mean = mean, sd = sd, nu = nu)$value)
	c(mean = mean, sd = sd, nu = nu)
	c(sd = sd, sigma = sigma)
	###
	
	# Test STD distribution function:
	# 5. Is the computed distribution the same as the integrated?	
	Probability1 = pstd(q = q, mean = mean, sd = sd, nu = nu)
	Probability2 = NULL
	for (i in 1:length(q)) Probability2 = c(Probability2, integrate(dstd, 
	  -Inf, q[i], mean = mean, sd = sd, nu = nu)$value)
	c(mean = mean, sd = sd, nu = nu)
	cbind(q, Probability1, Probability2)
	###
	
	
	# Test STD quantile function:	
	# 6. Is the quantile function the inverse of the distribution function?
	Probability = pstd(q = q, mean = mean, sd = sd, nu = nu)
	Quantile = qstd(p = Probability, mean = mean, sd = sd, nu = nu)
	c(mean = mean, sd = sd, nu = nu)
	cbind(q, Probability, Quantile)
	###
	
	
	# Test STD random variates:
	# 7. Have the random deviates proper mean and variance?
	r = rstd(n = n, mean = mean, sd = sd, nu = nu)
	c(mean = mean, sd = sd, nu = nu)
	c(mu = mean(r), sigma =  sqrt(var(r)))	
	# 8. Fits the true density to the histogram?
	x1 = x[1]*sd + mean
	x2 = x[length(x)]*sd + mean
	hist(r, n = 500, probability = TRUE, xlim = c(x1, x2))
	z = seq(x1, x2, length = 801)
	lines(z, dstd(z, mean = mean, sd = sd, nu = nu), col = "red")
	###
	

################################################################################
# SKEW-STUDENT-t SPECIFICATION TESTS


	# Test if the Skew Student-t has the desired properties:
	# 1. Fullfills the density the relation f(x|xi) = f(-x|1/xi)?
	# 2. Is the density normalized?
	# 3. Is the first moment equal to "mean"?
	# 4. Is the sqrt(variance) equal to "sd"?
	# 5. Is the computed distribution the same as the integrated?
	# 6. Is the quantile function the inverse of the distribution function?
	# 7. Have the random deviates proper mean and variance?
	# 8. Fits the true density to the histogram?
	###
	
	
	# Test SSTD Settings:
	par(mfrow = c(1, 1))
	mean = 0.5; sd = 1.5; nu = 2.5; xi = 1.5
	x = -4:4; q = x; n = 100000
	###
	
	
	# Test SSTD density function:
	# 1. Fullfills the density the relation f(x|xi) = f(-x|1/xi)?
	Density1 = dsstd( x, mean = 0, sd = 1, nu = nu, xi = xi)
	Density2 = dsstd(-x, mean = 0, sd = 1, nu = nu, xi = 1/xi)
	c(mean = mean, sd = sd, nu = nu, xi = xi)
	cbind(x, Density1, Density2)		
	# 2. Is the density normalized?
	f = function(x, mean, sd, nu, xi) { 
		dsstd(x, mean = mean, sd = sd, nu = nu, xi = xi) }
	Norm = integrate(f, -Inf, Inf, mean = mean, sd = sd, nu = nu, 
		xi = xi)$value
	c(mean = mean, sd = sd, nu = nu, xi = xi)
	c(norm = 1, Norm = Norm)		
	# 3. Is the first moment equal to "mean"?
	f = function(x, mean, sd, nu, xi) { 
		x * dsstd(x, mean = mean, sd = sd, nu = nu, xi = xi) }
	mu = integrate(f, -Inf, Inf, mean = mean, sd = sd, nu = nu, 
		xi = xi)$value
	c(mean = mean, sd = sd, nu = nu, xi = xi)
	c(mean = mean, mu = mu)	
	# 4. Is the variance equal to "sd"?
	f = function(x, mean, sd, nu, xi) { 
		(x-mean)^2 * dsstd(x, mean = mean, sd = sd, nu = nu, xi = xi) }
	sigma = sqrt(integrate(f, -Inf, Inf, mean = mean, sd = sd, 
		nu = nu, xi = xi)$value)
	c(sd = sd, sigma = sigma)
	
	
	# Test SSTD distribution function:
	# 5. Is the computed distribution the same as the integrated?	
	Probability1 = psstd(q = q, mean = mean, sd = sd, nu = nu, xi = xi)
	Probability2 = NULL
	for (i in 1:length(q)) Probability2 = c(Probability2, integrate(dsstd, 
		-Inf, q[i], mean = mean, sd = sd, nu = nu, xi = xi)$value)
	c(mean = mean, sd = sd, nu = nu, xi = xi)
	cbind(q, Probability1, Probability2)
	###
	

	# Test SSTD quantile function:	
	# 6. Is the quantile function the inverse of the distribution function?
	Probability = psstd(q = q, mean = mean, sd = sd, nu = nu, xi = xi)
	Quantile = qsstd(p = Probability, mean = mean, sd = sd, nu = nu, xi = xi)
	c(mean = mean, sd = sd, nu = nu, xi = xi)
	cbind(q, Probability, Quantile)
	###
	
	# Test SSTD random variates:
	# 7. Have the random deviates proper mean and variance?
	r = rsstd(n = n, mean = mean, sd = sd, nu = nu, xi = xi)
	c(mean = mean, sd = sd, nu = nu, xi = xi)
	c(mu = mean(r), sigma =  sqrt(var(r)))	
	# 8. Fits the true density to the histogram?
	x1 = x[1] * sd + mean
	x2 = x[length(x)] * sd + mean
	hist(r, n = 500, probability = TRUE, xlim = c(x1, x2))
	z = seq(x1, x2, length = 801)
	lines(z, dsstd(z, mean = mean, sd = sd, nu = nu, xi = xi), col = "red")


################################################################################
# SYMMETRIC-GED SPECIFICATION TESTS


	# Test if the Symmetric GED has the desired properties:
	# 1. Fullfills the density the relation f(x|xi) = f(-x|1/xi)?
	# 2. Is the density normalized?
	# 3. Is the first moment equal to "mean"?
	# 4. Is the sqrt(variance) equal to "sd"?
	# 5. Is the computed distribution the same as the integrated?
	# 6. Is the quantile function the inverse of the distribution function?
	# 7. Have the random deviates proper mean and variance?
	# 8. Fits the true density to the histogram?


	# Test GED Settings:
	par(mfrow = c(1, 1))
	mean = 0.5; sd = 1/2; nu = 2.5
	x = -4:4; q = x; n = 100000
	###
		
	
	# Test GED density function:
	# 1. Fullfills the density the relation f(x) = f(-x)?
	Density1 = dged( x, mean = 0, sd = 1, nu = nu)
	Density2 = dged(-x, mean = 0, sd = 1, nu = nu)
	c(mean = mean, sd = sd, nu = nu)
	cbind(x, Density1, Density2)
	# 2. Is the density normalized?
	f = function(x, mean, sd, nu) { 
	  dged(x, mean = mean, sd = sd, nu = nu) }
	Norm = integrate(f, -Inf, Inf, mean = mean, sd = sd, nu = nu)$value
	c(mean = mean, sd = sd, nu = nu)
	c(norm = 1, Norm = Norm)
	# 3. Is the first moment equal to "mean"?
	f = function(x, mean, sd, nu) { 
	  x * dged(x, mean = mean, sd = sd, nu = nu) }
	mu = integrate(f, -Inf, Inf, mean = mean, sd = sd, nu = nu)$value
	c(mean = mean, sd = sd, nu = nu)
	c(mean = mean, mu = mu)
	# 4. Is the variance equal to "sd"?
	f = function(x, mean, sd, nu) { 
	  (x-mean)^2 * dged(x, mean = mean, sd = sd, nu = nu) }
	sigma = sqrt(integrate(f, -Inf, Inf, mean = mean, sd = sd, nu = nu)$value)
	c(mean = mean, sd = sd, nu = nu)
	c(sd = sd, sigma = sigma)
	###
	
	
	# Test GED distribution function:
	# 5. Is the computed distribution the same as the integrated?	
	Probability1 = pged(q = q, mean = mean, sd = sd, nu = nu)
	Probability2 = NULL
	for (i in 1:length(q)) Probability2 = c(Probability2, integrate(dged, 
	  -Inf, q[i], mean = mean, sd = sd, nu = nu)$value)
	c(mean = mean, sd = sd, nu = nu)
	cbind(q, Probability1, Probability2)
	###
	
	# Test GED quantile function:	
	# 6. Is the quantile function the inverse of the distribution function?
	Probability = pged(q = q, mean = mean, sd = sd, nu = nu)
	Quantile = qged(p = Probability, mean = mean, sd = sd, nu = nu)
	c(mean = mean, sd = sd, nu = nu)
	cbind(q, Probability, Quantile)
	###
	
	
	# Test GED random variates:
	# 7. Have the random deviates proper mean and variance?
	r = rged(n = n, mean = mean, sd = sd, nu = nu)
	c(mean = mean, sd = sd, nu = nu)
	c(mu = mean(r), sigma =  sqrt(var(r)))
	# 8. Fits the true density to the histogram?
	x1 = 3 * x[1] * sd + mean
	x2 = 3 * x[length(x)] * sd + mean
	hist(r, n = 100, probability = TRUE, xlim = c(x1, x2))
	z = seq(x1, x2, length = 801)
	lines(z, dged(z, mean = mean, sd = sd, nu = nu), col = "red")
	###


################################################################################
# SKEW-GED SPECIFICATION TEST


	# Test if the Skew GED has the desired properties:
	# 1. Fullfills the density the relation f(x|xi) = f(-x|1/xi)?
	# 2. Is the density normalized?
	# 3. Is the first moment equal to "mean"?
	# 4. Is the sqrt(variance) equal to "sd"?
	# 5. Is the computed distribution the same as the integrated?
	# 6. Is the quantile function the inverse of the distribution function?
	# 7. Have the random deviates proper mean and variance?
	# 8. Fits the true density to the histogram?
	###

	
	# Test SGED Settings
	par(mfrow = c(1, 1))
	mean = 0.5; sd = 1/2; nu = 2.5; xi = 1.5
	x = -4:4; q = x; n = 100000
	###
	
		
	# Test SGED density function:
	# 1. Fullfills the density the relation f(x|xi) = f(-x|1/xi)?
	Density1 = dsged( x, mean = 0, sd = 1, nu = nu, xi = xi)
	Density2 = dsged(-x, mean = 0, sd = 1, nu = nu, xi = 1/xi)
	c(mean = mean, sd = sd, nu = nu, xi = xi)
	cbind(x, Density1, Density2)
	# 2. Is the density normalized?
	f = function(x, mean, sd, nu, xi) { 
	  dsged(x, mean = mean, sd = sd, nu = nu, xi = xi) }
	Norm = integrate(f, -Inf, Inf, mean = mean, sd = sd, nu = nu, 
	  xi = xi)$value
	c(mean = mean, sd = sd, nu = nu, xi = xi)
	c(norm = 1, Norm = Norm)
	###
	
	
	# 3. Is the first moment equal to "mean"?
	f = function(x, mean, sd, nu, xi) { 
	  x * dsged(x, mean = mean, sd = sd, nu = nu, xi = xi) }
	mu = integrate(f, -Inf, Inf, mean = mean, sd = sd, nu = nu, xi = xi)$value
	c(mean = mean, sd = sd, nu = nu, xi = xi)
	c(mean = mean, mu = mu)	
	# 4. Is the sqrt(variance) equal to "sd"?
	f = function(x, mean, sd, nu, xi) { 
	  (x-mean)^2 * dsged(x, mean = mean, sd = sd, nu = nu, xi = xi) }
	sigma = sqrt(integrate(f, -Inf, Inf, mean = mean, sd = sd, 
	  nu = nu, xi = xi)$value)
	c(mean = mean, sd = sd, nu = nu, xi = xi)
	c(sd = sd, sigma = sigma)
	###
	
	
	# Test SGED distribution function:
	# 5. Is the computed distribution the same as the integrated?	
	Probability1 = psged(q = q, mean = mean, sd = sd, nu = nu, xi = xi)
	Probability2 = NULL
	for (i in 1:length(q)) Probability2 = c(Probability2, integrate(dsged, 
	  -Inf, q[i], mean = mean, sd = sd, nu = nu, xi = xi)$value)
	c(mean = mean, sd = sd, nu = nu, xi = xi)
	cbind(q, Probability1, Probability2)
	###
	

	# Test SGED quantile function:	
	# 6. Is the quantile function the inverse of the distribution function?
	Probability = psged(q = q, mean = mean, sd = sd, nu = nu, xi = xi)
	Quantile = qsged(p = Probability, mean = mean, sd = sd, nu = nu, xi = xi)
	c(mean = mean, sd = sd, nu = nu, xi = xi)
	cbind(q, Probability, Quantile)
	###
	
	
	# Test SGED random variates:
	# 7. Have the random deviates proper mean and variance?
	r = rsged(n = n, mean = mean, sd = sd, nu = nu, xi = xi)
	c(mean = mean, sd = sd, nu = nu, xi = xi)
	c(mu = mean(r), sigma =  sqrt(var(r)))	
	# 8. Fits the true density to the histogram?
	x1 = 3 * x[1] * sd + mean
	x2 = 3 * x[length(x)] * sd + mean
	hist(r, n = 100, probability = TRUE, xlim = c(x1, x2))
	z = seq(x1, x2, length = 801)
	lines(z, dsged(z, mean = mean, sd = sd, nu = nu, xi = xi), col = "red")


################################################################################

