
# Exercise:

#	Compute the examples as given in Chapter 4 of the book of Zivot
# 	and Wang "Modelling Financial Time Series with S-Plus"

# CHAPTER 4: Unit Root Tests
#
#	4.1 Introduction
#	4.2 Tests for Nonstationary and Stationary
#	4.3 Autoregressive Unit Root Tests 
#   4.4 Stationarity Tests
#	4.5 References

	
################################################################################
# 4.2 - Testing for Nonstationary and Stationary


	# Simulate TS and DS data:
	# p. 106/107
	set.seed(201) # try other seeds ...
	e = rnorm(250)
	z.DS = cumsum(e)
	z.TS = arima.sim(model = list(ar = 0.75), n = 250, innov = e)
	y.DS = 5 + 0.1*seq(250) + z.DS
	y.TS = 5 + 0.1*seq(250) + z.TS
	# Plot - Figure 4.1:
	par(mfrow = c(2, 1), cex = 0.75)
	ts.plot(y.DS, y.TS, lty = c(1, 2))
	legend(25, max(y.TS), legend = c("I(1)", "I(0)"), lty = c(1, 2))
	###


################################################################################
# 4.3 - Autoregressive Unit Root Tests


	# Simulate Functions of Wiener processes:
	# p. 109
	wiener = function(nobs) {
		e = rnorm(nobs)
		y = cumsum(e)
		ym1 = y[1:(nobs-1)]
		intW2 = nobs^(-2) * sum(ym1^2)
		intWdW = nobs^(-1) * sum(ym1*e[2:nobs])
		ans = list(intW2 = intW2, intWdW = intWdW)
		ans 
	}
	wiener(250)
	###


	# Simulate DF Distributions and Plot Histograms:
	# p. 110
	set.seed(378)
	nobs = 1000
	nsim = 1000
	NB = rep(0, nsim)
	DF = rep(0, nsim)
	for (i in 1:nsim) {
		BN.moments = wiener(nobs)
		NB[i] = BN.moments$intWdW/BN.moments$intW2
		DF[i] = BN.moments$intWdW/sqrt(BN.moments$intW2) 
	}
	###
		

	# Plot histograms and QQ plots of densities:
	# p. 110
	# Plot - Figure 4.2:
	par(mfrow = c(2, 2), cex = 0.75)
	hist(DF, 
		main = "Simulated DF Distribution")
	hist(NB, 
		main = "Simulated Normalized Bias")
	# Plot Densities:
	plot(density(DF), type = "l", 
		main = "Density of Simulated DF", 
		xlab = "DF values", 
		ylab = "Density")
	plot(density(NB), type = "l", 
		main = "Density of Normalized Bias", 
		xlab = "Normalized Bias", 
		ylab = "Density")
	###
	
	
	# Compute the Empirical Quantiles:
	# p. 110/111
	quantile(DF, probs = c(0.01, 0.05, 0.1))
	quantile(NB, probs = c(0.01, 0.05, 0.1))
	###
	

	# Example: p-values and quantiles of DF and NB distributions 
	# using MacKinnon's programs
	args(qunitroot)
	# compute critical values appropriate for sample of size 100
	qunitroot(c(0.01, 0.05, 0.10), trend = "nc", statistic = "t", 
		n.sample = 100)
	qunitroot(c(0.01, 0.05, 0.10), trend = "nc", statistic = "n", 
		n.sample = 100)
	###	
		
	
	# Compute p-value based on the DF distribution:
	# p. 111
	args(punitroot)
	punitroot(-1.645, trend = "nc", statistic = "t")
	###
	
	
################################################################################
# 4.3.2 Trend Cases:
	

	# Simulate data from trend cases I, II 
	# p. 111
	set.seed(201) # try other seeds ...
	e = rnorm(250)
	y1.H0 = cumsum(e)
	y1.H1 = arima.sim(model = list(ar = 0.75), n = 250, innov = e)
	y2.H0 = 5 + y1.H0
	y2.H1 = 5 + y1.H1
	y3.H0 = 5 + 0.1*seq(250) + y1.H0
	y3.H1 = 5 + 0.1*seq(250) + y1.H1		
	# Plot - Figure 4.3:
	par(mfrow = c(3, 2), cex = 0.5)
	ts.plot(y2.H0, main = "Case I: I(1) data")
	ts.plot(y2.H1, main = "Case I: I(0) data")
	ts.plot(y3.H0, main = "Case II: I(1) data")
	ts.plot(y3.H1, main = "Case II: I(0) data")
	###


	# Generate p-values and quantiles for CASE I
	# p. 112/113
	# Constant Only
	qunitroot(c(0.01, 0.05, 0.10), trend = "c", statistic = "t", 
		n.sample = 100)
	qunitroot(c(0.01, 0.05, 0.10), trend = "c", statistic = "n", 
		n.sample = 100)
	punitroot(-1.645, trend = "c", statistic = "t", n.sample = 100)
	punitroot(-1.645, trend = "c", statistic = "n", n.sample = 100)
	###
	
	
	# Generate p-values and quantiles for CASE II
	# p. 113
	# Constant and Time Trend
	qunitroot(c(0.01, 0.05, 0.10), trend = "ct", statistic = "t", 
		n.sample = 100)
	qunitroot(c(0.01, 0.05, 0.10), trend = "ct", statistic = "n", 
		n.sample = 100)
	punitroot(-1.645, trend = "ct", statistic = "t", n.sample = 100)
	punitroot(-1.645, trend = "ct", statistic = "n", n.sample = 100)
	###
	

################################################################################
# 4.3.3 Dickey Fuller Unit Root Tests


    # The file "lexrates.dat.csv" contains spot and forward 
    # exchange rate data.
    # Source: Zivot, E. (2000). Cointegration and forward and spot exchange 
    # rate regressions. Journal of International Money and Finance, 
    # 19(6):785-812. 6:387-401. 
    dataPath = "src/library/urcdist/data/"
    lexrates.dat = read.timeSeries(
        paste(dataPath, "lexrates.dat.csv", sep = ""))
    lexrates.dat[1,]
    end(lexrates.dat)
    ###
	
	
	# Example 18:
	# p. 115-117
	# Testing for a unit root in exchange rate data using ADF tests
	uscn.spot = lexrates.dat[, "USCNS"]
	uscn.spot@title = "Log US/CN spot exchange rate"
	# Plot - Figure 4.4
	par(mfrow = c(3, 2), cex = 0.5)
	plot(uscn.spot, type = "l", col = "steelblue4",
		main = "Log of US/CN spot exchange rate")
	uscn.spot.vec = seriesData(uscn.spot)
	tmp = acf(uscn.spot.vec)
	uscn.diff = diff(uscn.spot.vec)
	plot(uscn.diff, type = "l", col = "steelblue4",
		main = "1st Diff of log US/CN FX Spot")
	tmp = acf(diff(uscn.spot.vec))
	###
	
		
	# Unit Root test:
	# p. 116/117
	# Show Arguments:
	args(unitroot)
	# ADF Test:
	# We have added a function adfTest 
	# > unitroot(uscn.spot, trend = "c", statistic = "t", 
	#	method = "adf", lags = 6)
	# Use:
	uscn.spot.vec = as.vector(uscn.spot)
	adft.out = adfTest(x = uscn.spot.vec, type = "c", lags = 6)
	class(adft.out)
	names(adft.out)
	adft.out
	summary(adft.out)
	###
	
	
	# Now with two lags:
	# p. 118
	adft.out = adfTest(uscn.spot.vec, type = "c", lags = 2)
	summary(adft.out)
	adfn.out = adfTest(uscn.spot.vec, type = "c", lags = 2, statistic = "n")
	adfn.out
	###


	# Example 19: 
	# p. 119/120
	# Testing for a unit root in log monthly stock prices:
	lnp = log(seriesData(singleIndex.dat[, 1]))
	# > lnp@title = "Log of S&P 500 Index"
	# Plot - Figure 4.5:
	par(mfrow = c(3, 2), cex = 0.5)
	plot(lnp, main = "Log of S&P 500 index")
	acf.plot(acf(lnp, plot=F))
	plot(diff(lnp), main = "First difference of log S&P 500 Index")
	acf.plot(acf(diff(lnp), plot = FALSE))
	###
	
	
	# Compute ADF Test:
	adft.out = unitroot(lnp, trend = "ct", lags = 4)
	summary(adft.out)
	###
	
	
################################################################################
# 4.4.1 Simulating the KPSS Distribution


	# Simulate asymptotic distributions for KPSS distribution
	wiener2 = function(nobs) {
		e = rnorm(nobs)
		# Create detrended errors
		e1 = e - mean(e)
		# > e2 = residuals(OLS(e ~ seq(1, nobs)))
		# ... use:
		e2 = residuals(lm(e ~ seq(1, nobs)))
		# Compute simulated Brownian Bridges
		y1 = cumsum(e1)
		y2 = cumsum(e2)
		intW2.1 = nobs^(-2) * sum(y1^2)
		intW2.2 = nobs^(-2) * sum(y2^2)
		ans = list(intW2.1 = intW2.1, intW2.2 = intW2.2)
		ans }
	wiener2(250)

		
	# Simulate KPSS distributions:
	nobs = 1000
	nsim = 1000
	KPSS1 = rep(0, nsim)
	KPSS2 = rep(0, nsim)
	# This can take some time ...
	for (i in 1:nsim) {
		BN.moments = wiener2(nobs)
		KPSS1[i] = BN.moments$intW2.1
		KPSS2[i] = BN.moments$intW2.2 }
	# Compute quantiles of distribution
	par (mfrow = c(2,2), cex= 0.75)
	hist(KPSS1, n = 100)
	hist(KPSS2, n = 100)
	quantile(KPSS1, probs = c(0.90, 0.925, 0.95, 0.975, 0.99))
	quantile(KPSS2, probs = c(0.90, 0.925, 0.95, 0.975, 0.99))
	###
	

################################################################################
# 4.4.2 Testing for Stationarity


	stationaryTest(x, trend = c("c", "ct"), bandwidth = NULL, na.rm = FALSE)

	ur.kpss(y, type = c("mu", "tau"), lags = c("short", "long", "nil"), 
		use.lag = NULL) 
    kpss.test(x, null = c("Level", "Trend"), lshort = TRUE) 
    
    kpss.test(uscn.spot.vec, null = "Level")
    

	# Stationary Test on Exchange Rates
	# > args(stationaryTest)
	# > kpss.out = stationaryTest(uscn.spot, trend = "c")
	# > class(kpss.out)
	# > kpss.out
	### 
	

################################################################################

	