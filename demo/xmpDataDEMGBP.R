

# Example:
#	Investigate the dem2gbp data set


	# Settings:
	par (mfrow = c(3, 2), cex = 0.7)
	
	
	# Load Data Set:
	data(dem2gbp)
	dem2gbp = dem2gbp[, 1]/100
	
	
	# Plot Log-Return Series:
	ts.plot(dem2gbp, xlab = "index", ylab = "log Return", col = "steelblue3",
		main = "DEMGBP FX Rate")
	grid()
	abline (h = 0, col = "grey")
	

	# Histogram Plot:
	hist(dem2gbp, col = "steelblue", probability = TRUE, 
		breaks = "FD", border = "white", xlim = c(-0.02, 0.02))
	abline (h = 0)
	x = seq(-0.02, 0.02, length = 201)
	lines(x, dnorm(x = x, mean = mean(dem2gbp), sd = sqrt(var(dem2gbp))))
	
	
	# Quantile - Quantile Plot:
    qqnorm(dem2gbp, xlab = "Normal Quantiles", ylab = "Empirical Quantiles",
    	col = "steelblue3", main = "Normal QQ-Plot", cex = 0.5) 
    qqline(dem2gbp)
    grid()
    
    
	# Partial Autocorrelation Plot:  
    pacf(dem2gbp)
  
    
    # ACF Volatilities Plot:
    acf(abs(dem2gbp))
   
	
	
	
	