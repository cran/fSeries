

# ------------------------------------------------------------------------------


# Settings:

	require(fBasics)
	par (mfrow = c(3, 2), cex = 0.5)
	
	
# Load Data Set:

	data(dem2gbp)
	demgbp = dem2gbp[, 1]/100
	mean(demgbp)
	skewness(demgbp)
	kurtosis(demgbp)
	
	
# Plot Log-Return Series:
	
	ts.plot(demgbp, xlab = "index", ylab = "log Return", col = "steelblue3")
	title(main = "DEM/GBP FX Rate")
	grid()
	abline (h = 0, col = "grey")
	

# Histogram Plot:

	hist(demgbp, n= 20, col = "steelblue", probability = TRUE, 
		border = "white", xlim = c(-0.02, 0.02) )
	abline (h = 0)
	x = seq(-0.02, 0.02, length = 201)
	lines(x, dnorm(x = x, mean = mean(demgbp), sd = sqrt(var(demgbp))))
	
# Quantile - Quantile Plot:
	
    qqnorm(demgbp, xlab = "Normal Quantiles", ylab = "Empirical Quantiles",
    	col = "steelblue3") 
    qqline(demgbp)
    grid()
    
# Partial Autocorrelation:
    
    acf(abs(demgbp))
   
		