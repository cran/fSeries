

# ------------------------------------------------------------------------------


# Settings:

	require(fBasics)
	par (mfrow = c(3, 2), cex = 0.5)
	
	
# Load Data Set:

	data(dem2gbp)
	demgbp = dem2gbp[, 1]/100
	
	
# 
	spec = garchSpec()
	spec
	
	