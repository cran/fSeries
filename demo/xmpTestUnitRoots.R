
#
# Example: 
#	Unit Root Test
#
# Description:
#	Here we show the usage of unit root tests using the wrapper functions
#	based on Bernhard Pfaff's contributed 'urca' package.
#
# Author:
#	(C) 2004, Diethelm Wuertz, GPL
#


################################################################################
# Load data:
	
	# Nelson Plotter GNP Data:
	data(nelsonplosser)
	GNP = na.omit(nelsonplosser[, "gnp.r"])
	par(mfrow = c(1, ,1))
	plot(1909:1970, GNP, xlab = "Year", main = "Nelson Plotter GNP Data")
	###
	

################################################################################
# Perform Tests:

	
	# ERS:
	ers = urersTest(GNP, type = "DF-GLS", model = "const", lag.max = 4)
	ers
	summary(ers)	
	###
	
	
	# KPSS:
	kpss = urkpssTest(log(GNP), type = "tau", lags = "short")
	kpss
	summary(kpss)	
	###
	
	
	# PP:
	pp = urppTest(GNP, type = "Z-tau", model = "trend", lags = "short")
	pp
	summary(pp)
	###
	

	# SP:
	sp = urspTest(GNP, type = "tau", pol.deg = 1, signif = 0.01)
	sp
	summary(sp)
	###
	
	
	# ZA:
	za = urzaTest(GNP, model = "both", lag = 2)
	za
	summary(za)
	###
	
	
################################################################################

	