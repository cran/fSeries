
#
# Example:
#   Unit root tests 
#
# Description:
#   This file shows the usage of unit root tests from Bernhard Pfaff's
#	contributed R package 'urca'. Note, you have to install the 'urca'
#   package and to load it calling "library(urca)" or "require(urca)".
#
# Author:
#   (C) 2004 Diethelm Wuertz, GPL
#


################################################################################
# Function Overview:

	# TABLE:
	#             ur.ers      ur.kpss     ur.pp       ur.sp          ur.za
	# ------------------------------------------------------------------------
	#  type       DF-GLS      mu          Z-alpha     tau            -
	#             P-test      tau         Z-tau       rho            -
	#  pol.deg    -           -           -           1:4            -
	#  signif                                         0.01,0.05,0.10 -
	#  model      constant    -           constant    -              intercept
	#             trend       -           trend       -              trend
	#             -           -           -           -              both
	#  lag.max    4           -           -           -              -
	#  lags       -           short       short       -              - 
	#             -           long        long        -              -
	#             -           nil         -           -              -
	#  use.lags   -           NULL        NULL        -              -
	# ------------------------------------------------------------------------


################################################################################
# Load Data:


	# Requirements:
	require(urca)
	###
	
	
	# GNP Data:
	data(nporg)
	gnp = na.omit(nporg[, "gnp.r"])
	par(mfrow = c(2, 2), cex = 0.7)
	gnp.log = log(gnp)
	###
	
	
	# Plot:
	plot(gnp, type = "b", main = "Nelson-Plosser Data")
	plot(diff(gnp), type = "b", main = "Differenced")
	abline(h = mean(diff(gnp)), lty = 3)
	# Plot Log Data:
	plot(gnp.log, type = "b", main = "Log: Nelson-Plosser Data")
	plot(diff(gnp.log), type = "b", main = "Differenced")
	abline(h = mean(diff(gnp.log)), lty = 3)
	###


################################################################################
# ERS - ERS Test


	# DF-GLS | constant:
	ers.gnp = ur.ers(diff(gnp), type = "DF-GLS", model = "constant", 
		lag.max = 4)
	print(ers.gnp)
	plot(ers.gnp)
	summary(ers.gnp)
	###
	
	
	# P-test | constant:
	ers.gnp = ur.ers(diff(gnp), type = "P-test", model = "constant", 
		lag.max = 4)
	print(ers.gnp)
	plot(ers.gnp)
	summary(ers.gnp)
	###
	
	
	# DF-GLS | trend:
	ers.gnp = ur.ers(gnp, type = "DF-GLS", model = "trend", 
		lag.max = 4)
	print(ers.gnp)
	plot(ers.gnp)
	summary(ers.gnp)
	###
	
	
	# P-test | trend:
	ers.gnp = ur.ers(gnp, type = "P-test", model = "trend", 
		lag.max = 4)
	print(ers.gnp)
	plot(ers.gnp)
	summary(ers.gnp)
	###


################################################################################
# KSPP - KSPP Test:
 

	# type - Type of deterministic part. 
	# lags - Maximum number of lags used for error term correction. 
	# use.lag - User specified number of lags. 

	# KPSS | tau |short:
	kpss.gnp = ur.kpss(gnp.log, type = "tau", lags = "short")
	print(kpss.gnp)
	plot(kpss.gnp)
	summary(kpss.gnp)
	###

	# KPSS | tau |long:
	kpss.gnp = ur.kpss(gnp.log, type = "tau", lags = "long")
	print(kpss.gnp)
	plot(kpss.gnp)
	summary(kpss.gnp)
	###
		
	# KPSS | tau |6 lags:
	kpss.gnp = ur.kpss(gnp.log, type = "tau", use.lag = 6)
	print(kpss.gnp)
	plot(kpss.gnp)
	summary(kpss.gnp)
	####
	
	# KPSS | mu | short:
	kpss.gnp = ur.kpss(diff(gnp.log), type = "mu", lags = "short")
	print(kpss.gnp)
	plot(kpss.gnp)
	summary(kpss.gnp)
	###
		
	# KPSS | mu | long:
	kpss.gnp = ur.kpss(diff(gnp.log), type = "mu", lags = "long")
	print(kpss.gnp)
	plot(kpss.gnp)
	summary(kpss.gnp)
	###
		
	# KPSS | mu | 6 lags:
	kpss.gnp = ur.kpss(diff(gnp.log), type = "mu", use.lag = 6)
	print(kpss.gnp)
	plot(kpss.gnp)
	summary(kpss.gnp)
	###


################################################################################
# PP - Phillips & Perron Unit Root Test


	# PP | Ztau | trend | short:
	pp.gnp = ur.pp(gnp, type = "Z-tau", model = "trend", lags = "short")
	print(pp.gnp)
	plot(pp.gnp)
	summary(pp.gnp)
	###
	
	# PP | Ztau | trend | long:
	pp.gnp = ur.pp(gnp, type = "Z-tau", model = "trend", lags = "long")
	print(pp.gnp)
	plot(pp.gnp)
	summary(pp.gnp)
	###
	
	# PP | Ztau | trend | 6 lags:
	pp.gnp = ur.pp(gnp, type = "Z-tau", model = "trend", use.lag = 6)
	print(pp.gnp)
	plot(pp.gnp)
	summary(pp.gnp)
	###
		
	# PP | Ztau | constant | short:
	pp.gnp = ur.pp(diff(gnp), type = "Z-tau", model = "constant", 
		lags = "short")
	print(pp.gnp)
	plot(pp.gnp)
	summary(pp.gnp)
	###
	
	# PP | Ztau | trend | long:
	pp.gnp = ur.pp(diff(gnp), type = "Z-tau", model = "constant", 
		lags = "long")
	print(pp.gnp)
	plot(pp.gnp)
	summary(pp.gnp)
	###
	
	# PP | Ztau | constant | 6 lags:
	pp.gnp = ur.pp(diff(gnp), type = "Z-tau", model = "constant", 
		use.lags = 6)
	print(pp.gnp)
	plot(pp.gnp)
	summary(pp.gnp)
	###
	
	# PP | Zalpha | trend:
	pp.gnp = ur.pp(gnp, type = "Z-alpha", model = "trend")
	print(pp.gnp)
	plot(pp.gnp)
	summary(pp.gnp)
	###
	
	# PP | Zalpha | constant:
	pp.gnp = ur.pp(diff(gnp), type = "Z-alpha", model = "constant")
	print(pp.gnp)
	plot(pp.gnp)
	summary(pp.gnp)
	###


################################################################################
# SP - Schmidt & Phillips Unit Root Test


	# SP:
	sp.gnp = ur.sp(gnp, type = "tau", pol.deg = 1, signif = 0.05) 
	print(sp.gnp)
	plot(sp.gnp)
	summary(sp.gnp)
    ###
    
	# SP:
	ur.sp(gnp, type = c("tau", "rho"), pol.deg = c(1, 2, 3, 4), 
	  signif = c(0.01, 0.05, 0.1)) 
	###
	

################################################################################
# ZA - Zivot and Andrews Unit Root Test


	# ZA | intercept | 2 lags:
	za.gnp = ur.za(gnp, model = "intercept", lag = 2)
	plot(za.gnp)
	summary(za.gnp)
	###
	
	
	# ZA | trend | 2 lags:
	za.gnp = ur.za(gnp, model = "trend", lag = 2)
	plot(za.gnp)
	summary(za.gnp)
	###
	
		
	# ZA | both | 2 lags:
	za.gnp = ur.za(gnp, model = "both", lag = 2)
	plot(za.gnp)
	summary(za.gnp)
	###
	
	
	# ZA:
	za.gnp = urzaTest(gnp, model = "intercept", lag = 2)
	summary(za.gnp)
	###
	
	# ZA:
	za.gnp = urzaTest(gnp, model = "trend", lag = 2)
	summary(za.gnp)
	###
	
	# ZA:
	za.gnp = urzaTest(gnp, model = "both", lag = 2)
	summary(za.gnp)
	###


################################################################################

