
#
# Example: 
#   Unit Root Test
#
# Description:
#   Here we show the usage of unit root tests using the wrapper functions
#   based on Bernhard Pfaff's contributed 'urca' package.
#
# Author:
#   (C) 2004, Diethelm Wuertz, GPL
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


	# ur.ers(y, type = c("DF-GLS", "P-test"), model = c("constant", "trend"), lag.max = 4)
	# ur.kpss(y, type = c("mu", "tau"), lags = c("short", "long", "nil"), use.lag = NULL)
	# ur.pp(x, type = c("Z-alpha", "Z-tau"), model = c("constant", "trend"), lags = c("short", "long"), use.lag = NULL)
	# ur.sp(y, type = c("tau", "rho"), pol.deg = c(1, 2, 3, 4), signif = c(0.01, 0.05, 0.1))
	# ur.za(y, model = c("intercept", "trend", "both"), lag)

    
    # ERS:
    urersTest(GNP, type = "DF-GLS", model = "const", lag.max = 4) 
    urersTest(GNP, type = "DF-GLS", model = "trend", lag.max = 4) 
    
    urersTest(GNP, type = "P-test", model = "const", lag.max = 4)   
    urersTest(GNP, type = "P-test", model = "trend", lag.max = 4)   
    ###
    
    
    # KPSS:
    urkpssTest(log(GNP), type = "tau", lags = "nil")  
    urkpssTest(log(GNP), type = "mu", lags = "nil") 
    
    urkpssTest(log(GNP), type = "tau", lags = "short")  
    urkpssTest(log(GNP), type = "mu", lags = "short") 
    
    urkpssTest(log(GNP), type = "tau", lags = "long")  
    urkpssTest(log(GNP), type = "mu", lags = "long")  
         
    urkpssTest(log(GNP), type = "tau", use.lag = 5)  
    urkpssTest(log(GNP), type = "mu", use.lag = 5)   
    ###
    
    
    # PP:
    urppTest(GNP, type = "Z-alpha", model = "constant", lags = "short")
    urppTest(GNP, type = "Z-tau", model = "constant", lags = "short")
    urppTest(GNP, type = "Z-alpha", model = "trend", lags = "short")
    urppTest(GNP, type = "Z-tau", model = "trend", lags = "short")
    
    urppTest(GNP, type = "Z-alpha", model = "constant", lags = "long")
    urppTest(GNP, type = "Z-tau", model = "constant", lags = "long")
    urppTest(GNP, type = "Z-alpha", model = "trend", lags = "long")
    urppTest(GNP, type = "Z-tau", model = "trend", lags = "long")
    
    urppTest(GNP, type = "Z-alpha", model = "constant", use.lag = 5)
    urppTest(GNP, type = "Z-tau", model = "constant", use.lag = 5)
    urppTest(GNP, type = "Z-alpha", model = "trend", use.lag = 5)
    urppTest(GNP, type = "Z-tau", model = "trend", use.lag = 5)
    ###
    

    # SP:
    urspTest(GNP, type = "tau", pol.deg = 1, signif = 0.01)
    urspTest(GNP, type = "tau", pol.deg = 2, signif = 0.01)
    urspTest(GNP, type = "tau", pol.deg = 3, signif = 0.01)
    urspTest(GNP, type = "tau", pol.deg = 4, signif = 0.01)
    
    urspTest(GNP, type = "rho", pol.deg = 1, signif = 0.01)
    urspTest(GNP, type = "rho", pol.deg = 2, signif = 0.01)
    urspTest(GNP, type = "rho", pol.deg = 3, signif = 0.01)
    urspTest(GNP, type = "rho", pol.deg = 4, signif = 0.01)
    ###
    
    
    # ZA:
	urzaTest(GNP, model = "intercept", lag = 2)
	urzaTest(GNP, model = "trend", lag = 2)
	urzaTest(GNP, model = "both", lag = 2)
    ###
    
    
################################################################################

 