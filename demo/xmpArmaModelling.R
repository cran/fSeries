
#
# Example: 
#   Analyse and Model a time series process using the ARIMA approach 
#   PART I:   ARIMA Analysis 
#   PART II:  another ARIMA Analysis 
#
# Description:
#	The first part of this example [formerly xmpArmaAnalysis] shows the 
#   different steps in an ARIMA Analysis:
#     	Step 1: Simulate ARIMA time series processes
#   	Step 2: Investigate Correlations
#   	Step 3: Compute Roots of AR / MA polynomials
#   	Step 4: Compute ACF of the true time series model
#   	Step 5: Identify the Order of an ARIMA process 
#	The second part of this example [formerly xmpArmaModelling] shows the 
#   different steps in an another ARIMA Analysis:
#      	Step 1: Model Selection
#      	Step 2: Parameter Estimation
#      	Step 3: Diagnostic Analysis
#      	Step 4: Forecasts
#
# Author:
#	(C) 2002, Diethelm Wuertz, GPL
#


################################################################################
# PART I: ARIMA Analysis and Modeling Process


	# PART I: Step 1 - Simulate ARIMA time series processes
   	# Settings:
    par(mfrow = c(4, 3), cex = 0.6)
    set.seed(471)
   	# AR(1) and AR(2) Examples:  
   	x = armaSim(model = list(ar = 0.5, d = 0, ma = 0), n = 1000)
      ts.plot(x, main = "AR(1): +0.5")
      acf(x, lag.max = 12)
      acf(x, type = "partial", lag.max = 12)
    x = armaSim(model = list(ar = -0.5, d = 0, ma = 0), n = 1000)
      ts.plot(x, main = "AR(1): -0.5")
      acf(x, lag.max = 12)
      acf(x, type = "partial", lag.max = 12)
    x = armaSim(model = list(ar = c(0.5, 0.3), d = 0, ma = 0), n = 1000)
      ts.plot(x, main = "AR(2): +0.5, +0.3")
      acf(x, lag.max = 12)
      acf(x, type = "partial", lag.max = 12)    
    x = armaSim(model=list(ar = c(-0.5, 0.3), d = 0, ma = 0), n = 1000)
      ts.plot(x, main = "AR(2): -0.5, +0.3")
      acf(x, lag.max = 12)
      acf(x, type = "partial",lag.max = 12)
    ###
        
    
	# PART I: Step 2 - Investigate Correlations
	# Settings:  
   	par(mfrow = c(4, 3), cex = 0.6)
   	set.seed(471)
    # MA and ARMA Examples:
    x = armaSim(n = 1000, model = list(ar = 0, d = 0, ma = 0.8))
      ts.plot(x, main = "MA(1): +0.8")
      acf(x, lag.max = 12)
      acf(x, type = "partial", lag.max = 12)
    x = armaSim(n = 1000, model = list(ar = 0, d = 0, ma = -0.8))
      ts.plot(x, main = "MA(1): -0.8")
      acf(x, lag.max = 12)
      acf(x, type = "partial", lag.max = 12)     
    x = armaSim(n = 1000, model = list(ar = 0.5, d = 0, ma = 0.8))
      ts.plot(x, main = "ARMA(1,1): +0.5, +0.8")
      acf(x, lag.max = 12)
      acf(x, type = "partial", lag.max = 12)   
    x = armaSim(n = 1000, model = list(ar = -0.5, d = 0, ma = 0.8))
      ts.plot(x, main = "ARMA(1,1): -0.5, +0.8")
      acf(x, lag.max = 12)
      acf(x, type = "partial", lag.max = 12)
	###

	
    # Part I: Step 3 - Compute Roots of AR / MA polynomials
	# Settings:
   	par(mfrow = c(2, 2))  
   	# Roots - Example 1: 
   	armaRoots(c( 1.6, 0.4))
   	armaRoots(c(-1.6, 0.4))   
    # Roots - Example 2:
    armaRoots(c( 0.5, -0.1))
    armaRoots(c(-0.5, -0.1))
    # Roots - Example 3:
    armaRoots(c( 0.5, -0.9,  0.1, 0.5))
    armaRoots(c(-0.5, -0.9, -0.1, 0.5))
	###
   
	
	# Part I: Step 4 - Compute ACF of the true time series model
	# Settings:
    par(mfcol = c(3, 2), cex = 0.6)   
    # First Model:
    model = list(ar = c(+0.5, 0.3), d = 0, ma = 0)
    x = armaSim(model, n = 1000)
    ts.plot(x, main = "AR(2): +0.5, 0.3")
    armaTrueacf(model, lag.max = 12)
    armaTrueacf(model, lag.max = 12, type = "partial")  
    # Second Model:
    model = list(ar = c(-0.5, 0.3), d = 0, ma = 0)
    x = armaSim(model, n = 1000)
    ts.plot(x, main = "AR(2): -0.5, 0.3")
    armaTrueacf(model, lag.max = 12)
    armaTrueacf(model, lag.max = 12, type = "partial")
    ###
    
    
	# Part I: Step 5 - Identify the Order of an ARIMA process   
    # Settings:
    par (mfrow = c(2, 1))
    set.seed(4711)
    # Simulate a Subset AR(5) Process:
    x = armaSim(model = list(ar = c(-0.4, 0.1, 0, 0, 0.1), 
      d = 0, ma = 0), n = 5000)
    # Plot:
    ts.plot(x, xlab="t")
    title(main = "Subset AR(5)") 
    # Indentify the model with the help of the PACF ...
    # PACF: Note, the PACF is consistent with a subset AR(5)  
    # model with nonzero parameters phi at positions 1, 2 and 5.
    acf(x, type = "partial", lag.max = 12)
    ###
   
  
################################################################################
# PART II: Another ARIMA Analysis and Modeling Process


	# PART II: Step 1 - Model Selection
    # Settings:
    par (mfrow = c(2, 1))
    set.seed(4732) 
    # Simulate time series:
    x = armaSim(model = list(ar = c(-0.4, 0.1, 0, 0, 0.1), 
      d = 0, ma = 0)) 
    # Estimate the Parameters:
    fit52 = armaFit(x ~ arima(5, 0, 2),
      fixed = c(NA, NA, NA, NA, NA, NA, NA, 0), transform.pars = FALSE)
    fit51 = armaFit(x ~ arima(5, 0, 2),
      fixed = c(NA, NA, NA, NA, NA,  0, NA, 0), transform.pars = FALSE)
    fit41 = armaFit(x ~ arima(5, 0, 2),
      fixed = c(NA, NA, NA,  0, NA,  0, NA, 0), transform.pars = FALSE)
    fit31 = armaFit(x ~ arima(5, 0, 2),
      fixed = c(NA, NA,  0,  0, NA,  0, NA, 0), transform.pars = FALSE)
    fit30 = armaFit(x ~ arima(5, 0, 2),
      fixed = c(NA, NA,  0,  0, NA,  0,  0, 0), transform.pars = FALSE)
    fit20 = armaFit(x ~ arima(5, 0, 2),
      fixed = c(NA,  0,  0,  0, NA,  0,  0, 0), transform.pars = FALSE)
    fit10 = armaFit(x ~ arima(5, 0, 2),
      fixed = c(NA,  0,  0,  0,  0,  0,  0, 0), transform.pars = FALSE)
    # Print and Plot Results:  
    p = c(5, 5, 4, 3, 3, 2, 1)
    q = c(2, 1, 1, 1, 0, 0, 0)
    sigma2 = c(fit52$sigma2, fit51$sigma2, fit41$sigma2, 
      fit31$sigma2, fit30$sigma2, fit20$sigma2, fit10$sigma2)
    aic = c(fit52$aic, fit51$aic, fit41$aic, fit31$aic, fit30$aic,
      fit20$aic, fit10$aic)
    plot(rev(p+q), rev(aic), xlab = "p+q", ylab = "aic", 
      main = "AIC Statistics")
    cbind(p, q, sigma2, aic)
	###
	
  
	# Part II: Step 2 - Parameter Estimation
    # Settings:
    par (mfrow = c(2, 1), cex = 0.7)
    set.seed(732)
    # Simulate a Subset AR(5) process ...
    x = armaSim(n = 2000, n.start = 100, rand.gen = rnorm,
      model = list(ar = c(-0.40, 0.10, 0, 0, 0.10), d = 0, ma = 0))
    ts.plot(x, main = "AR(5) Subset")
    # Run Maximum Loglikelihhood Estimation ...
    # Estimate Parameters:
    x = x-mean(x)
    fit = armaFit(x ~ arima(5, 0, 0), include.mean = FALSE, 
      fixed = c(NA, NA, 0, 0, NA))
    # Print Result:
    print(fit)
    # Simulate an ARMA(2,1) process ...
    x = armaSim(n = 2000, n.start = 100, rand.gen = rnorm,
      model = list(ar = c(0.4, -0.2), d = 0, ma = 0.1))
    ts.plot(x, main = "ARMA(2,1)")
    # Run Maximum Loglikelihhood Estimation ...
    # Estimate Parameters:
    x = x - mean(x)
    fit = armaFit(x ~ arma(2, 1), include.mean = FALSE)
    # Print Result:
    print(fit)
    ###
   
   
 	# Part II: Step 3 - Diagnostic Analysis
    # Settings:
    par(mfrow = c(2, 2), cex = 0.7)
    set.seed(732)  
    # Simulate a subset ARIMA process ...
    x = armaSim(model=list(ar=c(-0.4, 0.1, 0, 0, 0.1), d = 0, ma = 0), 
      n = 1000)
    # Perform Maximum Loglikelihhood Estimation ...
    fit = armaFit(x ~ arima(5, 0, 0), include.mean = TRUE, 
      fixed = c(NA, NA,  0,  0, NA, NA))
    print(fit) 
    # Compute diagnostics for ARIMA model ...
    summary(fit) 
    ###   

   
    # Part II: Step 4 - Forecasts   
    # Settings:
    set.seed(732)
    # Simulate a Subset AR(5) Process:
    x = armaSim(model = 
      list(ar = c(-0.4, 0.1, 0, 0, 0.1), d = 0, m = 0), n = 1000)  
    # Estimate the Parameters:
    fit = armaFit(x ~ arima(5, 0, 0), include.mean = TRUE,
      fixed = c(NA, NA,  0,  0, NA, NA), transform.pars = FALSE)
    print(fit)
    # Forecast 5 Steps ahead:
    predict(fit, n.ahead = 5)
    ###
    
    
################################################################################

