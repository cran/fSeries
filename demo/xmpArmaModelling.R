
#
# Example: 
#   Model a time series process using the ARIMA approach
#
# Description:
#   Part I: Model Selection
#   Part II: Parameter Estimation
#   Part III: Diagnostic Analysis
#   Part IV: Forecasts
#
# Author:
#	(C) 2002, Diethelm Wuertz, GPL
#


# ------------------------------------------------------------------------------


################################################################################
## Part I: Model Selection


   # Settings:
   par (mfrow = c(2, 1))
   set.seed(4732)
    
   # Simulate time series:
   x = armaSim(model = list(ar = c(-0.4, 0.1, 0, 0, 0.1), 
     d = 0, ma = 0))
    
   # Estimate the Parameters:
   fit52 = armaFit(x ~ arima(5, 0, 2),
     fixed=c(NA, NA, NA, NA, NA, NA, NA, 0), transform.pars = FALSE)
   fit51 = armaFit(x ~ arima(5, 0, 2),
     fixed=c(NA, NA, NA, NA, NA,  0, NA, 0), transform.pars = FALSE)
   fit41 = armaFit(x ~ arima(5, 0, 2),
     fixed=c(NA, NA, NA,  0, NA,  0, NA, 0), transform.pars = FALSE)
   fit31 = armaFit(x ~ arima(5, 0, 2),
     fixed=c(NA, NA,  0,  0, NA,  0, NA, 0), transform.pars = FALSE)
   fit30 = armaFit(x ~ arima(5, 0, 2),
     fixed=c(NA, NA,  0,  0, NA,  0,  0, 0), transform.pars = FALSE)
   fit20 = armaFit(x ~ arima(5, 0, 2),
     fixed=c(NA,  0,  0,  0, NA,  0,  0, 0), transform.pars = FALSE)
   fit10 = armaFit(x ~ arima(5, 0, 2),
     fixed=c(NA,  0,  0,  0,  0,  0,  0, 0), transform.pars = FALSE)

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

  
################################################################################
## Part II: Parameter Estimation


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
   Continue = readline("Press any key > ") 
    
   # Simulate an ARMA(2,1) process ...
   x = armaSim(n = 2000, n.start = 100, rand.gen = rnorm,
     model = list(ar = c(0.4, -0.2), d = 0, ma = 0.1))
   ts.plot(x, main = "ARMA(2,1)")

   # Run Maximum Loglikelihhood Estimation ...
   # Estimate Parameters:
   x = x-mean(x)
   fit = armaFit(x ~ arma(2, 1), include.mean = FALSE)
   # Print Result:
   print(fit)
   
   
################################################################################
## Part III: Diagnostic Analysis
   
   
   # Settings:
   par(mfrow = c(2, 2), cex = 0.7)
   set.seed(732)
    
   # Simulate a subset ARIMA process ...
   x = armaSim(model=list(ar=c(-0.4, 0.1, 0, 0, 0.1), d = 0, ma = 0), 
     n=1000)

   # Perform Maximum Loglikelihhood Estimation ...
   fit = armaFit(x ~ arima(5, 0, 0), include.mean = TRUE, 
     fixed = c(NA, NA,  0,  0, NA, NA))
   print(fit)
  
   # Compute diagnostics for ARIMA model ...
   summary(fit)    

   
################################################################################
## Part IV: Forecasts   
   
   
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
    
    
################################################################################

