
#
# Example: 
#   Analyse a time series process using the ARIMA approach 
#
# Description:
#   Part I: Simulate ARIMA time series processes
#   Part II: Investigate Correlations
#   PART III: Compute Roots of AR / MA polynomials
#   PART IV: Compute ACF of the true time series model
#   Part V: Identify the Order of an ARIMA process 
# Author:
#	(C) 2002, Diethelm Wuertz, GPL
#


# ------------------------------------------------------------------------------


################################################################################
## PART I: Simulate ARIMA time series processes


   # Settings:
    par(mfrow = c(4, 3), cex = 0.6)
    set.seed(471)
    
   # AR(1) and AR(2) Examples:  
   x = armaSim(model = list(ar = 0.5, d = 0, ma = 0), n = 1000)
     ts.plot(x, main = "AR(1): +0.5")
     acf(x, lag.max = 12)
     acf(x, type = "partial", lag.max = 12)
     Continue = readline("Press any key > ")  
   x = armaSim(model = list(ar = -0.5, d = 0, ma = 0), n = 1000)
     ts.plot(x, main = "AR(1): -0.5")
     acf(x, lag.max = 12)
     acf(x, type = "partial", lag.max = 12)
     Continue = readline("Press any key > ") 
   x = armaSim(model = list(ar = c(0.5, 0.3), d = 0, ma = 0), n = 1000)
     ts.plot(x, main = "AR(2): +0.5, +0.3")
     acf(x, lag.max = 12)
     acf(x, type = "partial", lag.max = 12)
     Continue = readline("Press any key > ")     
   x = armaSim(model=list(ar = c(-0.5, 0.3), d = 0, ma = 0), n = 1000)
     ts.plot(x, main = "AR(2): -0.5, +0.3")
     acf(x, lag.max = 12)
     acf(x, type = "partial",lag.max = 12)
    
     
################################################################################
## PART II: Investigate Correlations

  
   # Settings:  
   par(mfrow = c(4, 3), cex = 0.6)
   set.seed(471)
    
   # MA and ARMA Examples:
   x = armaSim(n = 1000, model = list(ar = 0, d = 0, ma = 0.8))
     ts.plot(x, main = "MA(1): +0.8")
     acf(x, lag.max = 12)
     acf(x, type = "partial", lag.max = 12)
     Continue = readline("Press any key > ")
   x = armaSim(n = 1000, model = list(ar = 0, d = 0, ma = -0.8))
     ts.plot(x, main = "MA(1): -0.8")
     acf(x, lag.max = 12)
     acf(x, type = "partial", lag.max = 12)
     Continue = readline("Press any key > ")      
   x = armaSim(n = 1000, model = list(ar = 0.5, d = 0, ma = 0.8))
     ts.plot(x, main = "ARMA(1,1): +0.5, +0.8")
     acf(x, lag.max = 12)
     acf(x, type = "partial", lag.max = 12)
     Continue = readline("Press any key > ")    
   x = armaSim(n = 1000, model = list(ar = -0.5, d = 0, ma = 0.8))
     ts.plot(x, main = "ARMA(1,1): -0.5, +0.8")
     acf(x, lag.max = 12)
     acf(x, type = "partial", lag.max = 12)


################################################################################
## Part III: Compute Roots of AR / MA polynomials
  
 
   # Settings:
   par(mfrow = c(2, 2))  
     
   # Roots - Example 1: 
   armaRoots(c( 1.6, 0.4))
   armaRoots(c(-1.6, 0.4))   
   Continue <- readline("Press any key > ")
   
   # Roots - Example 2:
   armaRoots(c( 0.5, -0.1))
   armaRoots(c(-0.5, -0.1))
   Continue <- readline("Press any key > ")
   
   # Roots - Example 3:
   armaRoots(c( 0.5, -0.9,  0.1, 0.5))
   armaRoots(c(-0.5, -0.9, -0.1, 0.5))

   
################################################################################
## Part IV: Compute ACF of the true time series model

   
   # Settings:
   par(mfcol = c(3, 2), cex = 0.6)
    
   # First Model:
   model = list(ar = c(+0.5, 0.3), d = 0, ma = 0)
   x = armaSim(model, n = 1000)
   ts.plot(x, main = "AR(2): +0.5, 0.3")
   armaTrueacf(model, lag.max = 12)
   armaTrueacf(model, lag.max = 12, type = "partial")
   Continue = readline("Press any key > ") 
    
   # Second Model:
   model = list(ar = c(-0.5, 0.3), d = 0, ma = 0)
   x = armaSim(model, n = 1000)
   ts.plot(x, main = "AR(2): -0.5, 0.3")
   armaTrueacf(model, lag.max = 12)
   armaTrueacf(model, lag.max = 12, type = "partial")
    
    
################################################################################
## Part V: Identify the Order of an ARIMA process   
 

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
   
  
################################################################################

