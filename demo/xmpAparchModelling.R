
#
# Example: 
#   Modelling asymmetric power arch processing
#
# Descriptions:
#   These examples show how to simulate and to analyse APARCH 
#   Models:
#   1 Aparch Model Simulations with several kind of innovations
#     formerly: xmpAparchInnovations
#   2 Aparch Model Simulations, continued
#     formerly: xmpAparchSimulation
#   3 Aparch Model Simulations and Parameter Estimations
#     formerly: xmpAparchEstimation
#   4 Aparch NYSE Composite Index Estimation
#     formerly: xmpAparchNYSERES
#
# Author:
#   (C) 2002, Diethelm Wuertz, GPL
#


################################################################################
# 1 Aparch Model Simulations with Several Kind of Innovations:


    # Settings:
    par(mfrow = c(3, 2), cex = 0.6)
    ###
  
    # Simulate GARCH(1,1) Normal Innovations
    set.seed(4711)
    x = aparchSim(model = list(omega = 1.02e-6, 
        alpha = 0.1, gamma = 0, alpha.lags = 1, beta = 0.6, 
        beta.lags = 1, delta = 2), innov = rnorm(2000), 
        start.innov = rnorm(100))   
    # Plot:
    ts.plot(x, ylim = c(-0.015, 0.015), 
        main = "GARCH(1,1) Normal Innovations")
    acf(abs(x), lag.max = 12)
    var(x)      
    ###
    
    # Simulate GARCH(1,1) Student t(4) Innovations
    set.seed(1858)
    x = aparchSim(model = list(omega = 0.45e-6, 
        alpha = 0.1, gamma = 0, alpha.lags = 1, beta = 0.6, 
        beta.lags = 1, delta = 2), innov = rt(2000, df = 4), 
        start.innov = rt(100, df = 4))      
    # Plot:
    ts.plot(x, ylim = c(-0.015, 0.015), 
        main = "GARCH(1,1) - Student t(4) Innovations")
    acf(abs(x), lag.max = 12)
    var(x)      
    ### 
    
    # Simulate GARCH(1,1) 1.9 Stable Innovations
    set.seed(1245)
    x = aparchSim(model = list(omega = 0.3e-6, 
        alpha = 0.1, gamma = 0, alpha.lags = 1, beta = 0.6, 
        beta.lags = 1, delta = 2), innov = rsymstb(2000, 
        alpha = 1.95), start.innov = rsymstb(100, alpha = 1.95))        
    # Plot:
    ts.plot(x, ylim = c(-0.015, 0.015), 
        main = "GARCH(1,1) - 1.95 Stable Innovations")
    acf(abs(x), lag.max = 12)
    var(x)
    ###
    
    
################################################################################
# 2 Aparch Model Simulations, continued:


    # Settings:
    par(mfrow = c(3, 2), cex = 0.6)

    # Bollerslev's GARCH(1,2) model
    # h(t)^2 = omega + alpha(i)*eps(t-i)^2 + beta(i)*h(t-i)^2 
    # Simulate:
    set.seed(4661)
    x = aparchSim(model = list(omega = 0.8e-5, 
      alpha = 0.1, gamma = 0, alpha.lags = 1, 
      beta = c(0.5, 0.2), beta.lags = c(1, 2), delta = 2), 
      innov = rnorm(5000), start.innov = rnorm(1000))
    # Plot:
    ts.plot(x, ylim = c(-0.04, 0.04), main = "GARCH(1,2) Model")
    acf(abs(x), lag.max = 20)
    var(x)
    ###

    # Taylor-Schwert GARCH(p,q) model
    # h(t) = omega + alpha(i)*|eps(t-i)| + beta(j)*h(t-j) 
    # Simulate:
    set.seed(5066)
    x = aparchSim(model=list(omega = 1.4e-3, 
      alpha = 0.1, gamma = 0, alpha.lags = 1, 
      beta = c(0.5, 0.2), beta.lags = c(1, 2), delta = 1), 
      innov = rnorm(5000), start.innov = rnorm(1000))
    # Plot:
    ts.plot(x, ylim = c(-0.04, 0.04), main = "TS GARCH(1,2) Model")
    acf(abs(x), lag.max = 20)
    var(x)
    ###
    
    # Glosten's et al. GJR(p,q) model
    # h(t)^2 = omega + alpha'(i)*eps(t-i)^2 + beta(j)*h(t-j)^2 + 
    #          gamma'(i)*S(i)*eps(t-i)
    #   S(i) = 1 if eps(t-i)>0 ; S(i) = 0 otherwise
    #   0 <= gamma(i) < 1
    # Simulate:
    set.seed(8531)
    x = aparchSim(model=list(omega = 0.8e-5, 
      alpha = 0.1, gamma = 0.2, alpha.lags = 1,
      beta = c(0.5, 0.2), beta.lags = c(1, 2), delta = 2), 
      innov = rnorm(5000), start.innov = rnorm(1000))
    # Plot:
    ts.plot(x, ylim = c(-0.04, 0.04), main = "GJR(1,2) Model")
    acf(abs(x), lag.max = 20)
    var(x)
    ###



################################################################################
# 3 Aparch Model Simulations and Parameter Estimations:


    # Settings:
    par(mfrow = c(3, 2), cex = 0.6)
    ###

    # Bollerslev's t-GARCH(1,1):
    # Simulate:
    set.seed(5066)
    x = aparchSim(model=list(omega = 1.0e-6, 
      alpha = 0.1, gamma = 0, alpha.lags = 1, 
      beta = 0.8, beta.lags = 1, delta = 2), 
      innov = rt(5000, df = 10), start.innov = rt(1000, df = 10))
    ts.plot(x, main = "GARCH(1,1) - Student-t")
    acf(abs(x), lag.max = 20)           
    # Estimate:
    aparchFit(x = x, 
      order = list(alpha.lags = 1, beta.lags = 1, delta = 2),
      opt = list(gamma = FALSE, delta = FALSE, disparm = TRUE),
      distribution = "t", disparm = 4, doprint = FALSE) 
    ###
    
	# Taylor-Schwert TS-GARCH(1,5[1,5]) subset model, delta=1:  
    # Simulate:
    set.seed(5066)
    x = aparchSim(model=list(omega = 1.0e-6, 
      alpha = 0.1, gamma = 0, alpha.lags = 1, 
      beta = c(0.5, 0.3), beta.lags = c(1, 5), delta = 1),
      innov = rnorm(5000), start.innov = rnorm(1000))
    ts.plot(x, main = "Subset TS GARCH(1,5) - Mormal")
    acf(abs(x), lag.max = 20)           
    # Estimate:
    aparchFit(x, 
      order = list(alpha.lags = 1, beta.lags = c(1, 5), delta = 1),
      opt = list(gamma=FALSE, delta=FALSE, disparm=FALSE),
      doprint = FALSE)        
    ###
        
	# Asymmetric GJR(1,1):                  
    # Simulate:
    set.seed(5066)
    x = aparchSim(model = list(omega = 1.0e-6, 
      alpha = 0.1, gamma = 0.2, alpha.lags = 1, 
      beta = 0.6, beta.lags = 1, delta = 2), 
      innov = rnorm(1000), start.innov = rnorm(100))
    ts.plot(x, main = "GJR(1,1) - Normal")
    acf(abs(x), lag.max = 20)           
    # Estimate:         
    aparchFit(x=x, 
      order = list(alpha.lags = 1, beta.lags = 1, delta = 2), 
      opt = list(gamma = TRUE, delta = FALSE, disparm = FALSE),
      doprint = FALSE)

    
################################################################################
# 4 Aparch NYSE Composite Index Estimation:


	# Settings: 
    par(mfrow = c(4, 2), cex = 0.6)
    set.seed(1257)          
    data(nyseres)
    s = 0.10
    ###
    
	# Model Data:
    nyseres = nyseres[abs(nyseres) < s]
    ts.plot(nyseres, ylim = c(-s, s), main = "NYSE- log Returns")
    loglik = 0
    variance = var(nyseres)
    skewness = skewness(nyseres)
    kurtosis = kurtosis(nyseres) 
    scaling = scalinglawPlot(nyseres, doplot = FALSE)$exponent
    ###
        
	# Model 1: Bollerslev GARCH(1,1):     
    # [1] -2.898508e+04  1.012818e-06  7.579743e-02  0.000000e+00  
    #      9.124155e-01  2.000000e+00  1.000000e+00
    fit = aparchFit(x = nyseres, 
      order = list(alpha.lags = 1, beta.lags = 1, delta = 2), 
      opt = list(gamma = FALSE, delta = FALSE, disparm = FALSE),
      distribution = c("norm", "t", "symstb"), disparm = c(1, 4, 1.9), 
      n.cond = NULL, doprint = FALSE, method = "Nelder-Mead") 
    x1 = aparchSim(model = list(omega = fit$omega,
      alpha = fit$alpha, gamma = fit$gamma, alpha.lags = 1, 
      beta = fit$beta, beta.lags = 1, delta = fit$delta), 
      innov = rnorm(length(nyseres)), start.innov = rnorm(5000))
    ts.plot(x = x1, ylim = c(-s, s), 
      main = paste("GARCH(1,1) - llh: ",
      as.character(floor(fit$value))))
    loglik = c(loglik, fit$value)
    variance = c(variance, var(x1))
    skewness = c(skewness, skewness(x1))
    kurtosis = c(kurtosis, kurtosis(x1))
    scaling = c(scaling, scalinglawPlot(x1, doplot = FALSE)$exponent)
    Continue = readline("Press any key > ") 
	###
	
	# Model 2: delta-GARCH(1,1):       
    # [1] -2.898861e+04  6.414045e-06  8.369417e-02  0.000000e+00  
    #      9.141348e-01  1.621929e+00  1.000000e+00
    fit = aparchFit(x = nyseres, 
      order = list(alpha.lags = 1, beta.lags = 1, delta = 2), 
      opt = list(gamma = FALSE, delta = TRUE, disparm = FALSE),
      distribution = c("norm", "t", "symstb"), disparm = c(1, 4, 1.9), 
      n.cond = NULL, doprint = FALSE, method = "Nelder-Mead")
    x2 = aparchSim(model = list(omega = fit$omega,
      alpha = fit$alpha, gamma = fit$gamma, alpha.lags = 1, 
      beta = fit$beta, beta.lags = 1, delta = fit$delta), 
      innov = rnorm(length(nyseres)), start.innov = rnorm(5000))
    ts.plot(x = x2, ylim = c(-s, s), 
      main = paste("Delta-GARCH(1,1) - llh: ",
      as.character(floor(fit$value))))
    loglik = c(loglik, fit$value)
    variance = c(variance, var(x2))
    skewness = c(skewness, skewness(x2))
    kurtosis = c(kurtosis, kurtosis(x2))
    scaling = c(scaling, scalinglawPlot(x2, doplot = FALSE)$exponent)   
    ### 
        
	# Model 3: asymmetric-GARCH(1,1): 
    # [1] -2.905463e+04  1.173995e-06  6.458562e-02  3.145372e-01  
    #      9.156142e-01  2.000000e+00  1.000000e+00
    fit = aparchFit(x = nyseres, 
      order = list(alpha.lags = 1, beta.lags = 1, delta = 2), 
      opt = list(gamma = TRUE, delta = FALSE, disparm = FALSE),
      distribution = c("norm", "t", "symstb"), disparm = c(1, 4, 1.9), 
      n.cond = NULL, doprint = FALSE, method = "Nelder-Mead") 
    x3 = aparchSim(model = list(omega = fit$omega,
      alpha = fit$alpha, gamma = fit$gamma, alpha.lags = 1, 
      beta = fit$beta, beta.lags = 1, delta = fit$delta), 
      innov = rnorm(length(nyseres)), start.innov = rnorm(5000))
    ts.plot(x = x3, ylim = c(-s,s), 
      main = paste("asymmetric-GARCH(1,1) - llh: ",
      as.character(floor(fit$value))))
    loglik = c(loglik, fit$value)
    variance = c(variance, var(x3))
    skewness = c(skewness, skewness(x3))
    kurtosis = c(kurtosis, kurtosis(x3))
    scaling = c(scaling, scalinglawPlot(x3, doplot = FALSE)$exponent)
    ###
                    
	# Model 4: delta-asymmetric-GARCH(1,1):    
    # [1] -2.907322e+04  3.534699e-05  7.570291e-02  4.236941e-01  
    #      9.196998e-01  1.316301e+00  1.000000e+00
    fit = aparchFit(x = nyseres, 
      order = list(alpha.lags = 1, beta.lags = 1, delta = 2), 
      opt = list(gamma = TRUE, delta = TRUE, disparm = FALSE),
      distribution = c("norm", "t", "symstb"), disparm = c(1, 4, 1.9), 
      n.cond = NULL, doprint = FALSE, method = "Nelder-Mead")
    x4 = aparchSim(model = list(omega = fit$omega,
      alpha = fit$alpha, gamma = fit$gamma, alpha.lags = 1, 
      beta = fit$beta, beta.lags = 1, delta = fit$delta), 
      innov = rnorm(length(nyseres)), start.innov = rnorm(5000))
    ts.plot(x = x4, ylim = c(-s, s), 
      main = paste("Delta-asym-GARCH(1,1) - llh: ",
      as.character(floor(fit$value))))
    loglik = c(loglik, fit$value)
    variance = c(variance, var(x4))
    skewness = c(skewness, skewness(x4))
    kurtosis = c(kurtosis, kurtosis(x4))
    scaling = c(scaling, scalinglawPlot(x4, doplot = FALSE)$exponent)
    ###
                    
	# Model 5: t-GARCH(1,1): 
    # [1] -2.921002e+04  5.354942e-07  3.990274e-02  0.000000e+00  
    #      9.341133e-01  2.000000e+00  7.138929e+00
    fit = aparchFit(x = nyseres, 
      order = list(alpha.lags = 1, beta.lags = 1, delta = 2), 
      opt = list(gamma = FALSE, delta = FALSE, disparm = TRUE),
      distribution = "t", disparm = 4, 
      n.cond = NULL, doprint = FALSE, method = "Nelder-Mead")
    x5 = aparchSim(model = list(omega = fit$omega,
      alpha = fit$alpha, gamma = fit$gamma, alpha.lags = 1, 
      beta = fit$beta, beta.lags = 1, delta = fit$delta), 
      innov = rt(length(nyseres), df = fit$disparm), 
      start.innov = rt(5000, df = fit$disparm))
    ts.plot(x = x5, ylim = c(-s,s), 
      main = paste("Student-t-GARCH(1,1) - llh: ",
      as.character(floor(fit$value))))
    loglik = c(loglik, fit$value)
    variance = c(variance, var(x5))
    skewness = c(skewness, skewness(x5))
    kurtosis = c(kurtosis, kurtosis(x5))
    scaling = c(scaling, scalinglawPlot(x5, doplot = FALSE)$exponent)
	###        
            
	# Model 6: asymmetric-t-GARCH(1,1): 
    # [1] -2.925041e+04  6.587438e-07  4.003464e-02  3.000886e-01  
    #      9.295025e-01  2.000000e+00  7.600815e+00
    fit = aparchFit(x = nyseres, 
      order = list(alpha.lags = 1, beta.lags = 1, delta = 2), 
      opt = list(gamma = TRUE, delta = FALSE, disparm = TRUE),
      distribution = "t", disparm = 4, 
      n.cond = NULL, doprint = FALSE, method = "BFGS")
    x6 = aparchSim(model = list(omega = fit$omega,
      alpha = fit$alpha, gamma = fit$gamma, alpha.lags = 1, 
      beta = fit$beta, beta.lags = 1, delta = fit$delta), 
      innov = rt(length(nyseres), df = fit$disparm), 
      start.innov = rt(5000, df = fit$disparm))
    ts.plot(x = x6, ylim = c(-s,s), 
      main = paste("asymmetric-t-GARCH(1,1) - llh: ",
      as.character(floor(fit$value))))
    loglik = c(loglik, fit$value)
    variance = c(variance, var(x6))
    skewness = c(skewness, skewness(x6))
    kurtosis = c(kurtosis, kurtosis(x6))
    scaling = c(scaling, scalinglawPlot(x6, doplot = FALSE)$exponent)
	###
	            
	# Model 7: delta-asymmetric-t-GARCH(1,1): 
    # [1] -2.926201e+04  1.435893e-04  5.902100e-02  4.701446e-01  
    #      9.268403e-01  9.959521e-01  6.568893e+00
    fit = aparchFit(x = nyseres, 
      order = list(alpha.lags = 1, beta.lags = 1, delta = 2), 
      opt = list(gamma = TRUE, delta = TRUE, disparm = TRUE),
      distribution = "t", disparm = 4, 
      n.cond = NULL, doprint = FALSE, method = "Nelder-Mead")
    x7 = aparchSim(model = list(omega = fit$omega,
      alpha = fit$alpha, gamma = fit$gamma, alpha.lags = 1, 
      beta = fit$beta, beta.lags = 1, delta = fit$delta), 
      innov = rt(length(nyseres), df = fit$disparm), 
      start.innov = rt(5000, df = fit$disparm))
    ts.plot(x = x7, ylim = c(-s,s), 
      main = paste("Delta-asym-t-GARCH(1,1) - llh: ",
      as.character(floor(fit$value))))
    loglik = c(loglik, fit$value)
    variance = c(variance, var(x7))
    skewness = c(skewness, skewness(x7))
    kurtosis = c(kurtosis, kurtosis(x7))
    scaling = c(scaling, scalinglawPlot(x7, doplot = FALSE)$exponent)
	###
	            
	# Summary:
    data.frame(cbind(loglik, variance, skewness, kurtosis, scaling))
    ###   
    
        
################################################################################

