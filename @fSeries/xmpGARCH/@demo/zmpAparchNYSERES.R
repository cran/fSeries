
#
# Example: Aparch NYSE Composite Index Estimation
#

# Settings:
		par(mfrow=c(4,2))
		set.seed(1257)			
	    data(nyseres)
	    s <- 0.10
	    nyseres <- nyseres[abs(nyseres)< s]
	    ts.plot(nyseres, ylim=c(-s,s), main="NYSE- log Returns")
	    loglik <- 0
	    variance <- var(nyseres)
	    skewness <- skewness(nyseres)
	    kurtosis <- kurtosis(nyseres) 
	    scaling <- scalinglaw(nyseres, doplot=FALSE)$exponent
	    Continue <- readline("Press any key > ") 
	    
# 1 - Bollerslev GARCH(1,1):  
		# [1] -2.898508e+04  1.012818e-06  7.579743e-02  0.000000e+00  
		#      9.124155e-01  2.000000e+00  1.000000e+00
		fit <- aparchFit(x=nyseres, 
			order=list(alpha.lags=1, beta.lags=1, delta=2), 
    		opt=list(gamma=FALSE, delta=FALSE, disparm=FALSE),
			distribution=c("norm", "t", "symstb"), disparm=c(1, 4, 1.9), 
			n.cond=NULL, doprint=FALSE, method="Nelder-Mead")	
		x1 <- aparchSim(model=list(omega=fit$omega,
  	    	alpha=fit$alpha, gamma=fit$gamma, alpha.lags=1, 
  	    	beta=fit$beta, beta.lags=1, delta=fit$delta), 
  	    	innov=rnorm(length(nyseres)), start.innov=rnorm(5000))
		ts.plot(x=x1, ylim=c(-s,s), main=paste("GARCH(1,1) - llh: ",
		 	as.character(floor(fit$value))))
		loglik <- c(loglik, fit$value)
		variance <- c(variance, var(x1))
		skewness <- c(skewness, skewness(x1))
	    kurtosis <- c(kurtosis, kurtosis(x1))
	    scaling <- c(scaling, scalinglaw(x1, doplot=FALSE)$exponent)
	    Continue <- readline("Press any key > ") 

# 2 - delta-GARCH(1,1): 
		# [1] -2.898861e+04  6.414045e-06  8.369417e-02  0.000000e+00  
		#      9.141348e-01  1.621929e+00  1.000000e+00
	    fit <- aparchFit(x=nyseres, 
			order=list(alpha.lags=1, beta.lags=1, delta=2), 
    		opt=list(gamma=FALSE, delta=TRUE, disparm=FALSE),
			distribution=c("norm", "t", "symstb"), disparm=c(1, 4, 1.9), 
			n.cond=NULL, doprint=FALSE, method="Nelder-Mead")
		x2 <- aparchSim(model=list(omega=fit$omega,
  	    	alpha=fit$alpha, gamma=fit$gamma, alpha.lags=1, 
  	    	beta=fit$beta, beta.lags=1, delta=fit$delta), 
  	    	innov=rnorm(length(nyseres)), start.innov=rnorm(5000))
		ts.plot(x=x2, ylim=c(-s,s), main=paste("Delta-GARCH(1,1) - llh: ",
			as.character(floor(fit$value))))
		loglik <- c(loglik, fit$value)
		variance <- c(variance, var(x2))
		skewness <- c(skewness, skewness(x2))
	    kurtosis <- c(kurtosis, kurtosis(x2))
	    scaling <- c(scaling, scalinglaw(x2, doplot=FALSE)$exponent)	
	    Continue <- readline("Press any key > ") 	

# 3 - asymmetric-GARCH(1,1): 
		# [1] -2.905463e+04  1.173995e-06  6.458562e-02  3.145372e-01  
		#      9.156142e-01  2.000000e+00  1.000000e+00
		fit <- aparchFit(x=nyseres, 
			order=list(alpha.lags=1, beta.lags=1, delta=2), 
    		opt=list(gamma=TRUE, delta=FALSE, disparm=FALSE),
			distribution=c("norm", "t", "symstb"), disparm=c(1, 4, 1.9), 
			n.cond=NULL, doprint=FALSE, method="Nelder-Mead")	
		x3 <- aparchSim(model=list(omega=fit$omega,
  	    	alpha=fit$alpha, gamma=fit$gamma, alpha.lags=1, 
  	    	beta=fit$beta, beta.lags=1, delta=fit$delta), 
  	    	innov=rnorm(length(nyseres)), start.innov=rnorm(5000))
		ts.plot(x=x3, ylim=c(-s,s), main=paste("asymmetric-GARCH(1,1) - llh: ",
			as.character(floor(fit$value))))
		loglik <- c(loglik, fit$value)
		variance <- c(variance, var(x3))
		skewness <- c(skewness, skewness(x3))
	    kurtosis <- c(kurtosis, kurtosis(x3))
	    scaling <- c(scaling, scalinglaw(x3, doplot=FALSE)$exponent)
	    Continue <- readline("Press any key > ") 
			
# 4 - delta-asymmetric-GARCH(1,1): 
		# [1] -2.907322e+04  3.534699e-05  7.570291e-02  4.236941e-01  
		#      9.196998e-01  1.316301e+00  1.000000e+00
	    fit <- aparchFit(x=nyseres, 
			order=list(alpha.lags=1, beta.lags=1, delta=2), 
    		opt=list(gamma=TRUE, delta=TRUE, disparm=FALSE),
			distribution=c("norm", "t", "symstb"), disparm=c(1, 4, 1.9), 
			n.cond=NULL, doprint=FALSE, method="Nelder-Mead")
		x4 <- aparchSim(model=list(omega=fit$omega,
  	    	alpha=fit$alpha, gamma=fit$gamma, alpha.lags=1, 
  	    	beta=fit$beta, beta.lags=1, delta=fit$delta), 
  	    	innov=rnorm(length(nyseres)), start.innov=rnorm(5000))
		ts.plot(x=x4, ylim=c(-s,s), main=paste("Delta-asym-GARCH(1,1) - llh: ",
			as.character(floor(fit$value))))
		loglik <- c(loglik, fit$value)
		variance <- c(variance, var(x4))
		skewness <- c(skewness, skewness(x4))
	    kurtosis <- c(kurtosis, kurtosis(x4))
	    scaling <- c(scaling, scalinglaw(x4, doplot=FALSE)$exponent)
	    Continue <- readline("Press any key > ") 
			
# 5 - t-GARCH(1,1): 
		# [1] -2.921002e+04  5.354942e-07  3.990274e-02  0.000000e+00  
		#      9.341133e-01  2.000000e+00  7.138929e+00
	    fit <- aparchFit(x=nyseres, 
			order=list(alpha.lags=1, beta.lags=1, delta=2), 
    		opt=list(gamma=FALSE, delta=FALSE, disparm=TRUE),
			distribution="t", disparm=4, 
			n.cond=NULL, doprint=FALSE, method="Nelder-Mead")
		x5 <- aparchSim(model=list(omega=fit$omega,
  	    	alpha=fit$alpha, gamma=fit$gamma, alpha.lags=1, 
  	    	beta=fit$beta, beta.lags=1, delta=fit$delta), 
  	    	innov=rt(length(nyseres), df=fit$disparm), 
  	    	start.innov=rt(5000, df=fit$disparm))
		ts.plot(x=x5, ylim=c(-s,s), main=paste("Student-t-GARCH(1,1) - llh: ",
			as.character(floor(fit$value))))
		loglik <- c(loglik, fit$value)
		variance <- c(variance, var(x5))
		skewness <- c(skewness, skewness(x5))
	    kurtosis <- c(kurtosis, kurtosis(x5))
	    scaling <- c(scaling, scalinglaw(x5, doplot=FALSE)$exponent)
	    Continue <- readline("Press any key > ") 
			
# 6 - asymmetric-t-GARCH(1,1): 
		# [1] -2.925041e+04  6.587438e-07  4.003464e-02  3.000886e-01  
		#      9.295025e-01  2.000000e+00  7.600815e+00
	    fit <- aparchFit(x=nyseres, 
			order=list(alpha.lags=1, beta.lags=1, delta=2), 
    		opt=list(gamma=TRUE, delta=FALSE, disparm=TRUE),
			distribution="t", disparm=4, 
			n.cond=NULL, doprint=FALSE, method="BFGS")
		x6 <- aparchSim(model=list(omega=fit$omega,
  	    	alpha=fit$alpha, gamma=fit$gamma, alpha.lags=1, 
  	    	beta=fit$beta, beta.lags=1, delta=fit$delta), 
  	    	innov=rt(length(nyseres), df=fit$disparm), 
  	    	start.innov=rt(5000, df=fit$disparm))
		ts.plot(x=x6, ylim=c(-s,s), main=paste("asymmetric-t-GARCH(1,1) - llh: ",
			as.character(floor(fit$value))))
		loglik <- c(loglik, fit$value)
		variance <- c(variance, var(x6))
		skewness <- c(skewness, skewness(x6))
	    kurtosis <- c(kurtosis, kurtosis(x6))
	    scaling <- c(scaling, scalinglaw(x6, doplot=FALSE)$exponent)
	    Continue <- readline("Press any key > ") 
			
# 7 - delta-asymmetric-t-GARCH(1,1): 
		# [1] -2.926201e+04  1.435893e-04  5.902100e-02  4.701446e-01  
		#      9.268403e-01  9.959521e-01  6.568893e+00
	    fit <- aparchFit(x=nyseres, 
			order=list(alpha.lags=1, beta.lags=1, delta=2), 
    		opt=list(gamma=TRUE, delta=TRUE, disparm=TRUE),
			distribution="t", disparm=4, 
			n.cond=NULL, doprint=FALSE, method="Nelder-Mead")
		x7 <- aparchSim(model=list(omega=fit$omega,
  	    	alpha=fit$alpha, gamma=fit$gamma, alpha.lags=1, 
  	    	beta=fit$beta, beta.lags=1, delta=fit$delta), 
  	    	innov=rt(length(nyseres), df=fit$disparm), 
  	    	start.innov=rt(5000, df=fit$disparm))
		ts.plot(x=x7, ylim=c(-s,s), main=paste("Delta-asym-t-GARCH(1,1) - llh: ",
			as.character(floor(fit$value))))
		loglik <- c(loglik, fit$value)
		variance <- c(variance, var(x7))
		skewness <- c(skewness, skewness(x7))
	    kurtosis <- c(kurtosis, kurtosis(x7))
	    scaling <- c(scaling, scalinglaw(x7, doplot=FALSE)$exponent)
	    Continue <- readline("Press any key > ") 
			
# Summary:
		data.frame(cbind(loglik, variance, skewness, kurtosis, scaling))
		
		