
library(garch)
require(tseries)
DEBUG = FALSE
data(dem2gbp)
x <- dem2gbp[,1]
	
	
# 	DEFAULT FITS:	
	
	# --------------------------------------------------
	f1 <- garchFit(x ~ garch(1, 1), trace=TRUE)
	f1b <- garchFit(x ~ garch(1, 1), trace=TRUE, ndigit=16, 
			gradtol=1e-9, steptol=1e-9, iterlim = 500)
	# Final Estimate:
    #    	omega       alpha1        beta1           mu 
 	# 0.010758811  0.153395918  0.805896005 -0.006209217
 	
 	# Coefficient(s):
    #   	  Estimate  Std. Error  t value Pr(>|t|)    
	# mu     -0.006209    0.008472   -0.733 0.463612    
	# omega   0.010759    0.003081    3.492 0.000479 ***
	# alpha1  0.153396    0.027775    5.523 3.34e-08 ***
	# beta1   0.805896    0.035767   22.532  < 2e-16 ***
 
	# --------------------------------------------------
	f2 <- garchFit(x ~ garch(1, 1), algorithm="optim", trace=TRUE)
	# Final Estimate:
    #       omega       alpha1        beta1           mu 
 	# 0.010755876  0.153533280  0.805813356 -0.006225622
    #       omega       alpha1        beta1           mu 
    # 0.010762279  0.153465344  0.805821265 -0.006213497
	#   BENCHMARK   
    # 0.01076      0.1531       0.8060      -0.00619       
    
    # Coefficient(s): 
    #         Estimate  Std. Error  t value Pr(>|t|)    
	# mu     -0.006226    0.008469   -0.735    0.462  
	# omega   0.010756    0.002632    4.086 4.38e-05 ***
	# alpha1  0.153533    0.025456    6.031 1.63e-09 ***
	# beta1   0.805813    0.031358   25.697  < 2e-16 *** 
	
	# --------------------------------------------------
	f2b <- garchOxFit(formula.var = ~ garch(1, 1), series=x)
	# Coefficient(s):
    #        	  Value Std.Error t.value
	# Cst(M)  -0.006144 0.0084372 -0.7282
	# Cst(V)   0.010770 0.0013244  8.1318
	# ARCH(1)  0.153380 0.0140020 10.9540
	# GARCH(1) 0.805850 0.0165750 48.6190

	# --------------------------------------------------
	f3 <- garch(x+0.00619041, c(1,1))
	# Coefficient(s):
    #		   Estimate Std.Error t.value
	# a0  	   0.010720 0.001283   8.355 
	# a1  	   0.153196 0.013787  11.111   
	# b1  	   0.806256 0.015961  50.513   


	
# 	ARCH(1):
	fit11 <- garchFit(x ~ arch(1))
	fit12 <- garch(x-fit1$coef[3], order=c(0, 1), grad="analytical")
	fit13 <- garch(x-fit1$coef[3], order=c(0, 1), grad="numerical")
	fit14 <- garchOxFit(formula.var = ~garch(1, 1), series=x)

	
# ------------------------------------------------------------------------------


# 	GARCH(1,1):
	coef.True = c(-0.00619041, 0.0107613, 0.153134, 0.805974)
	coef.Ox <- garchOxFit(formula.var = ~ garch(1, 1), series=x)$coef[,1]
	coef.tseries <- c(NA, garch(x-coef.True[1], order=c(1, 1))$coef)
	fit <- garchFit(x ~ garch(1,1))$coef
	coef.fGARCH <- c(fit[4], fit[1:3])
	round(cbind(coef.True, coef.Ox, coef.tseries, coef.fGARCH),6)
	
	fit <- garchFit(x ~ garch(1,1), doprint=TRUE, trace=TRUE,
		ndigit=16, gradtol=1e-12, steptol=1e-12, iterlim=1000)
	
	# THIS INVESTIGATION:
	# BENCHMARK    -0.00619       0.01076      0.1531       0.8060
	# fGARCH       -0.00621       0.01076      0.1534       0.8059
	# G@RCH-OX     -0.00614       0.01077      0.1534       0.8058
	# TSERIES     BENCHMARK       0.01072      0.1532       0.8063

    # FROM:
    # FCP          -0.00619       0.01076      0.1531       0.8060
    # G@RCH-OX     -0.00618       0.01076      0.1534       0.8059
    # Eviews       -0.00541       0.00958      0.1423       0.8213
    # PcGive       -0.00625       0.01076      0.1534       0.8059
    # TSP          -0.00619       0.01076      0.1531       0.8060
    # S-Plus       -0.00919       0.01170      0.1543       0.8003     
	
	# FROM:
	# E-VIEWS      -0.00540 -0.64 0.0096  8.01 0.143  11.09 0.821  53.83
	# GAUSS-FANPAC -0.00600 -0.75 0.0110  3.67 0.153   5.67 0.806  23.71
	# LIMDEP       -0.00619 -0.71 0.0108  3.44 0.153   5.61 0.806  26.73
	# MATLAB       -0.00619 -0.73 0.0108  8.13 0.153  10.96 0.806  48.67
	# MICROFIT     -0.00621 -0.73 0.0108  3.78 0.153   5.78 0.806  24.02
	# SAS          -0.00619 -0.74 0.0108  8.15 0.153  10.97 0.806  48.60
	# SHAZAM       -0.00613 -0.73 0.0107  5.58 0.154   7.91 0.806  36.96
	# RATS         -0.00625 -0.71 0.0108  3.76 0.153   5.79 0.806  23.93
	# TSP          -0.00619 -0.67 0.0108  1.66 0.153   2.86 0.806  11.11
	
	
# 	MORE GARCH(1,1) ...
	# GARCH(1,1) - with initial values:
	garchFit(x ~ garch(1, 1), init=c(0.0108, 0.153, 0.806, -0.00621))
	# GARCH(1,1) - mu fixed to mean(x):
	garchFit(x ~ garch(1, 1), fixed=c(F, F, F, T))
	# GARCH(1,1) - mu fixed to benchmark value -0.00619041:
	garchFit(x ~ garch(1, 1),
		init=c(var(x)/10, 0.1, 0.8, -0.00619041), fixed=c(F, F, F, T))

		
# 	GARCH(2,1):
	garchFit(x ~ garch(2, 1))

	
# ------------------------------------------------------------------------------

	
# 	GJR-ARCH(1,1)
	
	garchFit(x ~ gjrarch(1,1), doprint=TRUE, trace=TRUE,
		init=c(0.011231, 0.140777, 0.028321, 0.801357, -0.007906) )
    #    omega     alpha1     gamma1      beta1         mu  
 	# 0.011232   0.140773   0.028358   0.801348  -0.007936  
	
 	garchOxFit(formula.mean = ~ arma(0, 0), formula.var = ~ gjr(1, 1), 
  	  	series = x, trace = TRUE)
  	#    omega     alpha1     gamma1      beta1         mu  
	# 0.011231   0.140777   0.028321   0.801357  -0.007906


# ------------------------------------------------------------------------------


#	APARCH(1,1):
	garchFit(x ~ aparch(1,1), doprint=TRUE, trace=TRUE)
	#   omega    alpha1    gamma1     beta1        mu     delta  
    # 0.02318   0.17483   0.09533   0.79700  -0.00941   1.35396  
	#     llh -1102.616
	
    garchFit(x ~ aparch(1,1), doprint=TRUE, trace=TRUE,
		ndigit=16, gradtol=1e-12, steptol=1e-12, iterlim=1000)
	#   omega    alpha1    gamma1     beta1        mu     delta  
    # 0.02318   0.17484   0.09534   0.79699  -0.00941   1.35363  
    #     llh -1102.615
    
    garchFit(x ~ aparch(1,1), doprint=TRUE, trace=TRUE, 
		init = c(0.0243, 0.173, 0.0953, 0.797, -0.00961, 1.29) )		
	#   omega    alpha1    gamma1     beta1        mu     delta  
    # 0.02318   0.17484   0.09533   0.79700  -0.00941   1.35396  
	#     llh -1102.616

		
	garchOxFit(formula.mean = ~ arma(0, 0), formula.var = ~ aparch(1, 1), 
  	  	series = x, trace = TRUE)
	
	Cst(M)              -0.00961     -0.00941
	Cst(V)               0.0243       0.0232 
	ARCH(Alpha1)         0.173     	  0.175  
	GARCH(Beta1)         0.800        0.797
	APARCH(Gamma1)       0.100        0.0953
	APARCH(Delta)        1.29         1.35
                        -1101.849    -1102.616
                        
                        
                        
           llh      omega     alpha1     gamma1      beta1         mu      delta 
   -1102.98474    0.02430    0.17300    0.09530    0.79700   -0.00961    1.29000

   
# ------------------------------------------------------------------------------

"garchOxFit" <-
#  	  function(formula.mean = ~ arma(0, 0), formula.var = ~ garch(1, 1), 
#  	  series = x, cond.dist = c("gaussian", "t", "ged", "skewed-t"), 
#  	  arfima=FALSE, arch.in.mean=0, truncation=100, trace = TRUE)
