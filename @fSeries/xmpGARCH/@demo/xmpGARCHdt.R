
# 	SETTINGS *******************************************************************
	library(garch)
	require(tseries)
	DEBUG = FALSE
	data(dem2gbp)
	x <- dem2gbp[,1]
	

# COMPARE WITH *****************************************************************

	f0 = garchFit(formula.var = ~garch(1,1), trace=TRUE)
	f1 = garchFit(formula.var = ~garch(1,1), 
		 cond.dist="dt", dist.est=TRUE, dist.par=4, trace=TRUE)	
	f2 = garchFit(formula.var = ~ garch(1,1), 
		 cond.dist="dt", dist.est=FALSE, dist.par=4.3301, trace=TRUE)		 
	f3 = garchOxFit(formula.var = ~ garch(1, 1), series=x, cond.dist = "t") 

	
	#   		   llh     omega    alpha1     beta1        mu        df   		             
  	# ------------------------------------------------------------------
  	# f0:	  -1106.18  0.010741  0.015357  0.080585 -0.006128       Inf   
	# f1:      -989.37  0.002323  0.124711  0.884414  0.002198  4.124905 
    # f2:      -989.50  0.002345  0.122577  0.883219  0.002019  4.330119 FIXED
 	# fG@RCH   -989.73  0.002734  0.117290  0.882810  0.002170  4.330119
    # SPlus:	        0.001981  0.120837  0.887464  0.001175  4.163737  
 	
    
# ------------------------------------------------------------------------------


# garchOxFit <-
#  	  function(formula.mean = ~ arma(0, 0), formula.var = ~ garch(1, 1), 
#  	  series = x, cond.dist = c("gaussian", "t", "ged", "skewed-t"), 
#  	  arfima=FALSE, arch.in.mean=0, truncation=100, trace = TRUE)




coef1 = c( 0.002734, 0.117290, 0.882810, 0.002170, 4.330119 )
coef2 = c( 0.001393, 0.070060, 0.876073, 0.002186, 4.330119 )
x <<- dem2gbp[,1]


	garch11LLH = function(coef) {		
		omega = coef[1]
		alpha = coef[2]
		beta = coef[3]
		mu = coef[4]
		df = coef[5]		
		nt = length(x)
		h = rep(var(x), nt)	
		for (i in 2:nt) h[i] = omega + alpha*(x[i-1]-mu)^2 + beta*h[i-1]	   	 			   	   	
    	hh = sqrt(abs(h[2:nt]))*sqrt((df-2)/df)
    	llh = -sum(  log(  dt(  (x[2:nt]-mu)/hh, df  ) / hh ) )
		llh }
		
	garch11LLH(coef1)
	garch11LLH(coef2)
	
	coef1[2]+coef1[3]
	coef2[2]+coef2[3]


nu = 4
h = 1
z = seq(-5, 5, 0.1)

d = log(gamma((nu+1)/2)) - log(gamma(nu/2)) - 0.5*log(pi*(nu-2)) -
	0.5*log(h) -0.5 *(1+nu) *log(1+z*z/(nu-2))
	

plot(z,d)

