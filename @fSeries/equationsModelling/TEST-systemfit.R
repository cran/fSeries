## Examples from the 'systemfit' Package:
   data(kmenta)
   
################################################################################
# OLS Estimation:   

   	# Fit:
   	eqns = list(demand = q ~ p + d, supply = q ~ p + f + a )
   	fit.ols = systemfit(method = "OLS", eqns = eqns, data = kmenta)
   	class(fit.ols)

   	###
   	
   	# Methods ...
   	print(fit.ols)
   	plot(fit.ols # Error, ther is no plot method
   	summary(fit.ols) # same as print
   	###
   	
   	# More <stats> Methods ...
	coef(fit.ols)			# returns $b
	confint(fit.ols)        # for 0.95 returns kevels 0.025 and 0.975 !!
	fitted(fit.ols)
	residuals(fit.ols)
	vcov(fit.ols)
	
	# fEQNS Slots ...
	FITOLS@call
	FITOLS@data
	FITOLS@description
	FITOLS@formulas
	FITOLS@method
	FITOLS@title
	
	# "@fit" Slot ...
	FITOLS@fit$coef
    FITOLS@fit$confint
    FITOLS@fit$fitted
    FITOLS@fit$residuals
    FITOLS@fit$vcov
    
    
    FITOLS@fit$method
	FITOLS@fit$g
	FITOLS@fit$n
	FITOLS@fit$k
	FITOLS@fit$ki
	FITOLS@fit$df
  # FITOLS@fit$iter
	FITOLS@fit$b
  # FITOLS@fit$bt
	FITOLS@fit$se
	FITOLS@fit$t
	FITOLS@fit$p
	FITOLS@fit$bcov
  # FITOLS@fit$btcov
	FITOLS@fit$rcov
	FITOLS@fit$drcov
	FITOLS@fit$rcovest
	FITOLS@fit$rcor
	FITOLS@fit$olsr2
  # FITOLS@fit$mcelr2
	FITOLS@fit$y
	FITOLS@fit$x
  # FITOLS@fit$h				[only "2SLS" and "3SLS"]
	FITOLS@fit$data
  # FITOLS@fit$R.restr
	FITOLS@fit$q.restr
  # FITOLS@fit$TX
	FITOLS@fit$maxiter
	FITOLS@fit$tol
	FITOLS@fit$rcovformula
	FITOLS@fit$formula3sls
	FITOLS@fit$probdfsys
	FITOLS@fit$single.eq.sigma
	FITOLS@fit$solvetol
   
   
   
################################################################################   
# OLS Estimation with 2 Restrictions:
   
	Rrestr <- matrix(0, 2, 7)
   	qrestr <- matrix(0, 2, 1)
   	Rrestr[1,3] =  1
   	Rrestr[1,7] = -1
   	Rrestr[2,2] = -1
   	Rrestr[2,5] =  1
   	qrestr[2,1] =  0.5
   	FITOLS2 = eqnsFit(formulas, data = kmenta, R.restr = Rrestr, 
     	q.restr = qrestr)
   	FITOLS2
   
## Iterated SUR Estimation:
   FITSUR = eqnsFit(formulas, data = kmenta, method = "SUR", maxit = 100)
   FITSUR
   # Coefficients, Fitted Values, Residuals and Variance-Covariance Matrix:
   # Call by Method:
   coef(FITSUR)
   fitted(FITSUR)
   residuals(FITSUR)
   vcov(FITSUR)

## 2SLS Estimation:
   inst = ~ d + f + a
   FIT2SLS = eqnsFit(formulas, data = kmenta, method = "2SLS", inst = inst)
   FIT2SLS
   # Coefficients, Fitted Values, Residuals and Variance-Covariance Matrix:
   # Call by Slot:
   FIT2SLS@fit$coef
   FIT2SLS@fit$fitted
   FIT2SLS@fit$residuals
   FIT2SLS@fit$vcov

## 2SLS Estimation with Different Instruments in Each Equation:
   insts = list( ~ d + f, ~ d + f + a)
   FIT2SLS2 = eqnsFit(formulas, data = kmenta, method = "2SLS", inst = insts)
   FIT2SLS2

## 3SLS Estimation with GMM-3SLS Formula:
   instruments = ~ d + f + a
   FIT3SLS = eqnsFit(formulas, data = kmenta, method = "3SLS", 
   	 inst = instruments, formula3sls = "GMM")
   FIT3SLS