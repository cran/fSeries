

# ------------------------------------------------------------------------------

    
    ## Load "systemfit" package:
    require(fBasics)
    ###
     
    
    ## Data:
    data(kmenta)
    kmenta
    names(kmenta)
    ###
 
        
    ## Estimation:
    formula = list(demand = q ~ p + d, supply = q ~ p + f + a)
    fitOLS = systemFit(formula,  "OLS", data = kmenta )
    fitOLS
    ###



# ------------------------------------------------------------------------------


data(kmenta)  
      
spec = new("fSYSTEM", 
	formulas = list(demand = q ~ p + d, supply = q ~ p + f + a),
	data = kmenta,
	method = "SUR",
	fit = list(),
	title = "Equation System Estimation: Kmenta Data",
	description = date())

show(spec)

systemFit(spec)

show(spec)

   
    
# ------------------------------------------------------------------------------


# Examples from the book ...
    
# CHAPTER 10: Systems of Regression Equations
#	10.3 Linear Seemingly Unrelated Regressions 
#   10.4 Nonlinear Seemingly Unrelated Regressions


################################################################################
# Chapter 10.3 - Linear Seemingly Unrelated Regressions


# Example:  SUR estimation of exchange rate system


    # S-Plus: SUR
    #	SUR(formulas, data, subset, na.rm = F, start = NULL, end = NULL, 
    #     contrasts = NULL, df = 1, tol = 1e-06, iterate = F, trace = T, ...)
    # R: SUR - Wrapper


	# SUR estimation w/o iteration. 
	args(SUR)
	###
	
	
	# Load from File and Extract non NA Data Records:
	require(fBasics)
	args(read.timeSeries)
    # Where are the Data?
    dataPath = "library/fSeries/data/"
    ###
    
    
    # The file "surex1.ts.csv" contains ...
    # Data are ... downloadable and can be updated from Yahoo's web site.
    surex1.ts = read.timeSeries(paste(dataPath, "surex1.ts.csv", sep = ""))
    start(surex1.ts)
    end(surex1.ts)
    ###
    
	
    # Note, start sample in August 1978 to eliminate NAs for USJY:
    surex1.ts = cutSeries(surex1.ts, from = "1978-08-01", to = "1996-06-30")
    start(surex1.ts)
    end(surex1.ts)
	# Print:
	surex1.ts
	dim(seriesData(surex1.ts))
	###
	
	
	# Column Variables:
	colnames(seriesData(surex1.ts))
	# .FP.lag1 are one month forward premia
	# .diff are future returns on spot currency
	### 
	
	
	# Create List of Formulas for Regression:
	formula.list = list(
	  USCNS.diff ~ USCN.FP.lag1, 
	  USDMS.diff ~ USDM.FP.lag1,
	  USFRS.diff ~ USFR.FP.lag1,
	  USILS.diff ~ USIL.FP.lag1,
	  USJYS.diff ~ USJY.FP.lag1,
	  USUKS.diff ~ USUK.FP.lag1)
	###
		
		
	# SUR Estimation:
	sur.fit = SUR(formula.list, data = surex1.ts)
	class(sur.fit)
	###
	
	
	# Print the SUR Estimate:
	sur.fit
	# Note, the printed report from R is slightly different from
	# that one produced by S-Plus. The differences are:
	#   The time period for the input time series is not printed,
	#	Standard errors and t-values are printed additionally,
	#	SSR, MSE, and R-Squared measures are computed additionally 
	###
	
	
	# The summary method provides more detailed information ...
	summary(sur.fit)
	# Again the printed report from R is different from SPLus
	# Durbin-Watson Stat is not printed:
	#	dw <- colSums((diff(res))^2)/colSums(res^2)
	###
	
	
	# Compare Results from R and SPlus:
    # -----------------------------------------------------
	# Results from R:
    #              Estimate  Std. Error  t value  Pr(>|t|)   
	# (Intercept)   -0.0031     0.0012   -2.5911    0.0102   
	# USCN.FP.lag1  -1.6602     0.5886   -2.8207    0.0052  
	# (Intercept)    0.0006     0.0025    0.2545    0.7994 
	# USDM.FP.lag1   0.4847     0.2195    2.2085    0.0283   
	# (Intercept)    0.0013     0.0024    0.5415    0.5887    
	# USFR.FP.lag1   0.9995     0.2088    4.7867    0.0000 
	# (Intercept)   -0.0006     0.0029   -0.2110    0.8331  
	# USIL.FP.lag1   0.4589     0.3623    1.2664    0.2067  
	# (Intercept)    0.0078     0.0031    2.5302    0.0121
	# USJY.FP.lag1  -1.7748     0.6975   -2.5445    0.0117
	# (Intercept)   -0.0036     0.0027   -1.3325    0.1841  
	# USUK.FP.lag1  -1.3068     0.6342   -2.0605    0.0406
 	# ----------------------------------------------------
	# Results from S-Plus:
    #                 Value  Std. Error  t value  Pr(>|t|) 
	#  (Intercept)  -0.0031     0.0012   -2.5943    0.0101 
	# USCN.FP.lag1  -1.6626     0.5883   -2.8263    0.0052 
	#  (Intercept)   0.0006     0.0024    0.2384    0.8118  
	# USDM.FP.lag1   0.5096     0.2089    2.4392    0.0155  
	#  (Intercept)   0.0013     0.0024    0.5550    0.5795  
	# USFR.FP.lag1   1.0151     0.1993    5.0928    0.0000  
	#  (Intercept)  -0.0006     0.0028   -0.2071    0.8361 
	# USIL.FP.lag1   0.4617     0.3584    1.2883    0.1990 
	#  (Intercept)   0.0078     0.0031    2.5226    0.0124 
	# USJY.FP.lag1  -1.7642     0.6961   -2.5342    0.0120 
	#  (Intercept)  -0.0035     0.0027   -1.3256    0.1864 
	# USUK.FP.lag1  -1.2963     0.6317   -2.0519    0.0414 
	# ----------------------------------------------------
	

	# Plot:
	# plot(sur.fit)
	###
	
	
	# Compute Iterated Estimator:
	# (The above results are based on the non-iterated estimator)
	sur.fit2 = SUR(formula.list, data = surex1.ts, iterate = TRUE)
	sur.fit2
	# ... converged after 8 iterations
	###
	
	
	# Compare non-iterated and iterated SUR:
	cbind(coef(sur.fit),coef(sur.fit2))
	###

	
	# Residual Correlation Matrix:
	# > sd.vals = sqrt(diag(sur.fit$Sigma))
	# > cor.mat = sur.fit$Sigma/outer(sd.vals,sd.vals)
	# It's not necessary to do this, R's summary method for SUR
	# objects prints these matrices!
	summary(sur.fit)
	###
	
	
	# compare to OLS fits
	ols.uscn.fit = OLS(USCNS.diff ~ USCN.FP.lag1, data = surex1.ts)
	summary(ols.uscn.fit)
	###
	
	
	# Compute Wald Statistic.
	bigR = matrix(0, 6, 12)
	bigR[1,2] = bigR[2,4] = bigR[3,6] = bigR[4,8] = bigR[5,10] = bigR[6,12] = 1
	rr = rep(1,6)
	bHat = as.vector(coef(sur.fit))
	avar = bigR%*%vcov(sur.fit)%*%t(bigR)
	Wald = t((bigR%*%bHat-rr))%*%solve(avar)%*%(bigR%*%bHat-rr)
	Wald
	1 - pchisq(Wald, 6)
	# ... the data reject the unbiased hypothesis
	###
	
	
	# Compute LR statistic:
	# ... the restricted model must first be estimated
	formula.list = list(
		(USCNS.diff - USCN.FP.lag1) ~ 1, (USDMS.diff - USDM.FP.lag1) ~ 1,
		(USFRS.diff - USFR.FP.lag1) ~ 1, (USILS.diff - USIL.FP.lag1) ~ 1,
		(USJYS.diff - USJY.FP.lag1) ~ 1, (USUKS.diff - USUK.FP.lag1) ~ 1)
	sur.fit2r = SUR(formula.list, data = surex1.ts, iterate = TRUE)
	sur.fit2r
	###
	
	
	# compute LR statistic
	nobs = nrow(residuals(sur.fit2r))
	LR = nobs*(
		determinant(sur.fit2r$Sigma, logarithm = TRUE)$modulus -
		determinant(sur.fit2$Sigma, logarithm = TRUE)$modulus )
	as.numeric(LR)
	1 - pchisq(LR, 6)
	###
	

# ------------------------------------------------------------------------------
# Chapter 10.4 - Non-Linear Seemingly Unrelated Regressions

	
	# Example 60: Black CAPM Model
	# Estimating and Testing the CAPM
	###


	# Load Data Records:
	black.ts = read.timeSeries("src/library/fFinMetrics/data/black.ts.csv")
	###
	
	
	# Do LR test of CAPM by estimating restricted and unrestriced 
	# SUR model should be able to do this by specifying no constant 
	# in the regression with -1 in formula.
	colIds(black.ts)
	formula.list = list(
		BOISE ~ MARKET,
		CITCRP ~ MARKET,
		CONED ~ MARKET,
		CONTIL ~ MARKET,
		DATGEN ~ MARKET,
		DEC ~ MARKET,
		DELTA ~ MARKET,
		GENMIL ~ MARKET,
		GERBER ~ MARKET,
		IBM ~ MARKET,
		MOBIL ~ MARKET,
		PANAM ~ MARKET,
		PSNH ~ MARKET,
		TANDY ~ MARKET,
		TEXACO ~ MARKET,
		WEYER ~ MARKET)
	
	capm.fit = SUR(formula.list, data = black.ts)
	coef(capm.fit)
	
	# construct Wald statistic to test CAPM
	
	# fit restricted model for CAPM estimation
	restricted.formula.list = list(
		BOISE ~ MARKET-1,
		CITCRP ~ MARKET-1,
		CONED ~ MARKET-1,
		CONTIL ~ MARKET-1,
		DATGEN ~ MARKET-1,
		DEC ~ MARKET-1,
		DELTA ~ MARKET-1,
		GENMIL ~ MARKET-1,
		GERBER ~ MARKET-1,
		IBM ~ MARKET-1,
		MOBIL ~ MARKET-1,
		PANAM ~ MARKET-1,
		PSNH ~ MARKET-1,
		TANDY ~ MARKET-1,
		TEXACO ~ MARKET-1,
		WEYER ~ MARKET-1)
	
	capm.restricted.fit = SUR(restricted.formula.list, data = black.ts,
		iter = T)
	
	nobs = nrow(residuals(capm.restricted.fit))
	LR = nobs*(determinant(capm.fit$Sigma,log = TRUE)$modulus-
		determinant(capm.restricted.fit$Sigma,log = TRUE)$modulus)
	LR
	
	#
	# Nonlinear SUR models
	#
	args(NLSUR)
	formula.list = list(
		y1 ~ b10 + b11*x1^b,
		y2 ~ b20 + b21*x1^b,
		y3 ~ b30 + b31*x1^b)
	start.vals = c(0, 1 ,0, 1, 0, 1, 0.5)
	names(start.vals) = c("b10", "b11", "b20", "b21", "b30", "b31", "b")
	
	test.dat = mk.zero1[,1:4]
	formula.list = list(M.3~C1+C2*M.1^C3)
	nlsur.fit = NLSUR(formula.list,data=test.dat)
	
	# estimate Black form of CAPM using berndt data
	colIds(black.ts)
	# create formulas for restricted model
	formula.list = list(
		BOISE~(1-b1)*g + b1*MARKET,
		CITCRP~(1-b2)*g + b2*MARKET,
		CONED~(1-b3)*g + b3*MARKET,
		CONTIL~(1-b4)*g + b4*MARKET,
		DATGEN~(1-b5)*g + b5*MARKET,
		DEC~(1-b6)*g + b6*MARKET,
		DELTA~(1-b7)*g + b7*MARKET,
		GENMIL~(1-b8)*g + b8*MARKET,
		GERBER~(1-b9)*g + b9*MARKET,
		IBM~(1-b10)*g + b10*MARKET,
		MOBIL~(1-b11)*g + b11*MARKET,
		PANAM~(1-b12)*g + b12*MARKET,
		PSNH~(1-b13)*g + b13*MARKET,
		TANDY~(1-b14)*g + b14*MARKET,
		TEXACO~(1-b15)*g + b15*MARKET,
		WEYER~(1-b16)*g + b16*MARKET)
	
	start.vals = c(0, rep(1, 16))
	names(start.vals) = c("g", paste("b", 1:16, sep = ""))
	
	nlsur.fit = NLSUR(formula.list,data = black.ts,
		coef = start.vals, start = "Jan 1983", in.format = "%m %Y")
	names(nlsur.fit)
	nlsur.fit$message
	nlsur.fit
	 
	# compute standard error for b
	std.ers = sqrt(diag(vcov(nlsur.fit)))
	cbind(coef(nlsur.fit),std.ers)
	
	plot(nlsur.fit, ask=F)
	
	# unrestricted model
	colIds(black.ts)
	formula.list = list(
		BOISE~a1+b1*MARKET,
		CITCRP~a2+b2*MARKET,
		CONED~a3+b3*MARKET,
		CONTIL~a4+b4*MARKET,
		DATGEN~a5+b5*MARKET,
		DEC~a6+b6*MARKET,
		DELTA~a7+b7*MARKET,
		GENMIL~a8+b8*MARKET,
		GERBER~a9+b9*MARKET,
		IBM~a10+b10*MARKET,
		MOBIL~a11+b11*MARKET,
		PANAM~a12+b12*MARKET,
		PSNH~a13+b13*MARKET,
		TANDY~a14+b14*MARKET,
		TEXACO~a15+b15*MARKET,
		WEYER~a16+b16*MARKET)
	
	start.vals = c(rep(0, 16), rep(1, 16))
	names(start.vals) = 
		c(paste("a", 1:16, sep = ""), paste("b", 1:16, sep = ""))
	
	nlsur.fit2 = NLSUR(formula.list,data = black.ts,
		coef = start.vals, start = "Jan 1983", in.format = "%m %Y")
	
	# compute LR statistic
	nobs = nrow(residuals(nlsur.fit2))
	LR = nobs*(determinant(nlsur.fit$Sigma, log = TRUE)$modulus-
	determinant(nlsur.fit2$Sigma, log = TRUE)$modulus)
	as.numeric(LR)
	1 - pchisq(LR,16)
	
	#
	# Example: Estimating and testing Factor model a la 
	#	Burmiester and Elton see JBES 1988 article.
	#
	
	#
	# estimation of exchange rate system imposing cross equation restrictions
	#
	# NLSUR estimation of exchange rate system subject to cross equation restrictions
	# 
	formula.list = list(
		USCNS.diff ~ a1+g*USCN.FP.lag1,
		USDMS.diff ~ a2+g*USDM.FP.lag1,
		USFRS.diff ~ a3+g*USFR.FP.lag1,
		USILS.diff ~ a4+g*USIL.FP.lag1,
		USJYS.diff ~ a5+g*USJY.FP.lag1,
		USUKS.diff ~ a6+g*USUK.FP.lag1)
	
	start.vals = c(rep(0,6),1)
	names(start.vals) = c(paste("a",1:6,sep=""),"g")
	
	# NLSUR estimation. Note: start sample in 1978 to eliminate NA values
	# for USJY data
	args(NLSUR)
	nlsur.fit = NLSUR(formula.list, data = surex1.ts,
	start = "Aug 1978",in.format = "%m %Y")
	class(nlsur.fit)
	names(nlsur.fit)
	nlsur.fit
	
	# NLSUR estimation specifying starting values
	nlsur.fit = NLSUR(formula.list, data = surex1.ts, coef = start.vals,
	start="Aug 1978", in.format = "%m %Y")
	nlsur.fit
	
	# compute standard error for g
	sqrt(diag(vcov(nlsur.fit)))[7]
	
	# LR test of common parameter restriction
	nobs = nrow(residuals(nlsur.fit))
	LR = nobs*(determinant(nlsur.fit$Sigma, log = TRUE)$modulus
	-determinant(sur.fit2$Sigma, log = TRUE)$modulus)
	as.numeric(LR)
	1 - pchisq(LR, 6)
	

# ------------------------------------------------------------------------------

