
#
# Example:
#	Test file for the internal .urcval() function.
#
# Description:
#	This example tests the functionality of the internal function .urcval()
#	which is used for unit root testing. .urcval() computes the statistics 
#   from function based on J.G. McKinnon's Fortran program and data tables.
#	1 Compare Table 4.2 from Banerjee et al. with Internal urcval()
#	2 Compute p Values from "tseries" Package with those from JGM
#   3 Compute ADF Test for Chicken Data, compare with STATA/Splus
#   4 Test the function punitroot(), compare results with SPlus
#   5 Test the function qunitroot(), compare results with SPlus
#   6 Reproduce Figures 1 to 3 as plotted in the paper of JGM
#
# Author:
#	(C) 2004, Diethelm Wuertz, GPL
#


################################################################################
# 1 Compare Table 4.2 from Banerjee et al. with Internal .urcval()


	# In the book "Co-Integration, Error Correction, and the 
	#	Econometric Analysis of Non-Stationary Data" by A. Banaerjee
	#	et al. Chapter 4 lists in Table 4.2 the cumulative distribution
	#	for the regression types "nc", "c" and "ct". Compute the
	#	values from the function .urcval() and compare the results 
	#	with those given in Table 4.2.
	#
	# Table 4.2: Empirical Cumulative Distribution
	#	
	#	     0.01  .025  0.05   0.1   0.9  0.95  .975  0.99
	#	25  -2.66 -2.26 -1.95 -1.60  0.92  1.33  1.70  2.16
	#	50  -2.62 -2.25 -1.95 -1.60  0.91  1.31  1.66  2.08
	#	200 -2.60 -2.24 -1.95 -1.61  0.90  1.29  1.64  2.03
	#	250 -2.58 -2.23 -1.95 -1.61  0.89  1.29  1.63  2.01
	#	500 -2.58 -2.23 -1.95 -1.62  0.89  1.28  1.62  2.00
	#	Inf -2.58 -2.23 -1.95 -1.62  0.89  1.28  1.62  2.00
	#	
	#	     0.01  .025  0.05  0.10  0.90  0.95  .975  0.99
	#	25  -3.75 -3.33 -3.00 -2.63 -0.37  0.00  0.34  0.72
	#	50  -3.58 -3.22 -2.93 -2.60 -0.40 -0.03  0.29  0.66
	#	200 -3.51 -3.17 -2.89 -2.58 -0.42 -0.05  0.26  0.63
	#	250 -3.46 -3.14 -2.88 -2.57 -0.42 -0.06  0.24  0.62
	#	500 -3.44 -3.13 -2.87 -2.57 -0.43 -0.07  0.24  0.61
	#	Inf -3.43 -3.12 -2.86 -2.57 -0.44 -0.07  0.23  0.60
	#	
	#	     0.01  .025  0.05  0.10  0.90  0.95  .975  0.99
	#	25  -4.38 -3.95 -3.60 -3.24 -1.14 -0.80 -0.50 -0.15
	#	50  -4.15 -3.80 -3.50 -3.18 -1.19 -0.87 -0.58 -0.24
	#	200 -4.04 -3.73 -3.45 -3.15 -1.22 -0.90 -0.62 -0.28
	#	250 -3.99 -3.69 -3.43 -3.13 -1.23 -0.92 -0.64 -0.31
	#	500 -3.98 -3.68 -3.42 -3.13 -1.24 -0.93 -0.65 -0.32
	#	Inf -3.96 -3.66 -3.41 -3.12 -1.25 -0.94 -0.66 -0.33
	#	
	#	     0.01  0.02  0.05  0.10  0.90  0.95  0.97  0.99
	#	Inf -2.33 -1.96 -1.64 -1.28  1.28  1.64  1.96  2.33
		
	# Internal Function .urcval():
	args(.urcval)
	# Arguments: 
    #   arg - level of test (between .0001 and .9999) if nc = 1, 
    #       else test statistic if nc = 2
    #   nobs - sample size (0 for asymptotic)
    #   niv - number of integrated variables
    #   itt - 1 or 2 for tau or z test
    #   itv - 1, 2, 3, 4 for nc, c, ct, ctt
    #   nc - 1 or 2 for critical value or P value
    # Value:
    #   val - critical value if nc = 1 (returned by routine), 
    #       else P value if nc = 2 (returned by routine)
		
	# Compare:
	cat("\n\nTable 4.2: Empirical Cumulative Distribution\n\n")
	nobs.vec = c(25, 50, 200, 250, 500, 0)	
	x = c(0.01, 0.025, 0.05, 0.10, 0.90, 0.95, 0.975, 0.99)
	ECD = NULL
	for (itv in 1:3) {
	    ECD = rbind(ECD, x)
	    for ( nobs in nobs.vec ) {
			ans = .urcval(x, nobs = nobs, niv = 1, itt = 1, itv = itv, nc = 1)
			ECD = rbind(ECD, ans) 
		}
		ECD = round(ECD, 2) 
	}	
	NORM = matrix(round(qnorm(x), 2), byrow = TRUE, ncol = length(x))
	ECD = rbind(ECD, x, NORM)
	rowNames = c("", as.character(nobs.vec[1:5]), "Inf")
	rowNames = c(rep(rowNames, 3), "", "Inf")
	rownames(ECD) = rowNames
	colnames(ECD) = as.character(x)
	ECD[-1,] 
	###
	

################################################################################
# 2 Compute p.values from "tseries" Package with those from JGM


	# Compare p.values as used in Traplett's tseries Package with 
	# those from JGM:
	#	The ADF Test implemented in the adf.test() function from R's
	#	package 'tseries' interpolates the "p.values" as function of
	#	the statistics and the sample size. Compare for a sample of
	#	100 observations for the whole range of statistics the
	#	interpolated values from the Table 4.2 with those obtained
	#	from the function urcval().

    # Settings:
    n = 100
    STAT = (-80:20)/10 
    # Set Table:
    table = cbind(
    	c(4.38, 4.15, 4.04, 3.99, 3.98, 3.96), 
    	c(3.95, 3.80, 3.73, 3.69, 3.68, 3.66), 
    	c(3.60, 3.50, 3.45, 3.43, 3.42, 3.41), 
    	c(3.24, 3.18, 3.15, 3.13, 3.13, 3.12), 
    	c(1.14, 1.19, 1.22, 1.23, 1.24, 1.25), 
    	c(0.80, 0.87, 0.90, 0.92, 0.93, 0.94), 
    	c(0.50, 0.58, 0.62, 0.64, 0.65, 0.66), 
    	c(0.15, 0.24, 0.28, 0.31, 0.32, 0.33))
    table = -table
    tablen = dim(table)[2]
    tableT = c(25, 50, 100, 250, 500, 1.0e+05)
    tablep = c(0.01, 0.025, 0.050, 0.100, 0.900, 0.950, 0.975, 0.990)
    tableipl = numeric(tablen)
    ###
    
    # Interpolate:
    for (i in (1:tablen)) tableipl[i] = 
    	approx(x = tableT, y = table[, i], n, rule = 2)$y 	
   	p.ts = approx(tableipl, tablep, STAT, rule = 2)$y   	
   	###
   	
   	# Plot interpolatre data from Table 4.2:
   	par(mfrow = c(2, 2), cex = 0.75)
   	plot(STAT, log(p.ts), ylim = c(-12,0), type = "l", main = "ct")
   	# Add data from .urcval():
   	p.urcval = .urcval(STAT, nobs = n, niv = 1, itt = 1, itv = 3, nc = 2)
   	lines(STAT, log(p.urcval), col = "blue")  	
   	# Compute relative differences:
   	p.diff = abs(p.ts-p.urcval)/p.urcval
   	plot(STAT, log(p.diff), type = "l", main = paste("n =", n))
   	###
   	
   	# Print Results:
   	cbind(STAT, p.ts, p.urcval)
   	###
   		

################################################################################
# 3 Compute ADF Test for Chicken and Egg Data, compare with STATA/Splus


	# Compute the ADF Test for the "eggs" data and campare the results
	# with those obtained from the STATA and SPLUS functions:

	# Data:
    eggs = c(
    	3581, 3532, 3327, 3255, 3156, 3081, 3166, 3443, 3424, 3561, 3640,
		3840, 4456, 5000, 5366, 5154, 5130, 5077, 5032, 5148, 5404, 5322,
		5323, 5307, 5402, 5407, 5500, 5442, 5442, 5542, 5339, 5358, 5403,
		5345, 5435, 5474, 5540, 5836, 5777, 5629, 5704, 5806, 5742, 5502,
		5461, 5382, 5377, 5408, 5608, 5777, 5825, 5625, 5800, 5656)
	chic = c(
		468491, 449743, 436815, 444523, 433937, 389958, 403446, 423921, 
		389624, 418591, 438288, 422841, 476935, 542047, 582197, 516497, 
		523227, 467217, 499644, 430876, 456549, 430988, 426555, 398156, 
		396776, 390708, 383690, 391363, 374281, 387002, 369484, 366082, 
		377392, 375575, 382262, 394118, 393019, 428746, 425158, 422096, 
		433280, 421763, 404191, 408769, 394101, 379754, 378361, 386518, 
		396933, 400585, 392110, 384838, 378609, 364584)
	year = 1930:1983
	###
	
	# Plot:
    plot(chic, type = "l")  
    ### 
      
    # STATA - Results:
    statisticSTATA = c(nc = -0.712, c = -1.618, ct = -1.998)
    p.valueSTATA   = c(nc = NA,     c =  0.474, ct =  0.603)   
    # SPLUS - Results:
    if (FALSE) {
        # Splus Code:
        chic = as.vector(ts[, 1])
        unitroot(x = chic, trend = "nc", lags = 2, asymptotic = FALSE)
        unitroot(x = chic, trend = "c",  lags = 2, asymptotic = FALSE)
        unitroot(x = chic, trend = "ct", lags = 2, asymptotic = FALSE) }    
    statisticSPLUS = c(nc = -0.7122, c = -1.618,  ct = -1.998 )
    p.valueSPLUS   = c(nc =  0.4034, c =  0.4663, ct =  0.5886) 
    ###
  
    # Start Settings:
    ans1 = ans2 = NULL
    rowNames = c("nc", "c", "ct")   
    # Compute from adfTest:
    for (type in rowNames) {
        res1 = adfTest(x = chic, type = type, lags = 1)
        res2 = unitrootTest(x = chic, type = type, lags = 1)
        ans1 = rbind(ans1, c(res1$statistic, res1$p.value)) 
        ans2 = rbind(ans2, c(res2$statistic, res2$p.value)) 
    }        
    # Bind:
    ans = cbind(ans1, ans2)    
    # Add SPLUS results:
    ans = cbind(round(ans, digits = 3), statisticSPLUS, p.valueSPLUS)    
    # Add STATA results:
    ans = cbind(round(ans, digits = 3), statisticSTATA, p.valueSTATA)   
    # Add row and column names:
    rownames(ans) = rowNames
    colnames(ans) = c(
        "sADF", "pADF", "sJGM", "pJGM", "sSPLUS", "pSPLUS", "sSTATA", "pSTATA")
    ###
    
    # Print All Result:
    print(ans)
    ###
    
 
################################################################################
# 4 Test the function punitroot(), compare results with SPlus 


	# Compare results with those obtained from Splus 6.1:

	# Statistic: "t"
	round(punitroot(q=-6:2, n.sample = 100, trend = "nc", statistic = "t"), 4)
	# R:      0.0000 0.0000 0.0001 0.0030 0.0441 0.2829 0.6803 0.9155 0.9889
	# Splus:  0.0000 0.0000 0.0001 0.0030 0.0441 0.2829 0.6803 0.9155 0.9889
	round(punitroot(q=-6:2, n.sample = 100, trend = "c", statistic = "t"), 4)
	# R:      0.0000 0.0001 0.0021 0.0383 0.2865 0.7510 0.9558 0.9964 0.9999
	# Splus:  0.0000 0.0001 0.0021 0.0383 0.2865 0.7510 0.9558 0.9964 0.9999	
	round(punitroot(q=-6:2, n.sample = 100, trend = "ct", statistic = "t"), 4)
	# R:      0.0000 0.0004 0.0117 0.1375 0.5940 0.9387 0.9958 0.9999 1.0000
	# Splus:  0.0000 0.0004 0.0117 0.1375 0.5940 0.9387 0.9958 0.9999 1.0000	
	round(punitroot(q=-6:2, n.sample = 100, trend = "ctt", statistic = "t"), 4)
	# R:      0.0001 0.0020 0.0384 0.2979 0.8077 0.9847 0.9995 1.0000 1.0000
	# Splus:  0.0001 0.0020 0.0384 0.2979 0.8077 0.9847 0.9995 1.0000 1.0000	
	###
	
	# Statistic: "n"
	round(punitroot(q=-6:2, n.sample = 100, trend = "nc", statistic = "n"), 4)
	# R:      0.0872 0.1197 0.1654 0.2307 0.3262 0.4700 0.6803 0.9105 0.9884
	# Splus:  0.0872 0.1197 0.1654 0.2307 0.3262 0.4700 0.6803 0.9105 0.9884	
	round(punitroot(q=-6:2, n.sample = 100, trend = "c", statistic = "n"), 4)
	# R:      0.3375 0.4261 0.5318 0.6519 0.7759 0.8837 0.9558 0.9881 0.9975
	# Splus:  0.3375 0.4261 0.5318 0.6519 0.7759 0.8837 0.9558 0.9881 0.9975	
	round(punitroot(q=-6:2, n.sample = 100, trend = "ct", statistic = "n"), 4)
	# R:      0.7353 0.8144 0.8828 0.9349 0.9688 0.9874 0.9958 0.9988 0.9997
	# Splus:  0.7353 0.8144 0.8828 0.9349 0.9688 0.9874 0.9958 0.9988 0.9997	
	round(punitroot(q=-6:2, n.sample = 100, trend = "ctt", statistic = "n"), 4)
	# R:      0.9213 0.9542 0.9761 0.9889 0.9954 0.9984 0.9995 0.9998 1.0000
	# Splus:  0.9213 0.9542 0.9761 0.9889 0.9954 0.9984 0.9995 0.9998 1.0000
	###
	
	
################################################################################
# 5 Test the function qunitroot(), compare results with SPlus 


	# Compare results with those obtained from Splus 6.1
	# q-Values:
	q = c(0.15, 0.30, 0.45, 0.60, 0.75, 0.90)
	###
	 
	# Statistic: "t"
	round(qunitroot(q, n.sample = 100, trend = "nc", statistic = "t"), 4)
	# R:      -1.3979 -0.9576 -0.6118 -0.2328  0.2265  0.8967
	# Splus:  -1.3979 -0.9576 -0.6118 -0.2328  0.2265  0.8967 	
	round(qunitroot(q, n.sample = 100, trend = "c", statistic = "t"), 4)
	# R:      -2.3799 -1.9690 -1.6568 -1.3577 -1.0029 -0.4232    
	# Splus:  -2.3799 -1.9690 -1.6568 -1.3577 -1.0029 -0.4232	
	round(qunitroot(q, n.sample = 100, trend = "ct", statistic = "t"), 4)
	# R:      -2.9559 -2.5588 -2.2625 -1.9889 -1.6869 -1.2227     
	# Splus:  -2.9559 -2.5588 -2.2625 -1.9889 -1.6869 -1.2227	
	round(qunitroot(q, n.sample = 100, trend = "ctt", statistic = "t"), 4)
	# R:      -3.3927 -2.9955 -2.7002 -2.4292 -2.1354 -1.7073    
	# Splus:  -3.3928 -2.9955 -2.7002 -2.4291 -2.1354 -1.7073
	###
	
	# Statistic: "n"
	round(qunitroot(q, n.sample = 100, trend = "nc", statistic = "n"), 4)
	# R:      -4.2999 -2.2371 -1.1160 -0.3463  0.2795  0.9382    
	# Splus:  -4.2999 -2.2372 -1.1160 -0.3463  0.2795  0.9382	
	round(qunitroot(q, n.sample = 100, trend = "c", statistic = "n"), 4)
	# R:      -9.2935 -6.4935 -4.7588 -3.4215 -2.2137 -0.8166   
	# Splus:  -9.2934 -6.4936 -4.7588 -3.4215 -2.2137 -0.8165	
	round(qunitroot(q, n.sample = 100, trend = "ct", statistic = "n"), 4)
	# R:     -15.4278 -11.8711  -9.5254  -7.6163  -5.8205  -3.7052   
	# Splus: -15.4277 -11.8710  -9.5254  -7.6163  -5.8206  -3.7050	
	round(qunitroot(q, n.sample = 100, trend = "ctt", statistic = "n"), 4)
	# R:     -20.5330 -16.4756 -13.7175 -11.4106  -9.1839  -6.5128     
	# Splus: -20.5330 -16.4755 -13.7174 -11.4107  -9.1839  -6.5126
	###


################################################################################
# 6 Reproduce Figure 1 to 3 as plotted in the paper of J.G. MacKinnon


	# JGM: Figure 1 to 3
	figure1to3 = function(figure) {
		x = (-80:20)/10
		plot (x = c(-8, 2), y = c(0, 1), type = "n")
		for (i in 1:12) {
			y = .urcval(x, nobs = 0, niv = i, itt = 1, itv = figure, nc = 2)
			lines(x, y) 
		}
		if (figure == 3) {
			for (i in 1:12) {
				y = .urcval(x, nobs = 0, niv = i, itt = 1, itv = 4, nc = 2)
				lines(x, y, lty = 3) 
			} 
		}
		invisible() 
	}
	###
			
	# Plot:
	par(mfcol = c(3, 2), cex = 0.5)
	figure1to3(1)
	figure1to3(2)
	figure1to3(3)
	###


################################################################################
	
