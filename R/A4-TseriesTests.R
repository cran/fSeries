
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2, or (at your option)
# any later version.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.
#
# A copy of the GNU General Public License is available via WWW at
# http://www.gnu.org/copyleft/gpl.html.  You can also obtain it by
# writing to the Free Software Foundation, Inc., 59 Temple Place,
# Suite 330, Boston, MA  02111-1307  USA. 

# Copyrights (C)
# for this R-port: 
#   1999 - 2004, Diethelm Wuertz, GPL
#   Diethelm Wuertz <wuertz@itp.phys.ethz.ch>
#   info@rmetrics.org
#   www.rmetrics.org
# for the code accessed (or partly included) from other R-ports:
#   see R's copyright and license files
# for the code accessed (or partly included) from contributed R-ports
# and other sources
#   see Rmetrics's copyright file


################################################################################
# FUNCTION:					DESCRIPION
#  tsTest					 Time Series Test Suite
#   bdsTest					  Brock-Dechert-Scheinkman test for iid series
# FUNCTION:                 NORMALITY TESTS
#   jbTest			 		  Jarque-Bera Normality Test
# FUNCTION:                 NONLINEARITY TESTS
#   wnnTest				 	  White Neural Network Test for Nonlinearity
#   tnnTest			 		  Teraesvirta Neural Network Test for Nonlinearity
################################################################################
# BUILTIN - PACKAGE DESCRIPTION:
#  Package: tseries
#  Version: 0.9-21
#  Date: 2004-04-23
#  Title: Time series analysis and computational finance
#  Author: Compiled by Adrian Trapletti <a.trapletti@bluewin.ch>
#  Maintainer: Kurt Hornik <Kurt.Hornik@R-project.org>
#  Description: Package for time series analysis and computational finance
#  Depends: R (>= 1.9.0), quadprog
#  License: GPL (see file COPYING)
#  Packaged: Thu Apr 22 16:32:16 2004; hornik
################################################################################
# NOTES:
#  The runs.test is avbailable as dependency test in the fBasics Package
#    runs.test = function (x)
#  Most of the functions are BUILTIN from Adrian Trapletti's R package
#    tseries
################################################################################



tsTest = 
function(x, 
method = c("jb", "tnn", "wnn", "bds" ), ...)
{	# A function implemented by Diethelm Wuertz

	# Load Library:
    # require(tseries)
    
    # Transform:
	if (class(x) == "its") x = as.vector(x@.Data[, 1])
	x = as.vector(x)
	
	# Settings:
	method = method[1] 
	test = paste(method, "Test", sep = "")
	fun = match.fun(test)
	
	# Call:
    # cat("\nCall:\n")
    # cat(paste(deparse(match.call()), sep = "\n", collapse = "\n"), 
    #     "\n\n", sep = "")
    # cat("\nHypothesis Test:")    
    	
	# Return Result:
	fun(x = x, ...)
	
}


# ------------------------------------------------------------------------------


jbTest = 
function(x)
{	# A function implemented by Diethelm Wuertz

	# Description:
	#	Jarque-Bera Test
	
	# Notes:
	#	This function is a slightly modified copy of Adrian Trapletti's
	#	contributed function from his 'tseries' package.
	
	# FUNCTION:
	
	# Transform:
	x = as.vector(x)
	
	# Internal Function:
	# Source: tseries - Adrian Trapletti
	jarque.bera.test = 
	function(x) {
	    if (NCOL(x) > 1) stop("x is not a vector or univariate time series")
	    if (any(is.na(x))) stop("NAs in x")
	    DNAME = deparse(substitute(x))
	    n = length(x)
	    m1 = sum(x)/n
	    m2 = sum((x-m1)^2)/n
	    m3 = sum((x-m1)^3)/n
	    m4 = sum((x-m1)^4)/n
	    b1 = (m3/m2^(3/2))^2
	    b2 = (m4/m2^2)
	    STATISTIC = n*b1/6+n*(b2-3)^2/24
	    names(STATISTIC) = "X-squared"
	    PARAMETER = 2
	    names(PARAMETER) = "df"
	    PVAL = 1 - pchisq(STATISTIC,df = 2)
	    METHOD = "Jarque Bera Test"
	    structure(list(statistic = STATISTIC, parameter = PARAMETER,
	    	p.value = PVAL, method = METHOD, data.name = DNAME),
		class = "htest") }
	
	# Return Value:
	jarque.bera.test(x)
}


# ------------------------------------------------------------------------------


bdsTest = 
function(x, m = 3, eps = seq(0.5*sd(x),2*sd(x),length=4), trace = FALSE)
{	# A function implemented by Diethelm Wuertz

	# Description:
	#	Brock-Dechert-Scheinkman test for iid series
	
	# Notes:
	#	This function is a slightly modified copy of Adrian Trapletti's
	#	contributed function from his 'tseries' package.
	
	# Transform:
	x = as.vector(x)
	
	# FUNCTION:
	
	# Internal Function:
	# Source: tseries - Adrian Trapletti
	bds.test = 
	function(x, m = 3, eps = seq(0.5*sd(x),2*sd(x),length=4), trace = FALSE) {
	    if (NCOL(x) > 1) stop("x is not a vector or univariate time series")
	    if (any(is.na(x))) stop("NAs in x")
	    if (m < 2) stop("m is less than 2")
	    if (length(eps) == 0) stop("invalid eps")
	    if (any(eps <= 0)) stop("invalid eps")
	    DNAME = deparse(substitute(x))
	    n = length(x)
	    k = length(eps)
	    cc = double(m+1)
	    cstan = double(m+1)
	    STATISTIC = matrix(0,m-1,k)
	    for(i in (1:k)) {
	        res = .C("bdstest_main", as.integer(n), as.integer(m),
	        	as.vector(x, mode="double"), as.vector(cc),
	        	cstan = as.vector(cstan), as.double(eps[i]),
	            as.integer(trace), PACKAGE="fSeries")
	        STATISTIC[,i] = res$cstan[2:m+1] }
	    colnames(STATISTIC) = eps
	    rownames(STATISTIC) = 2:m
	    PVAL = 2 * pnorm(-abs(STATISTIC))
	    colnames(PVAL) = eps
	    rownames(PVAL) = 2:m
	    METHOD = "BDS Test"
	    PARAMETER = list(m = 2:m, eps = eps)
	    structure(list(statistic = STATISTIC, p.value = PVAL,
			method = METHOD, data.name = DNAME, parameter = PARAMETER), 
	        class = "bdstest") }
	
	# Internal Function:
	# Source: tseries - Adrian Trapletti
	print.bdstest = 
	function(x, digits = 4, ...) {
	    if (!inherits(x, "bdstest")) stop("method is only for bdstest objects")
	    cat("\n\t", x$method, "\n\n")
	    cat("data: ", x$data.name, "\n\n")
	    if (!is.null(x$parameter)) {
	        cat("Embedding dimension = ",
	            format(round(x$parameter$m, digits)), sep = " ", "\n\n")
	        cat("Epsilon for close points = ",
	            format(round(x$parameter$eps, digits)), sep = " ", "\n\n")}
	    if (!is.null(x$statistic)) {
	        colnames(x$statistic) <-
	            round(as.numeric(colnames(x$statistic)), digits)
	        colnames(x$statistic) = paste("[",colnames(x$statistic),"]")
	        rownames(x$statistic) <-
	            round(as.numeric(rownames(x$statistic)), digits)
	        rownames(x$statistic) = paste("[",rownames(x$statistic),"]")
	        cat("Standard Normal = \n")
	        print(round(x$statistic, digits))
	        cat("\n")}
	    if (!is.null(x$p.value)) {
	        colnames(x$p.value) <-
	            round(as.numeric(colnames(x$p.value)), digits)
	        colnames(x$p.value) = paste("[",colnames(x$p.value),"]")
	        rownames(x$p.value) <-
	            round(as.numeric(rownames(x$p.value)), digits)
	        rownames(x$p.value) = paste("[",rownames(x$p.value),"]")
	        cat("p-value = \n")
	        print(round(x$p.value, digits))
	        cat("\n") }
	    cat("\n")
	    invisible(x) }
	
	# Return Value:
	bds.test(x, m, eps, trace)
}
  

# ------------------------------------------------------------------------------

	
wnnTest = 
function(x, lag = 1, qstar = 2, q = 10, range = 4, 
type = c("Chisq", "F"), scale = TRUE, ...) 
{	# A function implemented by Diethelm Wuertz

	# Description:
	#	White's Neural Network Test for Nonlinearity
	
	# Notes:
	#	This function is a slightly modified copy of Adrian Trapletti's
	#	contributed function from his 'tseries' package.
	
	# FUNCTION:
	
	# Transform:
	x = as.vector(x)
	
	# Internal Function:
	# Source: tseries - Adrian Trapletti
	white.test = 
	function(x, ...) {
		UseMethod("white.test") }
	
	# Internal Function:
	# Source: tseries - Adrian Trapletti
	white.test.default = 
	function(x, y, qstar = 2, q = 10, range = 4,
	type = c("Chisq","F"), scale = TRUE, ...) {
	    DNAME = paste(deparse(substitute(x)), "and", deparse(substitute(y)))
	    x = as.matrix(x)
	    y = as.matrix(y)
	    if (any(is.na(x))) stop("NAs in x")
	    if (any(is.na(y))) stop("NAs in y")
	    nin = dim(x)[2]
	    t = dim(x)[1]
	    if (dim(x)[1] != dim(y)[1]) 
	        stop("number of rows of x and y must match")
	    if (dim(x)[1] <= 0) 
	        stop("no observations in x and y")
	    if (dim(y)[2] > 1)
	        stop("handles only univariate outputs")
	    if (!missing(type) && !is.na(pmatch(type, "chisq"))) {
	        warning(paste(
	        	"value `chisq' for `type' is deprecated,",
	          	"use `Chisq' instead"))
	        type = "Chisq" }
	    else
	        type = match.arg(type)
	    if (scale) {
	        x = scale(x)
	        y = scale(y) }
	    xnam = paste("x[,", 1:nin, "]", sep="")
	    fmla = as.formula(paste("y~",paste(xnam,collapse= "+")))
	    rr = lm(fmla)
	    u = residuals(rr)
	    ssr0 = sum(u^2)
	    max = range/2
	    gamma = matrix(runif ((nin+1)*q,-max,max),nin+1,q)
	    phantom = (1+exp(-(cbind(rep(1,t),x)%*%gamma)))^(-1)
	    phantomstar = as.matrix(prcomp(phantom,scale=TRUE)$x[,2:(qstar+1)])
	    xnam2 = paste("phantomstar[,", 1:qstar, "]", sep="")
	    xnam2 = paste(xnam2,collapse="+")
	    fmla = as.formula(paste("u~",paste(paste(xnam,collapse= "+"),
			xnam2,sep="+")))
	    rr = lm(fmla)
	    v = residuals(rr)
	    ssr = sum(v^2)
	    if (type == "Chisq") {
	        STAT = t*log(ssr0/ssr)
	        PVAL = 1-pchisq(STAT,qstar)
	        PARAMETER = qstar
	        names(STAT) = "X-squared"
	        names(PARAMETER) = "df" }
	    else if (type == "F") {
	        STAT = ((ssr0-ssr)/qstar)/(ssr/(t-qstar-nin))
	        PVAL = 1-pf(STAT,qstar,t-qstar-nin)
	        PARAMETER = c(qstar,t-qstar-nin)
	        names(STAT) = "F"
	        names(PARAMETER) = c("df1","df2") }
	    else
	        stop("invalid type")
	    ARG = c(qstar,q,range,scale)
	    names(ARG) = c("qstar","q","range","scale")
	    METHOD = "White Neural Network Test"
	    structure(list(statistic = STAT, parameter = PARAMETER,
			p.value = PVAL, method = METHOD, data.name = DNAME,
			arguments = ARG),
			class = "htest") }
	
	# Internal Function:
	# Source: tseries - Adrian Trapletti
	white.test.ts = 
	function(x, lag = 1, qstar = 2, q = 10, range = 4,
	type = c("Chisq","F"), scale = TRUE, ...) {
	    if (!is.ts(x)) stop("method is only for time series")
	    if (NCOL(x) > 1) stop("x is not a vector or univariate time series")
	    if (any(is.na(x))) stop("NAs in x")
	    if (lag < 1)  stop("minimum lag is 1")
	    if (!missing(type) && !is.na(pmatch(type, "chisq"))) {
	        warning(paste(
	        	"value `chisq' for `type' is deprecated,",
	     		"use `Chisq' instead"))
	        type = "Chisq" }
	    else
	        type = match.arg(type)
	    DNAME = deparse(substitute(x))
	    t = length(x)
	    if (scale) x = scale(x)
	    y = embed(x, lag+1)
	    xnam = paste("y[,", 2:(lag+1), "]", sep="")
	    fmla = as.formula(paste("y[,1]~",paste(xnam,collapse= "+")))
	    rr = lm(fmla)
	    u = residuals(rr)
	    ssr0 = sum(u^2)
	    max = range/2
	    gamma = matrix(runif ((lag+1)*q,-max,max),lag+1,q)
	    phantom = (1+exp(-(cbind(rep(1,t-lag),y[,2:(lag+1)])%*%gamma)))^(-1)
	    phantomstar = as.matrix(prcomp(phantom,scale=TRUE)$x[,2:(qstar+1)])
	    xnam2 = paste("phantomstar[,", 1:qstar, "]", sep="")
	    xnam2 = paste(xnam2, collapse="+")
	    fmla = as.formula(paste("u~",paste(paste(xnam,collapse= "+"),
			xnam2,sep="+")))
	    rr = lm(fmla)
	    v = residuals(rr)
	    ssr = sum(v^2)
	    if (type == "Chisq") {
	        STAT = t*log(ssr0/ssr)
	        PVAL = 1-pchisq(STAT,qstar)
	        PARAMETER = qstar
	        names(STAT) = "X-squared"
	        names(PARAMETER) = "df" } 
	    else if (type == "F") {
	        STAT = ((ssr0-ssr)/qstar)/(ssr/(t-lag-qstar))
	        PVAL = 1-pf(STAT,qstar,t-lag-qstar)
	        PARAMETER = c(qstar,t-lag-qstar)
	        names(STAT) = "F"
	        names(PARAMETER) = c("df1","df2") }
	    else
	        stop("invalid type")
	    ARG = c(lag,qstar,q,range,scale)
	    names(ARG) = c("lag","qstar","q","range","scale")
	    METHOD = "White Neural Network Test"
	    structure(list(statistic = STAT, parameter = PARAMETER,
			p.value = PVAL, method = METHOD, data.name = DNAME,
			arguments = ARG),
			class = "htest") }
			
	# UseMethod("white.test")
	white.test(x = as.ts(x), lag = lag, qstar = qstar, q = q, 
	 	range = range, type = type[1], scale = scale, ...) 
}



# ------------------------------------------------------------------------------


tnnTest = 
function(x, lag = 1, type = c("Chisq", "F"), scale = TRUE, ...) 
{	# A function implemented by Diethelm Wuertz

	# Description:
	#	Teraesvirta's Neural Network Test for Nonlinearity
	
	# Notes:
	#	This function is a slightly modified copy of Adrian Trapletti's
	#	contributed function from his 'tseries' package.
	
	# FUNCTION:
	
	# Transform:
	x = as.vector(x)
	
	# Internal Function:
	# Source: tseries - Adrian Trapletti
	terasvirta.test = 
	function(x, ...) {
		UseMethod("terasvirta.test") }
	
	# Internal Function:
	# Source: tseries - Adrian Trapletti
	terasvirta.test.default = 
	function(x, y, type = c("Chisq", "F"), scale = TRUE, ...) {
	    DNAME = paste(deparse(substitute(x)), "and", deparse(substitute(y)))
	    x = as.matrix(x)
	    y = as.matrix(y)
	    if (any(is.na(x))) stop("NAs in x")
	    if (any(is.na(y))) stop("NAs in y")
	    nin = dim(x)[2]
	    if (nin < 1) stop("invalid x")
	    t = dim(x)[1]
	    if (dim(x)[1] != dim(y)[1]) stop("number of rows of x and y must match")
	    if (dim(x)[1] <= 0) stop("no observations in x and y")
	    if (dim(y)[2] > 1) stop("handles only univariate outputs")
	    if (!missing(type) && !is.na(pmatch(type, "chisq"))) {
	        warning(paste(
	        	"value `chisq' for `type' is deprecated,",
	      		"use `Chisq' instead"))
	        type = "Chisq" }
	    else
	        type = match.arg(type)
	    if (scale) {
	        x = scale(x)
	        y = scale(y) }
	    xnam = paste("x[,", 1:nin, "]", sep="")
	    fmla = as.formula(paste("y~",paste(xnam,collapse= "+")))
	    rr = lm(fmla)
	    u = residuals(rr)
	    ssr0 = sum(u^2)
	    xnam2 = NULL
	    m = 0
	    for(i in (1:nin)) {
	        for(j in (i:nin)) {
	            xnam2 = c(xnam2,paste("I(x[,",i,"]*x[,",j,"])",sep=""))
	            m = m+1 } }
	    xnam2 = paste(xnam2,collapse="+")
	    xnam3 = NULL
	    for(i in (1:nin)) {
	        for(j in (i:nin)) {
	            for(k in (j:nin)) {
	                xnam3 = c(xnam3, paste("I(x[,", i, "]*x[,", j, "]*x[,",
	                	k ,"])", sep=""))
	                m = m+1 } } }
	    xnam3 = paste(xnam3,collapse="+")
	    fmla = as.formula(paste("u~",paste(paste(xnam,collapse= "+"),
	   		xnam2,xnam3,sep="+")))
	    rr = lm(fmla)
	    v = residuals(rr)
	    ssr = sum(v^2)
	    if (type == "Chisq") {
	        STAT = t*log(ssr0/ssr)
	        PVAL = 1-pchisq(STAT,m)
	        PARAMETER = m
	        names(STAT) = "X-squared"
	        names(PARAMETER) = "df" }
	    else if (type == "F") {
	        STAT = ((ssr0-ssr)/m)/(ssr/(t-nin-m))
	        PVAL = 1-pf(STAT,m,t-nin-m)
	        PARAMETER = c(m,t-nin-m)
	        names(STAT) = "F"
	        names(PARAMETER) = c("df1","df2") }
	    else
	        stop("invalid type")
	    METHOD = "Teraesvirta Neural Network Test"
	    ARG = scale
	    names(ARG) = "scale"
	    structure(list(statistic = STAT, parameter = PARAMETER,
	 		p.value = PVAL, method = METHOD, data.name = DNAME,
	     	arguments = ARG), class = "htest") }
		     	
	# Internal Function:
	# Source: tseries - Adrian Trapletti
	terasvirta.test.ts = 
	function(x, lag = 1, type = c("Chisq", "F"), scale = TRUE, ...) {
	    if (!is.ts(x)) stop("method is only for time series")
	    if (NCOL(x) > 1) stop("x is not a vector or univariate time series")
	    if (any(is.na(x))) stop("NAs in x")
	    if (lag < 1) stop("minimum lag is 1")
	    if (!missing(type) && !is.na(pmatch(type, "chisq"))) {
	        warning(paste(
	        	"value `chisq' for `type' is deprecated,",
	      		"use `Chisq' instead"))
	        type = "Chisq" }
	    else
	        type = match.arg(type)
	    DNAME = deparse(substitute(x))
	    t = length(x)
	    if (scale) x = scale(x)
	    y = embed(x, lag+1)
	    xnam = paste("y[,", 2:(lag+1), "]", sep="")
	    fmla = as.formula(paste("y[,1]~",paste(xnam,collapse= "+")))
	    rr = lm(fmla)
	    u = residuals(rr)
	    ssr0 = sum(u^2)
	    xnam2 = NULL
	    m = 0
	    for(i in (1:lag)) {
	        for(j in (i:lag)) {
	            xnam2 = c(xnam2,paste("I(y[,",i+1,"]*y[,",j+1,"])",sep=""))
	            m = m+1 } }
	    xnam2 = paste(xnam2,collapse="+")
	    xnam3 = NULL
	    for(i in (1:lag)) {
	        for(j in (i:lag)) {
	            for(k in (j:lag)) {
	                xnam3 = c(xnam3, paste("I(y[,", i+1, "]*y[,", j+1,
	                	"]*y[,", k+1, "])", sep=""))
	                m = m+1 } } }
	    xnam3 = paste(xnam3,collapse="+")
	    fmla = as.formula(paste("u~",paste(paste(xnam,collapse= "+"),
	    	xnam2,xnam3,sep="+")))
	    rr = lm(fmla)
	    v = residuals(rr)
	    ssr = sum(v^2)
	    if (type == "Chisq") {
	        STAT = t*log(ssr0/ssr)
	        PVAL = 1-pchisq(STAT,m)
	        PARAMETER = m
	        names(STAT) = "X-squared"
	        names(PARAMETER) = "df" }
	    else if (type == "F") {
	        STAT = ((ssr0-ssr)/m)/(ssr/(t-lag-m))
	        PVAL = 1-pf(STAT,m,t-lag-m)
	        PARAMETER = c(m,t-lag-m)
	        names(STAT) = "F"
	        names(PARAMETER) = c("df1","df2") }
	    else
	        stop("invalid type")
	    METHOD = "Teraesvirta Neural Network Test"
	    ARG = c(lag,scale)
	    names(ARG) = c("lag","scale")
	    structure(list(statistic = STAT, parameter = PARAMETER,
			p.value = PVAL, method = METHOD, data.name = DNAME,
	    	arguments = ARG), class = "htest") }
	    	
	# UseMethod("terasvirta.test")
	terasvirta.test(x = as.ts(x), lag = lag, type = type[1], 
	 	scale = scale, ...) 
}


# ******************************************************************************

