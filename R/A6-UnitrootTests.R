
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
#   Diethelm Wuertz <wuertz@itp.phys.ethz.ch>
# for the code accessed (or partly included) from other R-ports:
#   R: see R's copyright and license file
#   ts: collected by Brian Ripley. See SOURCES
#   tseries: Compiled by Adrian Trapletti <a.trapletti@bluewin.ch>
#   fracdiff: S original by Chris Fraley <fraley@stat.washington.edu>
#     R-port: by Fritz Leisch <leisch@ci.tu-wien.ac.at>
#     since 2003-12: Martin Maechler
#   lmtest: Torsten Hothorn <Torsten.Hothorn@rzmail.uni-erlangen.de>
#     Achim Zeileis <zeileis@ci.tuwien.ac.at>
#     David Mitchell
#   mda: S original by Trevor Hastie & Robert Tibshirani
#     R port by Friedrich Leisch, Kurt Hornik and Brian D. Ripley
#   mgcv: Simon Wood <simon@stats.gla.ac.uk>
#   modreg: Brian Ripley and the R Core Team
#   polspline: Charles Kooperberg <clk@fhcrc.org>
#   nnet: S original by Venables & Ripley. 
#     R port by Brian Ripley <ripley@stats.ox.ac.uk>
#       following earlier work by Kurt Hornik and Albrecht Gebhardt


################################################################################
# CLASSES AND METHODS:     DESCRIPTION:
#  fURTEST                  S4: Class for Hypothesis Tests
#  print.fURTEST            S3: Print Method for Hypothesis Test
#  summary.fURTEST          S3: Summary Method for Hypothesis Test
# FUNCTION:                PART I:
#  adfTest                  Augmented Dickey - Fuller Unit Root Test
#  unitrootTest             ADF Unit Root Test using McKinnons Test Statistics
#   .unitrootADF            Internal Function called by 'unitrootTest'
# FUNCTION:                PART II: TSERIES UNIT ROOT / COINTEGRATION TESTS
#  tsadfTest                Augmented Dickey - Fuller Unit Root Test
#  tsppTest                 Phillips - Perron Unit Root Test
#  tskpssTest               KPSS - Unit Root Test for Stationarity
#  tspoTest                 Phillips - Ouliaris Cointegration Test
# FUNCTION:                PART III: URCA UNIT ROOT / COINTEGRATION TESTS
#  urersTest                Elliott-Rothenberg-Stock test for unit roots
#  urkpssTest               KPSS unit root test for stationarity
#  urppTest                 Phillips-Perron test for unit roots
#  urspTest                 Schmidt-Phillips test for unit roots
#  urzaTest                 Zivot-Andrews test for unit roots
################################################################################


require(methods)


# ------------------------------------------------------------------------------


setClass("fURTEST", 
    representation(
        call = "call",
        data = "data.frame",
        data.name = "character",
        test = "list",     
        title = "character",
        description = "character") )
     
        
# ------------------------------------------------------------------------------
   

urTest = 
function(x, method = c("unitroot", "adf", "tsadf", "tskpss", "tspp", "urers", 
"urkpss", "urpp", "ursp", "urza"), title = NULL, description = NULL, ...)  
{   # A function implemented by Diethelm Wuertz
    
    # FUNCTION:
    
    # Match Function:
    funTest = match.fun(paste(method[1], Test, sep = ""))
    
    # Test:
    ans = funTest(x = x, ...)
    
    # Add:
    if (!is.null(title)) ans@title = as.character(title)
    if (!is.null(description)) ans@description = as.character(description)
    
    # Return Value:
    ans
    
}
     

# ******************************************************************************


print.fURTEST =
function(x, ...)
{   # A function implemented by Diethelm Wuertz

    # Source:
    #   This function copies code from base:print.htest
    
    # FUNCTION:
       
    # Test slot:
    ans = x
    x = x@test
    
    # Title:
    cat("\nTitle:\n", ans@title, "\n", sep = "")
    
    # Call:
    cat("\nCall:\n", deparse(ans@call), "\n", sep = "")
    
    # Data Name:
    # cat("\nData Name:\n", ans@data.name, "\n", sep = "")  
    
    # Test Results:
    cat("\nTest Results:\n", sep = "")
    
    if (class(x)[2] == "htest") {
        # Tests from tseries package:
        # Statistic:
        if (!is.null(x$statistic))
            cat(paste("  ", names(x$statistic), ": ", 
                format(round(x$statistic, 4)), sep = ""), "\n")
        # Parameter:
        if (!is.null(x$parameter))
            cat(paste("  ", names(x$parameter), ": ", 
                format(round(x$parameter, 3)), sep = ""), "\n") 
        # P-Value:
        if (!is.null(x$p.value))
            cat(paste("  p-value: ", format.pval(x$p.value, 
                digits = 4), sep = ""), "\n")   
        # Alternative Hypothesis:
        if (!is.null(x$alternative)) {
            cat("  Alternative Hypothesis: ")
            if (!is.null(x$null.value)) {
                if (length(x$null.value) == 1) {
                    alt.char <-
                        switch(x$alternative,
                            two.sided = "not equal to",
                            less = "less than",
                            greater = "greater than")
                    cat("true", names(x$null.value), "is", alt.char,
                    x$null.value, "\n")
                } else {
                    cat(x$alternative, "\nnull values:\n")
                    print(x$null.value, ...)
                }
            } else {
                cat(x$alternative, "\n")
            }
        }  
        # Confidence Interval:
        if (!is.null(x$conf.int)) {
            cat("  ")
            cat(format(100 * attr(x$conf.int, "conf.level")),
                "percent confidence interval:\n",
                format(c(x$conf.int[1], x$conf.int[2])), "\n")
        }  
        # Sample estimates:
        if (!is.null(x$estimate)) {
            cat("  Sample estimates:\n")
            print(x$estimate, ...)
        }
    } else {
        # Tests from urca package:
        # Extract parameters from test slot:
        test.name = as.character(ans@call)[1]
        if (test.name == "urersTest") {
            cat("  Type of Test:", attr(x, "type"))
            cat("\n  Model:", attr(x, "model"))
            cat("\n  Lags:", attr(x, "lag"))
            cat("\n  Critical Values:", attr(x, "cval"))
            cat(" for", colnames(attr(x, "cval")))
            cat("\n  Test Statistic:", attr(x, "teststat"))      
        }  
        if (test.name == "urkpssTest") {
            cat("  Type of Test:", attr(x, "type"))
            cat("\n  Lags:", attr(x, "lag"))
            cat("\n  Critical Values:", attr(x, "cval"))
            cat(" for", colnames(attr(x, "cval")))
            cat("\n  Test Statistic:", attr(x, "teststat"))      
        }  
        if (test.name == "urppTest") {
            cat("  Type of Test:", attr(x, "type"))
            cat("\n  Model:", attr(x, "model"))
            cat("\n  Lags:", attr(x, "lag"))
            cat("\n  Critical Values:", attr(x, "cval"))
            cat(" for", colnames(attr(x, "cval")))
            cat("\n  Test Statistic:", attr(x, "teststat"))  
            cat("\n  Auxiliary Z Statistic:", attr(x, "auxstat"))    
        }  
        if (test.name == "urspTest") {      
            cat("")
        }  
        if (test.name == "urzaTest") {
            cat("  Model:", attr(x, "model"))
            cat("\n  Lags:", attr(x, "lag"))
            cat("\n  Critical Values:", attr(x, "cval"))
            cat(" for", colnames(attr(x, "cval")))
            cat("\n  Test Statistic:", attr(x, "teststat"))  
            cat("\n  Break Point:", attr(x, "bpoint"))               
            cat("")
        }  
        cat("\n")  
    }

    # Description:
    cat("\nDescription:\n", ans@description, sep = "")   
    cat("\n\n")
    
    # Return Value:
    invisible(x)
}


# ------------------------------------------------------------------------------


# plot.fURTEST =
# function(x, y, ...)
# {  # A function implemented by Diethelm Wuertz

#    # FUNCTION:
    
#    # Return Value:
#    invisible(x)
# }


# ------------------------------------------------------------------------------


summary.fURTEST =
function(object, ...)
{   # A function implemented by Diethelm Wuertz

    # FUNCTION:
    
    # Print:
    print(object, ...)
    
    # Settings:
    # Test slot:
    ans = object
    object = object@test
    
    # Internal Function:
    
    # Summary:
    if (class(object)[2] == "htest") {
        return(invisible())
    } else {
        # Tests from urca package:
        # Extract parameters from test slot:
        test.name = as.character(ans@call)[1]
        if (test.name == "urersTest" || test.name == "urzaTest") {
            regression = attr(object, "testreg")
            regression$call = "Regression Details"
        }  
        print(regression)
    }
    
    # Return Value:
    invisible()
}
        

################################################################################
# PART I:


adfTest = 
function(x, type = c("nc", "c", "ct"), lags = 1)
{   # A function implemented by Diethelm Wuertz

    # Description:     
    #   Tests the null hypothesis of a unit root in y.

    # Arguments:
    #   x - numeric vector
    #   type - specifies the regression model to be estimatied and the 
    #       null hypothesis, "nc" no constant and no trend, "c" add 
    #       constant, "ct" add constant and trend.
    #   lags - specifies the number of lagged differences of x to be 
    #       included in the regression model. If 'lags' = h, a term 
    #       sum_{i=1}^{h-1} beta_i * diff(x)_(t-h) is added to the 
    #       regression equation. 

    # Value:      
    #   A list of class "htest" containing the following components:
    #   statistic - the value of the test statistic (t-statistic)
    #   parameter - the number of lags.
    #   p.value - the p-value of the test
    #   method - a character string indicating what type of test was performed
    #   data.name - a character string giving the name of the data y

    # Reference:   
    #   S. E. SAID and D. A. DICKEY (1984): Testing for Unit Roots in 
    #   Autoregressive-Moving Average Models of Unlag.diffnown Order. 
    #   Biometrika 71, 599–607.
    
    # Source:
    #   This function is an augmented version of Adrian Trapletti's
    #   function adf.test() which considers type "ct" only. We have added
    #   the types "c" and "nc" together with the appropriate statistics.

    # Check Arguments:
    if (ncol(as.matrix(x)) > 1) 
        stop("x is not a vector or univariate time series")
    if (any(is.na(x))) 
        stop("NAs in x")
    if (lags < 0) 
        stop("lags negative")
    
    # Settings:
    doprint = FALSE
    CALL = match.call()
    DNAME = deparse(substitute(x))  
    type = type[1]  
    x.name = deparse(substitute(x))
    lags = lags + 1
    y = diff(x)
    n = length(y)
    z = embed(y, lags)
    y.diff = z[, 1]
    y.lag.1 = x[lags:n]
    tt = lags:n
   
    # Regression:             
    if (lags > 1) {
        y.diff.lag = z[,2:lags]
        if (type == "nc"){
            res = lm(y.diff ~ y.lag.1 - 1 + y.diff.lag) }
        if (type == "c"){
            res = lm(y.diff ~ y.lag.1 + 1 +  y.diff.lag) }
        if (type == "ct") {
            res = lm(y.diff ~ y.lag.1 + 1 + tt + y.diff.lag) } }
    else {
        if (type == "nc") {
            res = lm(y.diff ~ y.lag.1 - 1) }
        if (type == "c"){
            res = lm(y.diff ~  y.lag.1 + 1) }
        if (type == "ct") {
            res = lm(y.diff ~ y.lag.1 + 1  + tt)  } }
    
    # Regression Summary:
    res.sum = summary(res)
    if (doprint) print(res.sum)
    
    # Statistic:
    if (type == "nc") coefNum = 1 else coefNum = 2
    STAT = res.sum$coefficients[coefNum, 1] / res.sum$coefficients[coefNum, 2]
        
    # Tables:
    if (type == "nc")
        table = cbind(
            c(-2.66, -2.26, -1.95, -1.60, +0.92, +1.33, +1.70, +2.16),
            c(-2.62, -2.25, -1.95, -1.61, +0.91, +1.31, +1.66, +2.08),
            c(-2.60, -2.24, -1.95, -1.61, +0.90, +1.29, +1.64, +2.03),
            c(-2.58, -2.23, -1.95, -1.62, +0.89, +1.29, +1.63, +2.01),
            c(-2.58, -2.23, -1.95, -1.62, +0.89, +1.28, +1.62, +2.00),
            c(-2.58, -2.23, -1.95, -1.62, +0.89, +1.28, +1.62, +2.00))
    if (type == "c")
        table = cbind(
            c(-3.75, -3.33, -3.00, -2.63, -0.37, +0.00, +0.34, +0.72),
            c(-3.58, -3.22, -2.93, -2.60, -0.40, -0.03, +0.29, +0.66),
            c(-3.51, -3.17, -2.89, -2.58, -0.42, -0.05, +0.26, +0.63),
            c(-3.46, -3.14, -2.88, -2.57, -0.42, -0.06, +0.24, +0.62),
            c(-3.44, -3.13, -2.87, -2.57, -0.43, -0.07, +0.24, +0.61),
            c(-3.43, -3.12, -2.86, -2.57, -0.44, -0.07, +0.23, +0.60))
    if (type == "ct")
        table = cbind(
            c(-4.38, -3.95, -3.60, -3.24, -1.14, -0.80, -0.50, -0.15),
            c(-4.15, -3.80, -3.50, -3.18, -1.19, -0.87, -0.58, -0.24),
            c(-4.04, -3.73, -3.45, -3.15, -1.22, -0.90, -0.62, -0.28),
            c(-3.99, -3.69, -3.43, -3.13, -1.23, -0.92, -0.64, -0.31),
            c(-3.98, -3.68, -3.42, -3.13, -1.24, -0.93, -0.65, -0.32),
            c(-3.96, -3.66, -3.41, -3.12, -1.25, -0.94, -0.66, -0.33))
            
    # Critical values:  
    table = t(table)
    tablen = dim(table)[2]
    tableT = c(25, 50, 100, 250, 500, 1e+05)
    tablep = c(0.01, 0.025, 0.05, 0.1, 0.9, 0.95, 0.975, 0.99)
    tableipl = numeric(tablen)
    for (i in (1:tablen)) 
    	tableipl[i] = approx(tableT, table[, i], n, rule = 2)$y   
    PVAL = approx(tableipl, tablep, STAT, rule = 2)$y 
    if (is.na(approx(tableipl, tablep, STAT, rule = 1)$y)) {
        if (PVAL == min(tablep)) {
            warning("p-value smaller than printed p-value") 
        } else {
            warning("p-value greater than printed p-value") 
        }
	}
	        
    # Names:
    PARAMETER = lags - 1
    names(PARAMETER) = "Lag order"
    METHOD = "Augmented Dickey-Fuller Test"
    names(STAT) = "Dickey-Fuller"
   
    # Return Value:
    test = list(
        statistic = STAT, 
        parameter = PARAMETER, 
        p.value = PVAL, 
        method = METHOD, 
        data.name = DNAME)   
    class(test) = c("list", "htest")  
        
    # Add:
    title = test$method
    description = date()
    
    # Return Value:
    new("fURTEST", 
        call = CALL,
        data = as.data.frame(x),
        data.name = test$DNAME,
        test = test, 
        title = as.character(title),
        description = as.character(description)
        )   
}

  
# ------------------------------------------------------------------------------


.unitrootADF = 
function(x, trend, statistic, lags)
{   # A function implemented by Diethelm Wuertz

    # Description:     
    #   Tests the null hypothesis of a unit root in y.

    # Arguments:
    #   x - numeric vector
    #   trend - specifies the regression model to be estimatied and the 
    #       null hypothesis, "nc" no constant and no trend, "c" add 
    #       constant, "ct" add constant and trend.
    #   lags - specifies the number of lagged differences of x to be 
    #       included in the regression model. If 'lags' = h, a term 
    #       sum_{i=1}^{h-1} beta_i * diff(x)_(t-h) is added to the 
    #       regression equation. 

    # Value:      
    #   A list with class "htest" containing the following components:
    #   statistic - the value of the test statistic (t-statistic)
    #   parameter - the number of lags.
    #   p.value - the p-value of the test
    #   method - a character string indicating what "trend" type of 
    #       the test was performed
    #   data.name - a character string giving the name of the data y

    # Reference:   
    #   Said S.E., Dickey D.A. (1984): Testing for Unit Roots in 
    #   Autoregressive-Moving Average Models of Unlag.diffnown Order. 
    #   Biometrika 71, 599–-607.
    
    # Source:
    #   This function is an augmented version of Adrian Trapletti's
    #   function adf.test() which considers trend "ct" only. We have added
    #   the trend types "c" and "nc" together with the appropriate statistics.

    # FUNCTION:
    
    # Check Arguments:
    if (ncol(as.matrix(x)) > 1) 
        stop("x is not a vector or univariate time series")
    if (any(is.na(x))) 
        stop("NAs in x")
    if (lags < 0) 
        stop("lags negative")
    
    # Settings:
    DNAME = deparse(substitute(x))   
    x.name = deparse(substitute(x))
    lags = lags + 1
    y = diff(x)
    n = length(y)
    z = embed(y, lags)
    y.diff = z[, 1]
    y.lag.1 = x[lags:n]
    tt = lags:n
   
    # Regression:   
    type = trend[1]          
    if (lags > 1) {
        y.diff.lag = z[,2:lags]
        if (type == "nc"){
            res = lm(y.diff ~ y.lag.1 - 1 + y.diff.lag) }
        if (type == "c"){
            res = lm(y.diff ~ y.lag.1 + 1 +  y.diff.lag) }
        if (type == "ct") {
            res = lm(y.diff ~ y.lag.1 + 1 + tt + y.diff.lag) } 
        if (type == "ctt") {
            res = lm(y.diff ~ y.lag.1 + 1 + tt + tt^2 + y.diff.lag) } }
    else {
        if (type == "nc") {
            res = lm(y.diff ~ y.lag.1 - 1) }
        if (type == "c"){
            res = lm(y.diff ~  y.lag.1 + 1) }
        if (type == "ct") {
            res = lm(y.diff ~ y.lag.1 + 1  + tt)  } 
        if (type == "ctt") {
            res = lm(y.diff ~ y.lag.1 + 1 + tt + tt^2) } }
    
    # Statistic:
    res.sum = summary(res)
    if (type == "nc") coefNum = 1 else coefNum = 2
    STAT = res.sum$coefficients[coefNum, 1] / res.sum$coefficients[coefNum, 2]
        
    if (type == "nc") { itv = 1 }
    if (type == "c")  { itv = 2 }
    if (type == "ct") { itv = 3 }
    if (type == "ctt"){ itv = 4 }
    
    # P Value:
    if (statistic == "t") itt = 1
    if (statistic == "n") itt = 2
    PVAL = .urcval(arg = STAT, nobs = n, niv = 1, itt = itt, itv = itv, nc = 2)
            
    # Names:
    PARAMETER = lags - 1
    names(PARAMETER) = "Lag order"
    METHOD = "Augmented Dickey-Fuller Test"
    names(STAT) = "Dickey-Fuller"
   
    # Return Value:
    list(
        statistic = STAT, 
        parameter = PARAMETER, 
        p.value = PVAL, 
        method = METHOD, 
        data.name = DNAME, 
        regression = res.sum, 
        class = "htest")     
}


# ------------------------------------------------------------------------------


unitrootTest =
function(x, trend = c("nc", "c", "ct"), statistic = c("t", "n"), 
method = "adf", lags = 1)
{   # A function implemented by Diethelm Wuertz

    # Function:
    
    # Convert to numeric vector:
    if (class(x) == "timeSeries") x = x@Data
    
    # Settings:
    CALL = match.call()
    
    # Test:
    test = .unitrootADF(x = x, trend = trend[1], statistic = statistic[1], 
        lags = lags) 
    class(test) = c("list", "htest")  
        
    # Add:
    title = test$method
    description = date()
    
    # Return Value:
    new("fURTEST", 
        call = CALL,
        data = as.data.frame(x),
        data.name = test$DNAME,
        test = test, 
        title = as.character(title),
        description = as.character(description)
        )   
}

    

################################################################################
# PART II:

#
# PACKAGE DESCRIPTION:
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
#


tsadfTest = 
function(x, alternative = c("stationary", "explosive"), 
k = trunc((length(x)-1)^(1/3)))
{   # A function implemented by Diethelm Wuertz

    # Description:
    #   Augmented Dickey - Fuller Unit Root Test
    
    # FUNCTION:
    
    # Call:  
    CALL = match.call()
    alternative = alternative[1]
    
    # Internal Function:
    # Source: tseries - Adrian Trapletti
    .adf.test = function(x, alternative, k) {
        if (NCOL(x) > 1) stop("x is not a vector or univariate time series")
        if (any(is.na(x))) stop("NAs in x")
        if (k < 0) stop("k negative")
        # alternative = match.arg(alternative)
        DNAME = deparse(substitute(x))
        k = k+1
        y = diff(x)
        n = length(y)
        z = embed(y, k)
        yt = z[,1]
        xt1 = x[k:n]
        tt = k:n
        if (k > 1) {
            yt1 = z[,2:k] 
            res = lm(yt ~ xt1 + 1 + tt + yt1) 
        } else {
            res = lm(yt ~ xt1 + 1 + tt) }
        res.sum = summary(res)
        STAT = res.sum$coefficients[2,1] / res.sum$coefficients[2,2]
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
        tableT = c(25, 50, 100, 250, 500, 100000)
        tablep = c(0.01, 0.025, 0.05, 0.10, 0.90, 0.95, 0.975, 0.99)
        tableipl = numeric(tablen)
        for(i in (1:tablen))
            tableipl[i] = approx(tableT, table[, i], n, rule=2)$y
        interpol = approx(tableipl, tablep, STAT, rule=2)$y
        if (is.na(approx(tableipl, tablep, STAT, rule=1)$y))
            if (interpol == min(tablep))
                warning("p-value smaller than printed p-value")
            else
                warning("p-value greater than printed p-value")
        if (alternative == "stationary")
            PVAL = interpol
        else if (alternative == "explosive")
            PVAL = 1 - interpol
        else 
            stop("irregular alternative") 
        PARAMETER = k-1
        METHOD = "Augmented Dickey-Fuller Test"
        names(STAT) = "Dickey-Fuller"
        names(PARAMETER) = "Lag order"
        structure(list(statistic = STAT, parameter = PARAMETER,
            alternative = alternative, p.value = PVAL, method = METHOD,
            data.name = DNAME), 
            class = "htest") 
        }
    
    # Object:
    test = .adf.test(x = x, alternative = alternative[1], k = k)
    class(test) = c("list", "htest")
    
    # Add:
    title = test$method
    description = date()
    
    # Return Value:
    new("fURTEST", 
        call = CALL,
        data = as.data.frame(x),
        data.name = test$data.name,
        test = test, 
        title = as.character(title),
        description = as.character(description)
        )   
}


# ------------------------------------------------------------------------------

    
tsppTest = 
function(x, alternative = c("stationary", "explosive"),
type = c("Z(alpha)", "Z(t_alpha)"), lshort = TRUE)
{   # A function implemented by Diethelm Wuertz

    # Description:
    #   Phillips - Perron Unit Root Test
    
    # FUNCTION:
    
    # Call:  
    CALL = match.call()
    alternative = alternative[1]
    type = type[1]
    
    # Internal Function:
    # Source: tseries - Adrian Trapletti
    .pp.test = function(x, alternative, type, lshort) {
        if (NCOL(x) > 1) stop("x is not a vector or univariate time series")
        # type = match.arg(type)
        # alternative = match.arg(alternative)
        DNAME = deparse(substitute(x))
        z = embed(x, 2)
        yt = z[,1]
        yt1 = z[,2]
        n = length(yt)
        tt = (1:n)-n/2
        res = lm(yt ~ 1 + tt + yt1)
        if (res$rank < 3) stop("Singularities in regression")
        res.sum = summary(res)
        u = residuals(res)
        ssqru = sum(u^2)/n
        if (lshort) 
            l = trunc(4*(n/100)^0.25)
        else 
            l = trunc(12*(n/100)^0.25)
        ssqrtl = .C("R_pp_sum", as.vector(u, mode = "double"),
            as.integer(n), as.integer(l), ssqrtl = as.double(ssqru),
            PACKAGE = "fSeries")$ssqrtl
        n2 = n^2
        trm1 = n2*(n2-1)*sum(yt1^2)/12
        trm2 = n*sum(yt1*(1:n))^2
        trm3 = n*(n+1)*sum(yt1*(1:n))*sum(yt1)
        trm4 = (n*(n+1)*(2*n+1)*sum(yt1)^2)/6
        Dx = trm1-trm2+trm3-trm4
        if (type == "Z(alpha)") {
            alpha = res.sum$coefficients[3,1]
            STAT = n*(alpha-1)-(n^6)/(24*Dx)*(ssqrtl-ssqru)
            table = cbind(
                c(22.5, 25.7, 27.4, 28.4, 28.9, 29.5),
                c(19.9, 22.4, 23.6, 24.4, 24.8, 25.1),
                c(17.9, 19.8, 20.7, 21.3, 21.5, 21.8),
                c(15.6, 16.8, 17.5, 18.0, 18.1, 18.3),
                c(3.66, 3.71, 3.74, 3.75, 3.76, 3.77),
                c(2.51, 2.60, 2.62, 2.64, 2.65, 2.66),
                c(1.53, 1.66, 1.73, 1.78, 1.78, 1.79),
                c(0.43, 0.65, 0.75, 0.82, 0.84, 0.87)) 
        } else if (type == "Z(t_alpha)") {
            tstat <- (res.sum$coefficients[3,1]-1)/res.sum$coefficients[3,2]
            STAT = sqrt(ssqru)/sqrt(ssqrtl)*tstat-(n^3) /
                (4*sqrt(3)*sqrt(Dx)*sqrt(ssqrtl))*(ssqrtl-ssqru)
            table = cbind(
                c(4.38, 4.15, 4.04, 3.99, 3.98, 3.96),
                c(3.95, 3.80, 3.73, 3.69, 3.68, 3.66),
                c(3.60, 3.50, 3.45, 3.43, 3.42, 3.41),
                c(3.24, 3.18, 3.15, 3.13, 3.13, 3.12),
                c(1.14, 1.19, 1.22, 1.23, 1.24, 1.25),
                c(0.80, 0.87, 0.90, 0.92, 0.93, 0.94),
                c(0.50, 0.58, 0.62, 0.64, 0.65, 0.66),
                c(0.15, 0.24, 0.28, 0.31, 0.32, 0.33)) 
        } else {
            stop("Irregular type") }
        table = -table
        tablen = dim(table)[2]
        tableT = c(25, 50, 100, 250, 500, 100000)
        tablep = c(0.01, 0.025, 0.05, 0.10, 0.90, 0.95, 0.975, 0.99)
        tableipl = numeric(tablen)
        for (i in (1:tablen))
            tableipl[i] = approx(tableT, table[, i], n, rule=2)$y
        interpol = approx(tableipl, tablep, STAT, rule=2)$y
        if (is.na(approx(tableipl, tablep, STAT, rule=1)$y))
            if (interpol == min(tablep)) {
                warning("p-value smaller than printed p-value")
            } else {
                warning("p-value greater than printed p-value") }
        if (alternative == "stationary") {
            PVAL = interpol
        } else if (alternative == "explosive") {
            PVAL = 1 - interpol
        } else {
            stop("irregular alternative") }
        PARAMETER = l
        METHOD = "Phillips-Perron Unit Root Test"
        names(STAT) = paste("Dickey-Fuller", type)
        names(PARAMETER) = "Truncation lag parameter"
        structure(list(statistic = STAT, parameter = PARAMETER,
            alternative = alternative, p.value = PVAL, method = METHOD,
            data.name = DNAME), class = "htest") }
    
    # Object:
    test = .pp.test(x = x, alternative = alternative, type = type,
        lshort = lshort)
    class(test) = c("list", "htest")
    
    # Add:
    title = test$method
    description = description = date()
    
    # Return Value:
    new("fURTEST", 
        call = CALL,
        data = as.data.frame(x),
        data.name = test$data.name,
        test = test, 
        title = as.character(title),
        description = as.character(description)
        )   
}


# ------------------------------------------------------------------------------


tspoTest = 
function(x, demean = TRUE, lshort = TRUE)
{   # A function implemented by Diethelm Wuertz

    # Description:
    #   Phillips - Ouliaris Cointegration Test
    
    # FUNCTION:
    
    # Call:  
    CALL = match.call()
    
    # Internal Function:
    # Source: tseries - Adrian Trapletti
    .po.test = function(x, demean, lshort) {
        if (NCOL(x) <= 1) stop("x is not a matrix or multivariate time series")
        DNAME = deparse(substitute(x))
        x = as.matrix(x)
        dimx = ncol(x)
        if (dimx > 6) stop("no critical values for this dimension")
        if (demean) 
            res = lm(x[,1]~x[,-1])
        else 
            res = lm(x[,1]~x[,-1]-1)
        z = embed(residuals(res), 2)
        ut = z[,1]
        ut1 = z[,2]
        n = length(ut)
        res = lm(ut ~ ut1 - 1)
        if (res$rank < 1) stop("Singularities in regression")
        res.sum = summary(res)
        k = residuals(res)
        ssqrk = sum(k^2)/n
        if (lshort) 
            l = trunc(n/100)
        else 
            l = trunc(n/30)
        ssqrtl = .C("R_pp_sum", as.vector(k, mode="double"),
            as.integer(n), as.integer(l), ssqrtl=as.double(ssqrk),
            PACKAGE="fSeries")$ssqrtl
        alpha = res.sum$coefficients[1,1]
        STAT = n*(alpha-1)-0.5*n^2*(ssqrtl-ssqrk)/(sum(ut1^2))
        if (demean) {
            table = cbind(
                c(28.32, 34.17, 41.13, 47.51, 52.17),
                c(23.81, 29.74, 35.71, 41.64, 46.53),
                c(20.49, 26.09, 32.06, 37.15, 41.94),
                c(18.48, 23.87, 29.51, 34.71, 39.11),
                c(17.04, 22.19, 27.58, 32.74, 37.01),
                c(15.93, 21.04, 26.23, 31.15, 35.48),
                c(14.91, 19.95, 25.05, 29.88, 34.20)) 
        } else {
            table = cbind(
                c(22.83, 29.27, 36.16, 42.87, 48.52),
                c(18.89, 25.21, 31.54, 37.48, 42.55),
                c(15.64, 21.48, 27.85, 33.48, 38.09),
                c(13.81, 19.61, 25.52, 30.93, 35.51),
                c(12.54, 18.18, 23.92, 28.85, 33.80),
                c(11.57, 17.01, 22.62, 27.40, 32.27),
                c(10.74, 16.02, 21.53, 26.17, 30.90)) 
        }
        table = -table
        tablep = c(0.01, 0.025, 0.05, 0.075, 0.10, 0.125, 0.15)
        PVAL = approx(table[dimx-1,], tablep, STAT, rule=2)$y   
        if (is.na(approx(table[dimx-1, ], tablep, STAT, rule = 1)$y))
            if (PVAL == min(tablep))
                warning("p-value smaller than printed p-value")
            else
                warning("p-value greater than printed p-value") 
        PARAMETER = l
        METHOD = "Phillips-Ouliaris Cointegration Test"
        if (demean) names(STAT) = "Phillips-Ouliaris demeaned"
        else names(STAT) = "Phillips-Ouliaris standard"
        names(PARAMETER) = "Truncation lag parameter"
        structure(list(statistic = STAT, parameter = PARAMETER,
            p.value = PVAL, method = METHOD, data.name = DNAME),
            class = "htest") }
               
    # Object:
    test = .po.test(x = x, demean = demean, lshort = lshort)
    class(test) = c("list", "htest")
    
    # Add:
    title = test$method
    description = date()
    
    # Return Value:
    new("fURTEST", 
        call = CALL,
        data = as.data.frame(x),
        data.name = test$data.name,
        test = test, 
        title = as.character(title),
        description = as.character(description)
        )   
}


# ------------------------------------------------------------------------------


tskpssTest = 
function(x, nullhyp = c("level", "trend"), lshort = TRUE)
{   # A function implemented by Diethelm Wuertz

    # Description:
    #   KPSS - Unit Root Test for Stationarity
    
    # FUNCTION:
    
    # Call:  
    CALL = match.call()
    nullhyp = nullhyp[1]
    
    # Internal Function:
    # Source: tseries - Adrian Trapletti
    .kpss.test = function(x, nullhyp, lshort) {
        if (NCOL(x) > 1) 
        	stop("x is not a vector or univariate time series")
        DNAME = deparse(substitute(x))
        n = length(x)
        if (nullhyp == "trend") {
            t = 1:n
            e = residuals(lm(x ~ t))
            table = c(0.216, 0.176, 0.146, 0.119) 
        } else if (nullhyp == "level") {
            e = residuals(lm(x ~ 1))
            table = c(0.739, 0.574, 0.463, 0.347) 
        }
        tablep = c(0.01, 0.025, 0.05, 0.10)
        s = cumsum(e)
        eta = sum(s^2)/(n^2)
        s2 = sum(e^2)/n
        if (lshort) {
          	l = trunc(3*sqrt(n)/13)
    	} else {
	    	l = trunc(10*sqrt(n)/14)
    	}
        s2 = .C("R_pp_sum",
             as.vector(e, mode = "double"),
             as.integer(n),
             as.integer(l),
             s2 = as.double(s2),
             PACKAGE = "fSeries")$s2
        STAT = eta/s2
        PVAL = approx(table, tablep, STAT, rule = 2)$y
        if (is.na(approx(table, tablep, STAT, rule = 1)$y)) {
            if (PVAL == min(tablep)) {
                warning("p-value smaller than printed p-value")
        	} else {
                warning("p-value greater than printed p-value") 
        	}
    	}
        PARAMETER = l
        METHOD = paste("KPSS Test for", nullhyp, "Stationarity")
        names(STAT) = paste("KPSS", nullhyp)
        names(PARAMETER) = "Truncation lag parameter"
        ans = structure(list(statistic = STAT, parameter = PARAMETER,
            p.value = PVAL, method = METHOD, data.name = DNAME),
            class = "htest") 
        # Return Value:
        ans
    }

    # Object:
    test = .kpss.test(x = x, nullhyp = nullhyp, lshort = lshort)
    class(test) = c("list", "htest")

    # Add:
    title = test$method
    description = date()
    
    # Return Value:
    new("fURTEST", 
        call = CALL,
        data = as.data.frame(x),
        data.name = test$data.name,
        test = test, 
        title = as.character(title),
        description = as.character(description)
        )   
}


################################################################################
# PART III:


#
# Package: urca
# Version: 0.5-4
# Date: 2004-05-13
# Title: Unit root and cointegration tests for time series data
# Author: Bernhard Pfaff <bernhard.pfaff@drkw.com>
# Maintainer: Bernhard Pfaff <bernhard.pfaff@drkw.com>
# Depends: R (>= 1.8.1), methods, nlme
# Description: Unit root and cointegration tests encountered in applied
#        econometric analysis are implemented.
# License: GPL version 2 or newer
# URL: http://www.r-project.org
# Packaged: Mon May 17 23:03:59 2004; bp
# Built: R 1.9.0; ; 2004-05-18 21:58:08; windows
#

#
# hTest Value:
#   VARABLE:        OUTPUT:
#   method
#   data.name       Name of the data set
#   statistic
#   parameter
#   p.value         p.value =
#   alternative     alternative hypothesis
#   null.value
#   method
#   conf.int        percent confidence interval:
#   estimate        sample estimates:
#


################################################################################
# FUNCTION:            DESCRIPTOION:
#  urersTest            Elliott-Rothenberg-Stock test for unit roots
#  urkpssTest           KPSS unit root test for stationarity
#  urppTest             Phillips-Perron test for unit roots
#  urspTest             Schmidt-Phillips test for unit roots
#  urzaTest             Zivot-Andrews test for unit roots
################################################################################


urersTest =
function(x, type = c("DF-GLS", "P-test"), model = c("constant", "trend"),
lag.max = 4)
{   # A function implemented by Diethelm Wuertz
    
    # Description:  
    #   Elliott-Rothenberg-Stock test for unit roots   

    # NoteS:
    #   Requires "urca" which is not part of this distribution
    
    # FUNCTION:
    
    # Settings:
    CALL = match.call()
    method = paste(model[1], "model  and test type", type[1])
    data.name = as.character(deparse(substitute(x)))
    
    # Test:
    test = .urers(x, type = type[1], model = model[1], lag.max = lag.max)    
    class(test) = c("list", "htest")
    
    # Add:
    title = "Elliott-Rothenberg-Stock Unit Root Test"
    description = date()
    
    # Return Value:
    new("fURTEST", 
        call = CALL,
        data = as.data.frame(x),
        data.name = data.name, 
        test = test, 
        title = as.character(title),
        description = as.character(description)
        )   
}
    

# ------------------------------------------------------------------------------

    
urkpssTest =
function(x, type = c("mu", "tau"), lags = c("short", "long", "nil"),
use.lag = NULL)
{   # A function implemented by Diethelm Wuertz
    
    # Description:  
    #   KPSS unit root test for stationarity

    # Note:
    #   Requires "urca" which is not part of this distribution
    
    # FUNCTION:
    
    # Settings:
    CALL = match.call()
    data.name = as.character(deparse(substitute(x)))
    
    # Test:
    test = .urkpss(x, type = type[1], lags = lags[1], use.lag = use.lag)
    class(test) = c("list", "ur.kpss")
  
    # Add:
    title = "KPSS Unit Root Test"
    description = date
    
    # Return Value:
    new("fURTEST", 
        call = CALL,
        data = as.data.frame(x),
        data.name = data.name, 
        test = test, 
        title = as.character(title),
        description = as.character(description)
        )   
}   


# ------------------------------------------------------------------------------


urppTest =
function(x, type = c("Z-alpha", "Z-tau"), model = c("constant", "trend"),
lags = c("short", "long"))
{   # A function implemented by Diethelm Wuertz
    
    # Description:  
    #   Phillips-Perron test for unit roots

    # Note:
    #   Requires "urca" which is not part of this distribution
    
    # FUNCTION:
    
    # Settings:
    CALL = match.call()
    data.name = as.character(deparse(substitute(x)))
    use.lag = NULL
    if (is.numeric(lags)) use.lag = lags
    
    # Test:
    test = .urpp(x, type = type[1], model = model[1], lags = lags[1], 
        use.lag = use.lag)
    class(test) = c("list", "ur.pp")
    
    # Add:
    title = "Phillips-Perron Unit Root Test"
    description = date()
    
    # Return Value:
    new("fURTEST", 
        call = CALL,
        data = as.data.frame(x),
        data.name = data.name, 
        test = test, 
        title = as.character(title),
        description = as.character(description)
        )   
}   
    

# ------------------------------------------------------------------------------


urspTest =
function(x, type = c("tau", "rho"), pol.deg = c(1, 2, 3, 4),
signif = c(0.01, 0.05, 0.10))
{   # A function implemented by Diethelm Wuertz
    
    # Description:  
    #   Schmidt-Phillips test for unit roots

    # Note:
    #   Requires "urca" which is not part of this distribution
    
    # FUNCTION:
    
    # Settings:
    CALL = match.call()
    data.name = as.character(deparse(substitute(x)))
    
    # Test:
    test = .ursp(x, type = type[1], pol.deg = pol.deg, signif = signif)
    class(test) = c("list", "ur.sp")
    
    # Add:
    title = "Schmidt-Phillips Unit Root Test"
    description = date()
    
    # Return Value:
    new("fURTEST", 
        call = CALL,
        data = as.data.frame(x),
        data.name = data.name, 
        test = test, 
        title = as.character(title),
        description = as.character(description)
        )   
}   


# ------------------------------------------------------------------------------


urzaTest =
function(x, model = c("intercept", "trend", "both"), lag = 2)
{   # A function implemented by Diethelm Wuertz
    
    # Description:  
    #   Zivot-Andrews test for unit roots
    
    # Note:
    #   Requires "urca" which is not part of this distribution

    # FUNCTION:
    
    # Settings:
    CALL = match.call()
    data.name = as.character(deparse(substitute(x)))
    
    # Test:
    test = .urza(x, model = model[1], lag = lag)
    class(test) = c("list", "ur.za")
    
    # Add:
    title = "Zivot-Andrews Unit Root Test"
    description = date()
    
    # Return Value:
    new("fURTEST", 
        call = CALL,
        data = as.data.frame(x),
        data.name = data.name, 
        test = test, 
        title = as.character(title),
        description = as.character(description)
        )   
}   


################################################################################
# URCA BUILTIN:


.urers = 
function(y, type = c("DF-GLS", "P-test"), model = c("constant", "trend"), 
lag.max = 4) {
    
    # Description:
    # Elliott, Rothenberg and Stock Unit Root Test
    
    # Author:
    #   Bernhard Pfaff
    #   modified by Diethelm Wuertz
    
    # Source:
    #   Contributed R package "urca" Version 0.5-5, GPL
    
    # FUNCTION:
    
    # Settings:
    type = type[1]
    model = model[1]
    lag.max = as.integer(lag.max)
    if (lag.max < 0) {
        warning("\nlag.max bust be greater or equal to one and integer.")
        warning("\nSetting lag.max = 4")
        lag.max = 4
    }
    
    lag.max = lag.max + 1
    idx = 2:lag.max
    y = na.omit(as.vector(y))
    nobs = length(y)
    
    if (nobs < 50) {
        rowsel = 1
    } else if (nobs < 100) {
        rowsel = 2
    } else if (nobs <= 200) {
        rowsel = 3
    } else if (nobs > 200) {
        rowsel = 4
    }
        
    if (model == "constant") {
        ahat = 1 - 7.0/nobs
        ya = c(y[1], y[2:nobs]-ahat*y[1:(nobs-1)])
        za1 = c(1, rep(1-ahat, nobs-1))
        yd.reg = summary(lm(ya ~ -1 + za1))
        yd = y - coef(yd.reg)[1]
    } else if (model == "trend") {
        ahat = 1 - 13.5/nobs
        ya = c(y[1], y[2:nobs]-ahat*y[1:(nobs-1)])
        za1 = c(1, rep(1-ahat, nobs-1))
        trd = 1:nobs
        za2 = c(1, trd[2:nobs]-ahat*trd[1:(nobs-1)])
        yd.reg = summary(lm(ya ~ -1 + za1 + za2))
        yd = y - coef(yd.reg)[1] - coef(yd.reg)[2]*trd 
    }
    
    what = function(x, z = y) {
        z.l =  z[1:(nobs-1)]
        z.diff = diff(z)
        z.dlags = embed(diff(z), x)[, -1]
        data.what = data.frame(cbind(
            z.diff[-(1:(x-1))], 
            z.l[-(1:(x-1))], z.dlags))
        bic = BIC(lm(data.what))
        return(bic)
    }
    
    if (type == "P-test") {
        cvals.ptest = array(c(
            1.87, 1.95, 1.91, 1.99, 2.97, 3.11, 3.17, 3.26, 3.91, 
            4.17, 4.33, 4.48, 4.22, 4.26, 4.05, 3.96, 5.72, 5.64, 
            5.66, 5.62, 6.77, 6.79, 6.86, 6.89), c(4, 3, 2))
        res = residuals(yd.reg)
        if (model == "constant") {
            null.res = c(0, diff(y))
            cvals = as.matrix(t(cvals.ptest[rowsel, , 1]))
            model = "with intercept"
        } else if (model == "trend") {
            null.res = c(0, diff(y))
            null.res = null.res - mean(null.res)
            cvals = as.matrix(t(cvals.ptest[rowsel, , 2]))
            model = "with intercept and trend"
        }
        sig.null = sum(null.res^2)
        sig.res = sum(res^2)
        if (lag.max > 1) {
            bic = sapply(idx, what, z = y)
            BIC.opt = which.min(bic) + 1
            y.l =  y[1:(nobs-1)]
            y.diff = diff(y)
            y.dlags = embed(diff(y), BIC.opt)[, -1]
            data.what = data.frame(cbind(
                y.diff[-(1:(BIC.opt-1))], y.l[-(1:(BIC.opt-1))], y.dlags))
            what.reg = summary(lm(data.what))
            npar = nrow(what.reg$coef)
            sumlc = sum(what.reg$coef[3:npar,1])
            lag.max = BIC.opt-1
        } else if (lag.max <= 1) {
            y.diff = diff(y)
            y.l = y[1:(nobs-1)]
            what.reg = summary(lm(y.diff ~ y.l))
            sumlc = 0
            lag.max = lag.max-1
        }
        what.sq = what.reg$sigma^2/(1-sumlc)^2
        teststat = (sig.res - ahat*sig.null)/what.sq
        test.reg = NULL
    } else if (type == "DF-GLS") {
        if (model=="constant") {
            cvals = as.matrix(t(c(-2.5658-1.960/nobs-10.04/(nobs**2),
                -1.9393-0.398/nobs,-1.6156-0.181/nobs)))
            model = "with intercept"
        } else if (model == "trend") {
            cvals.dfgls.tau = matrix(-1*c(3.77, 3.58, 3.46, 3.48, 3.19, 
                3.03, 2.93, 2.89, 2.89, 2.74, 2.64, 2.57), nrow = 4, ncol = 3)
            cvals = as.matrix(t(cvals.dfgls.tau[rowsel,]))
            model = "with intercept and trend"
        }
        yd.l =  yd[1:(nobs-1)]
        yd.diff = diff(yd)
        if (lag.max > 1) {
            yd.dlags = embed(diff(yd), lag.max)[, -1]
            data.dfgls = data.frame(cbind(yd.diff[-(1:(lag.max-1))], 
                yd.l[-(1:(lag.max-1))], yd.dlags))
            colnames(data.dfgls) = c("yd.diff", "yd.lag", 
                paste("yd.diff.lag", 1:(lag.max-1), sep=""))
            dfgls.form = formula(paste("yd.diff ~ -1 + ", 
                paste(colnames(data.dfgls)[-1], collapse=" + ")))
        } else if (lag.max <= 1) {
          data.dfgls = data.frame(cbind(yd.diff, yd.l))
          colnames(data.dfgls) = c("yd.diff", "yd.lag")
          dfgls.form = formula("yd.diff ~ -1 + yd.lag")
        }
        dfgls.reg = summary(lm(dfgls.form, data=data.dfgls))
        teststat = coef(dfgls.reg)[1,3]
        test.reg = dfgls.reg
        lag.max = lag.max-1
    }
    
    # Add Names:
    colnames(cvals) = c("1%", "5%", "10%")
    rownames(cvals) = c("critical values")
    
    # Result:
    ans = list(test = "ur.ers", y = y, yd = yd, type = type, 
        model = model, lag = as.integer(lag.max), cval = 
        round(cvals, 2), teststat = teststat, testreg = test.reg, 
        test.name = "Elliot, Rothenberg \& Stock",
        
        statistic = teststat, 
        parameter = c(type, model), 
        method = "Elliot, Rothenberg \& Stock"
        )
        
    # Return Value:
    class(ans) = "htest"
    ans
}


# ------------------------------------------------------------------------------


.urkpss = 
function(y, type = c("mu", "tau"), lags = c("short", "long", "nil"), 
use.lag = NULL) {
    
    # Description:
    #   KPSS Unit Root Test
    
    # Author:
    #   Bernhard Pfaff
    #   modified by Diethelm Wuertz
    
    # Source:
    #   Contributed R package "urca" Version 0.5-5, GPL
    
    # FUNCTION:
    
    y = na.omit(as.vector(y))
    n = length(y)
    type = match.arg(type)
    lags = match.arg(lags)
    
    if (!(is.null(use.lag))) {
        lmax = as.integer(use.lag)
        if (lmax < 0) {
            warning("\nuse.lag has to be positive and integer; lags='short' used.")
            lmax = trunc(4*(n/100)^0.25)
        }
    } else if (lags == "short") {
        lmax = trunc(4*(n/100)^0.25)
    } else if (lags == "long") {
        lmax = trunc(12*(n/100)^0.25)
    } else if (lags == "nil") {
        lmax = 0
    }
    
    if (type == "mu") {
        cval = as.matrix(t(c(0.347, 0.463, 0.574, 0.739)))
        colnames(cval) = c("10%", "5%", "2.5%", "1%")
        rownames(cval) = "critical values"
        res = y - mean(y)
    } else if (type == "tau") {
        cval = as.matrix(t(c(0.119, 0.146, 0.176, 0.216)))
        colnames(cval) = c("10%", "5%", "2.5%", "1%")
        rownames(cval) = "critical values"
        trend = 1:n
        res = residuals(lm(y ~ trend))
    }
    
    S = cumsum(res)
    nominator = sum(S^2)/n^2
    s2 = sum(res^2)/n
    
    if (lmax == 0) {
        denominator = s2
    } else {
        index = 1:lmax
        x.cov = sapply(index, function(x) t(res[-c(1:x)]) %*%
            res[-c((n-x+1):n)])
        bartlett = 1-index/(lmax+1)
        denominator = s2 + 2/n*t(bartlett)%*%x.cov
    }
    
    # Result:
    teststat = nominator/denominator
    ans = list(test = "ur.kpss", y = y, type = type, lag = as.integer(lmax), 
        teststat=as.numeric(teststat), cval = cval, res = res, 
        test.name = "KPSS") 
        
    # Return Value:
    class(ans) = "htest"
    ans
}


# ------------------------------------------------------------------------------


.urpp = 
function(x, type = c("Z-alpha", "Z-tau"), model = c("constant", "trend"), 
lags = c("short", "long"), use.lag = NULL) {
    
    # Description:
    #   Phillips-Perron Unit Root Test
    
    # Author:
    #   Bernhard Pfaff
    #   modified by Diethelm Wuertz
    
    # Source:
    #   Contributed R package "urca" Version 0.5-5, GPL
    
    # FUNCTION:
        
    x = na.omit(as.vector(x))
    n = length(x)
    y = x[-1]
    y.l1 = x[-n]
    n = n-1
    lags = match.arg(lags)
    model = match.arg(model)
    type = match.arg(type)
    if (!(is.null(use.lag))) {
        lmax = as.integer(use.lag)
        if (lmax < 0) {
            warning("\nuse.lag has to be positive and integer")
            warning("\nlags = 'short' used.")
            lmax = trunc(4*(n/100)^0.25)}
    } else if (lags == "short") {
        lmax = trunc(4*(n/100)^0.25)
    } else if (lags == "long") {
        lmax = trunc(12*(n/100)^0.25)
    }
    if (model == "trend") {
        cval = as.matrix(t(c(-3.9638-8.353/n-47.44/(n^2), 
            -3.4126-4.039/n-17.83/(n^2), -3.1279-2.418/n-7.58/(n^2))))
        colnames(cval) = c("1%", "5%", "10%")
        rownames(cval) = "critical values"
        model = "with intercept and trend"
        trend = (1:n) - n/2
        test.reg = summary(lm(y ~ y.l1 + trend))
        res = residuals(test.reg)
        my.tstat = coef(test.reg)[1, 3]
        beta.tstat = coef(test.reg)[3, 3]
        res = residuals(test.reg)
        s = 1/n*(sum(res^2))
        myybar = (1/n^2)*sum((y-mean(y))^2)
        myy = (1/n^2)*sum(y^2)
        mty = (n^(-5/2))*(t(1:n)%*%y)
        my = (n^(-3/2))*sum(y)
        idx = 1:lmax
        coprods = sapply(idx, 
            function(l) t(res[-c(1:l)])%*%(res[-c((n-l+1):n)]))
        weights = 1 - idx/(lmax+1)
        sig = s + (2/n)*(t(weights)%*%coprods)
        lambda = 0.5*(sig-s)
        lambda.prime = lambda/sig
        M = (1-n^(-2))*myy - 12*mty^2 + 12*(1 + 1/n)*mty*my - 
            (4 + 6/n + 2/n^2)*my^2
        my.stat = sqrt(s/sig)*my.tstat - lambda.prime*sqrt(sig)*my /
            (sqrt(M)*sqrt((M+my^2)))
        beta.stat = sqrt(s/sig)*beta.tstat - 
            lambda.prime*sqrt(sig)*(0.5*my - mty)/(sqrt(M/12)*sqrt(myybar))
        aux.stat = as.matrix(c(round(my.stat, 4), round(beta.stat, 4)))
        rownames(aux.stat) = c("Z-tau-mu", "Z-tau-beta")
        colnames(aux.stat) = "aux. Z statistics"
        if (type=="Z-tau") {
              tstat = (coef(test.reg)[2, 1] - 1)/coef(test.reg)[2, 2]
              teststat = sqrt(s/sig)*tstat-lambda.prime*sqrt(sig)/sqrt(M)
        } else if (type=="Z-alpha") {
              alpha = coef(test.reg)[2, 1]
              teststat = n*(alpha-1)-lambda/M
        }
    } else if (model=="constant") {
        cval = as.matrix(t(c(-3.4335-5.999/n-29.25/(n^2), 
            -2.8621-2.738/n-8.36/(n^2), -2.5671-1.438/n-4.48/(n^2))))
        colnames(cval) = c("1%", "5%", "10%")
        rownames(cval) = "critical values"
        model = "with intercept"
        test.reg = summary(lm(y ~ y.l1))
        my.tstat = coef(test.reg)[1, 3]
        res = residuals(test.reg)
        s = 1/n*(sum(res^2))
        myybar = (1/n^2)*sum((y-mean(y))^2)
        myy = (1/n^2)*sum(y^2)
        my = (n^(-3/2))*sum(y)
        idx = 1:lmax
        coprods = sapply(idx, function(l) t(res[-c(1:l)]) %*%
            (res[-c((n-l+1):n)]))
        weights = 1 - idx/(lmax+1)
        sig = s + (2/n)*(t(weights)%*%coprods)
        lambda = 0.5*(sig-s)
        lambda.prime = lambda/sig
        my.stat = sqrt(s/sig)*my.tstat + lambda.prime*sqrt(sig)*my /
            (sqrt(myy)*sqrt(myybar))
        aux.stat = as.matrix(round(my.stat, 4))
        rownames(aux.stat) = "Z-tau-mu"
        colnames(aux.stat) = "aux. Z statistics"
        if (type=="Z-tau") {
              tstat = (coef(test.reg)[2, 1] - 1)/coef(test.reg)[2, 2]
              teststat = sqrt(s/sig)*tstat-lambda.prime*sqrt(sig)/sqrt(myybar)
        } else if (type=="Z-alpha") {
              alpha = coef(test.reg)[2, 1]
              teststat = n*(alpha-1)-lambda/myybar
        }
    }
    
    # Result:
    ans = list(test = "ur.pp", y = y, type = type, model = model, 
        lag = as.integer(lmax), cval = cval, teststat = as.numeric(teststat), 
        testreg = test.reg, auxstat = aux.stat, res = res, 
        test.name = "Phillips-Perron")
        
    # Return Value:
    class(ans) = "htest"
    ans
}


# ------------------------------------------------------------------------------


.ursp = 
function(y, type = c("tau", "rho"), pol.deg = c(1, 2, 3, 4), 
signif = c(0.01, 0.05, 0.1)) {
    
    # Description:
    #   Schmidt-Phillips Unit Root Test
    
    # Author:
    #   Bernhard Pfaff
    #   modified by Diethelm Wuertz
    
    # Source:
    #   Contributed R package "urca" Version 0.5-5, GPL
    
    # FUNCTION:
    
    y = na.omit(as.vector(y))
    type = match.arg(type)
    signif.val = c(0.01, 0.05, 0.1)
    
    if (!(signif %in% signif.val)) {
        warning("\nPlease, provide as signif one of c(0.01, 0.05, 0.1)")
        warning("\nsignif = 0.01 used")
        signif = 0.01
    }
    
    if (!(pol.deg %in% c(1:4))) {
        warning("\nPlease, provide as polynomial degree one of c(1, 2, 3, 4)")
        warning("\npol.deg = 1 used")
        pol.deg = 1
    }
    
    n = length(y)
    lag = trunc(12*(n/100)^0.25)
    idx = 1:lag
    trend1 = 1:n
    y.diff = diff(y)
    
    if (pol.deg == 1) {
        yi.hat = (y[n]-y[1])/(n-1)
        phi.y = y[1]-yi.hat
        S.hat = y - phi.y - yi.hat*trend1
        S.hat.l1 = S.hat[-n]
        test.reg = summary(lm(y.diff ~ 1 + S.hat.l1))
        sp.data = data.frame(cbind(y[2:n],  y[1:(n-1)], trend1[2:n]))
        colnames(sp.data) = c("y", "y.lagged", "trend.exp1")
        corr.reg = summary(lm(sp.data))
        res = residuals(corr.reg)
        sig.eps = (1/n)*sum(res^2)
        coprods = sapply(idx, 
            function(x) t(res[-c(1:x)])%*%(res[-c((n-x):(n-1))]))
        weights = (2*(lag-idx)/lag)^2
        sig = sig.eps + (2/n)*(t(weights)%*%coprods)
        omega2.hat = sig.eps/sig
    } else if (pol.deg == 2) {
        trend2 = trend1^2
        S.hat = c(0, cumsum(residuals(
            summary(lm(y.diff ~ trend1[2:n])))))
        test.reg = summary(lm(y.diff ~ S.hat[-n] + trend1[2:n]))
        sp.data = data.frame(cbind(y[2:n],  y[1:(n-1)], 
            trend1[2:n], trend2[2:n]))
        colnames(sp.data) = c("y", "y.lagged", "trend.exp1", "trend.exp2")
        corr.reg = summary(lm(sp.data))
        res = residuals(corr.reg)
        sig.eps = (1/n)*sum(res^2)
        coprods = sapply(idx, 
            function(x) t(res[-c(1:x)])%*%(res[-c((n-x):(n-1))]))
        weights = (2*(lag-idx)/lag)^2
        sig = sig.eps + (2/n)*(t(weights)%*%coprods)
        omega2.hat = sig.eps/sig
    } else if (pol.deg==3) {
        trend2 = trend1^2
        trend3 = trend1^3    
        S.hat = c(0, cumsum(residuals(
            summary(lm(y.diff ~ trend1[2:n] + trend2[2:n])))))
        test.reg = summary(lm(y.diff ~ S.hat[-n] + trend1[2:n] + trend2[2:n]))
        sp.data = data.frame(cbind(y[2:n],  y[1:(n-1)], trend1[2:n], 
            trend2[2:n], trend3[2:n]))
        colnames(sp.data) = c("y", "y.lagged", "trend.exp1", "trend.exp2", 
        	"trend.exp3")
        corr.reg = summary(lm(sp.data))
        res = residuals(corr.reg)
        sig.eps = (1/n)*sum(res^2)
        coprods = sapply(idx, 
            function(x) t(res[-c(1:x)])%*%(res[-c((n-x):(n-1))]))
        weights = (2*(lag-idx)/lag)^2
        sig = sig.eps + (2/n)*(t(weights)%*%coprods)
        omega2.hat = sig.eps/sig
    } else if (pol.deg==4) {
        trend2 = trend1^2
        trend3 = trend1^3
        trend4 = trend1^4
        S.hat = c(0, cumsum(residuals(
            summary(lm(y.diff ~ trend1[2:n] + trend2[2:n] + trend3[2:n])))))
        test.reg = summary(lm(y.diff ~ S.hat[-n] + trend1[2:n] + 
            trend2[2:n] + trend3[2:n]))
        sp.data = data.frame(cbind(y[2:n],  y[1:(n-1)], trend1[2:n], 
            trend2[2:n], trend3[2:n], trend4[2:n]))
        colnames(sp.data) = c("y", "y.lagged", "trend.exp1", "trend.exp2", 
            "trend.exp3", "trend.exp4")
        corr.reg = summary(lm(sp.data))
        res = residuals(corr.reg)
        sig.eps = (1/n)*sum(res^2)
        coprods = sapply(idx, function(x) t(res[-c(1:x)]) %*% 
            (res[-c((n-x):(n-1))]))
        weights = (2*(lag-idx)/lag)^2
        sig = sig.eps + (2/n)*(t(weights)%*%coprods)
        omega2.hat = sig.eps/sig
    }
  
    if (type == "rho") {
        rho = n*coef(test.reg)[2,1]
        teststat = rho/omega2.hat
        cval = .Spcv(obs = n, type = "rho", pol.deg = pol.deg, 
            signif = signif)
    } else if (type == "tau") {
        tau = coef(test.reg)[2,3]
        teststat = tau/sqrt(omega2.hat)
        cval = .Spcv(obs = n, type = "tau", pol.deg = pol.deg, 
            signif = signif)
    }
  
    # Result:
    ans = list(test = "ur.sp", y = y, type = type, polynomial = 
        as.integer(pol.deg), teststat = as.numeric(teststat), 
        cval = cval, signif = signif, res = res, testreg = corr.reg, 
        test.name = "Schmidt-Phillips")
        
    # Return Value:
    class(ans) = "htest"
    ans
}


# ------------------------------------------------------------------------------


.Spcv = 
function(obs, type = c("tau", "rho"), pol.deg = c(1, 2, 3, 4), 
signif=c(0.01, 0.025, 0.05, 0.1)) 
{
    
    # Function for critical values of ur.sp
    
    obs = as.integer(obs)
    obs.ranges = c(25, 50, 100, 200, 500, 1000, 1.0e30)
    type = match.arg(type)
    dim.1 = which(obs.ranges >= obs, arr.ind = TRUE)[1]
    signif.val = c(0.01, 0.05, 0.1)
    if (!(signif %in% signif.val)) {
        warning("\nPlease, provide as signif one of c(0.01, 0.05, 0.1)")
        warning("\nsignif = 0.01 used")
        signif = 0.01
    }
    dim.2 = which(signif==signif.val, arr.ind=TRUE)
    if (!(pol.deg %in% c(1:4))) {
        warning("\nPlease, provide as polynomial degree one of c(1, 2, 3, 4)")
        warning("\npol.deg = 1 used")
        pol.deg = 1
    }
    dim.3 = pol.deg
    if (type == "tau") {
        cvs.tau = -1*c(3.9, 3.73, 3.63, 3.61, 3.59, 3.58, 3.56, 3.18, 3.11, 
            3.06, 3.04, 3.04, 3.02, 3.02, 2.85, 2.8, 2.77, 2.76, 2.76, 2.75, 
            2.75, 4.52, 4.28, 4.16, 4.12, 4.08, 4.06, 4.06, 3.78, 3.77, 
            3.65, 3.6, 3.55, 3.55, 3.53, 3.52, 3.41, 3.34, 3.31, 3.28, 
            3.26, 3.26, 3.26, 5.07, 4.73, 4.59, 4.53, 4.5, 4.49, 4.44, 
            4.26, 4.08, 4.03, 3.99, 3.96, 3.95, 3.93, 3.89, 3.77, 3.72, 
            3.69, 3.68, 3.68, 3.67, 5.57, 5.13, 4.99, 4.9, 4.85, 4.83, 
            4.81, 4.7, 4.47, 4.39, 4.33, 4.31, 4.31, 4.29, 4.3, 4.15, 4.1, 
            4.06, 4.03, 4.03, 4.01)
        cv.array = array(cvs.tau, dim = c(7, 3, 4), dimnames = 
            c("obs", "signif", "pol.deg"))
        cval = cv.array[dim.1, dim.2, dim.3]
    } else if (type == "rho") {
        cvs.rho = -1*c(20.4, 22.8, 23.8, 24.8, 25.3, 25.3, 25.2, 15.7, 
            17.0, 17.5, 17.9, 18.1, 18.1, 18.1, 13.4, 14.3, 14.6, 14.9, 
            15.0, 15.0, 15.0, 24.6, 28.4, 30.4, 31.8, 32.4, 32.5, 32.6, 
            20.1, 22.4, 23.7, 24.2, 24.8, 24.6, 24.7, 17.8, 19.5, 20.4, 
            20.7, 21.0, 21.1, 21.1, 28.1, 33.1, 36.3, 38.0, 39.1, 39.5, 
            39.7, 23.8, 27.0, 29.1, 30.1, 30.6, 30.8, 30.6, 21.5, 24.0, 
            25.4, 26.1, 26.6, 26.7, 26.7, 31.0, 37.4, 41.8, 44.0, 45.3, 
            45.7, 45.8, 26.9, 31.2, 34.0, 35.2, 36.2, 36.6, 36.4, 24.7,
            28.1, 30.2, 31.2, 31.8, 32.0, 31.9)
        cv.array = array(cvs.rho, dim = c(7, 3, 4), dimnames = 
            c("obs", "signif", "pol.deg"))
        cval = cv.array[dim.1, dim.2, dim.3]
    }
    
    # Return Value:
    return(cval)
}


# ------------------------------------------------------------------------------


.urza = 
function(y, model=c("intercept", "trend", "both"), lag) {
    
    # Description:
    #   Zivot-Andrews Unit Root Test
    
    # Author:
    #   Bernhard Pfaff
    #   modified by Diethelm Wuertz
    
    # Source:
    #   Contributed R package "urca" Version 0.5-5, GPL
    
    # FUNCTION:
    
    y = na.omit(as.vector(y))
    n = length(y)
    model = match.arg(model)
    lag = as.integer(lag)
    
    if (length(lag) > 1 || lag < 1) {
        warning("\nPlease, specify max. number of lags for differenced series")
        warning("\nas positive integer; lags = 1 is now used.")
        lags = 1
    }
        
    datmat = matrix(NA, n, lag + 3)
    if (n < ncol(datmat) + 2) {
        stop("\nInsufficient number of obeservations.")
    }
    
    idx = 1:(n-1)
    trend = seq(1, n)
    datmat[,1] = y
    datmat[,2] = c(NA, y)[1:n]
    datmat[,3] = trend
    
    for (i in 1:lag) {
        datmat[ , i + 3] = c(rep(NA, i + 1), diff(y))[1:n]}
        datmat = as.data.frame(datmat)
        colnames(datmat) = c("y", "y.l1", "trend", paste("y.dl", 1:lag, 
            sep = ""))
        if (model=="intercept") {
            roll = 
            function(z) {
                du = c(rep(0, z), rep(1, (n-z)))
                rollmat = cbind(datmat, du)
                roll.reg = coef(summary(lm(rollmat)))
                (roll.reg[2,1]-1.0)/roll.reg[2,2]
            }
            roll.stat = sapply(idx, roll)
            cval = c(-5.34, -4.8, -4.58)
            bpoint = which.min(roll.stat)
            du = c(rep(0, bpoint), rep(1, (n-bpoint)))
            testmat = cbind(datmat, du)
            test.reg = summary(lm(testmat)) 
        } else if (model=="trend") {
            roll = 
            function(z) {
                dt = c(rep(0, z), 1:(n-z))
                rollmat = cbind(datmat, dt)
                roll.reg = coef(summary(lm(rollmat)))
                (roll.reg[2,1]-1.0)/roll.reg[2,2]
            }
            roll.stat = sapply(idx, roll)
            cval = c(-4.93, -4.42, -4.11)
            bpoint = which.min(roll.stat)
            dt = c(rep(0, bpoint), 1:(n-bpoint))
            testmat = cbind(datmat, dt)
            test.reg = summary(lm(testmat)) 
        } else if (model=="both") {
            test.reg = summary(lm(datmat))
            roll = 
            function(z) {
                du = c(rep(0, z), rep(1, (n-z)))
                dt = c(rep(0, z), 1:(n-z))
                rollmat = cbind(datmat, du, dt)
                roll.reg = coef(summary(lm(rollmat)))
                (roll.reg[2,1]-1.0)/roll.reg[2,2]
            }
            roll.stat = sapply(idx, roll)
            cval = c(-5.57, -5.08, -4.82)
            bpoint = which.min(roll.stat)
            du = c(rep(0, bpoint), rep(1, (n-bpoint)))
            dt = c(rep(0, bpoint), 1:(n-bpoint))
            testmat = cbind(datmat, du, dt)
            test.reg = summary(lm(testmat)) 
        }
    
    # Result:
    teststat = roll.stat[bpoint]
    ans = list(test = "ur.za", y = y, model = model, lag = lag, 
        teststat = teststat, cval = cval, bpoint = bpoint, 
        tstats = roll.stat, res = test.reg$residuals, testreg = test.reg, 
        test.name = "Zivot-Andrews")
        
    # Return Value:
    class(ans) = "htest"
    ans
}


# ------------------------------------------------------------------------------


.plot.urers = 
function(x) {
  
    if (is.null(x@testreg)) {
        stop("No plot method for P-test available")
    }
    oldpar = par(no.readonly = TRUE)
    on.exit(par(oldpar))
    par(mfrow = c(1, 1))
    layout(matrix(c(1, 2, 3, 1, 2, 4), 3 , 2))
    suppressWarnings(plot.ts(diff(x@yd)[-c(1:x@lag)], 
        main = "Diagram of fit for test regression", 
        sub = paste("detrending ", x@model, " and ", 
        x@lag, " lagged differences used in test regression",  
        sep=""), ylab = "Actual and fitted values", xlab = ""))
    lines(diff(x@yd)[-c(1:x@lag)] - resid(x@testreg), col = "seagreen")
    plot.ts(resid(x@testreg), main = "Residuals", ylab = "", xlab = "")
    abline(h = 0, col = "red")
    acf(resid(x@testreg), main = "Autocorrelations of Residuals")
    pacf(resid(x@testreg), main = "Partial Autocorrelations of Residuals")
}


# ------------------------------------------------------------------------------


.plot.urkpss = 
function(x) {
    
    oldpar = par(no.readonly = TRUE)
    on.exit(par(oldpar))
    par(mfrow = c(1,1))
    layout(matrix(c(1, 2, 1, 3), 2 , 2))
    plot.ts(x@res, main = paste("Residuals from test regression of type:", 
        x@type, " with", x@lag, "lags", sep = " "), ylab = "residuals", 
        xlab = "")
    abline(h = 0, col = "red")
    acf(x@res, main = "Autocorrelations of Residuals")
    pacf(x@res, main = "Partial Autocorrelations of Residuals")
}


# ------------------------------------------------------------------------------


.plot.urpp = 
function(x) {
  
    oldpar = par(no.readonly = TRUE)
    on.exit(par(oldpar))
    par(mfrow=c(1,1))
    layout(matrix(c(1, 2, 3, 1, 2, 4), 3 , 2))
    plot.ts(x@y[-1], main =
        paste("Diagram of fit for model", x@model, sep = " "), 
        ylab="Actual and fitted values", xlab = "")
    lines(x@y - x@res, col="seagreen")
    plot.ts(x@res, main = "Residuals", ylab = "", xlab = "")
    abline(h = 0, col = "red")
    acf(x@res, main = "Autocorrelations of Residuals")
    pacf(x@res, main = "Partial Autocorrelations of Residuals")
}


# ------------------------------------------------------------------------------


.plot.ursp =  
function(x) {
    
    oldpar = par(no.readonly = TRUE)
    on.exit(par(oldpar))
    par(mfrow=c(1,1))
    layout(matrix(c(1, 2, 3, 1, 2, 4), 3 , 2))
    plot.ts(x@y[-1], main = 
        paste("Diagram of fit for model with polynomial degree of ", 
        x@polynomial, sep="") , ylab = "Actual and fitted values", 
        xlab = "")
    lines(x@y[-1] - x@res, col="seagreen")
    plot.ts(x@res, main = "Residuals", ylab = "", xlab = "")
    abline(h = 0, col = "red")
    acf(x@res, main = "Autocorrelations of Residuals")
    pacf(x@res, main = "Partial Autocorrelations of Residuals")
}


# ------------------------------------------------------------------------------


.plot.urza =  
function(x) {
    
    oldpar = par(no.readonly = TRUE)
    on.exit(par(oldpar))
    par(mfrow=c(1,1))
    yvals = sort(c(x@cval, x@tstats))
    n = length(x@y)
    xvals = pretty(1:n)
    plot.ts(x@tstats, main = "Zivot and Andrews Unit Root Test", 
        ylab = "t-statistics for lagged endogenous variable", 
        ylim = c(min(yvals), max(yvals)))
    abline(h = x@cval, col = c("red", "blue", "seagreen"))
    if (x@teststat < x@cval[3]) {
        abline(v = x@bpoint, col = "red", lty = 2)
    }
    mtext(paste("Model type:", x@model, sep = " "), side = 1, line = 4)
    legend(x = n, y = max(yvals), c("1% c.v.", "2.5% c.v.", "5% c.v."), 
        col = c("red", "blue", "seagreen"), xjust = 1, yjust = 1, 
        lty = 1, horiz = TRUE, cex = 0.66, bty = "n")
}
  

################################################################################

