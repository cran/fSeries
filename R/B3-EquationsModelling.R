
# This library is free software; you can redistribute it and/or
# Modify it under the terms of the GNU Library General Public
# License as published by the Free Software Foundation; either
# Version 2 of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Library General Public License for more details.
#
# You should have received a copy of the GNU Library General
# Public License along with this library; if not, write to the
# Free Foundation, Inc., 59 Temple Place, Suite 330, Boston,
# MA  02111-1307  USA

# Copyrights (C)
# For this R-port: 
#   1999 - 2004, Diethelm Wuertz, GPL
#   Diethelm Wuertz <wuertz@itp.phys.ethz.ch>
#   info@rmetrics.org
#   www.rmetrics.org
# For the code accessed (or partly included) from other R-ports:
#   see R's copyright and license files
# For the code accessed (or partly included) from contributed R-ports
# And other sources
#   see Rmetrics's copyright file


################################################################################
# FUNCTION:             REGRESSION MODELLING:
#  eqnsFit               Wrapper Function for "systemfit" and "sem" Models:
#  * OLS                  Ordinary Least Squares
#  * WLS                  Weighted Least Squares
#  * SUR                  Seemingly Unrelated Regressions
#  * 2SLS                 Two-Stage Least Squares
#  * W2SLS                Weighted Two-Stage Least Squares
#  * 3SLS                 Three-Stage Least Squares
#  * W3SLS                Weighted Three-Stage Least Squares
# S3-METHODS:           DESCRIPTION:
#  print.fEQNS           S3: Print method for an object of class fEQNS  
#  plot.fEQNS            S3: Plot method for an object of class fEQNS
#  summary.fEQNS         S3: Summary method for an object of class fEQNS
# S3-METHODS:           DESCRIPTION:
#  coef.fEQNS            S3: Method for coefficients
#  fitted.fEQNS          S3: Method for fitted values
#  residuals.fEQNS       S3: Method for residual values
#  vcov.fEQNS            S3: Method for variance-covariance Matrix
# S-PLUS LIKE:          WRAPPER:
#  SUR                   SUR Wrapper
# S3-METHODS:           DESCRIPTION:
#  predict.fEQNS         S3: Pfrediction method for an object of class fEQNS
################################################################################
# REQUIRED PACLKAGE - PACKAGE DESCRIPTION:
# Package: systemfit
# Version: 0.7-2
# Date: 2004/08/19
# Title: Simultaneous Equation Estimation Package
# Author: Jeff D. Hamann <jeff.hamann@forestinformatics.com> and
#   Arne Henningsen <ahenningsen@agric-econ.uni-kiel.de>
# Maintainer: Jeff D. Hamann <jeff.hamann@forestinformatics.com>
# Depends: R (>= 1.8.0)
# Description: This package contains functions for fitting simultaneous
#   systems of linear and nonlinear equations using Ordinary Least
#   Squares (OLS), Weighted Least Squares (WLS), Seemingly Unrelated
#   Regressions (SUR), Two-Stage Least Squares (2SLS), Weighted
#   Two-Stage Least Squares (W2SLS), Three-Stage Least Squares (3SLS),
#   and Weighted Three-Stage Least Squares (W3SLS).
# License: GPL version 2 or newer
# URL: http://www.r-project.org, http://www.forestinformatics.com, 
#   http://www.arne-henningsen.de
################################################################################


require(methods)


# ------------------------------------------------------------------------------


setClass("fEQNS", 
    representation(
        call = "call",
        formulas = "list",
        data = "data.frame",
        method = "character",
        fit = "list",   
        title = "character",
        description = "character"))
        
      
# ------------------------------------------------------------------------------


eqnsFit = 
function(formulas, data = list(), 
method = c("OLS", "WLS", "SUR", "2SLS", "W2SLS", "3SLS", "W3SLS"), 
title = NULL, description = NULL, ...)
{   # A function implemented by Diethelm Wuertz

    # Description:
    #   Common function call for several system equation models.
    
    # Arguments:
    #   formulas - the list of formulas describing the system of 
    #       equations
    #   data - the input data set in form of a 'data.frame' or 
    #       'timeSeries' object
    #   method - a character string describing the desired method, 
    #       one of: "OLS", "WLS", "SUR", "2SLS", "W2SLS", "3SLS", 
    #       or "W3SLS".
    #   title - a character string which allows for a project title
    #   description - a character string which allows for a brief 
    #       description
    #   ... - additional optional arguments to be passed to the 
    #       underlying function 'systemfit' 
    
    # Value:
    # The function 'eqnaFit' returns an object of class "fEQNS"
    # with the following slots:
    #   @call - the matched function call
    #   @data - the input data in form of a 'data.frame' or a 
    #       'timeSeries' object
    #   @description - a character string which allows for a brief 
    #       project description
    #   @method - the character string describing the desired method
    #   @formulas - the list of formulas describing the system of 
    #       equations
    #   @title - a character string which allows for a project title
    #   @fit - a summary of the  results as a list returned from the 
    #       underlying functions from the 'systemfit' package.  

    # FUNCTION:
    
    # Fit:
    fit = .systemfit(method = method[1], eqns = formulas, 
        data = as.data.frame(data), ...)    
    class(fit) = c("list", "systemfit")
    
    # Internal Function:
    .confint.systemfit = 
    function(object, parm = NULL, level = 0.95, ...)  {
        a = (1 - level) / 2
        a = c(a, 1 - a)
        pct = paste(round(100 * a, 1), "%")
        ci = array(NA, dim = c(length(object$b), 2),
                dimnames = list(names(object$b), pct))
        j = 1
        for (i in 1:object$g) {
            object$eq[[i]]$dfSys = object$df
            object$eq[[i]]$probdfsys = object$probdfsys
            ci[j:(j+object$eq[[i]]$k-1),] = 
                .confint.systemfit.equation(object$eq[[i]])
            j = j + object$eq[[i]]$k
        }
        ci }

    # Internal Function:
    .confint.systemfit.equation = 
    function(object, parm = NULL, level = 0.95, ...) {
        a = (1 - level) / 2
        a = c(a, 1 - a)
        pct = paste(round(100 * a, 1), "%")
        ci = array(NA, dim = c(length(object$b), 2),
                dimnames = list(names(object$b), pct))
        if (object$ probdfsys) {
            fac = qt(a, object$dfSys)
        } else {
            fac = qt(a, object$df)
        }
        ci[] = object$b + object$se %o% fac
        ci 
    }
    
    # Title:
    if (is.null(title)) title = paste(method[1], "Estimation")
    
    # Description:
    if (is.null(description)) description = as.character(date())
    
    # Result:
    ans = new("fEQNS",     
        call = as.call(match.call()),
        formulas = formulas, 
        data = as.data.frame(data),
        method = as.character(method[1]), 
        fit = fit,
        title = as.character(title), 
        description = as.character(description))
        
    # Add:
    ans@fit$coef = coef(ans)
    ans@fit$fitted = fitted(ans)
    ans@fit$residuals = residuals(ans)
    ans@fit$vcov = vcov(ans)
    
    # Return Value:
    ans
}

    
# ------------------------------------------------------------------------------    


print.fEQNS = 
function(x, ...)
{   # A function implemented by Diethelm Wuertz

    # Description:
    #   S3 Print Method for an object of class "fEQNS"
    
    # Arguments:
    #   x - an object of class "fEQNS"
    
    # FUNCTION:
    
    # Title:
    cat("\nTitle:\n")
    cat(as.character(x@title), "\n")
    
    # Formulas:   
    cat("\nFormula List:\n")
    for (i in 1: length(x@formulas)) 
        cat(as.character(x@formulas[i]), "\n") 
        
    # Internal Function:
    .print.systemfit = 
    function(x, digits = 6, ...) {
        # Settings:
        save.digits = unlist(options(digits = digits))
        on.exit(options(digits = save.digits))
        table = NULL
        labels = NULL
        # Print:
        cat("\n")
        cat("Systemfit Results \n")
        cat("Method: ")
        if (!is.null(x$iter)) if (x$iter>1) cat("iterated ")
        cat(paste(x$method, "\n\n"))
        if (!is.null(x$iter)) {
            if (x$iter>1) {
                if (x$iter<x$maxiter) {
                    cat(paste("convergence achieved after",
                        x$iter, "iterations\n\n"))
                } else {
                    cat(paste("warning: convergence not achieved after", 
                        x$iter, "iterations\n\n")) } } }
        for (i in 1:x$g) {
            row = NULL
            row = cbind(round(x$eq[[i]]$n, digits), round(x$eq[[i]]$df, digits),
                round(x$eq[[i]]$ssr, digits), round(x$eq[[i]]$mse, digits),
                round(x$eq[[i]]$rmse, digits), round(x$eq[[i]]$r2, digits),
                round(x$eq[[i]]$adjr2, digits))
            table = rbind(table, row)
            labels = rbind(labels, x$eq[[i]]$eqnlabel) }
        rownames(table) = c(labels)
        colnames(table) = c("N","DF", "SSR", "MSE", "RMSE", "R2", "Adj R2")
        print.matrix(table, quote = FALSE, right = TRUE)
        cat("\n")
        # The Covariance Matrix of the Residuals Used for Estimation
        if (!is.null(x$rcovest)) {
            cat("The Covariance Matrix of the Residuals Used for Estimation\n")
            rcov = x$rcovest
            rownames(rcov) = labels
            colnames(rcov) = labels
            print(rcov)
            cat("\n")
            if (min(eigen(x$rcov, only.values = TRUE)$values) < 0) {
                cat("warning: ")
                cat("Covariance Matrix is NOT positive semidefinit!\n\n") } }
        # The Covariance Matrix of the Residuals:
        cat("The Covariance Matrix of the Residuals\n")
        rcov = x$rcov
        rownames(rcov) = labels
        colnames(rcov) = labels
        print(rcov)
        cat("\n")  
        # The Correlations of the Residuals:
        cat("The Correlations of the Residuals\n")
        rcor = x$rcor
        rownames(rcor) = labels
        colnames(rcor) = labels
        print(rcor)
        cat("\n")
        # The Determinant of the Residual Covariance Matrix:
        cat("The Determinant of the Residual Covariance Matrix: ", 
            x$drcov, "\n")  
        # OLS R-Squared Value of the System:
        cat("OLS R-squared Value of the System: ", x$olsr2, "\n")
        # McElroy's R-Squared Value for the System:
        if (!is.null(x$mcelr2)) 
            cat("McElroy's R-Squared Value for the System: ", x$mcelr2, "\n")  
        # Now print the individual equations: x$eq[[i]]
        for (i in 1:x$g) 
            .print.systemfit.equation(x$eq[[i]], digits)  
        # Return value:
        invisible(x) 
    }
    
    # Internal Function:
    .print.systemfit.equation = 
    function(x, digits = 6, ...) {
        # Print Equation:
        cat("\n")
        cat(x$method, "Estimates for", x$eqnlabel, " (Equation", x$i, ")")
        cat("\nModel Formula: ")
        print(x$formula)
        if (!is.null(x$inst)) {
            cat("Instruments: ")
            print(x$inst)
        }
        cat("\n")
        Signif = symnum(x$p, corr = FALSE, na = FALSE, cutpoints = 
            c(0, 0.001, 0.01, 0.05, 0.1, 1),
            symbols = c("***", "**", "*", "." ," "))
         table = cbind(round(x$b, digits), round(x$se, digits),
            round(x$t, digits), round(x$p,  digits), Signif)
        rownames(table) = names(x$b)
        colnames(table) = c("Estimate","Std. Error","t value","Pr(>|t|)", "")
        print.matrix(table, quote = FALSE, right = TRUE)
        cat("---\nSignif. codes: ", attr(Signif,"legend"), "\n")
        cat(paste("\nResidual Standard Error:", round(x$s, digits),
            "on", x$df, "degrees of freedom\n"))
        cat(paste("Number of observations:", round(x$n, digits),
            "Degrees of Freedom:", round(x$df, digits), "\n"))
        cat(paste("SSR:", round(x$ssr, digits), "MSE:", round(x$mse, 
            digits), "Root MSE:", round(x$rmse,  digits), "\n"))
        cat(paste("Multiple R-Squared:", round(x$r2, digits),
            "Adjusted R-Squared:", round(x$adjr2, digits), "\n\n")) 
        invisible()
    }       
       
    # Printing:
    .print.systemfit(x@fit, ...)

    # Footer:
    if (x@description != "") {
        cat("Description:\n")
        cat(as.character(x@description), "\n") }
    
    # Return Value:
    invisible(x)
}


# ------------------------------------------------------------------------------


plot.fEQNS = 
function(x, ...)
{   # A function implemented by Diethelm Wuertz

    invisible(x)
}


# ------------------------------------------------------------------------------


summary.fEQNS = 
function(object, ...)
{   # A function implemented by Diethelm Wuertz

    # Description:
    #   S3 Summary Method for an object of class "fEQNS"
    
    # Arguments:
    #   object - an object of class "fEQNS"
    
    # FUNCTION:
    
    # Print:
    print(object, ...)
    
    # Return Value:
    invisible(object)  
}


# ******************************************************************************


coef.fEQNS = 
function(object, ...) 
{
    # Description:
    #   Returns the fitted values
    
    # FUNCTION:
    
    # From Slot:
    coef = object@fit$b
    coef = as.data.frame(coef)
   
    # Return Value:
    coef
}


# ------------------------------------------------------------------------------


fitted.fEQNS = 
function(object, ...) 
{
    # Description:
    #   Returns the fitted values
    
    # FUNCTION:
    
    # From Slot:
    x = object@fit
    
    # Fitted values:
    fitted = array(NA, c(length(x$eq[[1]]$fitted), x$g))
    colnames(fitted) = as.character(1:ncol(fitted))
    for (i in 1:x$g)  {
        fitted[, i] = x$eq[[i]]$fitted
        colnames(fitted)[i] = paste("eq", as.character(i), sep = "")
    }
    
    # Return Value:
    fitted
}


# ------------------------------------------------------------------------------


residuals.fEQNS = 
function(object, ...) 
{
    # Description:
    #   Returns all residuals
    
    # FUNCTION:
    
    # From Slot:
    x = object@fit
    
    # Residuals:
    residuals = data.frame(eq1 = x$eq[[1]]$residuals)
    if (x$g > 1) {
        for (i in 2:x$g) {
            residuals = cbind(residuals, new = x$eq[[i]]$residuals)
                names(residuals)[i] = paste("eq", as.character(i), sep = "")
        }
    }
    
    # Return Value:
    residuals
}


# ------------------------------------------------------------------------------


vcov.fEQNS = 
function(object, ...) 
{
    # Description:
    #   Returns the variance-covariance matrix of the coefficients
    
    # FUNCTION:
    
    # Return Value:
    object@fit$bcov
}


# ******************************************************************************


SUR = 
function(formulas, data = list(), ...)
{
    # Fit:
    ans = eqnsFit(formulas = formulas, data = data, method = "SUR", ...)
        
    # Return Value:
    ans
}


# ------------------------------------------------------------------------------


predict.fEQNS = 
function (object, newdata = object@data, se.fit = FALSE, se.pred = FALSE, 
interval = "none", ci = 0.95, ...) 
{   
    
    # Predict:
    ans = .predict.systemfit(object = object@fit, data = newdata,
        se.fit = se.fit, se.pred = se.pred, interval = interval, 
        level = ci, ...)
    
    # Return Value:
    ans
}


################################################################################
# BUILTIN


# This library is free software; you can redistribute it and/or
# Modify it under the terms of the GNU Library General Public
# License as published by the Free Software Foundation; either
# Version 2 of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Library General Public License for more details.
#
# You should have received a copy of the GNU Library General
# Public License along with this library; if not, write to the
# Free Foundation, Inc., 59 Temple Place, Suite 330, Boston,
# MA  02111-1307  USA


# PART I: linear modelling
# PART II: nonlinear modelling


# ------------------------------------------------------------------------------


# $Id: systemfit.R,v 1.24 2004/04/17 18:35:50 henningsena Exp $
# Simultaneous Equation Estimation for R
# Copyright 2002-2004 
#   Jeff D. Hamann <jeff.hamann@forestinformatics.com>
#   Arne Henningsen <http://www.arne-henningsen.de>
#
# This file is part of the nlsystemfit library for R and related languages.
# It is made available under the terms of the GNU General Public License, 
# Version 2, or at your option, any later version, incorporated herein by 
# Reference.
#
# This program is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or 
# FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for 
# More details.
#
# You should have received a copy of the GNU General Public License along 
# with this program; if not, write to the Free Software Foundation, Inc., 
# 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA


# ------------------------------------------------------------------------------


.systemfit = function(method, eqns, 
eqnlabels = c(as.character(1:length(eqns))), inst = NULL, data = list(), 
R.restr = NULL, q.restr = matrix(0,max(nrow(R.restr),0),1), TX = NULL, 
maxiter = 1, tol = 1.0e-5, rcovformula = 1, formula3sls = "GLS", 
probdfsys = !(is.null(R.restr) & is.null(TX)),
single.eq.sigma=(is.null(R.restr) & is.null(TX)), 
solvetol = .Machine$double.eps, saveMemory = FALSE)
{

    # FUNCTION:
    
    # Some Tests:
    if (!(method == "OLS" | method == "WLS" | method == "SUR" | 
        method == "2SLS" | method == "W2SLS" | method == "3SLS")) {
        stop("Method must be 'OLS', 'WLS', 'SUR', '2SLS', 'W2SLS' or '3SLS'")
    }
    if ((method == "2SLS" | method == "W2SLS" | 
        method == "3SLS") & is.null(inst)) {
        stop("Methods '2SLS', 'W2SLS' and '3SLS' need instruments!")
    }

    # Settings:
    # Results to be returned
    results = list()       
    # Results for the individual equations
    results$eq = list()    
    # Results of the ith Equation
    resulti = list()       
    # Residuals equation wise
    residi = list()        
    # Number of Iterations
    iter = NULL            
    # Number of Equations
    G = length(eqns)       
    # Endogenous variables equation wise
    y = list()             
    # Stacked endogenous variables
    Y = matrix(0, 0, 1)    
    # Regressors equation-wise
    x = list()             
    # Stacked matrices of all regressors (unrestricted)
    X = matrix(0, 0, 0)    
    # Number of observations in Each Equation:
    n = array(0, c(G)) 
    # Number of (unrestricted) Coefficients/Regressors in Each Equation:    
    k = array(0, c(G))     
    # List of the instruments for Each Equation:
    instl = list() 
    # Sum of Squared Residuals of Each Equation:       
    ssr = array(0, c(G))   
    # Mean square error (residuals) of Each Equation:
    mse = array(0, c(G))   
    # Root of mse
    rmse = array(0, c(G))  
    # R-squared value
    r2 = array(0, c(G)) 
    # Adjusted R-squared value   
    adjr2 = array(0, c(G)) 
    # Names of Regressors
    xnames = NULL 
    
    # Call:
    call = match.call() 
    # Without ... -expansion:
    m0 = match.call(expand.dots = FALSE) 
    
    # Arguments for Model Matrices:
    temp = c("", "data", "weights", "subset", "na.action")              
    
    # Positions of temp-arguments:
    m0 = m0[match(temp, names(m0), nomatch = 0)]
    m0[[1]] = as.name("model.frame")
            
    # Find Matrices for Individual Models:
    for (i in 1:G) {
        m = m0
        Terms = terms(eqns[[i]], data = data)
        m$formula = Terms
        m = eval(m, parent.frame())
        weights = model.extract(m, "weights")
        y[[i]] = model.extract(m, "response")
        x[[i]] = model.matrix(Terms, m)
        Y = c(Y,y[[i]])
        X = rbind(cbind(X, matrix(0, nrow(X), ncol(x[[i]]))),
            cbind(matrix(0, nrow(x[[i]]), ncol(X)), x[[i]]))
        n[i] = length(y[[i]])
        k[i] = ncol(x[[i]])
        for (j in 1:k[i]) {
            xnames = c(xnames, paste("eq",as.character(i),colnames(x[[i]])[j]))
        }
    }

    # Total Number of Observations:
    N = sum(n)
    # Total Number of (unrestricted) Coefficients/Regressors:
    K = sum(k) 
    # Total Number of Linear Independent Coefficients:
    Ki = K    
    # Total Number of linear Independent Coefficients in Each Equation:  
    ki = k        
   
    if (!is.null(TX)) {
        XU = X
        X = XU %*% TX
        Ki = Ki - (nrow(TX) - ncol(TX))
        for (i in 1:G) {
            ki[i] = ncol(X)
            for (j in 1: ncol(X)) {
                if (sum(X[(1+sum(n[1:i])-n[i]):(sum(n[1:i])),j]^2) == 0) 
                    ki[i] = ki[i]-1
            }
        }
    }
   
    if (!is.null(R.restr)) {
        Ki = Ki - nrow(R.restr)
        if (is.null(TX)) {
            for (j in 1:nrow(R.restr)) {
                for (i in 1:G) {  
                    # Search for restrictions that are NOT cross-equation:
                    if (sum(R.restr[j, 
                        (1+sum(k[1:i])-k[i]):(sum(k[1:i]))]^2) == 
                        sum(R.restr[j,]^2)) {
                    ki[i] = ki[i]-1
                }
                }
            }
        }
    }
    
    # Degress of Freedom of Each Equation:
    df = n - ki    

    # Only for OLS, WLS and SUR Estimation:
    if (method == "OLS" | method == "WLS" | method == "SUR") {
        if (is.null(R.restr)) {
            b = solve(crossprod(X), crossprod(X, Y), tol = solvetol)
        } else {
            W = rbind(cbind(t(X) %*% X, t(R.restr)),
                  cbind(R.restr, matrix(0, nrow(R.restr), nrow(R.restr))))
            V = rbind(t(X) %*% Y , q.restr)
            b = (solve(W, tol = solvetol) %*% V)[1:ncol(X)]
        }
    }

    # Only for OLS Estimation:
    if (method == "OLS") {
        resids = Y - X %*% b                  
        for (i in 1:G) 
            residi[[i]] = resids[(1+sum(n[1:i])-n[i]):(sum(n[1:i]))]
        if (single.eq.sigma) {
            # Residual Covariance Matrix:
            rcov = matrix(0, G, G)   
            for (i in 1:G) 
                rcov[i,i] = sum(residi[[i]]*residi[[i]])/df[i]
            # Omega Inverse:
            Oinv = solve(rcov, tol = solvetol) %x% diag(1,n[1],n[1])      
            if (is.null(R.restr)) {
                # Coefficient Covariance Matrix:
                bcov = solve(t(X) %*% Oinv %*% X, tol = solvetol)
            } else {
                W = rbind(cbind(t(X) %*% Oinv %*% X, t(R.restr)),
                    cbind(R.restr, matrix(0, nrow(R.restr), nrow(R.restr))))
                bcov = solve(W, tol = solvetol)[1:ncol(X),1:ncol(X)]
            }
        } else {
            # Sigma Squared:
            s2 = sum(resids^2)/(N-Ki)       
            if (is.null(R.restr)) {
                # Coefficient Covariance Matrix:
                bcov = s2 * solve(crossprod(X), tol = solvetol)
            } else {
                # Coefficient Covariance Matrix:
                bcov = s2 * solve(W, tol = solvetol)[1:ncol(X),1:ncol(X)]                   
            }
        }
    }

    # Only for WLS estimation:
    if (method == "WLS") {
        # Coefficients of Previous Step:
        bl = b  
        # Difference of Coefficients Between this and Previous Step:  
        bdif = b   
        # Residual Covariance Matrix:
        rcov = matrix(0, G, G)      
        iter = 0
        while ((sum(bdif^2)/sum(bl^2))^0.5>tol & iter < maxiter) {
            iter = iter+1
            # Coefficients of Previous Step:
            bl = b   
            # Residuals             
            resids = Y - X %*% b     
            for (i in 1:G) 
                residi[[i]] = resids[(1+sum(n[1:i])-n[i]):(sum(n[1:i]))]
            for (i in 1:G)  {
                if (rcovformula == 0) {
                    rcov[i,i] = sum(residi[[i]]*residi[[i]])/n[i]
                } else {
                    rcov[i,i] = sum(residi[[i]]*residi[[i]])/df[i]
                }
            }
            # Omega Inverse (= weight. matrix):
            Oinv = solve(rcov, tol = solvetol) %x% diag(1,n[1],n[1])
            if (is.null(R.restr)) {
                # Coefficients
                b = solve(t(X) %*% Oinv %*% X, tol = solvetol) %*% 
                    t(X) %*% Oinv %*%Y
            } else {
                W = rbind(cbind(t(X) %*% Oinv %*% X, t(R.restr)),
                    cbind(R.restr, matrix(0, nrow(R.restr), nrow(R.restr))))
                V = rbind(t(X) %*% Oinv %*% Y , q.restr)
                Winv = solve(W, tol = solvetol)
                # Restricted Coefficients:
                b = (Winv %*% V)[1:ncol(X)]     
            }
            # Difference of Coefficients Between this and Previous Step:
            bdif = b-bl 
        }
    
        if (is.null(R.restr)) {
            # Final Step Coefficient Covariance Matrix:
            bcov = solve(t(X) %*% Oinv %*% X, tol = solvetol)
        } else {
            # Coefficient Covariance Matrix:
            bcov = Winv[1:ncol(X),1:ncol(X)]     
        }
    
        # Residuals:
        resids = Y - X %*% b                        
        for (i in 1:G) 
            residi[[i]] = resids[(1+sum(n[1:i])-n[i]):(sum(n[1:i]))]
    }

    # Only for SUR Estimation:
    if (method == "SUR") {
        # Coefficients of Previous Step:
        bl = b    
        # Difference of Coefficients Between this and Previous Step:
        bdif = b    
        # Residual Covariance Matrix:
        rcov = matrix(0, G, G)   
        iter = 0
        while ((sum(bdif^2)/sum(bl^2))^0.5>tol & iter < maxiter) {
            iter = iter + 1
            # Coefficients of Previous Step:
            bl = b                           
            # Residuals:
            resids = Y - X %*% b                     
            for (i in 1:G) 
                residi[[i]] = resids[(1+sum(n[1:i]) - n[i]):(sum(n[1:i]))]
            for (i in 1:G) {
                for (j in 1:G) {
                    if (rcovformula == 0 | rcovformula == 1) {
                        rcov[i,j] = sum(residi[[i]]*residi[[j]]) / (
                            sqrt((n[i]-rcovformula*ki[i])*(n[j] - 
                                rcovformula*ki[j])))
                    } else {
                        rcov[i,j] = sum(residi[[i]]*residi[[j]])/(
                            n[i]-ki[i]-ki[j] + sum(diag(
                            solve(crossprod(x[[i]]), tol = solvetol) %*%
                            crossprod(x[[i]], x[[j]]) %*%
                            solve(crossprod(x[[j]]), tol = solvetol) %*%
                            crossprod(x[[j]], x[[i]]))))
                    }
                }
            }
            # Omega Inverse (= Weighting Matrix):
            Oinv = solve(rcov, tol = solvetol) %x% diag(1,n[1],n[1])   
            if (is.null(R.restr)) {
                # Coefficients:
                b = solve(t(X) %*% Oinv %*% X, tol = solvetol) %*% 
                    t(X) %*% Oinv %*%Y
            } else {
                W = rbind(cbind(t(X) %*% Oinv %*% X, t(R.restr)),
                    cbind(R.restr, matrix(0, nrow(R.restr), nrow(R.restr))))
                V = rbind(t(X) %*% Oinv %*% Y , q.restr)
                Winv = solve(W, tol = solvetol)
                # Restricted Coefficients:
                b = (Winv %*% V)[1:ncol(X)]     
            }
            # Difference of Coefficients Between this and Previous Step:
            bdif = b - bl  
        }
        if (is.null(R.restr)) {
            # Final Step Coefficient Covariance Matrix:
            bcov = solve(t(X) %*% Oinv %*% X, tol = solvetol)         
        } else {
            # Coefficient Covariance Matrix:
            bcov = Winv[1:ncol(X),1:ncol(X)]     
        }
        # Residuals:
        resids = Y - X %*% b                        
        for (i in 1:G) 
            residi[[i]] = resids[(1+sum(n[1:i])-n[i]):(sum(n[1:i]))]
    }

    # Only for 2SLS, W2SLS and 3SLS Estimation:
    if (method == "2SLS" | method == "W2SLS" | method == "3SLS") {
        for (i in 1:G) {
            if (is.list(inst)) {
                instl[[i]] = inst[[i]]
            } else {
                instl[[i]] = inst
            }
        }
        # Fitted X values:
        Xf = array(0,c(0,ncol(X)))       
        # Stacked matrices of All Instruments:
        H = matrix(0, 0, 0)           
        h = list()
        for (i in 1:G) {
            # Regressors of the i-th Equation (including zeros)
            Xi = X[(1+sum(n[1:i])-n[i]):(sum(n[1:i])),]
            # h[[i]] = model.matrix(instl[[i]])
            # The following lines have to be substituted for the previous
            # Line due to changes in the data handling.
            # Code provided by Ott Toomet
            m = m0
            Terms = terms(instl[[i]], data = data)
            m$formula = Terms
            m = eval(m, parent.frame())
            h[[i]] = model.matrix(Terms, m)
            if (nrow(h[[i]]) != nrow(Xi)) {
                stop(paste(
                    "The instruments and the regressors of Equation", 
                    as.character(i),
                    "have different numbers of observations."))
            }
            # Extract instrument matrix:
            # 'fitted' X-values:
            Xf = rbind(Xf, h[[i]] %*% solve(crossprod(h[[i]]), tol = solvetol)
                  %*% crossprod(h[[i]], Xi))       
            H = rbind(cbind(H, matrix(0, nrow(H), ncol(h[[i]]))),
                cbind(matrix(0, nrow(h[[i]]), ncol(H)), h[[i]]))
        }
        if (is.null(R.restr)) {
            # 2nd Stage Coefficients:
            b = solve(crossprod(Xf), crossprod(Xf, Y), tol = solvetol)     
        } else {
            W = rbind(cbind(crossprod(Xf), t(R.restr)),
                cbind(R.restr, matrix(0, nrow(R.restr), nrow(R.restr))))
            V = rbind(t(Xf) %*% Y , q.restr)
            # Restricted Coefficients:
            b = (solve(W, tol = solvetol) %*% V)[1:ncol(X)] 
        }
        b2 = b
    }

    # Only for 2SLS Estimation:
    if (method == "2SLS") {
        # Residuals:
        resids = Y - X %*% b                        
        for (i in 1:G) 
            residi[[i]] = resids[(1+sum(n[1:i])-n[i]):(sum(n[1:i]))]
        if (single.eq.sigma) {
            # Residual Covariance Matrix:
            rcov = matrix(0, G, G)   
            for (i in 1:G) 
                rcov[i,i] = sum(residi[[i]]*residi[[i]])/(df[i])
            if (is.null(R.restr)) {
                # bcov = solve(t(Xf) %*% Oinv %*% Xf, tol = solvetol)
                # The following 2 lines are substituted for the previous 
                #   line to increase
                # Speed (suggested by Ott Toomet)
                Xfs = Xf*rep(1/diag(rcov), n)
                # Coefficient Covariance Matrix:
                bcov = solve(t(Xfs) %*% Xf, tol = solvetol)                      
            } else {
                # Omega Inverse:
                Oinv = solve(rcov, tol = solvetol) %x% diag(1,n[1],n[1]) 
                W = rbind(cbind(t(Xf) %*% Oinv %*% Xf, t(R.restr)),
                    cbind(R.restr, matrix(0, nrow(R.restr), nrow(R.restr))))
                bcov = solve(W, tol = solvetol)[1:ncol(X),1:ncol(X)]
            }
        } else {
            # Sigma Squared:
            s2 = sum(resids^2)/(N-Ki) 
            if (is.null(R.restr)) {
                # Coefficient Covariance Matrix:
                bcov = s2 * solve(crossprod(Xf), tol = solvetol)            
            } else {
                # Coeff. Covariance Matrix:
                bcov = s2 * solve(W, tol = solvetol)[1:ncol(X),1:ncol(X)]
            }
        }
    }

    # Only for W2SLS Estimation:
    if (method == "W2SLS") {
        # Coefficients of Previous Step:
        bl = b   
        # Difference of Coefficients Between this and Previous Step:
        bdif = b  
        # Residual Covariance Matrix: 
        rcov = matrix(0, G, G) 
        iter = 0
        while ((sum(bdif^2)/sum(bl^2))^0.5>tol & iter < maxiter) {
            iter = iter + 1
            # Coefficients of Previous Step:
            bl = b  
            # Residuals                         
            resids = Y - X %*% b                
            for (i in 1:G) 
                residi[[i]] = resids[(1+sum(n[1:i])-n[i]):(sum(n[1:i]))]
            for (i in 1:G)  {
                if (rcovformula == 0) {
                    rcov[i,i] = sum(residi[[i]]*residi[[i]])/n[i]
                } else {
                    rcov[i,i] = sum(residi[[i]]*residi[[i]])/df[i]
                }
            }
            # Omega Inverse(= weight. matrix):
            Oinv = solve(rcov, tol = solvetol) %x% diag(1,n[1],n[1])        
            if (is.null(R.restr)) {
                # (unrestr.) coeffic.
                b = solve(t(Xf) %*% Oinv %*% Xf, tol = solvetol) %*% 
                    t(Xf) %*% Oinv %*% Y
            } else {
                W = rbind(cbind(t(Xf) %*% Oinv %*% Xf, t(R.restr)),
                    cbind(R.restr, matrix(0, nrow(R.restr), nrow(R.restr))))
                V = rbind(t(Xf) %*% Oinv %*% Y , q.restr)
                Winv = solve(W, tol = solvetol)
                # Restricted Coefficients:
                b = (Winv %*% V)[1:ncol(X)]    
            }
            # Difference of Coefficients Between this and Previous Step:
            bdif = b - bl 
        }
        if (is.null(R.restr)) {
            # Coefficient Covariance Matrix:
            bcov = solve(t(Xf) %*% Oinv %*% Xf, tol = solvetol) 
        } else {
            # Coefficient Covariance Matrix:
            bcov = Winv[1:ncol(X),1:ncol(X)]     
        }
        # Residuals:
        resids = Y - X %*% b                        
        for (i in 1:G) 
            residi[[i]] = resids[(1+sum(n[1:i])-n[i]):(sum(n[1:i]))]
    }

    # Only for 3SLS Estimation:
     if (method == "3SLS") {
        # Coefficients of Previous Step
        bl = b  
        # Difference of Coefficients Between this and Previous Step:
        bdif = b  
        # Residual Covariance Matrix:
        rcov = matrix(0, G, G)  
        iter = 0
        while ((sum(bdif^2)/sum(bl^2))^0.5>tol & iter < maxiter) {
            iter = iter+1
            # Coefficients of Previous Step:
            bl = b                           
            # Residuals:
            resids = Y-X%*%b                     
            for (i in 1:G) 
                residi[[i]] = resids[(1+sum(n[1:i])-n[i]):(sum(n[1:i]))]
            for (i in 1:G) {
                for (j in 1:G) {
                    if (rcovformula == 0 | rcovformula==1) {
                        rcov[i,j] = sum(residi[[i]]*residi[[j]]) / 
                            (sqrt((n[i]-rcovformula*ki[i])*(n[j] - 
                            rcovformula*ki[j])))
                    } else {
                        rcov[i,j] = sum(residi[[i]]*residi[[j]])/(
                            n[i]-ki[i]-ki[j] + sum(diag(
                            solve(crossprod(x[[i]]), tol = solvetol) %*%
                            crossprod(x[[i]], x[[j]]) %*%
                            solve(crossprod(x[[j]]), tol = solvetol) %*%
                            crossprod(x[[j]], x[[i]]))))
                    }
                }
            }
            # Omega Inverse (= Weighting Matrix)
            Oinv = solve(rcov, tol = solvetol) %x% diag(1, n[1], n[1])
            # GLS:
            if (formula3sls == "GLS") {
                # (unrestr.) coeffic.
                if (is.null(R.restr)) {
                    b = solve(t(Xf) %*% Oinv %*% Xf, tol = solvetol) %*% 
                        t(Xf) %*% Oinv %*% Y
                } else {
                    W = rbind(cbind(t(Xf) %*% Oinv %*% Xf, t(R.restr)),
                        cbind(R.restr, matrix(0, nrow(R.restr), nrow(R.restr))))
                    V = rbind(t(Xf) %*% Oinv %*% Y , q.restr)
                    Winv = solve(W, tol = solvetol)
                    # Restricted Coefficients:
                    b = (Winv %*% V)[1:ncol(X)]     
                }
            }
            # IV:
            if (formula3sls == "IV") {
                if (is.null(R.restr)) {
                    # (unrestr.) coeffic.
                    b = solve(t(Xf) %*% Oinv %*% X, tol = solvetol) %*% 
                    t(Xf) %*% Oinv %*% Y                    
                } else {
                    W = rbind(cbind(t(Xf) %*% Oinv %*% X, t(R.restr)),
                        cbind(R.restr, matrix(0, nrow(R.restr), 
                        nrow(R.restr))))
                    V = rbind(t(Xf) %*% Oinv %*% Y , q.restr)
                    Winv = solve(W, tol = solvetol)
                    # Restricted Coefficients:
                    b = (Winv %*% V)[1:ncol(X)]     
                }
            }
            # GMM:
            if (formula3sls == "GMM") {
                if (is.null(R.restr)) {
                    #(unrestr.) coeffic.
                    b = solve(t(X) %*% H %*% solve(t(H) %*% (rcov %x% 
                        diag(1,n[1],n[1])) %*% H, tol = solvetol) %*% t(H) %*% 
                        X, tol = solvetol) %*% t(X) %*% H %*% solve(t(H) %*% 
                        (rcov %x% diag(1,n[1],n[1])) %*% H, tol = solvetol) %*% 
                        t(H) %*% Y  
                } else {
                    W = rbind(cbind(t(X) %*% H %*% solve(t(H) %*% 
                        (rcov %x% diag(1,n[1],n[1])) %*% H, tol = solvetol) %*% 
                        t(H) %*% X, t(R.restr)), cbind(R.restr, matrix(0, 
                        nrow(R.restr), nrow(R.restr))))
                    V = rbind(t(X) %*% H %*% solve(t(H) %*% (rcov %x% 
                        diag(1,n[1],n[1])) %*% H, tol = solvetol) %*% t(H) %*% 
                        Y , q.restr)
                    Winv = solve(W, tol = solvetol)
                    # Restricted Coefficients:
                    b = (Winv %*% V)[1:ncol(X)]     
                }
            }
            # Schmidt:
            if (formula3sls == "Schmidt") {
                if (is.null(R.restr)) {
                    # (unrestr.) coeffic.
                    b = solve(t(Xf) %*% Oinv %*% Xf, tol = solvetol) %*% 
                        (t(Xf) %*% Oinv %*% H %*% solve(crossprod(H), 
                        tol = solvetol) %*% crossprod(H, Y))                
                } else {
                    W = rbind(cbind(t(Xf) %*% Oinv %*% XF, t(R.restr)),
                        cbind(R.restr, matrix(0, nrow(R.restr), 
                        nrow(R.restr))))
                    V = rbind(t(Xf) %*% Oinv %*% H %*% solve(crossprod(H), 
                        tol = solvetol) %*% crossprod(H, Y), q.restr)
                    Winv = solve(W, tol = solvetol)
                    # Restricted Coefficients:
                    b = (Winv %*% V)[1:ncol(X)]     
                }
            }
            # Eviews:
            if (formula3sls == "EViews") {
                if (is.null(R.restr)) {
                    # (unrestr.) coeffic.
                    b = b2 + solve(t(Xf) %*% Oinv %*% Xf, tol = solvetol) %*% 
                        (t(Xf) %*% Oinv %*% (Y -  X %*% b2))   
                } else {
                    W = rbind(cbind(t(Xf) %*% Oinv %*% Xf, t(R.restr)),
                        cbind(R.restr, matrix(0, nrow(R.restr), 
                        nrow(R.restr))))
                    V = rbind(t(Xf) %*% Oinv %*% (Y -  X %*% b2) , q.restr)
                    Winv = solve(W, tol = solvetol)
                    # Restricted Coefficients:
                    b = b2 + (Winv %*% V)[1:ncol(X)]     
                }
            }           
            # Difference of Coefficients Between this and Previous Step
            bdif = b - bl           
        }
            
        
        if (formula3sls == "GLS") {
            if (is.null(R.restr)) {
                # Coefficient Covariance Matrix:
                bcov = solve(t(Xf) %*% Oinv %*% Xf, tol = solvetol)  
            } else {
                # Coefficient Covariance Matrix:
                bcov = Winv[1:ncol(X),1:ncol(X)] 
            }
        }
            
        if (formula3sls == "IV") {
            if (is.null(R.restr)) {
                # Final Step Coefficient Covariance Matrix:
                bcov = solve(t(Xf) %*% Oinv %*% X, tol = solvetol)
            
            } else {
                # Coefficient Covariance Matrix:
                bcov = Winv[1:ncol(X),1:ncol(X)] 
            }
        }
            
        if (formula3sls == "GMM") {
            if (is.null(R.restr)) {
                # Final Step Coefficient Covariance Matrix:
                bcov = solve(t(Xf) %*% Oinv %*% Xf, tol = solvetol)
            } else {
                # Coefficient Covariance Matrix:
                bcov = Winv[1:ncol(X),1:ncol(X)] 
            }
        }
        
        if (formula3sls == "Schmidt") {
            if (is.null(R.restr)) {
                # Final Step Coefficient Covariance Matrix:
                bcov = solve(t(Xf) %*% Oinv %*% Xf, tol = solvetol)
            } else {
                # Coefficient Covariance Matrix:
                bcov = Winv[1:ncol(X),1:ncol(X)]
            }
        }
        
        if (formula3sls == "EViews") {
            if (is.null(R.restr)) {
                # Final Step Coefficient Covariance Matrix:
                bcov = solve(t(Xf) %*% Oinv %*% Xf, tol = solvetol)
            } else {
                W = rbind(cbind(t(Xf) %*% Oinv %*% Xf, t(R.restr)),
                    cbind(R.restr, matrix(0,K-Ki,K-Ki)))
                V = rbind(t(Xf) %*% Oinv %*% Y , q.restr)
                # Coefficient Covariance Matrix:
                bcov = solve(W, tol = solvetol)[1:ncol(X),1:ncol(X)]
            }
        }
        
        # Residuals:
        resids = Y - X %*% b                        
        for (i in 1:G) 
            residi[[i]] = resids[(1+sum(n[1:i]) - n[i]):(sum(n[1:i]))]
    }

    # For all Estimation Methods:
    
    # Fitted endogenous values:
    fitted = X %*% b                              
    bt = NULL
    btcov = NULL
    if (!is.null(TX)) {
        bt = b
        b = TX %*% bt
        btcov = bcov
        bcov = TX %*% btcov %*% t(TX)
    }
    
    # Standard Errors of all estimated Coefficients
    se = diag(bcov)^0.5 
                          
    # T-values of all estimated Coefficients
    t = b/se                                 
    if (probdfsys) {
        # p-values of all estimated Coefficients
        prob = 2*(1-pt(abs(t), N - Ki))            
    } else {
        # p-values of all estimated Coefficients
        prob = matrix(0, 0, 1)                    
    }

    # Equation Wise Results:
    for (i in 1:G) {
        # Estimated Coefficients of Equation i
        bi = b[(1+sum(k[1:i])-k[i]):(sum(k[1:i]))]   
        # Std. errors of est. param. of Equation i
        sei = c(se[(1+sum(k[1:i])-k[i]):(sum(k[1:i]))])
        # T-values of estim. param. of Equation i
        ti = c(t[(1+sum(k[1:i])-k[i]):(sum(k[1:i]))])           
        # Covariance Matrix of Estimated Coefficients of Equation i
        bcovi = bcov[(1+sum(k[1:i])-k[i]):(sum(k[1:i])),
            (1+sum(k[1:i])-k[i]):(sum(k[1:i]))]         
        bi = array(bi,c(k[i],1))
        rownames(bi) = colnames(x[[i]])
        attr(bi,"names") = colnames(x[[i]])

        if (probdfsys) {
            # p-values of estim. param. of Equation i:
            probi = c(prob[(1+sum(k[1:i])-k[i]):(sum(k[1:i]))])
        } else {
            # p-values of estim. param. of Equation i:
            probi = 2*(1 - pt(abs(ti), df[i]))
            # p-values of all estimated Coefficients     
            prob = c(prob,probi) 
        }
    
        # Sum of Squared Residuals:
        ssr = sum(residi[[i]]^2)   
                              
        # Estimated Variance of Residuals:
        mse = ssr/df[i]          
                            
        # Estimated Standard Error of Residuals:
        rmse = sqrt(mse)                                
        r2 = 1 - ssr/(t(y[[i]])%*%y[[i]]-n[i]*mean(y[[i]])^2)
        adjr2 = 1 - ((n[i]-1)/df[i])*(1-r2)
        fittedi = fitted[(1+sum(n[1:i])-n[i]):(sum(n[1:i]))]
        # Datai = model.frame(eqns[[i]])
        #   the following lines have to be substituted for the previous
        #   line due to changes in the data handling.
        #   code provided by Ott Toomet
        m = m0
        Terms = terms(eqns[[i]], data = data)
        m$formula = Terms
        m = eval(m, parent.frame())
        datai = model.frame(Terms, m)
    
        if (method == "2SLS" | method == "3SLS") {
            # Datai = cbind(datai, model.frame(instl[[i]]))
            # The following lines have to be substituted for the previous
            # Line due to changes in the data handling.
            # Code provided by Ott Toomet
            m = m0
            Terms = terms(instl[[i]], data = data)
            m$formula = Terms
            m = eval(m, parent.frame())
            datai = cbind(datai, model.frame(Terms, m))
        }
    
        if (i == 1) {
            # Dataframe for All Data Used for Estimation:
            alldata = datai                    
        } else {
            # Dataframe for All Data Used for Estimation:
            alldata = cbind(alldata, datai)  
        }

        # Build the "return" Structure for the Equation:
        resulti$method = method
        # Equation Number:
        resulti$i = i               
        resulti$eqnlabel = eqnlabels[[i]]
        resulti$formula = eqns[[i]]
        # Number of Observations:
        resulti$n = n[i]            
        # Number of Coefficients/Regressors:
        resulti$k = k[i]            
        # Number of Linear Independent Coefficients:
        resulti$ki = ki[i]           
        # Degrees of freedom of Residuals:
        resulti$df = df[i]           
        # Degrees of Freedom of Residuals of the Whole System:
        resulti$dfSys = N - Ki           
        resulti$probdfsys = probdfsys     
        # Estimated Coefficients:
        resulti$b = c(bi)  
        # Standard Errors of Estimated Coefficients:       
        resulti$se = c(sei)   
        # T-values of Estimated Coefficients:     
        resulti$t = c(ti)    
        # p-values of Estimated Coefficients:     
        resulti$p = c(probi)  
        # Covariance Matrix of Estimated Coefficients:    
        resulti$bcov = bcovi    
        # Vector of Endogenous Variables:       
        resulti$y = y[[i]] 
        # Matrix of Regressors:         
        resulti$x = x[[i]]   
        # Data Frame of this Equation (incl. instruments)       
        resulti$data = datai  
        # Fitted Values:        
        resulti$fitted = fittedi
        # Residuals:      
        resulti$residuals = residi[[i]]
        # Sum of Squared Errors/Residuals:    
        resulti$ssr = ssr             
        # Estimated Variance of the Residuals (Mean Squared Error):
        resulti$mse = mse 
        # The Same (sigma hat squared):            
        resulti$s2 = mse             
        # Estimated Standard Error of the Residuals:
        resulti$rmse = rmse   
        # The Same (sigma hat):        
        resulti$s = rmse            
        # R-sqared Value:
        resulti$r2 = r2    
        # Adjusted R-Squared Value:         
        resulti$adjr2 = adjr2           
        # Matrix of Instrumental Variables:
        if (method == "2SLS" | method == "3SLS") {
          resulti$inst = instl[[i]]
          resulti$h = h[[i]]          
        }
        class(resulti) = "systemfit.equation"
        results$eq[[i]] = resulti
    }

    # Results of the Total System:
    
    # OLS system R2
    # olsr2 = 1 - t(resids) %*% resids / (t(Y) %*% (diag(1, G, G) %x% 
    #   (diag(1, n[1], n[1]) - rep(1, n[1]) %*% t(rep(1,n[1])) / n[1])) %*% Y)
    # The following lines are substituted for the previous 2 lines to 
    #   increase Speed (idea suggested by Ott Toomet)
   
    # Compute Mean of Y by Equations:
    meanY = numeric(length(Y)) 
    for (i in 1:G) {
        meanY[(1+sum(n[1:i])-n[i]):(sum(n[1:i]))] = 
            mean(Y[(1+sum(n[1:i])-n[i]):(sum(n[1:i]))])
    }
    
    # OLS System R2:
    olsr2 = 1 - t(resids) %*% resids / sum((Y - meanY)^2)
                                         
    # Residual Covariance Matrix used for Estimation
    if (method == "SUR" | method == "3SLS") {
        rcovest = rcov                   
    }
    
    # Residual Covariance Matrix:
    rcov = matrix(0, G, G)    
    rcor = matrix(0, G, G)                        
    for (i in 1:G) { 
        for (j in 1:G) {
            rcor[i,j] = sum(residi[[i]] * residi[[j]]) / (
                sqrt(sum(residi[[i]]^2) * sum(residi[[j]]^2)))
            if (rcovformula == 0 | rcovformula==1) {
                rcov[i,j] = sum(residi[[i]]*residi[[j]])/(
                    sqrt((n[i]-rcovformula*ki[i])*(n[j]-rcovformula*ki[j])))
            } else {
                rcov[i,j] = sum(residi[[i]]*residi[[j]])/(
                    n[i]-ki[i]-ki[j] + sum(diag(
                    solve(crossprod(x[[i]]), tol = solvetol) %*%
                    crossprod(x[[i]], x[[j]]) %*%
                    solve(crossprod(x[[j]]), tol = solvetol) %*%
                    crossprod(x[[j]], x[[i]]))))
            }
        } 
    }
    
    drcov = det(rcov, tol = solvetol)
  
    if (!saveMemory) {
        # McElroy's (1977a) R2:
        mcelr2 = 1 - (t(resids) %*% (solve(rcov, tol = solvetol) %x%
            diag(1, n[1],n[1])) %*% resids) / (t(Y) %*% (solve(rcov, 
            tol = solvetol) %x% (diag(1,n[1],n[1]) - rep(1,n[1]) %*%
            t(rep(1,n[1])) / n[1])) %*% Y)   
    } else {
        mcelr2 = NA
    }

    b = c(b)
    names(b) = xnames
    se = c(se)
    names(se) = xnames
    t = c(t)
    names(t) = xnames
    prob = c(prob)
    names(prob) = xnames

    # Build the "return" Structure for the Whole System:
    results$method = method
    # Number of Equations:
    results$g = G              
    # Total number of Observations:
    results$n = N  
    # Total Number of Coefficients:            
    results$k = K        
    # Total number of Linear Independent Coefficients:      
    results$ki = Ki        
    # Degrees of Freedom of the Whole System:     
    results$df = N - Ki   
    # All Estimated Coefficients      
    results$b = b  
    # Transformed vector of Estimated Coefficients            
    results$bt = bt  
    # Standard Errors of Estimated Coefficients:           
    results$se = se  
    # T-values of Estimated Coefficients:           
    results$t = t 
    # p-values of Estimated Coefficients:             
    results$p = prob           
    # Coefficients Covariance Matrix:
    results$bcov = bcov   
    # Covariance Matrix for Transformed Coeff. Vector:       
    results$btcov = btcov      
    # Residual Covariance Matrix:  
    results$rcov = rcov   
    # Determinant of Residual Covariance Matrix:        
    results$drcov = drcov   
    # Residual Correlation Matrix:       
    results$rcor = rcor  
    # R-squared Value of the Equation System:        
    results$olsr2 = olsr2        
    # Residual Correlation Matrix: 
    results$iter = iter  
    # Vector of All (Stacked) Endogenous Variables:         
    results$y = y 
    # Matrix of All (Diagonally Stacked) Regressors:             
    results$x = X    
    # Vector of All (Stacked) Residuals:         
    results$resids = resids 
    # Data Frame for All Data Used in the System:        
    results$data = alldata        
    # Matrix of All (Diagonally Stacked) instr. Variables:
    if (method == "2SLS" | method == "3SLS") {
        results$h = H            
    }
    if (method == "SUR" | method == "3SLS") {
        # Residual Covariance Matrix Used for Estimation:
        results$rcovest = rcovest      
        # McElroy's R-Squared Value for the Equation System:
        results$mcelr2 = mcelr2       
    }
  
    # Continue ...
    results$R.restr = R.restr
    results$q.restr = q.restr
    results$TX = TX
    results$maxiter = maxiter
    results$tol = tol
    results$rcovformula = rcovformula
    results$formula3sls = formula3sls
    results$probdfsys = probdfsys
    results$single.eq.sigma = single.eq.sigma
    results$solvetol = solvetol
    class(results) = "systemfit"

    # Return Value:
    results
}


# ******************************************************************************


.predict.systemfit = 
function(object, data = object$data, se.fit = FALSE, se.pred = FALSE,
interval = "none", level = 0.95, ...) 
{
    # Description:
    #   Calculate predicted values, its Standard Errors and the 
    #   prediction intervals

    # FUNCTION:
    
    # Settings:
    attach(data)
    on.exit(detach(data))
    predicted = data.frame(obs = seq(nrow(data)))
    colnames(predicted) = as.character(1:ncol(predicted))
    g = object$g
    n = array(NA,c(g))
    eqns = list()
    
    # Regressors equation-wise:
    x = list() 
                  
    # Stacked matrices of all regressors (unrestricted):
    X = matrix(0, 0, 0)    
    for (i in 1:g)  {
        eqns[[i]] = object$eq[[i]]$formula
        x[[i]] =  model.matrix(eqns[[i]])
        X =  rbind(cbind(X, matrix(0, nrow(X), ncol(x[[i]]))),
            cbind(matrix(0, nrow(x[[i]]), ncol(X)), x[[i]]))
        n[i] =  nrow(x[[i]])
    }
   
    Y = X %*% object$b
    if (object$method == "SUR" | object$method == "3SLS") {
        if (se.fit | interval == "confidence") {
            ycovc = X %*% object$bcov %*% t(X)
        }
        if (se.pred | interval == "prediction") {
            ycovp = X %*% object$bcov %*% t(X) + 
                object$rcov %x% diag(1,n[1],n[1])
        }
    }
    
    for (i in 1:g) {
        # Fitted values
        Yi = Y[(1+sum(n[1:i])-n[i]):sum(n[1:i]),]
        predicted = cbind(predicted, Yi)
        names(predicted)[length(predicted)] = 
            paste("eq", as.character(i), ".pred", sep = "")
        # Calculate variance covariance matrices
        if (se.fit | interval == "confidence") {
            if (object$method == "SUR" | object$method == "3SLS") {
                ycovci = ycovc[(1 + sum(n[1:i]) - n[i]) : sum(n[1:i]),
                    (1 + sum(n[1:i]) - n[i]) : sum(n[1:i])]
            } else {
                ycovci = x[[i]] %*% object$eq[[i]]$bcov %*% t(x[[i]])
            }
        }
        if (se.pred | interval == "prediction") {
            if (object$method == "SUR" | object$method == "3SLS") {
                ycovpi = ycov[(1 + sum(n[1:i]) - n[i]) : sum(n[1:i]),
                    (1 + sum(n[1:i]) - n[i]) : sum(n[1:i])]
            } else {
                ycovpi = x[[i]] %*% object$eq[[i]]$bcov %*% t(x[[i]]) +
                    object$eq[[i]]$s2
            }
        }
        
        # Standard Errors of fitted values:
        if (se.fit) {
            if (nrow(data) == 1) {
                predicted = cbind(predicted, sqrt(ycovci))
            } else {
                predicted = cbind(predicted, sqrt(diag(ycovci)))
            }
            names(predicted)[length(predicted)] = 
                paste("eq", as.character(i), ".se.fit", sep = "")
        }
      
        # Standard Errors of prediction:
        if (se.pred) {
            if (nrow(data) == 1) {
                predicted = cbind(predicted, sqrt(ycovpi))
            } else {
                predicted = cbind(predicted, sqrt(diag(ycovpi)))
            }
            names(predicted)[length(predicted)] = 
                paste("eq", as.character(i), ".se.pred", sep = "")
        }

        # Confidence intervals:
        if (interval == "confidence") {
            if (object$probdfsys) {
                tval = qt(1 - (1 - level)/2, object$df)
            } else {
                tval = qt(1 - (1 - level)/2, object$eq[[i]]$df)
            }
            if (nrow(data) == 1) {
                seci = sqrt(ycovci)
            } else {
                seci = sqrt(diag(ycovci))
            }
            predicted = cbind(predicted, Yi - (tval * seci))
            names(predicted)[length(predicted)] <-
                paste("eq", as.character(i), ".lwr", sep = "")
            predicted = cbind(predicted, Yi + (tval * seci))
            names(predicted)[length(predicted)] <-
                paste("eq", as.character(i), ".upr", sep = "")
        }
      
        # Prediction intervals:
        if (interval == "prediction") {
            if (object$probdfsys) {
                tval = qt(1 - (1- level)/2, object$df)
            } else {
                tval = qt(1 - (1- level)/2, object$eq[[i]]$df)
            }
            if (nrow(data) == 1) {
                sepi = sqrt(ycovpi)
            } else {
                sepi = sqrt(diag(ycovpi))
            }
            predicted = cbind(predicted, Yi - (tval * sepi))
            names(predicted)[length(predicted)] <-
                paste("eq", as.character(i), ".lwr", sep = "")
            predicted = cbind(predicted, Yi + (tval * sepi))
            names(predicted)[length(predicted)] <-
                paste("eq", as.character(i), ".upr", sep = "")
        }
    }
   
    # Return Value:
    predicted[2: length(predicted)]
}


# ******************************************************************************


.hausmanTest = 
function(li.results, fi.results) 
{
    # Description:
    
    # FUNCTION:
    
    # This function returns test statistic for
    # The hausman test which.... i forget, but people want to see it...
    # From the sas docs
    # given 2 estimators, b0 abd b1, where under the null hypothesis,
    # both are consistent, but only b0 is asympt. efficient and
    # under the alter. hypo only b1 is consistent, so the statistic (m) is
    
    # The m-statistic is then distributed with k degrees of freedom, where k
    # is the rank of the matrix .A generalized inverse is used, as
    # Recommended by Hausman (1982).
    
    # you need to fix this up to return the test statistic, df, and the p value
    
    # Man is this wrong...

    # Build the variance-Covariance Matrix
    # For the full information and the limited information
    # Matricies
    
    ficovb = NULL
    licovb = NULL
    lib = NULL
    fib = NULL

    # Build the final large matrix...
    for (i in 1:li.results$g) {
        fitr = NULL
        litr = NULL
        ## get the dimensions of the current matrix
        for (j in 1:li.results$g) {
            if (i == j) {
                litr = cbind(litr, li.results$eq[[i]]$bcov)
            } else {
                ## bind the zero matrix to the row
                di = dim(li.results$eq[[i]]$bcov)[1]
                dj = dim(li.results$eq[[j]]$bcov)[1]
                litr = cbind(litr, matrix(0, di, dj))
            }
        }
        licovb = rbind(licovb, litr)
        # Now add the rows of the parameter estimates
        # To the big_beta matrix to compute the differences
        lib = c(lib, li.results$eq[[i]]$b)
        fib = c(fib, fi.results$eq[[i]]$b)
    }
    vli = licovb
    vfi = fi.results$bcov
    q = fib - lib
    
    hausman = t(q) %*% solve(vli - vfi) %*% q
    
    # Return Value:
    hausman
}


# ------------------------------------------------------------------------------


.lrTest = 
function(resultc, resultu) 
{
    # Description:
    #   Likelihood Ratio Test
    
    # FUNCTION:
    
    lrtest = list()
    if (resultc$method == "SUR" & resultu$method == "SUR") {
        n = resultu$eq[[1]]$n
        lrtest$df = resultu$ki - resultc$ki
        if (resultc$rcovformula != resultu$rcovformula) {
          stop(paste("both estimations must use the same formula to calculate",
                       "the residual Covariance Matrix!"))
        }
        if (resultc$rcovformula == 0) {
          lrtest$lr = n * (log(resultc$drcov) - log(resultu$drcov))
        } else {
          residc = array(resultc$resids,c(n,resultc$g))
          residu = array(resultu$resids,c(n,resultu$g))
          lrtest$lr = n * (log(det((t(residc) %*% residc))) -
                             log(det((t(residu) %*% residu))))
        }
        lrtest$p = 1-pchisq(lrtest$lr, lrtest$df)
    }
    
    # Return Value:
    lrtest
}


################################################################################
# DON'T USE - NOT YET READY !!


#   $Id: nlsystemfit.r,v 1.13 2004/08/18 17:33:14 hamannj Exp $ 
#
#            Simultaneous Nonlinear Least Squares for R
#
# Copyright 2003-2004 Jeff D. Hamann <jeff.hamann@forestinformatics.com>
#
# This file is part of the nlsystemfit library for R and related languages.
# It is made available under the terms of the GNU General Public
# License, version 2, or at your option, any later version,
# incorporated herein by reference.
#
# This program is distributed in the hope that it will be
# useful, but WITHOUT ANY WARRANTY; without even the implied
# warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
# PURPOSE.  See the GNU General Public License for more
# Details.
#
# You should have received a copy of the GNU General Public
# License along with this program; if not, write to the Free
# Software Foundation, Inc., 59 Temple Place - Suite 330, Boston,
# MA 02111-1307, USA

# uses the Dennis + Schnabel Minimizer which is the one utilized by R's nlm()


# ------------------------------------------------------------------------------


# This function is the "driver" function for the minimization...
.knls = 
function(theta, eqns, data, fitmethod = "OLS", parmnames, instr = NULL, 
S = NULL) 
{

	r = matrix()               # Residuals equation wise
	r = NULL
	
	gmm.resids = matrix()
	gmm.resids = NULL
	
	residi = list()               # Residuals equation wise
	lhs = list()
	rhs = list()
	neqs = length(eqns)
	nobs = dim(data)[[1]]            # Number of nonmissing observations  
	
	## GMM specific variables, in this case... g = 2, k = 3
	#  V = matrix(0, g*k, g*k)          # V is a 6x6 matrix
	moments = list()
	mn = array()
	
	moments = NULL
	mn = NULL
	lhs = NULL
	rhs = NULL
	residi = NULL
	
	## get the values of the parameters
	for (i in 1:length(parmnames)) {
		name = names(parmnames)[i]
		val = theta[i]
		storage.mode(val) =  "double"
		assign(name, val)
	}                          
    
  	## build the residual vector...
  	for (i in 1:length(eqns)) {
	   lhs[[i]] = as.matrix(eval(as.formula(eqns[[i]])[[2]]))
	   rhs[[i]] = as.matrix(eval(as.formula(eqns[[i]])[[3]]))
	   residi[[i]] = lhs[[i]] - rhs[[i]]
	   r = rbind(r, as.matrix(residi[[i]]))    
	   if (fitmethod == "GMM") {
	     	gmm.resids = cbind(gmm.resids, as.matrix(residi[[i]]))
	   }   
 	}

	# These are the objective functions for the various fitting methods
	if (fitmethod == "OLS") {
		obj = crossprod(r)
	}
	if (fitmethod == "2SLS") {
		# W is premultiplied == (diag(neqs) %x% W)
		# obj = (t(r) %*% S %*% r)
		obj = crossprod(t(crossprod(r,S)),r)
	}
	if (fitmethod == "SUR") {
		# S is premultiplied == (qr.solve(S) %x% diag(nobs))
		# obj = (t(r) %*% S %*% r)
		obj = crossprod(t(crossprod(r,S)),r)
	}
	if (fitmethod == "3SLS") {
		# S is premultiplied == (qr.solve(S) %x% W)
		# obj = (t(r) %*% S %*% r)
		obj = crossprod(t(crossprod(r,S)),r)
	}
	if (fitmethod == "GMM") {
		# This just can't be correct... or can it...
		# S is a gx x gk matrix
		# g = number of eqns, k = number of inst variables 
		z = as.matrix(model.frame(instr))
		for (t in 1:nobs) {
		  	moments = rbind(moments, gmm.resids[t,] %x% z[t,])
		}
		# Number of Equations
		g = length(eqns)                 
		k = dim(as.matrix(model.frame(inst)))[[2]]
		gk = g*k
		for (i in 1:gk) {
		  	mn = rbind(mn, mean(moments[,i]))
		}
		# obj = (t(nobs*mn) %*% S %*% (nobs*mn)) / nobs
		# obj = (t(mn) %*% S %*% (mn))
		obj = crossprod(t(crossprod(mn,S)),mn)
	}
  
  ## it would be nice to place the gradient and/or hessian attributes...
  ## how can I make this work???
  ## Attr(obj, "gradient") = "hi mom"
  ## Attr(obj, "hessian") = hessian...

  obj
}


# ------------------------------------------------------------------------------


.nlsystemfit = 
function(method = "OLS", eqns, startvals,
eqnlabels = c(as.character(1:length(eqns))), inst = NULL, data = list(),
solvtol=.Machine$double.eps, pl = 0, maxiter = 1000) 
{
  
	attach(data)
	
	## Some tests
	if (!(method == "OLS" | method == "SUR" | method == "2SLS" | method == "3SLS" | method == "GMM")){
		stop("The method must be 'OLS', 'SUR', '2SLS', or '3SLS'")}
	if ((method == "2SLS" | method == "3SLS" | method == "GMM") & is.null(inst)) {
		stop("The methods '2SLS', '3SLS' and 'GMM' need instruments!")}
  
	lhs = list()
	rhs = list()
	derivs = list()
	
	results = list()               # Results to be returned
	results$eq = list()            # Results for the individual equations
	resulti = list()               # Results of the ith Equation
	residi = list()               # Residuals equation wise
	iter = NULL                 # Number of iterations
	G = length(eqns)       # Number of Equations
	n = array(0, c(G))      # Number of observations in Each Equation
	k = array(0, c(G))     # Number of (unrestricted) coefficients/regressors in Each Equation
	df = array(0, c(G))     # Degrees of freedom in Each Equation
	instl = list()               # List of the instruments for Each Equation
	ssr = array(0, c(G))      # Sum of Squared Residuals of Each Equation
	mse = array(0, c(G))      # Mean square error (residuals) of Each Equation
	rmse = array(0, c(G))      # Root of mse
	r2 = array(0, c(G))      # R-squared value
	adjr2 = array(0, c(G))      # Adjusted R-squared value
	nobs = dim(data)[[1]]
	S = matrix(0, G, G)               # Covariance Matrix of the residuals  
	X = array()
	x = list()
  
	resids = array()
	resids = NULL
	
	if (method == "OLS") {
		est = nlm(.knls, startvals, gradtol = solvtol, typsize = 
			abs(startvals), print.level = pl, iterlim = maxiter, 
			steptol = solvtol, eqns = eqns, data = data, fitmethod = 
			method, parmnames = startvals)
	}
	
  	if (method == "2SLS") {
    	# just fit and part out the return structure 
	    z = as.matrix(model.frame(inst))
	    Wt = z %*% qr.solve(crossprod(z), tol=solvtol) %*% t(z)
	    W2 = diag(length(eqns)) %x% Wt        
	    est = nlm(.knls, startvals, gradtol = solvtol, typsize = 
			abs(startvals), print.level = pl, iterlim = maxiter, 
			steptol = solvtol, eqns = eqns, data = data, fitmethod = 
			method, parmnames = startvals, S = W2)
  
	}
  	if (method == "SUR" || method == "3SLS" || method == "GMM") {
    	# Fit ols/2sls, build the cov matrix for Estimation and refit
    	if (method == "SUR") {
      		estols = nlm(.knls, startvals,
                    gradtol=solvtol,typsize=abs(startvals),print.level=pl,iterlim=maxiter,steptol=solvtol,
                    eqns=eqns,
                    data=data,
                    fitmethod="OLS",
                    parmnames=startvals)
    }
    if (method == "3SLS" || method == "GMM") {
      z = as.matrix(model.frame(inst))
      W = z %*% qr.solve(crossprod(z), tol=solvtol) %*% t(z)
      W2 = (diag(length(eqns)) %x% W)
      estols = nlm(.knls, startvals,
                    gradtol=solvtol,typsize=abs(startvals),print.level=pl,iterlim=maxiter,steptol=solvtol,
                    eqns=eqns,
                    data=data,
                    fitmethod="2SLS",
                    parmnames=startvals,
                    instr=inst,
                    S=W2)
    }
        
    ## build the S matrix
    names(estols$estimate) = names(startvals)
    for (i in 1:length(estols$estimate)) {
      name = names(estols$estimate)[i]
      val = estols$estimate[i]
      storage.mode(val) =  "double"
      assign(name, val)
    }                          
    
    ## get the rank for the eqns, compute the first-stage
    ## Cov matrix to finish the SUR and 3SLS methods
    for (i in 1:G) {
      lhs[[i]] = as.matrix(eval(as.formula(eqns[[i]])[[2]]))
      rhs[[i]] = as.matrix(eval(as.formula(eqns[[i]])[[3]]))
      residi[[i]] = lhs[[i]] - rhs[[i]]
      derivs[[i]] = deriv(as.formula(eqns[[i]]), names(startvals))
      ## Computing the jacobian to get the rank to get the number of variables...
      jacobian = attr(eval(derivs[[i]]), "gradient")
      n[i] =  length(lhs[[i]])
      k[i] = qr(jacobian)$rank
      df[i] = n[i] - k[i]
    }
    
    ## Covariance Matrix of the residuals from the first stage...
    Sols = matrix(0, G, G)
    rcovformula = 1
    for (i in 1:G) {
      for (j in 1:G) {
        Sols[i,j] = sum(residi[[i]]*residi[[j]])/(
                                                   sqrt((n[i]-rcovformula*k[i])*(n[j]-rcovformula*k[j])))
      }
    }

    if (method == "SUR") {
      Solsinv = qr.solve(Sols, tol=solvtol) %x% diag(nobs)      
      est = nlm(.knls,estols$estimate,
                 gradtol=solvtol,typsize=abs(startvals),print.level=pl,iterlim=maxiter,steptol=solvtol,
                 eqns=eqns, data=data, fitmethod=method, parmnames=startvals, S=Solsinv)
    }
    if (method == "3SLS") {
      z = as.matrix(model.frame(inst))
      W = z %*% qr.solve(crossprod(z), tol=solvtol) %*% t(z)
      Solsinv = qr.solve(Sols, tol=solvtol) %x% W      
      est = nlm(.knls, estols$estimate,
                 gradtol=solvtol,typsize=abs(startvals),print.level=pl,iterlim=maxiter,steptol=solvtol,
                 eqns=eqns, data=data, fitmethod=method, parmnames=startvals, S=Solsinv, instr=z)
    }
    if (method == "GMM") {
      resids = NULL
      for (i in 1:G) {
        resids = cbind(resids, residi[[i]])
      }
      z = as.matrix(model.frame(inst))
      moments = list()
      moments = NULL
      for (t in 1:nobs) {
        moments = rbind(moments, resids[t,] %x% z[t,])
      }
      v2sls = qr.solve(var(moments), tol=solvtol)
      est = nlm(.knls,estols$estimate,
                 gradtol=solvtol,typsize=abs(startvals),print.level=pl,iterlim=maxiter,steptol=solvtol,
                 eqns=eqns, data=data, fitmethod="GMM", parmnames=startvals, S=v2sls, instr=inst)
    }
  }

	## Done with the fitting...
	## Now, part out the results from the nlm function
	## To rebuild the equations and return object
	## get the parameters for each of the equations and 

  
	## Evaluate the residuals for eqn
	## get the values of the final parameters
		names(est$estimate) = names(startvals)
		for (i in 1:length(est$estimate)) {
		name = names(est$estimate)[i]
		### I wonder if I need to clear out the variables before assigning them for good measure...
		assign(name, NULL)
		val = est$estimate[i]
		storage.mode(val) =  "double"
		assign(name, val)
	} 

	## get the rank for the eqns, compute the first-stage
	## Cov matrix to finish the SUR and 3SLS methods
	X = NULL
	results$resids = array()
	results$resids = NULL

	## you're working on parsing out the parameters and the estiamtes for the return structure...
	for (i in 1:G) {
		lhs[[i]] = as.matrix(eval(as.formula(eqns[[i]])[[2]]))
		rhs[[i]] = as.matrix(eval(as.formula(eqns[[i]])[[3]]))
		residi[[i]] = lhs[[i]] - rhs[[i]]
		derivs[[i]] = deriv(as.formula(eqns[[i]]), names(startvals))
		## Computing the jacobian to get the rank to get the number of variables...
		jacobian = attr(eval(derivs[[i]]), "gradient")
		n[i] =  length(lhs[[i]])
		k[i] = qr(jacobian)$rank
		df[i] = n[i] - k[i]
		ssr[i] = crossprod(residi[[i]])
		mse[i] = ssr[i] / (n[i] - k[i])
		rmse[i] = sqrt(mse[i])
		X = rbind(X, jacobian)
		results$resids = cbind(results$resids, as.matrix(residi[[i]]))
	}
  
    # Compute the final Covariance Matrix
    # you really should use the code below to handle weights...
    rcovformula = 1
    for (i in 1:G) {
        for (j in 1:G) {
            S[i,j] = sum(residi[[i]]*residi[[j]])/(
                sqrt((n[i]-rcovformula*k[i])*(n[j]-rcovformula*k[j])))
        }
    }


  
    # Get the Variance-Covariance Matrix
  	if (method == "OLS") {
		SI = diag(diag(qr.solve(S, tol=solvtol))) %x% diag(nrow(data))
	    covb = qr.solve(t(X) %*% SI %*% X, tol=solvtol)
  	}
  	if (method == "2SLS") {
	    Z = model.matrix(inst)
	    W = Z %*% qr.solve(crossprod(Z), tol=solvtol) %*% t(Z)
	    SW = diag(diag(qr.solve(S, tol=solvtol))) %x% W
	    covb = qr.solve(t(X) %*% SW %*% X, tol=solvtol)
  	}
  	if (method == "SUR") {
	    SI = qr.solve(S, tol=solvtol) %x% diag(nrow(data))
	    covb = qr.solve(t(X) %*% SI %*% X, tol=solvtol)
  	}
  	if (method == "3SLS") {
	    Z = model.matrix(inst)
	    W = Z %*% qr.solve(crossprod(Z), tol=solvtol) %*% t(Z)
	    SW = qr.solve(S, tol=solvtol) %x% W
	    covb = qr.solve(t(X) %*% SW %*% X, tol=solvtol)
  	}
  	if (method == "GMM") {
		# print("obtaining GMM(SE)")
		z = as.matrix(model.frame(inst))
		moments = list()
		moments = NULL
		resids = NULL
		for (i in 1:G) {
		  	resids = cbind(resids, residi[[i]])
		}
		for (t in 1:nobs) {
		  	moments = rbind(moments, resids[t,] %x% z[t,])
		}
		# print(var(moments))
    	Vinv = qr.solve(var(moments), tol=solvtol)
		# print(Vinv)
    	Y = diag(length(eqns)) %x% t(z)
		# print("covb now...")
		# print(dim(Y))
		# print(dim(X))
        covb = qr.solve(t(Y %*% X) %*% Vinv %*% (Y %*% X ), tol=solvtol)
		# print(covb)
	}

  
  
    ## bind the Standard Errors to the parameter estimate matrix
    se2 = sqrt(diag(covb))
    t.val = est$estimate / se2
    prob = 2*(1 - pt(abs(t.val), sum(n) - sum(k))) ### you better check this...
    
    results$method = method
    results$n = sum(n)
    results$k = sum(k)
    results$b = est$estimate
    results$se = se2
    results$t = t.val
    results$p = prob

    # Build the results structure...
    for (i in 1:G) {
        # You may be able to shrink this up a little and write the 
        # values directly to the return structure...
        eqn.terms = vector()
        eqn.est = vector()    
        eqn.se = vector()    
        jacob = attr(eval(deriv(as.formula(eqns[[i]]), 
        	names(startvals))), "gradient")
        for (v in 1:length(est$estimate)) {
          	if (qr(jacob[,v])$rank > 0) {
            	eqn.terms = rbind(eqn.terms, name = names(est$estimate)[v])
            	eqn.est = rbind(eqn.est, est$estimate[v])
            	eqn.se = rbind(eqn.se, se2[v])
          	}
        }


		## build the "return" structure for the equations
	    resulti$method = method
	    resulti$i = i               # Equation number
	    resulti$eqnlabel = eqnlabels[[i]]
	    resulti$formula = eqns[[i]]
	    resulti$b = as.vector(eqn.est)
	    names(resulti$b) = eqn.terms
	    resulti$se = eqn.se 
	    resulti$t = resulti$b / resulti$se 
	    resulti$p = 2*(1-pt(abs(resulti$t), n[i] - k[i]))
	    resulti$n = n[i]            # Number of observations
	    resulti$k = k[i]            # Number of coefficients/regressors
	    resulti$df = df[i]           # Degrees of freedom of Residuals    
	    resulti$predicted = rhs[[i]]           # predicted values
	    resulti$residuals = residi[[i]]     # Residuals
	    resulti$ssr = ssr[i]             # Sum of squared errors/residuals
	    resulti$mse = mse[i]             # Estimated Variance of the residuals (mean squared error)
	    resulti$s2 = mse[i]             #        the same (sigma hat squared)
	    resulti$rmse = rmse[i]            # Estimated Standard Error of the residuals
	    resulti$s = rmse[i]            #        the same (sigma hat)    
	
		# you'll need these to compute the correlations...
		# print(paste("eqn ", i))

	    resulti$covb = covb[(1+sum(k[1:i])-k[i]):(sum(k[1:i])),
	    	(1+sum(k[1:i])-k[i]):(sum(k[1:i]))]
	
		# resulti$x = model.frame(as.formula(eqns[[i]])[[3]])
		# print(resulti$x)    
		# print(model.frame(eval(eqns[[i]])))
	    
	    # Fix this to allow for multiple instruments?
	    if (method == "2SLS" | method == "3SLS" | method == "GMM") {
			resulti$inst = inst
			# resulti$inst = inst[[i]]
			# resulti$inst = instl[[i]]
			# Matrix of instrumental variables
			#  Resulti$h = h[[i]]          
	    }

		resulti$r2 = 1 - ssr[i] / ((crossprod(lhs[[i]])) - mean(lhs[[i]])^2 * nobs)
		resulti$adjr2 = 1 - ((n[i]-1)/df[i])*(1-resulti$r2)
		
		class(resulti)   = "nlsystemfit.equation"
		results$eq[[i]] = resulti
  	}

	results$solvtol = solvtol
	results$covb = covb
	results$rcov = S
	results$rcor = cor(results$resids)
	# Det(rcov, tol = solvetol)
	results$drcov = det(results$rcov)          
  
  	if (method == "2SLS" || method == "3SLS") {
    	##      results$h = H            # Matrix of all (diagonally stacked) instrumental variables
  	}
  	
  	if (method == "SUR" || method == "3SLS" || method == "GMM") {
    	results$rcovest = Sols      # Residual Covariance Matrix used for Estimation
    	##results$mcelr2 = mcelr2       # McElroy's R-squared value for the equation system
  	}
  
	## build the "return" structure for the Whole System
	results$method = method
	results$g = G              # Number of Equations
	results$nlmest = est
	
	class(results) = "nlsystemfit.system"
	
	detach(data)
	results
}


# ------------------------------------------------------------------------------


## print the (summary) results that belong to the Whole System
.summary.nlsystemfit.system = 
function(object,...) 
{
  summary.nlsystemfit.system = object
  summary.nlsystemfit.system
}


# ------------------------------------------------------------------------------


## print the results that belong to the Whole System
.print.nlsystemfit.system = 
function(x, digits=6,...) 
{
    object = x

    save.digits = unlist(options(digits = digits))
    on.exit(options(digits=save.digits))

    table = NULL
    labels = NULL

    cat("\n")
    cat("nlsystemfit results \n")
    cat("method: ")

    cat(paste(object$method, "\n\n"))

    cat(paste("convergence achieved after",
        object$nlmest$iterations,"iterations\n"))
    cat(paste("nlsystemfit objective function value:",
        object$nlmest$minimum,"\n\n"))

    for (i in 1:object$g) {
        row = NULL
        row = cbind(round(object$eq[[i]]$n,digits),
            round(object$eq[[i]]$df, digits),
            round(object$eq[[i]]$ssr, digits),
            round(object$eq[[i]]$mse, digits),
            round(object$eq[[i]]$rmse, digits),
            round(object$eq[[i]]$r2, digits),
            round(object$eq[[i]]$adjr2, digits))
        table = rbind(table, row)
        labels = rbind(labels, object$eq[[i]]$eqnlabel)
    }
    
    rownames(table) = c(labels)
    colnames(table) = c("N","DF", "SSR", "MSE", "RMSE", "R2", "Adj R2")

    print.matrix(table, quote = FALSE, right = TRUE)
    cat("\n")

    # Check this code before release...
    if (!is.null(object$rcovest)) {
        cat("The Covariance Matrix of the residuals used for Estimation\n")
        rcov = object$rcovest
        rownames(rcov) = labels
        colnames(rcov) = labels
        print(rcov)
        cat("\n")
        if (min(eigen(object$rcovest, only.values=TRUE)$values) < 0) {
            cat("warning: Covariance Matrix is NOT positive semidefinit!\n")
            cat("\n")
        }
    }

    cat("The Covariance Matrix of the residuals\n")
    rcov = object$rcov
    rownames(rcov) = labels
    colnames(rcov) = labels
    print(rcov)
    cat("\n")
    
    cat("The correlations of the residuals\n")
    rcor = object$rcor
    rownames(rcor) = labels
    colnames(rcor) = labels
    print(rcor)
    cat("\n")
    
    cat("The determinant of the residual Covariance Matrix: ")
    cat(object$drcov)
    cat("\n")

    # Now print the individual equations
    for (i in 1:object$g) {
        print(object$eq[[i]], digits)
    }
}


# ------------------------------------------------------------------------------


## print the (summary) results for a single equation
.summary.nlsystemfit.equation = 
function(object,...) 
{
    summary.nlsystemfit.equation = object
    summary.nlsystemfit.equation
}


# ------------------------------------------------------------------------------


## print the results for a single equation
.print.nlsystemfit.equation = 
function(x, digits=6, ...) 
{
    object = x

    save.digits = unlist(options(digits=digits))
    on.exit(options(digits=save.digits))

    cat("\n")
    cat(paste(object$method, " estimates for ", object$eqnlabel, 
        " (equation ", object$i, ")\n", sep = ""))

    cat("Model Formula: ")
    print(object$formula)
    if (!is.null(object$inst)) {
        cat("Instruments: ")
        print(object$inst)
    }
    
    cat("\n")

    Signif = symnum(object$p, corr = FALSE, na = FALSE,
        cutpoints = c(0,  .001,.01,.05, .1, 1),
        symbols = c("***","**","*","."," "))

    table = cbind(round(object$b,  digits),
        round(object$se, digits),
        round(object$t,  digits),
        round(object$p,  digits),
        Signif)

    rownames(table) = names(object$b)
    colnames(table) = c("Estimate","Std. Error","t value","Pr(>|t|)","")

    print.matrix(table, quote = FALSE, right = TRUE)
    cat("---\nSignif. codes: ",attr(Signif,"legend"),"\n")

    # S ist the variance, isn't it???
    cat(paste("\nResidual Standard Error:", round(object$s, digits),  
        "on", object$df, "degrees of freedom\n"))

    cat(paste("Number of observations:", round(object$n, digits),
        "Degrees of Freedom:", round(object$df, digits),"\n"))

    cat(paste("SSR:",     round(object$ssr,    digits),
        "MSE:", round(object$mse, digits),
        "Root MSE:",   round(object$rmse,  digits), "\n"))

    cat(paste("Multiple R-Squared:", round(object$r2,    digits),
        "Adjusted R-Squared:", round(object$adjr2, digits),
        "\n"))
        
    cat("\n")
}


################################################################################

