
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Library General Public
# License as published by the Free Software Foundation; either
# version 2 of the License, or (at your option) any later version.
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
# FUNCTION:             REGRESSION MODELLING:
#  regFit                Wrapper Function for regression Models
#  * lm                   Linear Regression Model
#  * glm                  Generalized Linear Model
#  * gam                  Generalized Additive Model
#  * ppr                  Projection Pursuit Regression Model
#  * mars                 Multivariate Adaptive Regression Spline
#  * polymars             Polytochomous MARS
#  * nnet                 Feedforward Neural Network Model
# S3-METHODS:           DESCRIPTION:
#  print                 Prints results from a regression model fit     
#  plot                  Plots fit and diagnostics for a regression model
#  summary               Summarizes fit and diagnostics for a regression model
#  predict               Predicts values from a fitted regression model
#  fitted.values         Returns fitted values from a fitted regression model
#  residulals            Returns residuals from a fitted regression model
################################################################################
# BUILTIN - PACKAGE DESCRIPTION:
#  Package: mda
#  Version: 0.2-23
#  Author: S original by Trevor Hastie & Robert Tibshirani.  R port by
#    Friedrich Leisch, Kurt Hornik and Brian D. Ripley.
#  Maintainer: Kurt Hornik <Kurt.Hornik@R-project.org>
#  Description: Mixture and flexible discriminant analysis, multivariate
#    additive regression splines (MARS), BRUTO, ...
#  Title: Mixture and flexible discriminant analysis
#  Depends: class, R (>= 1.5.0)
#  License: GPL version 2
#  Packaged: Sat Jan 31 13:31:19 2004; hornik
################################################################################
# BUILTIN - PACKAGE DESCRIPTION:
#  Package: polspline
#  Version: 1.0.5
#  Date: 2004-04-22
#  Title: Polynomial spline routines
#  Author: Charles Kooperberg <clk@fhcrc.org>
#  Maintainer: Charles Kooperberg <clk@fhcrc.org>
#  Depends: R
#  Description: Routines for the polynomial spline fitting routines
#    hazard regression, hazard estimation with flexible tails, logspline,
#    lspec, polyclass, and polymars, by C. Kooperberg and co-authors
#  License: GPL version 2 or newer
#  Packaged: Thu Apr 22 13:59:50 2004; hornik
################################################################################


################################################################################
# MODEL:        PACKAGE     print  plot  summary  print   predict
#                                                 summary
#   lm          base        x      x     x        x       x
#   glm         base        x      -     x        x       x
#   gam         mgcv        x      x     x        x       x
#   ppr         modreg      x      x     x        x       x
#   mars*       mda         -      -     -        -       x 
#   polymars*   polspline   -      x     x        -       x
#   nnet        nnet        x      -     x        x       x
#
#   *BUILTIN
################################################################################


require(methods)


# ------------------------------------------------------------------------------


setClass("fREG", 
    representation(
        call = "call",
        formula = "formula",
        family = "character",
        data = "data.frame",
        method = "character",
        fit = "list",
        title = "character",
        description = "character")  
)
    

# ------------------------------------------------------------------------------


regFit = 
function (formula, family = gaussian(), data = list(), 
method = c("LM", "GLM", "GAM", "PPR", "MARS", "POLYMARS", "NNET"), 
nterms = NA, size = NA, title = "", description = "", ...) 
{   # A function implemented by Diethelm Wuertz

    # Description:
    #   Common function call for several selected regression models.
    
    # Details:
    #   This is a wrapper function for the following regrssion models:
    #   LM          Linear Regression Modelling
    #   GLM         Generalized Linear Modelling
    #   GAM         Generalized Additive Modelling
    #   PPR         Projection Pursuit Regression
    #   MARS        Multivariate Adaptive Regression Splines
    #   POLYMARS    Polytochomous MARS Modeling
    #   NNET        Feedforward Neural Net
    
    # Notes:
    #   Available Methods are
    #   "print", "plot", "summary", and "predict" method
    #   "residuals" and "fitted.values" method
    
    # FUNCTION:
    
    # Settings:
    method = method[1]
    
    # Title:
    if (title == "") {
        if (method == "LM") title = "Linear Regression Modelling"
        if (method == "GLM") title = "Generalized Linear Modelling"
        if (method == "GAM") title = "Generalized Additive Modelling"
        if (method == "PPR") title = "Projection Pursuit Regression"
        if (method == "MARS") title = "Multivariate Adaptive Regression Splines"
        if (method == "POLYMARS") title = "Polytochomous MARS Modeling"
        if (method == "NNET") title = "Feedforward Neural Network Modelling" }  
    
    # Internal Function: lm
    lmFit = function (formula, data, ...)  { 
        # From: R-package: base
        # Fit:
        fit = lm(formula = formula, data = data, ...)   
        # Return Value:
        fit$family = c("", "")
        fit$parameters = as.vector(fit$coefficients)
        fit }
        
    # Internal Function: glm
    glmFit = function (formula, family, data, ...) {
        # From R-package: base
        # Fit:
        fit = glm(formula = formula, family = family, data = data, ...)         
        # Return Value:
        fit$family = c(family$family, family$link)
        fit$parameters = as.vector(fit$coefficients)
        fit }

    # Internal Function: gam
    gamFit = function (formula, family, data, ...) {
        # From R-package: mgcv  
        # Fit:
        fit = gam(formula = formula, family = family, data = data, ...)
        # Return Value:
        fit$family = c(family$family, family$link)
        fit$parameters = as.vector(fit$coefficients)
        fit }

    # Internal Function: ppr
    pprFit = function (formula, data, nterms, ...) {
        # From R package: modreg
        if (is.na(nterms)) stop("Argument nterms must be specified")
        # Fit:
        fit = ppr(formula = formula, data = data, nterms = nterms, ...)     
        # Return Value:
        fit$family = c("", "")
        fit$parameters = c(as.vector(fit$alpha), as.vector(fit$beta))   
        fit }

    # Internal Function: mars
    # BuiltIn Uses: BImars()
    marsFit = function(formula, data, ...) {
        # From R-package: mda
        # Settings:
        m = match.call(expand = FALSE)
        m$contrasts = m$... = NULL
        m[[1]] = as.name("model.frame")
        m = eval(m, parent.frame())
        na.act = attr(m, "na.action")
        Terms = attr(m, "terms")
        attr(Terms, "intercept") = 0
        X = model.matrix(Terms, m, contrasts)
        Y = model.extract(m, response)
        w = model.extract(m, weights)
        if (length(w) ==  0) w = rep(1, nrow(X))
        # Fit:
        fit = BImars(X, Y, w, ...)
        fit$terms = Terms   
        # Return Value:
        fit$family = c("", "")
        fit$parameters = as.vector(fit$coefficients)
        fit }
        
    # Internal Function: polymars
    # BuiltIn Uses: BIpolymars
    polymarsFit = function (formula, data, ...) {
        # From R-Package: polspline
        # Fit:
        pmars.formula =
            function(formula, data = sys.parent(), gcv = 4.0, additive = FALSE,
            knot.space = 3, tolerance = 1e-5, verbose = FALSE, ...) {
            # A function implemented by Diethelm Wuertz
            # Function:
            m = match.call(expand = FALSE)
            m$contrasts = m$... = NULL
            m[[1]] = as.name("model.frame")
            m = eval(m, parent.frame())
            na.act = attr(m, "na.action")
            Terms = attr(m, "terms")
            attr(Terms, "intercept") = 0
            X = model.matrix(Terms, m, contrasts)
            Y = model.extract(m, response)
            # fit is of class "polymars"
            fit = BIpolymars(responses = Y, predictors = X, gcv = gcv,
                additive = additive, knot.space = knot.space, 
                tolerance = tolerance, verbose = verbose, ...)
            fit$terms = Terms
            fit$contrasts = contrasts
            fit$fitted.values = Y - fit$residuals
            # Return Value:
            fit }
        fit = pmars.formula(formula = formula, data = data, ...)
        # Return Value:
        fit$family = c("", "")
        fit$parameters = as.vector(fit$model$coefs)
        fit }
        
    # Internal Function: nnet
    nnetFit = function (formula, data, size, ...) {
        # From R-package: nnet
        if (is.na(size)) stop("Argument size must be specified")
        # Fit:
        fit = nnet.formula(formula = formula, data = data, size = size, ...)    
        # Return Value:
        fit$formula = formula
        fit$family = c("", "")
        fit$parameters = as.vector(fit$wts)
        fit }   
    
    # Fit:
    fit = NULL
    
    # Linear Modelling: [base:lm]
    if (method ==  "LM") 
        fit = lmFit(formula = formula, data = data, ...)
        
    # Generalized Linear Modelling: [base:glm]
    if (method ==  "GLM")
        fit = glmFit(formula = formula, family = family, data = data, ...)
        
    # Generalized Additive Modelling: [mgcv:gam]
    if (method ==  "GAM") 
        fit = gamFit(formula = formula, family = family, data = data, ...)
        
    # Projection Pursuit Regression: [modreg:ppr]
    if (method ==  "PPR") 
        fit = pprFit(formula = formula, data = data, nterms = nterms, ...)
        
    # MARS Regression: [mda:mars]
    if (method ==  "MARS") 
        fit = marsFit(formula = formula, data = data, ...)
        
    # POLYMARS Regression: [polspline:polymars]
    if (method ==  "POLYMARS") 
        fit = polymarsFit(formula = formula, data = data, ...)
        
    # Neural Network Regression: [nnet:nnet]
    if (method ==  "NNET") 
        fit = nnetFit(formula = formula, data = data, size = size, ...)
        
    # Add to Fit:
    object.family = fit$family
    fit$call = match.call() 
    fit$family = family
    fit$residuals = as.vector(fit$residuals)    
    fit$fitted.values = as.vector(fit$fitted.values)
    class(fit) = "list"
    
    # Return Value:
    new("fREG",     
        call = as.call(match.call()),
        formula = as.formula(formula), 
        family = as.character(object.family),
        data = as.data.frame(data),
        method = as.character(method), 
        fit = fit,
        title = as.character(title), 
        description = as.character(description) )
}


# ------------------------------------------------------------------------------


print.fREG =
function(x, ...)
{   # A function implemented by Diethelm Wuertz

    # Description:
    #   Print method for Regression Modelling, an object of class "fREG"
    
    # FUNCTION:
    
    # Settings:
    object = x
    
    # Title:
    cat("\nTitle:\n")
    cat(as.character(object@title), "\n")
    
    # Call:
    cat("\nCall:\n")
    cat(paste(deparse(object@call), sep = "\n", collapse = "\n"), 
        "\n", sep = "") 
        
    # Formula:
    cat("\nFormula:\n")
    cat(as.character(object@formula), "\n")
    
    # Formula:
    if (object@family[1] != "" && object@family[2] != "") {     
        cat("\nFamily:\n")
        cat(as.character(object@family), "\n") }
    
    # Digits:
    digits = max(4, getOption("digits") - 4)
        
    # Model Parameters:
    cat("\nModel Parameters:\n")        
        
        # Regression Model LM:
        if (object@method == "LM") {
            print.default(format(object@fit$coef, digits = digits), 
                print.gap = 2, quote = FALSE) }
        
        # Regression Model GLM:
        if (object@method == "GLM") {
            if (length(object@fit$coef)) {
            if (is.character(co = object@fit$contrasts)) 
                cat("  [contrasts: ", apply(cbind(names(co), co), 
                    1, paste, collapse = "="), "]")
                # cat(":\n")
                print.default(format(object@fit$coefficients,
                    digits = digits), print.gap = 2, quote = FALSE)}
            else { cat("No coefficients\n\n") } }   
        
        # Regression Model GAM:
        if (object@method == "GAM") {
            print.default(format(object@fit$coef, digits = digits), 
                print.gap = 2, quote = FALSE) }     
        
        # Regression Model PPR:
        if (object@method == "PPR") {
            cat("-- Projection Direction Vectors --\n")
            print(object@fit$alpha)
            cat("-- Coefficients of Ridge Terms --\n")
            print(object@fit$beta) }    
        
        # Regression Model MARS:
        if (object@method == "MARS") {      
            Parameters = round(object@fit$parameters, digits = digits)      
            print(data.frame(Parameters)) }             
        
        # Regression Model POLYMARS:
        if (object@method == "POLYMARS") {
            print(object@fit$model) }  
        
        # Regression Model NNET:
        if (object@method == "NNET") {
            cat("   a ",object@fit$n[1], "-", object@fit$n[2], "-", 
                object@fit$n[3], " network", " with ", 
                length(object@fit$wts), " weights\n", sep="")
            cat("   options were -")
            tconn = diff(object@fit$nconn)
            if (tconn[length(tconn)] > object@fit$n[2]+1) 
                cat(" skip-layer connections ")
            if (object@fit$nunits > object@fit$nsunits && 
                !object@fit$softmax) 
                cat(" linear output units ")
            if (object@fit$entropy) 
                cat(" entropy fitting ")
            if (object@fit$softmax) 
                cat(" softmax modelling ")
            if (object@fit$decay[1] > 0) 
                cat(" decay=", object@fit$decay[1], sep="")
            cat("\n")
            Weights = object@fit$wts
            print(Weights) } 
        
    # Residual Variance:
    # cat("\nResidual Variance:\n", var(object@fit$residuals))
    cat("\n")
        
    # Return Value:
    invisible(object)
}               


# ------------------------------------------------------------------------------

    
plot.fREG =
function(x, y, ...)
{   # A function imlemented by Diethelm Wuertz

    # Description:
    #   Plot method for Regression Modelling, an object of class "fREG"
    
    # FUNCTION:
    
    # Settings:
    r = x@fit$residuals
    v = x@fit$fitted.values
    
    # Residuals Plot:
    ts.plot(as.ts(r), xlab = "Index", ylab = "Residuals", 
        main = "Residual Series")
    abline(h = 0, col = "steelblue3")
    
    # Quantile Plot:
    rs = (r - mean(r))/sqrt(var(r))
    span = 5
    lim = c(-span, span)
    qqnorm(rs[abs(rs) < span], xlim = lim, ylim = lim, 
        ylab = "Standardized Residuals")
    qqline(rs, col = "steelblue3")
    
    # Fitted Values vs. Residuals Plot:
    plot(v, r, xlab = "Fitted Values", ylab = "Residuals", 
            main = "Fitted vs. Residual Values")
        
    # Return Value:
    invisible(x)
}       


# ------------------------------------------------------------------------------

    
summary.fREG =
function(object, ...)
{   # A function imlemented by Diethelm Wuertz

    # Description:
    #   Summary method for Regression Modelling, an object of class "fREG"
    
    # FUNCTION:
    
    # Digits:
    digits = max(4, getOption("digits") - 4)
    
    # Print all from print Method:
    print(object)
    
    # Add Residual Variance:
    cat("Residual Variance:\n", var(object@fit$residuals))
    cat("\n\n")

    # Internal Function: fResiduals
    fResiduals = function(x, digits) {
        cat("Non-Weighted Residuals:\n")
        names = c("Min", "1Q", "Median", "3Q", "Max")
        rq = structure(quantile(x), names = names)
        print(rq, digits = digits) 
        names = c("Variance", "StDev", "Skewness", "Kurtosis")
        skewness = sum((x - mean(x))^3/sqrt(var(x))^3)/length(x)
        kurtosis = sum((x - mean(x))^4/var(x)^2)/length(x) - 3
        rq = structure(c(var(x), sqrt(var(x)), skewness, kurtosis), 
            names = names)
        print(rq, digits = digits) 
        print("done")
        cat("\n") 
        invisible() }
        
    # Internal Function: print.summary.LM
    print.summary.LM = function (x, ...) {
        digits = max(4, getOption("digits") - 4)
        symbolic.cor = x$symbolic.cor
        signif.stars = getOption("show.signif.stars")
        # cat("\nCall:\n")
        # cat(paste(deparse(x$call), sep = "\n", collapse = "\n"), 
        #   "\n\n", sep = "")
        resid = x$residuals
        df = x$df; rdf = df[2]
        cat(if (!is.null(x$w) && diff(range(x$w))) 
            "Weighted ", "Residuals:\n", sep = "")
        if (rdf > 5) {
            nam = c("Min", "1Q", "Median", "3Q", "Max")
            rq = if (length(dim(resid)) == 2) 
                structure(apply(t(resid), 1, quantile), 
                    dimnames = list(nam, dimnames(resid)[[2]]))
            else structure(quantile(resid), names = nam)
            print(rq, digits = digits, ...) }
        else if (rdf > 0) {
            print(resid, digits = digits, ...)}
        else {
            cat("ALL", df[1], "residuals are 0: no residual ",
                "degrees of freedom!\n") }
        if (length(x$aliased) == 0) {
            cat("\nNo Coefficients\n") }
        else {
            if (nsingular<-df[3] - df[1]) 
                cat("\nCoefficients: (", nsingular, " not defined ",
                    "because of singularities)\n", sep = "")
            else cat("\nCoefficients:\n")
            coefs = x$coefficients
            if (!is.null(aliased = x$aliased) && any(aliased)) {
                cn = names(aliased)
                coefs = matrix(NA, length(aliased), 4, dimnames = 
                    list(cn, colnames(coefs)))
                coefs[!aliased, ] = x$coefficients }
            printCoefmat(coefs, digits = digits, signif.stars = 
                signif.stars, na.print = "NA", ...) }
        cat("\nResidual standard error:", format(signif(x$sigma, 
            digits)), "on", rdf, "degrees of freedom\n")
        if (!is.null(x$fstatistic)) {
            cat("Multiple R-Squared:", formatC(x$r.squared, 
                digits = digits))
            cat(",  Adjusted R-squared:", formatC(x$adj.r.squared, 
                digits = digits), "\nF-statistic:", 
                formatC(x$fstatistic[1], digits = digits), "on", 
                x$fstatistic[2], "and", x$fstatistic[3], 
                "DF,  p-value:", format.pval(pf(x$fstatistic[1], 
                x$fstatistic[2], x$fstatistic[3], lower.tail = FALSE), 
                digits = digits), "\n") }
        correl = x$correlation
        if (!is.null(correl)) {
            p = NCOL(correl)
            if (p > 1) {
                cat("\nCorrelation of Coefficients:\n")
                if (is.logical(symbolic.cor) && symbolic.cor) {
                    print(symnum(correl, abbr.col = NULL)) }
                else {
                    correl = format(round(correl, 2), nsmall = 2, 
                      digits = digits)
                    correl[!lower.tri(correl)] = ""
                    print(correl[-1, -p, drop = FALSE], quote = FALSE) }} }
        cat("\n")
        invisible() }
        
    # Internal Function: print.summary.GLM
    print.summary.GLM = function (x, ...) {
        digits = max(4, getOption("digits") - 4)
        symbolic.cor = x$symbolic.cor
        signif.stars = getOption("show.signif.stars")
        #cat("\nCall:\n")
        #cat(paste(deparse(x$call), sep = "\n", collapse = "\n"), 
        #   "\n\n", sep = "")
        cat("Deviance Residuals: \n")
        if (x$df.residual > 5) {
            x$deviance.resid = quantile(x$deviance.resid, na.rm = TRUE)
            names(x$deviance.resid) = c("Min", "1Q", "Median", "3Q", 
                "Max") }
        print.default(x$deviance.resid, digits = digits, na = "", 
            print.gap = 2)
        if (length(x$aliased) == 0) {
            cat("\nNo Coefficients\n") }
        else {
            if (!is.null(df = x$df) && (nsingular = df[3] - df[1])) 
                cat("\nCoefficients: (", nsingular, " not defined ",
                "because of singularities)\n", sep = "")
            else cat("\nCoefficients:\n")
            coefs = x$coefficients
            if (!is.null(aliased = x$aliased) && any(aliased)) {
                cn = names(aliased)
                coefs = matrix(NA, length(aliased), 4, dimnames = 
                    list(cn, colnames(coefs)))
                coefs[!aliased, ] = x$coefficients }
            printCoefmat(coefs, digits = digits, signif.stars = 
                signif.stars, na.print = "NA", ...) }
        cat("\n(Dispersion parameter for ", x$family$family, 
            " family taken to be ", format(x$dispersion), ")\n\n", 
            apply(cbind(paste(format.char(c("Null", "Residual"), 
            width = 8, flag = ""), "deviance:"), 
            format(unlist(x[c("null.deviance", "deviance")]), 
            digits = max(5, digits + 1)), " on", 
            format(unlist(x[c("df.null", "df.residual")])), 
            " degrees of freedom\n"), 1, paste, collapse = " "), 
            "AIC: ", format(x$aic, digits = max(4, digits + 1)), 
            "\n\n", "Number of Fisher Scoring iterations: ", 
            x$iter, "\n", sep = "")
        correl = x$correlation
        if (!is.null(correl)) {
            p = NCOL(correl)
            if (p > 1) {
                cat("\nCorrelation of Coefficients:\n")
                if (is.logical(symbolic.cor) && symbolic.cor) {
                    print(symnum(correl, abbr.col = NULL)) }
                else {
                    correl = format(round(correl, 2), nsmall = 2, 
                      digits = digits)
                    correl[!lower.tri(correl)] = ""
                    print(correl[-1, -p, drop = FALSE], quote = FALSE) }}}
        cat("\n")
        invisible() }

    # Internal Function: print.summary.GAM
    print.summary.GAM = function(x, ...) { 
        if (length(x$p.coeff) > 0) {
            cat("Parametric coefficients:\n")
            width = max(nchar(names(x$p.coeff)))
            cat(rep(" ",width), "   Estimate  std. err.    t ratio",
                "    Pr(>|t|)\n", sep = "")
            for (i in 1:length(x$p.coeff))
                cat(formatC(names(x$p.coeff)[i], width = width), " ",
                    formatC(x$p.coeff[i], width=10, digits=5), " ",
                    formatC(x$se[i], width = 10, digits = 4), " ",
                    formatC(x$p.t[i], width = 10, digits = 4), "    ",
                    format.pval(x$p.pv[i]), "\n", sep="") }
        cat("\n")
        if (x$m > 0) { 
            cat("Approximate significance of smooth terms:\n")
            width = max(nchar(names(x$chi.sq)))
            cat(rep(" ",width), "        edf       chi.sq     ",
                "p-value\n", sep = "")
            for (i in 1:x$m)
                cat(formatC(names(x$chi.sq)[i], width = width), " ", 
                    formatC(x$edf[i], width = 10, digits = 4), "   ",
                    formatC(x$chi.sq[i], width = 10, digits = 5), "     ",
                    format.pval(x$s.pv[i]), "\n", sep = "") }
        cat("\nR-sq.(adj) = ", formatC(x$r.sq, digits = 3, width = 5),
            "   Deviance explained = ", formatC(x$dev.expl*100, 
            digits = 3, width = 4), "%", sep = "")
        if (is.null(x$ubre)) 
            cat("\nGCV score = ", formatC(x$gcv, digits = 5), " ", 
                sep = "")
        else 
            cat("\nUBRE score = ", formatC(x$ubre, digits = 5), 
                sep = "")
        cat("  Scale est. = ", formatC(x$scale, digits = 5, 
            width = 8, flag = "-"), "  n = ", x$n, "\n", sep = "") 
        invisible() }           

    # Fit:
    fit = object@fit 
            
    # Regression Model: LM
    if (object@method == "LM") {
        class(fit) = "lm"
        ans = summary.lm(object = fit, ...)
        print.summary.LM(x = ans, ...) }    
    
    # Regression Model: GLM
    if (object@method == "GLM") {
        class(fit) = c("glm", "lm")
        ans = summary.glm(object = fit, ...)
        print.summary.GLM(x = ans, ...) }   
    
    # Regression Model: GAM
    if (object@method == "GAM") {
        class(fit) = "gam"
        ans = summary.gam(object = fit, ...)
        print.summary.GAM(x = ans, ...) }   
    
    # Regression Model: GAM
    if (object@method == "GAM") {
        class(fit) = "gam"
        ans = summary.gam(object = fit, ...)
        print.summary.GAM(x = ans, ...) }
        
    # Regression Model: PPR
    if (object@method == "PPR") {
        # This is what print.ppr produces.
        mu = fit$mu; ml = fit$ml
        cat("Goodness of fit:\n")
        gof = fit$gofn; names(gof) = paste(1:ml, "terms")
        print(format(gof[mu:ml], ...), quote = FALSE)
        # This is what summary.ppr produces.
        if (any(fit$edf > 0)) {
            cat("\nEquivalent df for ridge terms:\n")
            edf = fit$edf
            names(edf) = paste("term", 1:fit$mu)
            print(round(edf, 2), ...)} }        
                
    # Regression Model: MARS
    if (object@method == "MARS") {
        # Use the print Method 
        } 
            
    # Regression Model: POLYMARS
    if (object@method == "POLYMARS") {
        # This is what summary.polymars produces.
        # There is no print.summary.polymars.
        cat("Model Fitting:\n")
        print(fit$fitting)
        if(fit$responses != 1)
            cat("\nResponses:", fit$responses, "\n")
        if(!is.null(fit$Rsquared))
            cat("\nRsquared:",round(fit$Rsquared, 4),"\n") 
        cat("\n") }
        
    # Regression Model: NNET
    if (object@method == "NNET") {
        # Use the print Method
        }       

    # Return Value:
    invisible() 
}


# ******************************************************************************


predict.fREG =
function(object, newdata, type = "response", ...)
{   # A function implemented by Diethelm Wuertz

    # Description:
    #   Predict method for Regression Modelling, an object of class "fREG"
    
    # FUNCTION:
    
    # Fit:
    fit = object@fit
    
    # Regression Model: LM
    if (object@method == "LM") {
        class(fit) = "lm"
        ans = predict(object = fit, newdata = newdata, 
            se.fit = TRUE, type = type, ...) }
            
    # Regression Model: GLM
    if (object@method == "GLM") {
        class(fit) = c("glm", "lm")
        ans = predict(object = fit, newdata = newdata, 
            se.fit = TRUE, type = type, ...) }
            
    # Regression Model: GAM
    if (object@method == "GAM") {
        class(fit) = "gam"
        ans = predict(object = fit, newdata = newdata, 
            se.fit = TRUE, type = type, ...) }  
                
    # Regression Model: PPR
    if (object@method == "PPR") {
        class(fit) = "ppr"
        ans = NULL
        ans$fit = predict(object = fit, newdata = newdata, ...) }
            
    # Regression Model: MARS
    # BuitIn Uses: predictBImars()
    if (object@method == "MARS") {
        ans = predictBImars(object = fit, newdata, ...) }
        
    # Regression Model: POLYMARS
    # BuitIn Uses: predictBIpolymars()
    if (object@method == "POLYMARS") {
        # type not used
        predict.POLYMARS = function (object, newdata, ...) {
            # Settings:
            classify = FALSE
            intercept = TRUE    
            # Continue:
            class(object) = "polymars"
            newdata = as.data.frame(newdata)
            RowNames = row.names(newdata)
            Terms = delete.response(object$terms)
            m = model.frame(Terms, newdata, na.action = na.omit)
            keep = match(row.names(m), RowNames)
            X = model.matrix(Terms, m, contrasts = object$contrasts)    
            # Predict:
            result = predictBIpolymars(object = object, x = X,
                classify = classify, intercept = intercept, ...)
            # Return Value:
            as.vector(result) } 
        ans = NULL
        ans$fit = predictBIpolymars(object = fit, newdata, ...) }    
        
    # Regression Model: NNET
    if (object@method == "NNET") {
        if (type == "response") type = "raw"
        class(fit) = c("nnet.formula", "nnet")
        ans = NULL
        ans$fit = predict(object = fit, newdata, type = type, ...) }        
            
    # Return Value:
    ans$fit = as.vector(ans$fit)
    ans
}


# ******************************************************************************

        
fitted.values.fREG = 
function(object, ...)
{   # A function imlemented by Diethelm Wuertz

    # Description:
    #   Fitted values method for Regression Modelling, an object of 
    #   class "fREG".
    
    # FUNCTION:
    
    # Fitted Values:
    ans = as.vector(object@fit$fitted.values) 
            
    # Return Value:
    ans
}
        

# ------------------------------------------------------------------------------


                        
residuals.fREG = 
function(object, ...)
{   # A function imlemented by Diethelm Wuertz

    # Description:
    #   Residuals method for Regression Modelling, an object of 
    #   class "fREG".
    
    # FUNCTION:
    
    # Residuals:
    ans = as.vector(object@fit$residuals) 
            
    # Return Value:
    ans
}


################################################################################


# Builtin mars:
# This is a partly copy from R's mda package Version 0.2-23


# Functions:
#   BImars
#   predictBImars


# ------------------------------------------------------------------------------


BImars = 
function (x, y, w = rep(1, nrow(x)), wp, degree = 1, nk = max(21, 
2 * ncol(x) + 1), penalty = 2, thresh = 0.001, prune = TRUE, 
trace.mars = FALSE, forward.step = TRUE, prevfit = NULL, ...) 
{
    # Internal Function:
    model.matrix.BImars <-
    function (object, x, which = object$selected.terms, full = FALSE, ...) {
        if (missing(x)) return(object$x)
        x = as.matrix(x)
        dd = dim(x)
        n = dd[1]
        p = dd[2]
        nterms = length(which)
        dir = object$factor
        cut = object$cuts
        if (full) {
            bx = matrix(0, nrow = n, ncol = object$lenb)
            bx[, 1] = 1 }
        else bx = matrix(1, nrow = n, ncol = nterms)
        which = which[-1]
        for (i in seq(along = which)) {
            j = which[i]
            if (all(dir[j, ] == 0)) { stop("error in factor or which") }
            temp1 = 1
            for (k in 1:p) {
                if (dir[j, k] != 0) {
                    temp2 = dir[j, k] * (x[, k] - cut[j, k])
                    temp1 = temp1 * temp2 * (temp2 > 0) } }
            if (full) bx[, j] = temp1
            else bx[, i + 1] = temp1}
        bx }
        
    # MARS:
    this.call = match.call()
    if ((nk%%2) != 1) 
        nk = nk - 1
    x = as.matrix(x)
    np = dim(x)
    n = np[1]
    p = np[2]
    y = as.matrix(y)
    nclass = ncol(y)
    if (is.null(np)) {
        np = c(length(x), 1)
        x = as.matrix(x) }
    if (forward.step) {
        interms = 1
        lenb = nk
        bx = matrix(rep(0, nrow(x) * nk), nrow = n)
        res = matrix(rep(0, nrow(x) * ncol(y)), nrow = n)
        fullin = rep(0, nk)
        cuts = NULL
        factor = NULL }
    else {
        bx = model.matrix.BImars(prevfit, x, full = TRUE)
        interms = ncol(bx)
        lenb = prevfit$lenb
        o = prevfit$all.terms
        fullin = rep(0, ncol(bx))
        fullin[o] = 1
        res = prevfit$res
        factor = prevfit$factor
        cuts = prevfit$cuts
        if (missing(penalty)) 
            penalty = prevfit$penalty
        degree = prevfit$degree
        nk = lenb
        thresh = prevfit$thresh }
    if (missing(penalty) & (degree > 1)) 
        penalty = 3
    if (!missing(wp)) {
        if (any(wp <= 0)) 
            stop("wp should all be positive")
        wp = sqrt(wp/sum(wp))
        y = y * outer(rep(1, n), wp) }
    else wp = NULL
    tagx = x
    storage.mode(tagx) = "integer"
    for (j in 1:p) {
        tagx[, j] = order(x[, j]) }
    bestin = rep(0, nk)
    flag = matrix(rep(0, nk * p), nrow = nk, ncol = p)
    if (is.null(cuts)) 
        cuts = matrix(rep(0, nk * p), nrow = nk, ncol = p)
    if (is.null(factor)) {
        dir = matrix(rep(0, nk * p), nrow = nk, ncol = p) }
    else {
        dir = factor }
    alpha = rep(0, nclass)
    beta = matrix(rep(0, nk * nclass), nrow = nk)
    bestgcv = 0
    storage.mode(y) = "double"
    storage.mode(x) = "double"
    storage.mode(bx) = "double"
    storage.mode(flag) = "integer"
    storage.mode(cuts) = "double"
    storage.mode(dir) = "double"
    storage.mode(res) = "double"
    storage.mode(beta) = "double"
    lenscrat = 1 + n + 2 * n * nk + 4 * nk * nk + 3 * nk + 3 * 
        nk * nclass + 3 * nclass + 28 * n + 51
    junk = .Fortran("marss", as.integer(n), as.integer(n), as.integer(p), 
        as.integer(nclass), as.matrix(y), as.matrix(x), as.double(w), 
        as.matrix(tagx), as.integer(degree), as.integer(nk), 
        as.double(penalty), as.double(thresh), as.logical(forward.step), 
        as.integer(interms), as.logical(prune), bx = as.matrix(bx), 
        fullin = as.integer(fullin), lenb = as.integer(lenb), 
        bestgcv = as.double(bestgcv), bestin = as.integer(bestin), 
        flag = as.matrix(flag), cuts = as.matrix(cuts), dir = as.matrix(dir), 
        res = as.matrix(res), alpha = as.double(alpha), beta = as.matrix(beta), 
        double(lenscrat), integer(4 * nk), trace.mars, PACKAGE = "fSeries")
    lenb = junk$lenb
    all.terms = seq(lenb)[junk$fullin[1:lenb] == 1]
    selected.terms = seq(lenb)[junk$bestin[1:lenb] == 1]
    coefficients = junk$beta[seq(selected.terms), , drop = FALSE]
    residuals = junk$res
    fitted.values = y - residuals
    if (!is.null(wp)) {
        TT = outer(rep(1, n), wp)
        residuals = residuals/TT
        fitted.values = fitted.values/TT
        coefficients = coefficients/outer(rep(1, length(selected.terms)), 
            wp) }
    dir = junk$dir[seq(lenb), , drop = FALSE]
    dimnames(dir) = list(NULL, dimnames(x)[[2]])
    cutss = junk$cuts[seq(lenb), , drop = FALSE]
    x = junk$bx[, selected.terms, drop = FALSE]
    
    # Return Value:
    structure(list(call = this.call, all.terms = all.terms, 
        selected.terms = selected.terms, 
        penalty = penalty, degree = degree, nk = nk, thresh = thresh, 
        gcv = junk$bestgcv, factor = dir, cuts = cutss, 
        residuals = residuals, 
        fitted.values = fitted.values, lenb = junk$lenb, 
        coefficients = coefficients, x = x), 
        class = "mars")
}
 

# ------------------------------------------------------------------------------


predictBImars = 
function (object, newdata, ...) 
{   
    # Internal Function:
    model.matrix.BImars = 
    function (object, x, which = object$selected.terms, full = FALSE, ...) {
        if (missing(x)) return(object$x)
        x = as.matrix(x)
        dd = dim(x)
        n = dd[1]
        p = dd[2]
        nterms = length(which)
        dir = object$factor
        cut = object$cuts
        if (full) {
            bx = matrix(0, nrow = n, ncol = object$lenb)
            bx[, 1] = 1 }
        else bx = matrix(1, nrow = n, ncol = nterms)
        which = which[-1]
        for (i in seq(along = which)) {
            j = which[i]
            if (all(dir[j, ] == 0)) { stop("error in factor or which") }
            temp1 = 1
            for (k in 1:p) {
                if (dir[j, k] != 0) {
                    temp2 = dir[j, k] * (x[, k] - cut[j, k])
                    temp1 = temp1 * temp2 * (temp2 > 0) } }
            if (full) bx[, j] = temp1
            else bx[, i + 1] = temp1}
        bx }
        
    # Remove Response from Data:
    newdata = model.matrix.BImars(delete.response(terms(object)), newdata) 
    
    # Predict:
    result = model.matrix.BImars(object = object, x = newdata, 
        which = object$selected.terms, full = FALSE) %*% 
        object$coefficients
        
    # Return Value:
    list(fit = as.vector(result)) 
}


################################################################################




