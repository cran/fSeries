
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Library General Public
# License as published by the Free Software Foundation; either
# version 2 of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Library General Public License for more details.
#
# You should have received a copy of the GNU Library General
# Public License along with this library; if not, write to the
# Free Foundation, Inc., 59 Temple Place, Suite 330, Boston,
# MA  02111-1307  USA

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
# FUNCTION:               DESCRIPTION:
#  armaSim                 Simulates an AIMA Time Series Process
#   fARMA                   Class Representation for "fARMA" Objects
#  armaFit                 Fits Model Parameters for ARMA Time Series Process
# METHODS:                DESCRIPTION:
#  predict.fARMA           S3: Predicts from an ARMA Time Series Process 
#  print.fARMA             S3: Prints a Fitted ARMA Time Series Object
#  plot.fARMA              S3: Plots Stylized Facts of a Fitted ARMA Object
#  summary.fARMA           S3: Analyzes a Fitted ARMA Time Series Object
#  fitted.values.fARMA     S3: Returns Fitted Values from a Fitted ARMA Object
#  residuals.fARMA         S3: Returns Residuals from a Fitted ARMA Object
# FUNCTION:               DESCRIPTION:
#  armaTrueacf             True ARMA Autocorrelation Function
#  armaRoots               Roots of the ARMA Characteristic Polynomial
################################################################################
# BUILTIN - PACKAGE DESCRIPTION:
#  Package: fracdiff
#  Version: 1.1-1
#  Title: Fractionally differenced ARIMA (p,d,q) models
#  Date: 2004-01-12
#  Author: S original by Chris Fraley <fraley@stat.washington.edu>.
#    R port by Fritz Leisch <leisch@ci.tu-wien.ac.at>;
#    since 2003-12: Martin Maechler
#  Maintainer: Martin Maechler <maechler@stat.math.ethz.ch>
#  Description: Maximum likelihood estimation of the parameters of a 
#    fractionally differenced ARIMA(p,d,q) model (Haslett and Raftery, 
#    Appl.Statistics, 1989).
#  License: GPL version 2 or later
#  Packaged: Mon Jan 12 11:22:27 2004; maechler
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


armaSim = 
function(model = list(ar = c(0.5, -0.5), d = 0, ma = 0.1), n = 100,
innov = NULL, n.start = 100, start.innov = NULL, rand.gen = rnorm, ...) 
{   # A function implemented by Diethelm Wuertz

    # Description:
    #   Simulates an ARIMA Time Series Process
    
    # Note:
    #   Splus-Like argument list
    
    # FUNCTION:
    
    # Simulate:
    if (!is.list(model)) stop("model must be list")
    if (is.null(innov)) innov = rand.gen(n, ...)
    n = length(innov) 
    if (is.null(start.innov)) start.innov = rand.gen(n, ...) 
    n.start = length(start.innov)

    # AR PART:
    p = length(model$ar)
    if (p == 1 && model$ar == 0) p = 0
    if (p) { minroots = min(Mod(polyroot(c(1, -model$ar))))
        if (minroots <= 1) stop("ar part of model is not stationary") }
    
    # MA PART:
    q = length(model$ma)
    if (q == 1 && model$ar == 0) q = 0
    if (n.start < p + q) stop("burn-in must be as long as ar + ma")
    
    # DIFFERENCING:
    ## if (model$d < 0) stop("d must be positive ") 
    dd = length(model$d)    
    if (dd) { 
        # FRACDIFF if "dd" is a non-integer value:
        d = model$d
        if (d != round(d) ) { TSMODEL = "FRACDIFF" }
        else { TSMODEL = "ARIMA" } }
    else {
        d = 0 
        TSMODEL = "ARIMA" } 
    
    # ARMA:
    if (TSMODEL == "ARIMA") {
        x = ts(c(start.innov, innov), start = 1 - n.start) 
        if (length(model$ma)) x = filter(x, c(1, model$ma), sides = 1)
        if (length(model$ar)) x = filter(x, model$ar, method = "recursive")
        x = x[-(1:n.start)]
        if (d > 0) x = diffinv(x, differences = d) }
        
    if (TSMODEL == "FRACDIFF") {
        if (p == 0) ar = 0
        if (q == 0) ma = 0
        mu = 0
        # Use Fortran Routine from R's contributed fracdiff package:
        # This is a BUILTIN function ...
        x = .Fortran("fdsim", as.integer(n), as.integer(p), as.integer(q), 
            as.double(ar), as.double(ma), as.double(d), as.double(mu), 
            as.double(rnorm(n + q)), x = double(n + q), 
            as.double(.Machine$double.xmin), as.double(.Machine$double.xmax), 
            as.double(.Machine$double.neg.eps), as.double(.Machine$double.eps), 
            PACKAGE = "fSeries")$x[1:n] }
               
    # Return Value:
    as.ts(x)
}


# ******************************************************************************


require(methods)


# ------------------------------------------------------------------------------


setClass("fARMA", 
    representation(
        call = "call",
        formula = "formula",
        method = "character",
        parameter = "list",
        data = "data.frame",
        fit = "list",
        residuals = "numeric",
        fitted.values = "numeric",
        title = "character",
        description = "character")  
)


# ------------------------------------------------------------------------------


armaFit = 
function(
formula = x ~ arima(2, 0, 1),  
method = c("CSS-ML", "ML", "CSS", "yw", "burg", "ols", "mle"), 
include.mean = TRUE, fixed = NULL, fracdiff.M = 100, fracdiff.h = -1, 
title = "", description = "", ...)
{   # A function implemented by Diethelm Wuertz

    # Description:
    #   Fits Model Parameters for an ARMA Time Series Process
    
    # Notes:
    #   Valid formulas are:
    #       "ar", "arma", "arima", "fracdiff"
    
    # Example:
    #   x = armaSim(); fit = armaFit(x ~ arima(2, 0, 1)); fit
    
    # FUNCTION:
    
    # Transform x:
    ## x <<- as.ts(as.vector(x))
    
    # Call:
    fit = NULL
    call = match.call()
    M = fracdiff.M
    h = fracdiff.h
    
    # Check for Formula length:
    m = length(formula)
    if (m != 3) stop("Formula misspecified")
    
    # Get Series:
    ts = eval(formula[[2]], + sys.parent())
    
    # Allow for univariate 'timeSeries' Objects:
    # Added 2004-09-04 DW
    if (class(ts) == "timeSeries") ts = as.vector(ts)

    # Check for Method:
    # ar.method       = c("yw", "burg", "ols", "mle")
    # arma.method     = c("CSS")
    # arima.method    = c("CSS-ML", "ML", "CSS")    
    # fracdiff.method = NA
    method = method[1]
    
    # Which Model?
    regexpr("\\(", as.character(formula[3]))
    end = regexpr("\\(", as.character(formula[3]))-1
    tsmodel =  substr(as.character(formula[3]), 1, end)
    
    # Valid Model?
    valid = FALSE
    if (tsmodel == "ar" ) valid = TRUE
    if (tsmodel == "arma") valid = TRUE
    if (tsmodel == "arima") valid = TRUE
    if (tsmodel == "fracdiff") valid = TRUE
    if (!valid) stop("Invalid Formula Specification")
    
    # Internal Function: ar
    if (tsmodel == "ar") {
    .arFit <<- function (x, order, include.mean, fixed = NULL,
        method = c("yw", "burg", "ols", "mle"), M = NULL, h = NULL, ...) {
        # Fit:
        call = match.call()
        method = method[1]
        fit = ar(x = x, aic = FALSE, order.max = order, method = method, 
            demean = include.mean, ...)  
        # Add and Modify:
        fit$call = call
        fit$tstitle = paste("AR(", 
            as.character(order), ") with method: ", method, sep = "")
        fit$order = order
        # Residuals:
        fit$residuals = fit$resid
        fit$fitted.values = x - fit$resid
        fit$sigma2 = fit$var.pred 
        # Coefficients:
        fit$coef = fit$ar
        names(fit$coef) = c(paste("ar", 1:order, sep=""))
        if (include.mean) {
            coeff = c(fit$coef, fit$x.mean)
            names(coeff) = c(names(fit$coef), "intercept") 
            fit$coef = coeff} 
        if (method == "ols") { 
            fit$se.coef = fit$asy.se.coef
            n = sqrt(length(as.vector(fit$se.coef)))
            fit$var.coef = matrix(rep(NA, times = n*n), ncol = n) }
        else { 
            fit$var.coef = fit$asy.var.coef
            fit$se.coef = sqrt(diag(fit$asy.var.coef))  
            if (include.mean) {        
                m = dim(fit$asy.var.coef)[1] + 1
                var.coef = matrix(rep(NA, times = m*m), m, m)
                for ( i in 1:(m-1) ) { 
                    for( j in 1:(m-1) ) {
                        var.coef[i,j] = fit$var.coef[i,j] } }
                fit$var.coef = var.coef
                fit$se.coef = c(fit$se.coef, NA) } }
        fit$x = x
        # Return Value:
        fit } }
 
    # Internal Function: arma
    # Use: tseries::arma
    # BUILTIN: Function - tseries::arma
    arma = function(x, order = c(1, 1), lag = NULL, coef = NULL,
    include.intercept = TRUE, series = NULL, qr.tol = 1e-07, ...) {
        seqN = function(N) { 
            if (0==length(N)) NULL else if (N<=0) NULL else seq(N)}
        err = function(coef) {
            u = double(n)
            u[seqN(max.order)] = 0
            u = .C("arma", as.vector(x, mode = "double"),
                u = as.vector(u), as.vector(coef, mode = "double"),
                as.integer(lag$ar), as.integer(lag$ma), as.integer(ar.l),
                as.integer(ma.l), as.integer(max.order), as.integer(n),
                as.integer(include.intercept), PACKAGE="fSeries")$u
            return(sum(u^2)) }
        resid = function(coef) {
            u = double(n)
            u[seqN(max.order)] = 0
            u = .C("arma", as.vector(x, mode = "double"), u = as.vector(u),
                as.vector(coef, mode = "double"), as.integer(lag$ar),
                as.integer(lag$ma), as.integer(ar.l), as.integer(ma.l),
                as.integer(max.order), as.integer(n),
                as.integer(include.intercept), PACKAGE="fSeries")$u
            return(u) }   
        arma.init = function() {
            k = round(1.1*log(n))
            e = na.omit(drop(ar.ols(x, order.max = k, aic = FALSE,
                demean = FALSE, intercept = include.intercept)$resid))
            ee = embed(e, max.order+1)
            xx = embed(x[-(1:k)], max.order+1)
            if (include.intercept == TRUE) {
                if (is.null(lag$ar)) 
                    coef = lm(xx[,1]~ee[,lag$ma+1])$coef
                else if (is.null(lag$ma))
                    coef = lm(xx[,1]~xx[,lag$ar+1])$coef
                else coef = lm(xx[,1]~xx[,lag$ar+1]+ee[,lag$ma+1])$coef
                coef = c(coef[-1], coef[1]) } 
            else {
                if (is.null(lag$ar))
                    coef = lm(xx[,1]~ee[,lag$ma+1]-1)$coef
                else if (is.null(lag$ma))
                    coef = lm(xx[,1]~xx[,lag$ar+1]-1)$coef
                else coef = lm(xx[,1]~xx[,lag$ar+1]+ee[,lag$ma+1]-1)$coef }
            return(coef) }
        if (!is.null(order) & !is.null(lag)) warning("order is ignored")
        if (is.null(order) & is.null(lag)) stop("order or lag must be given")
        if (is.null(lag) & !is.null(order))
            lag = list(ar=seqN(order[1]), ma=seqN(order[2]))
        lag$ar = unique(lag$ar)
        lag$ma = unique(lag$ma)
        max.order = max(unlist(lag),0)
        ar.l = length(lag$ar)
        ma.l = length(lag$ma)
        if (NCOL(x) > 1)
            stop("x is not a vector or univariate time series")
        if (is.null(series)) series = deparse(substitute(x))
        ists = is.ts(x)
        x = as.ts(x)
        xfreq = frequency(x)
        if (any(is.na(x))) stop("NAs in x")
        if (ists) xtsp = tsp(x)
        n = length(x)
        if (!is.null(unlist(lag)))
            if (min(unlist(lag)) < 1 | max(unlist(lag)) > (n-1))
                stop("invalid lag")
        ncoef = length(unlist(lag))+as.numeric(include.intercept)
        if (is.null(coef)) {
            if (!is.null(unlist(lag))) coef = arma.init()
            else coef = 0 }
        if (length(coef) != ncoef) stop("invalid coef")
        md = optim(coef, err, gr=NULL, hessian=TRUE, ...)
        coef = md$par
        rank = qr(md$hessian, qr.tol)$rank
        if (rank != ncoef) {
            se = rep(NA, ncoef)
            cat("Warning: singular Hessian\n") }
        else {
            di = diag(2*md$value/n*solve(md$hessian))
            if (any(di < 0)) cat("Warning: Hessian negative-semidefinite\n")
            se = sqrt(di) }
        e = resid(coef)
        e[seqN(max.order)] = NA
        f = x-e
        if (ists) {
            attr(e, "tsp") = xtsp
            attr(e, "class") = "ts"
            attr(f, "tsp") = xtsp
            attr(f, "class") = "ts" }
        nam.ar = if (!is.null(lag$ar))
            paste("ar", lag$ar, sep = "")
        else NULL
        nam.ma = if (!is.null(lag$ma))
            paste("ma", lag$ma, sep = "")
        else NULL
        nam.int = if (include.intercept) "intercept" else NULL
        nam.coef = c(nam.ar, nam.ma, nam.int)
        names(coef) = nam.coef
        names(se) = nam.coef
        arma = list(coef = coef, css = md$value, n.used = n,
            residuals = e, fitted.values = f, series = series,
            frequency = xfreq, call = match.call(), asy.se.coef = se,
            lag = lag, convergence = md$convergence,
            include.intercept = include.intercept)
        class(arma) = "arma"
        return(arma) }
    # Continue:
    if (tsmodel == "arma") {
    .armaFit <<- function(x, order, include.mean, fixed, 
        method = NULL, M = NULL, h = NULL, ...){
        # Fit:
        call = match.call()
        method = "CSS" # fix
        fit = arma(x = x, order = order, include.intercept = include.mean, ...)
        # Note, Residuals and Fitted Values are returned from "arma"       
        # Add and Modify:
        # # fit$sigma2 = var(na.remove(fit$residuals))
        fit$sigma2 = var(as.vector(na.omit(fit$residuals)))
        fit$tstitle = paste("ARMA(", 
            as.character(order[1]), ",", 
            as.character(order[2]), ") with method: ", method[1], sep = "")
        fit$order = order
        fit$se.coef = fit$asy.se.coef
        fit$x = x
        # Return Value:
        fit$call = call
        fit } }
    
    # Internal Function: arima
    # Use: stats ...
    if (tsmodel == "arima") {
    .arimaFit <<- function (x, order, include.mean, fixed,  
        method = c("CSS-ML", "ML", "CSS"), M = NULL, h = NULL, ...) {
        # Fit:
        call = match.call()
        fit = arima(x = x, order = order, method=  method[1], 
            include.mean = include.mean, fixed = fixed, ...) 
        # Added:
        fit$tstitle = paste("ARIMA(", 
            as.character(order[1]), ",", as.character(order[2]), ",",
            as.character(order[3]), ") with method: ", method[1], sep="")
        fit$x = x   
        fit$fitted.values = fit$x - fit$residuals
        fit$se.coef = sqrt(diag(fit$var.coef))  
        # Return Value:
        fit$call = call
        fit } }
    
    # Internal Function: 
    # Use: fracdiff ...   
    if (tsmodel == "fracdiff") {
        # Internal Function: 
        BIfracdiff = function(x, nar = 0, nma = 0, ar = rep(NA, max(nar, 1)), 
        ma = rep(NA, max(nma, 1)), dtol = NULL, drange = c(0, 0.5), h = -1, 
        M = 100) {    
            # A Builtin Copy from R's fracdiff Package 
            # Arguments:
            #   x      - time series for the ARIMA model
            #   nar    - number of autoregressive parameters
            #   nma    - number of moving average parameters
            #   ar     - initial autoregressive parameters
            #   ma     - initial moving average parameters
            #   dtol   - desired accurcay for d, by default (and if 
            #            negative), (4th root of machine precision)
            #            is used.  dtol will be changed internally if 
            #            necessary
            #   drange - interval over which the likelihood function is 
            #            to be maximized as a function of d
            #   h      - finite difference interval
            #   M      - number of terms in the likelihood approximation
            #           (see Haslett and Raftery 1989) 
            # FRACDIFF:
            if (any(is.na(x)))
                stop("missing values not allowed in time series")
            if (is.matrix(x) && ncol(x) > 2)
                stop("multivariate time series not allowed")
            n = length(x)
            npq = nar + nma
            npq1 = npq + 1
            lwork = max(npq + 2 * (n + M), 3 * n + (n + 6) * npq + 
                npq %/% 2 + 1, (3 + 2 * npq1) * npq1 + 1)
            ar[is.na(ar)] = 0
            ma[is.na(ma)] = 0
            if (is.null(dtol)) dtol = .Machine$double.eps^0.25 # ~ 1.22e-4
            ## if dtol < 0: the fortran code will choose defaults
            result = .Fortran("fracdf", as.double(x), as.integer(n), 
                as.integer(M), as.integer(nar), as.integer(nma), 
                dtol = as.double(dtol), drange = as.double(drange),
                hood = double(1), d = double(1), ar = as.double(ar), 
                ma = as.double(ma), w = double(lwork), as.integer(lwork), 
                info = integer(1), .Machine$double.xmin, 
                .Machine$double.xmax, .Machine$double.neg.eps,
                .Machine$double.eps, PACKAGE = "fSeries")
            if (result$info) switch(result$info,
                stop("insufficient workspace"),
                stop("error in gamma function"),
                stop("invalid MINPACK input"),
                warning("warning in gamma function"),
                warning("optimization failure"),
                warning("optimization limit reached"))
            hess = .Fortran("fdhpq",
                 as.double(x), hess = double(npq1 * npq1), as.integer(npq1),
                 result$w, PACKAGE = "fSeries")$hess
            temp = .Fortran("fdcov", as.double(x), as.double(result$d),
                 h = as.double(if (missing(h)) -1 else h), hd = double(npq1),
                 cov = hess, as.integer(npq1), cor = hess, as.integer(npq1), 
                 se = double(npq1), result$w, info = integer(1), 
                 PACKAGE = "fSeries")
            if (temp$info) switch(temp$info,
                 warning("warning in gamma function"),
                 warning("singular Hessian"),
                 warning("unable to compute correlation matrix"),
                 stop("error in gamma function"))
            if (npq == 0) {
                result$ar = NULL
                result$ma = NULL }
            nam = "d"
            if (nar) nam = c(nam, paste("ar", 1:nar, sep = ""))
            if (nma) nam = c(nam, paste("ma", 1:nma, sep = ""))
            hess = matrix(hess, nrow = npq1, ncol = npq1, 
                dimnames = list(nam, nam))
            hess[1, ] = temp$hd
            hess[row(hess) > col(hess)] = hess[row(hess) < col(hess)]
            se.ok = temp$info != 0 || temp$info < 3
            list(log.likelihood = result$hood,
                d = result$d, ar = result$ar, ma = result$ma,
                covariance.dpq = array(temp$cov, c(npq1, npq1), 
                list(nam, nam)), stderror.dpq = if (se.ok) temp$se, # else NULL
                correlation.dpq = 
                    if (se.ok) array(temp$cor, c(npq1, npq1)), # else NULL
                h = temp$h, d.tol = result$dtol, M = M, hessian.dpq = hess)}
        # Next:
        .resFRACDIFF <<- function (object) {
            n = 0:object$M
            w = gamma(-object$d+n)/(gamma(-object$d)*gamma(n+1)) 
            filter(object$x, w, sides=1) }
        # Next:
        .fracdiffFit <<- function (x, order, include.mean, fixed, 
            method = "FRACDIFF", M = 100, h = -1) {
            # Fit:
            M = filter
            call = match.call()
            fit = BIfracdiff (x = x, nar = order[1], nma = order[2], 
                ar = rep(NA, max(order[1], 1)), ma = rep(NA, max(order[2], 1)), 
                drange = c(0, 0.5), M = M)
            # Added:
            fit$tstitle = paste("FRACDIFF(", as.character(order[1]), ",", 
                as.character(order[2]), ") with method: ", method[1], sep="")
            fit$x = x   
                fit$coef = c(fit$d, fit$ar, fit$ma)
            namesCoef = "d"
            if (order[1] > 0) {
                names.ar = c(paste("ar", 1:order[1], sep=""))
                namesCoef = c(namesCoef, names.ar) }
            if (order[2] > 0) {
                names.ma = c(paste("ma", 1:order[2], sep=""))
                namesCoef = c(namesCoef, names.ma) }
            names(fit$coef) = namesCoef
            fit$var.coef = fit$correlation.dpq  
            fit$fitted.values = .resFRACDIFF(fit)
            fit$residuals = x - fit$fitted.values
            fit$se.coef = fit$stderror.dpq    
            fit$fracdiff = c(M, h)  
            # Return Value:
            fit$call = call
            fit } }
                
    # Which Function?
    fun = match.fun(paste(".", tsmodel, "Fit", sep=""))
    
    # Which Order?
    start = regexpr("\\(", as.character(formula[3]))+1
    end   = regexpr("\\)", as.character(formula[3]))-1
    order = substr(as.character(formula[3]), start, end)
    if (tsmodel == "ar") {
        order = as.integer(order) }
    if (tsmodel == "arma") {
        pos = regexpr(",", order)
        p = as.integer(substr(order, 1, pos-1))
        q = as.integer(substr(order, pos+1, nchar(order)))
        order = c(p, q)}
    if (tsmodel == "arima") {
        pos = regexpr(",", order)   
        p = as.integer(substr(order, 1, pos-1))
        order = substr(order, pos+2, nchar(order))
        d = as.integer(substr(order, 1, pos-1))
        q = as.integer(substr(order, pos+1, nchar(order)))
        order = c(p, d, q)}
    if (tsmodel == "fracdiff") {
        pos = regexpr(",", order)
        p = as.integer(substr(order, 1, pos-1))
        q = as.integer(substr(order, pos+1, nchar(order)))
        order = c(p, q)}
    
    # Fit:
    filter = 100
    fit = fun(x = ts, order = order, include.mean = include.mean, 
        method = method[1], fixed = fixed, M = M, h = h, ...)  
    # "ols" specific:
    if (method == "ols") {
        se.coef = unlist(fit$se.coef)
        if (include.mean){
            ols.mean = se.coef[1]
            fit$se.coef = c(se.coef[-1], ols.mean) } } 
    fit$call = call
    fit$tsmodel = tsmodel
    fit$class = "fARMA"
    class(fit) = "list"
    
    # Return Value:
    new("fARMA",     
        call = as.call(match.call()),
        formula = as.formula(formula), 
        method = as.character(method),
        parameter = list(include.mean = include.mean, fixed = fixed, 
            fracdiff.M = fracdiff.M, fracdiff.h = fracdiff.h),
        data = as.data.frame(x),
        fit = fit,
        residuals = as.vector(fit$residuals),
        fitted.values = as.vector(fit$fitted.values),
        title = as.character(title), 
        description = as.character(description) )
}


# ******************************************************************************


predict.fARMA = 
function (object, n.ahead = 10, n.back = 50, conf = c(80, 95), 
doplot = TRUE, doprint = TRUE, ...) 
{   # A function implemented by Diethelm Wuertz

    # Description:
    #   Predicts from an ARMA Time Series Process
    
    # FUNCTION:
    
    # Object
    object = object@fit
    class(object) = object$class
    TS = "stats"
    
    # FRACDIFF:
    if (object$tsmodel == "fracdiff") {
        warning("Prediction for FRACDIFF not yet implemented")
        return(NA) }
    
    # Internal Function:
    predict.Arima = function (object, n.ahead = 1, newxreg = NULL, 
        se.fit = TRUE, ...) {
        myNCOL = function(x) if (is.null(x)) 0 else NCOL(x)
        rsd = object$residuals
        xr = object$call$xreg
        xreg = if (!is.null(xr)) eval.parent(xr) else NULL
        ncxreg = myNCOL(xreg)
        if (myNCOL(newxreg) != ncxreg)
            stop("xreg and newxreg have different numbers of columns")
        class(xreg) = NULL
        xtsp = tsp(rsd)
        n = length(rsd)
        arma = object$arma
        coefs = object$coef
        narma = sum(arma[1:4])
        if (length(coefs) > narma) {
            if (names(coefs)[narma + 1] == "intercept") {
                xreg = cbind(intercept = rep(1, n), xreg)
                newxreg = cbind(intercept = rep(1, n.ahead), newxreg)
                ncxreg = ncxreg + 1}
            xm = drop(as.matrix(newxreg) %*% coefs[-(1:narma)])}
        else xm = 0
        if (arma[2] > 0) {
            ma = coefs[arma[1] + 1:arma[2]]
            if (any(Mod(polyroot(c(1, ma))) < 1))
                warning("ma part of model is not invertible")}
        if (arma[4] > 0) {
            ma = coefs[sum(arma[1:3]) + 1:arma[4]]
            if (any(Mod(polyroot(c(1, ma))) < 1))
                warning("seasonal ma part of model is not invertible")}
        z = KalmanForecast(n.ahead, object$mod)
        pred = ts(z[[1]] + xm, start = xtsp[2] + deltat(rsd),
               frequency = xtsp[3])
        if (se.fit) {
            se = ts(sqrt(z[[2]] * object$sigma2),
                 start = xtsp[2] + deltat(rsd),
                 frequency = xtsp[3])
            return(pred, se)}
        else return(pred) }

    # Internal Function:
    predictTS <<- function (object, n.ahead, se.fit, ...) {
        
        # Predict "ar":
        if (object$tsmodel == "ar") {
            # Internal Function "predict.ar":
            # This patches missing predict.ar from R's ts package.
            predict.ar = function(object, newdata, n.ahead = 1, 
            se.fit=TRUE, ...) { 
                if (missing(newdata)) {
                    newdata = eval.parent(parse(text=object$series))
                    if (!is.null(nas = object$call$na.action))
                        newdata = eval.parent(call(nas, newdata)) }
                nser = NCOL(newdata)
                ar = object$ar
                p = object$order
                st = tsp(as.ts(newdata))[2]
                dt = deltat(newdata)
                xfreq = frequency(newdata)
                tsp(newdata) = NULL
                class(newdata) = NULL
                if (NCOL(ar) != nser)
                    stop("number of series in fit and newdata do not match")
                n = NROW(newdata)
                if (nser > 1) {
                    if (is.null(object$x.intercept)) xint = rep(0, nser)
                    else xint = object$x.intercept
                    x = rbind(sweep(newdata, 2, object$x.mean),
                        matrix(rep(0, nser), n.ahead, nser, byrow = TRUE))
                    if (p > 0) {
                        for(i in 1:n.ahead) {
                            x[n+i,] = ar[1,,] %*% x[n+i-1,] + xint
                            if (p > 1) for(j in 2:p)
                                x[n+i,] = x[n+i,] + ar[j,,] %*% x[n+i-j,] }
                        pred = x[n+(1:n.ahead), ] } 
                    else {
                        pred = matrix(xint, n.ahead, nser, byrow=TRUE) }
                    pred = pred + matrix(object$x.mean, n.ahead, nser, 
                        byrow = TRUE)
                    colnames(pred) = colnames(object$var.pred)
                    if (se.fit) {
                        warning(
                        "se.fit not yet implemented for multivariate models")
                        se = matrix(NA, n.ahead, nser) } } 
                else {
                    if (is.null(object$x.intercept)) xint = 0
                    else xint = object$x.intercept
                    x = c(newdata - object$x.mean, rep(0, n.ahead))
                    if (p > 0) {
                        for(i in 1:n.ahead) {
                            x[n+i] = sum(ar * x[n+i - (1:p)]) + xint }
                        pred = x[n+(1:n.ahead)]
                        if (se.fit) {
                            npsi = n.ahead - 1
                            psi = .C("artoma",
                                    as.integer(object$order), as.double(ar),
                                    psi = double(npsi+object$order+1),
                                    as.integer(npsi), PACKAGE = TS)$psi[1:npsi]
                            vars = cumsum(c(1, psi^2))
                            se = sqrt(object$var.pred*vars)[1:n.ahead] } }
                    else {
                        pred = rep(xint, n.ahead)
                        if (se.fit) se = rep(sqrt(object$var.pred), n.ahead) }
                    pred = pred + rep(object$x.mean, n.ahead) }
                pred = ts(pred, start = st + dt, frequency=xfreq)
                if (se.fit) se = ts(se, start = st + dt, frequency=xfreq)
                if (se.fit) return(pred, se) else return(pred) } 
                # predict.ar done
            
             # Continue:
            result = predict.ar(object = object, newdata = object$x, 
                n.ahead=n.ahead, se.fit = se.fit, ...) 
            #result 
            }  # end of "ar" if     
        
        # Predict "arma":
        if (object$tsmodel == "arma") {
            object$arma = c(
                object$order[1], 
                object$order[2], 0, 0, 1, 0, 0)
            object$mod = makeARIMA(
                phi = object$coef[1:object$order[1]], 
                theta = object$coef[(object$order[1]+1):sum(object$order)], 
                Delta = numeric(), 
                kappa = 1e6)
            if (!exists("xreg")) xreg = NULL
            if (!exists("newxreg")) newxreg = NULL
            result = predict.Arima(object = object, n.ahead = n.ahead, 
                newxreg = newxreg, se.fit = se.fit, xreg = xreg, ...) 
            }
        
        # Predict "arima"
        if (object$tsmodel == "arima") {
            if (!exists("xreg")) xreg = NULL
            if (!exists("newxreg")) newxreg = NULL
            result = predict.Arima(object = object, n.ahead = n.ahead, 
                newxreg = newxreg, se.fit = se.fit, xreg = xreg, ...) 
            } 
            
        # Predict "fracdiff" 
        # Not yet implemented ...
        
        # Return Value: 
        result }
    
    # Prediction:
    options(warn = -1)
    pred = predictTS(object = object, n.ahead = n.ahead, se.fit = TRUE)

    # if (doprint) print(pred)
    nint = length(conf)
    upper = lower = matrix(NA, ncol = nint, nrow = length(pred$pred))
    for (i in 1:nint) {
        qq = qnorm(0.5 * (1 + conf[i]/100))
        lower[, i] = pred$pred - qq * pred$se
        upper[, i] = pred$pred + qq * pred$se}    
    colnames(lower) = colnames(upper) = paste(conf, "%", sep = "")
        
    # Colors:
    shadecols = switch(1 + (length(conf) > 1), 7, length(conf):1)
    shadepalette = heat.colors(length(conf))
    col = 1
   
    # Data:  
    data = as.ts(object$x)
    freq = frequency(data)
    start = start(data)
    n = length(data)
    #upper = as.matrix(upper)
    #lower = as.matrix(lower)  
    
    # Plot History:   
    if (doplot) {
        pred.mean = pred$pred
        npred = length(pred.mean)
        ylim = range( c(data[(n-n.back+1):n], pred.mean), na.rm = TRUE)
        ylim = range(ylim, lower, upper, na.rm = TRUE)   
        ylab = paste("Series: ", object$series)
        plot(ts(c(data[(n-n.back+1):n], pred.mean[1], rep(NA, npred-1)), 
            end = tsp(data)[2] + npred/freq, f = freq), ylim = ylim, 
            ylab = ylab)
        title(main = paste(object$tstitle)) }
         
    # Confidence Intervals:
    xx = tsp(data)[2] + (1:npred)/freq
    idx = rev(order(conf))
    if (nint > 1) palette(shadepalette)     
    for (i in 1:nint) { polygon(c(xx, rev(xx)), c(lower[, idx[i]], 
        rev(upper[, idx[i]])), col = shadecols[i], border = FALSE) }
    palette("default")
    
    # Mean:
    lines(ts(pred.mean, start=tsp(data)[2]+1/freq, f = freq), lty = 1, 
        col = 4)
   
    # Printout:
    nconf = length(conf)
    out = pred.mean
    upper = as.matrix(upper)
    lower = as.matrix(lower)
    names = "Forecast"
    for (i in nconf:1) {
        out = cbind(out, lower[, i])
        names = c(names, paste("Low", conf[i])) }
    out = cbind(out, pred.mean)
    names = c(names, "Forecast")
    for (i in 1:nconf) {
        out = cbind(out, upper[, i])
        names = c(names, paste("High", conf[i])) }
    out = round(out, digits=4)[,2:(2*nconf+2)]
    colnames(out) = names[2:(2*nconf+2)]
    if (doprint) print(out)
    options(warn = 0)
    
    # Return Value:
    result = list(pred = pred$pred, se = pred$se)
    invisible(result)
}


# ------------------------------------------------------------------------------


print.fARMA = 
function(x, ...)
{   # A function implemented by Diethelm Wuertz

    # Description:
    #   Prints a Fitted ARMA Time Series Object
    
    # FUNCTION:
    
    # Call:
    # object$call, object$coef
    object = x@fit
    
    # Call:
    cat("\nCall:\n")
    cat(paste(deparse(object$call), sep = "\n", collapse = "\n"), 
        "\n", sep = "")
      
    # Model: 
    cat("\nModel:\n", object$tstitle, "\n", sep = "")
    
    # Coefficients:
    cat("\nCoefficient(s):\n")
    digits = max(4, getOption("digits") - 4) 
    print.default(format(object$coef, digits = digits), print.gap = 2, 
        quote = FALSE)
        
    # Return Value
    cat("\n")
    invisible()
}


# ------------------------------------------------------------------------------


summary.fARMA = 
function (object, doplot = TRUE, ...) 
{   # A function implemented by Diethelm Wuertz

    # Description:
    #   Analyzes a Fitted ARMA Time Series Object
    
    # FUNCTION:
        
    # Initialize:
    ans = NULL
    
    # Fit Call and Model:
    x = object
    object = x@fit
    ans$call = object$call
    ans$tsmodel = object$tstitle
    
    # Calculate Residuals and Variance:
    # ans$residuals = na.remove(object$residuals)
    ans$residuals = as.vector(na.omit(object$residuals))
    if (length(ans$residuals) == 0) { 
        ans$var = 0 }
    if (length(ans$residuals) > 0) { 
        ans$var = var(ans$residuals) }
    ans$sigma2 = object$sigma2
    
    # Generate Coef Matrix:
    tval = object$coef/object$se.coef
    prob = 2 * (1 - pnorm(abs(tval)))
    ans$coefmat = cbind(object$coef, object$se.coef, tval, prob)
    dimnames(ans$coefmat) = list(names(object$coef), 
        c(" Estimate", " Std. Error", " t value", "Pr(>|t|)"))
   
    # More Parameters: aic, etc ...
    if (object$tsmodel == "ar") {
        ans$aic = (object$n.used * (1 + log(2 * pi)) + object$n.used * 
            log(ans$var) + 2 * length(object$coef)) }
    if (object$tsmodel == "arma") {
        ans$aic = (object$n.used * (1 + log(2 * pi)) + object$n.used * 
            log(ans$var) + 2 * length(object$coef))
        ans$css = object$css }
    if (object$tsmodel == "arima") {
        ans$aic = object$aic
        ans$loglik = object$loglik }
    if (object$tsmodel == "fracdiff") {
        doplot = FALSE }
    
    # Print:
    print.fARMA(x)
    
    # Residuals:
    digits = max(4, getOption("digits") - 4)
    if (length(object$residuals) > 2) {
        cat("Residuals:\n")
        rq = structure(quantile(ans$residuals), 
            names = c("Min", "1Q", "Median", "3Q", "Max"))
        print(rq, digits = digits)
        skewness = sum((ans$residuals - mean(ans$residuals))^3 /
            sqrt(var(ans$residuals))^3)/length(ans$residuals)
        kurtosis = sum((ans$residuals - mean(ans$residuals))^4 /
            var(ans$residuals)^2)/length(ans$residuals) - 3 
        stats = structure(c(skewness, kurtosis), 
            names = c("Skewness", "Kurtosis"))
        print(stats, digits = digits) }
    
    # Coef Matrix:
    cat("\nCoefficient(s):\n")
    signif.stars = getOption("show.signif.stars")
    printCoefmat(ans$coefmat, digits = digits, 
        signif.stars = signif.stars, ...)
    
    # Fit:
    cat("\n")
    if (x@fit$tsmodel == "ar") {
        cat("sigma^2 estimated as:       ", 
            format(object$var, digits = digits), "\n")
        cat("AIC Criterion:              ", 
            format(round(object$aic, 2)), "\n") }
    if (x@fit$tsmodel == "arma") {
        cat("sigma^2 estimated as:       ", 
            format(object$sigma2, digits = digits), "\n")
        cat("Conditional Sum-of-Squares: ", 
            format(round(object$css, digits=2)), "\n")
        ## cat("AIC Criterion:              ", 
        ##    format(round(object$aic, digits=2)), "\n") 
        }  
    if (x@fit$tsmodel == "arima") {
        cm = object$call$method
        if (is.null(cm) || cm != "CSS")
            cat(
              "sigma^2 estimated as: ", format(object$sigma2, digits = digits),
            "\nlog likelihood:       ", format(round(object$loglik, 2)),
            "\nAIC Criterion:        ", format(round(object$aic, 2)), 
            "\n", sep = "")
        else
            cat(
              "sigma^2 estimated as: ", format(object$sigma2, digits = digits),
            "\npart log likelihood:  ", format(round(object$loglik,2)),
            "\n", sep = "") }
       
    # Doplot:
    if (doplot) plot.fARMA(x, ...)
    
    # Return Value:
    cat("\n")
    invisible(ans)
}


# ------------------------------------------------------------------------------


plot.fARMA = 
function(x, gof.lag = 10, reference.grid = TRUE, col = "steelblue4", ...)
{   # A function implemented by Diethelm Wuertz

    # Description:
    #   Plots stylized facts of a fitted ARMA object
    
    # FUNCTION:
    
    # Standardized Residuals:
    object = x@fit
    rs = as.vector(na.omit(object$residuals))
    stdres = rs/sqrt(object$sigma2)
    plot(stdres, type = "h", main = "Standardized Residuals", 
        ylab = "Residuals", col = col)
    if (reference.grid) grid()
    abline(h = 0, col = "grey")
    
    # Plot ACF:
    acf(object$residuals, plot = TRUE, main = "ACF of Residuals", 
        na.action = na.pass)
    if (reference.grid) grid()
    
    # QQ Plot of Residuals:
    qqnorm(stdres, xlab = "Normal Quantiles", ylab = "Residual Quantiles", 
        main = "QQ-Plot of Residuals", col = col)
    qqline(stdres, col = "grey")
    if (reference.grid) grid()
    
    # LB Statistic:
    nlag = gof.lag
    pval = numeric(nlag)
    for (i in 1:nlag) pval[i] = Box.test(rs, i, type = "Ljung-Box")$p.value
    plot(1:nlag, pval, xlab = "lag", ylab = "p value", ylim = c(0, 1), 
        main = "Ljung-Box p-values")
    abline(h = 0.05, lty = 2, col = "grey")
    if (reference.grid) grid()

    # Return Value:
    invisible()
}


# ******************************************************************************


fitted.values.fARMA = 
function(object, ...)
{   # A function implemented by Diethelm Wuertz

    # Description:
    #   Returns Fitted Values from a Fitted ARMA Object
    
    # FUNCTION:
    
    # Return Value:
    return(object@fitted.values)
}


# ------------------------------------------------------------------------------


residuals.fARMA = 
function(object, ...)
{   # A function implemented by Diethelm Wuertz

    # Description:
    #   Returns Residuals from a Fitted ARMA Object
    
    # FUNCTION:
    
    # Return Value:
    return(object@fit$residuals)
}


################################################################################
# FUNCTION:                 DESCRIPTION:
#  armaTrueacf               True ARMA Autocorrelation Function
#  armaRoots                 Roots of the ARMA Characteristic Polynomial
################################################################################


armaTrueacf = 
function(model, lag.max = 20, type = "correlation", doplot = TRUE)
{   # A function implemented by Diethelm Wuertz

    # Description:
    #   A synonyme to ARMAacf

    # Notes:
    #   A synonyme for arma.tacf under R. See R's .First.lib.
    #   Implemented from ARMAacf
    
    # FUNCTION:
    
    # Settings:
    lag = 0:lag.max
    result = NA
    if (type=="partial" || type=="p" || type=="both" || type=="b") {
        main = ylab = "True PACF"
        lag = 1:lag.max
        pacf = ARMAacf(model$ar, model$ma, lag.max=lag.max, 
            pacf=TRUE)
        result = data.frame(cbind(lag, pacf))
        if (doplot) {
            plot(x=lag, y=pacf, type = "n", xlab = "Lag", 
                ylab = ylab, main = main, 
                ylim = c(min(c(pacf, 0)), 1) )
            lines(x = lag, y = pacf, type = "h")
            abline(h = 0)}}
    if (type == "correlation" || type == "c" || type == "both" || type=="b") {
        main = ylab = "True ACF"
        lag = 0:lag.max
        acf = ARMAacf(model$ar, model$ma, lag.max=lag.max, 
            pacf=FALSE)
        result = data.frame(cbind(lag, acf))
        if (doplot) {
            plot(x=lag, y = acf, type = "n", xlab = "Lag", 
                ylab = ylab, main = main, 
                ylim = c(min(c(acf, 0)), 1) )
            lines(x=lag, y=acf, type = "h")
            abline(h = 0) } }   
            
    # Return Value:
    result
}


# ------------------------------------------------------------------------------


armaRoots = 
function(coefficients, n.plot = 400, digits = 4, ...)
{   # A function implemented by Diethelm Wuertz

    # Description:
    #   Calculates the roots of a characteristc polynomial

    # FUNCTION:
    
    # Algorithm:
    root = polyroot(c(1, -coefficients))
    real.root = Re(root)
    im.root = Im(root)
    xrange = range(real.root)
    xrange = c(xrange[1] - 1.2*abs(xrange[1]), 
        xrange[2]+1.2 * abs(xrange[2]))
    xplot = seq(xrange[1], xrange[2], length = n.plot)
    fpoly = 1
    for(i in 1:length(coefficients)) {
        fpoly = fpoly - xplot^i * coefficients[i] }
    plot(xplot, fpoly, type = "l", xlab = "B", ylab = "Function", ...)
    title(main = "Polynomial Function vs. B")
    abline(h = 0)
    distance = sqrt(real.root^2 + im.root^2)
    root.mat = cbind(round(real.root, digits = digits),
        round(im.root, digits = digits), 
        round(distance, digits = digits))
    dimnames(root.mat) = list(1:nrow(root.mat), 
        c("re", "im", "dist"))
    size.limit = max(abs(real.root), 1.5, abs(im.root))
    plot(root, xlim = c( - size.limit, size.limit),
        ylim = c( - size.limit, size.limit), 
        xlab = "", ylab = "", ...)
    x = (2*pi/360)*(0:360)
    # symbols(0, 0, circles = 1, add = TRUE, inches = FALSE, col = 6)
    lines(sin(x), cos(x))
    abline(h = 0)
    abline(v = 0)
    title("Roots and Unit Circle", 
        xlab = "Real Part", ylab = "Imaginary Part")
    result = root.mat
        
    # Return Value:
    data.frame(result)
}


# -----------------------------------------------------------------------------


