
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
# FUNCTION:                 GARCH TIME SERIES MODELLING:
#  garchSim                  Simulates a GARCH Time Series Process
#  garchFit                  Fits Model Parameters for a GARCH Process
# METHODS:                  DESCRIPTION:
#  print.fGARCH              Prints a Fitted GARCH Time Series Object
#  summary.fGARCH            Analyzes a Fitted GARCH Time Series Object
#  print.summary.fGARCH      Prints Summary Report for a Fitted GARCH Object
#  plot.fGARCH               Plots Stylized Facts of a Fitted GARCH Object
#  fitted.values.fGARCH      Returns Fitted Values from a Fitted GARCH Object
#  residuals.fGARCH          Returns Residuals from a Fitted GARCH Object
# FUNCTION:                 GARCH TIME SERIES MODELLING:
#  aparchSim                 Simulates a APARCH Time Series Process
#  aparchFit                 Fits Model Parameters for a APARCH Process
# METHODS:                  DESCRIPTION:
#  print.fAPARCH             Prints a Fitted APARCH Time Series Object
#  summary.fAPARCH           Analyzes a Fitted APARCH Time Series Object
#  print.summary.fAPARCH     Prints Summary Report for a Fitted APARCH Object
#  NA: fitted.values.fAPARCH Returns Fitted Values from a Fitted APARCH Object
#  NA: residuals.fAPARCH     Returns Residuals from a Fitted APARCH Object
################################################################################


fGARCH = 
function(x, ...)
{   # A function implemented by D. Wuertz

    UseMethod("fGARCH")
}


# ******************************************************************************


garchSim = 
function(model = list(omega = 1e-6, alpha = 0.1, beta = 0.8, mu = 0), 
n = 100, innov = NULL, n.start = 100, start.innov = NULL, rand.gen = rnorm, 
...)
{   # A function implemented by Diethelm Wuertz 
    
    # Description:
    #   Generates a GARCH(p,q) process. It is possible to generate
    #   ARCH(p) processes as GARCH(p,0) processes.
    
    # Arguments:
    #   model       - a list to specify the GARCH(p,q) model with the
    #                 following elements:
    #    omega      - Constant variance parameter, a numeric value.
    #    alpha      - Alpha parameter which iterates the squared time 
    #                 series values "x^2", a numeric value or a numeric 
    #                 vector of length "p".
    #    beta       - Beta parameter which iterates the conditional
    #                 variances "h", a numeric value or a numeric 
    #                 vector of length "q".
    #   mu          - Constant mean value.
    #   n           - Number of time series points to be generated.
    #   innov       - The innovations for the generation of the time
    #               - series "x", a vector of length n, If "innov" is 
    #                 set to "NULL", the innovations will be generated 
    #                 with the help of the random number generator
    #                 "rand.gen=rnorm".
    #   n.start     - Number of start values of innovations for the 
    #                 generation of the time series "x".
    #   start.innov - The innovations for the generation of the time
    #               - series "x", a vector of length n. If "start.innov"  
    #                 is set to "NULL", the innovations will be generated 
    #                 with the help of the random number generator
    #                 "rand.gen=rnorm".
    #   rand.gen    - a function which is called to generate the innovations. 
    #                 Usually, "rand.gen" will be a random number generator. 
    #   ...         - Optional Arguments passed to the function "rand.gen".
    
    # Note:
    #   To generate the subclass of ARCH models, use GARCH(p,0).
    #   This is a Splus like function call.
    
    # FUNCTION:
    
    # Doesn't, replace the three following three lines ...
    # If Missing from Model List - Add to List:
    # if (!exists("model$alpha")) model$alpha = 0
    # if (!exists("model$beta")) model$beta = 0
    # if (!exists("model$mu")) model$mu = 0
    
    # with ...
    if (is.null(model$alpha)) model$alpha = 0
    if (is.null(model$beta)) model$beta = 0
    if (is.null(model$mu)) model$mu = 0
    
    # Innovations:
    max.order = max(length(model$alpha), length(model$beta))
    if (n.start < max.order)
        stop("n.start must be greater or equal max(alpha,beta)")   
    if (is.null(start.innov)) start.innov = rand.gen(n.start, ...)  
    if (is.null(innov)) innov = rand.gen(n, ...)
    
    # Setting Start Values and Vectors:
    h = x = z = c(start.innov, innov)  
        
    # Initialization:
    for (i in 1:max.order) {
        h[i] = model$omega/(1-sum(model$alpha)-sum(model$beta))
        x[i] = sqrt(h[i]) * z[i] + model$mu }
        
    # Iteration:
    n.alpha = length(model$alpha)
    n.beta = length(model$beta)
    for (i in (max.order+1):(n.start+n)) {     
        h[i] = model$omega +    
            sum(model$alpha*x[i-(1:n.alpha)]^2) +
            sum(model$beta*h[i-(1:n.beta)]) 
        x[i] = sqrt(h[i]) * z[i] + model$mu }
        
    # Return Value:
    as.ts(x[-(1:n.start)])
}


# ------------------------------------------------------------------------------


garchFit = 
function (x, order = c(1, 1), ...)
{   # A function imlemented by Diethelm Wuertz

    # Description:
    #   Fits the parameters of a GARCH Model 
    
    # Details:
    #   Calls -
    #     garch (x, order = c(1, 1), coef = NULL, itmax = 200, eps = NULL, 
    #       grad = c("analytical", "numerical"), series = NULL, trace = TRUE, 
    #       ...) 
    #   from contributed R-package [tseries:garch]
    
    # Notes:
    #   Use A. Trapletti's [tseries:garch]
    
    # FUNCTION:
    
    # Internsal Function:
    # A copy from Adrian Trapletti's contribute R package, tseries::garch
    BIgarch = function (x, order = c(1, 1), coef = NULL, itmax = 200, 
    eps = NULL, grad = c("analytical", "numerical"), series = NULL, 
    trace = TRUE, ...) {
        if (NCOL(x) > 1) stop("x is not a vector or univariate time series")
        if (!is.vector(order)) stop("order is not a vector")
        grad = match.arg(grad)
        switch(grad,
            analytical = (agrad = TRUE),
            numerical = (agrad = FALSE))
        if (is.null(series)) series = deparse(substitute(x))
        ists = is.ts(x)
        x = as.ts(x)
        xfreq = frequency(x)
        if (any(is.na(x))) stop("NAs in x")
        if (ists) xtsp = tsp(x)
        x = as.matrix(x)
        n = nrow(x)
        e = double(n)
        ncoef = order[1]+order[2]+1
        hess = matrix(0.0, ncoef, ncoef)
        small = 0.05
        if (is.null(coef)) 
            coef = c(var(x)*(1.0-small*(ncoef-1)),rep(small,ncoef-1))
        if (!is.vector(coef)) 
            stop("coef is not a vector")
        if (ncoef != length(coef)) 
            stop("incorrect length of coef")
        if (is.null(eps)) 
            eps = .Machine$double.eps
        nlikeli = 1.0e+10
        fit = .C("fit_garch", as.vector(x, mode = "double"),
            as.integer(n), coef = as.vector(coef, mode = "double"),
            as.integer(order[1]), as.integer(order[2]), as.integer(itmax),
            as.double(eps), nlikeli = as.double(nlikeli), as.integer(agrad),
            as.integer(trace), PACKAGE="fSeries")
        pred = .C("pred_garch", as.vector(x, mode = "double"),
            e = as.vector(e, mode = "double"), as.integer(n),
            as.vector(fit$coef, mode = "double"), as.integer(order[1]),
            as.integer(order[2]), as.integer(FALSE), PACKAGE = "fSeries")
        com.hess = .C("ophess_garch", as.vector(x, mode = "double"),
            as.integer(n), as.vector(fit$coef, mode = "double"),
            hess = as.matrix(hess), as.integer(order[1]),
            as.integer(order[2]), PACKAGE="fSeries")
        rank = qr(com.hess$hess, ...)$rank
        if (rank != ncoef) {
            se.garch = rep(NA, ncoef)
            cat("Warning: singular information\n") 
        } else {
            se.garch = sqrt(diag(solve(com.hess$hess)))
        }
        sigt = sqrt(pred$e)
        sigt[1:max(order[1],order[2])] = rep(NA, max(order[1],order[2]))
        f = cbind(sigt,-sigt)
        colnames(f) = c("sigt","-sigt")
        e = as.vector(x)/sigt  
        if (ists) {
            attr(e, "tsp") =  attr(f, "tsp") = xtsp
            attr(e, "class") = attr(f, "class") = "ts" 
        }
        names(order) = c("p", "q")
        coef = fit$coef
        nam.coef = "a0"
        if (order[2] > 0)
            nam.coef = c(nam.coef, paste("a", seq(order[2]), sep = ""))
        if (order[1] > 0)
            nam.coef = c(nam.coef, paste("b", seq(order[1]), sep = ""))
        names(coef) = nam.coef
        names(se.garch) = nam.coef
        
        # Return Value:
        garch = list(order = order, coef = coef, n.likeli = fit$nlikeli,
            n.used = n, residuals = e, fitted.values = f, series = series, 
            frequency = xfreq, call = match.call(), asy.se.coef = se.garch)
        return(garch) 
    }
    
    # Call:
    call = match.call()
    
    # Fit:
    sink("@sink@")
    fit = BIgarch(x = x, order = order, ...)
    sink()
    unlink("@sink@")
    
    # Add to result:
    names(fit$coef)[1] = "omega"
    fit$call = call
    fit$model = paste("GARCH(", 
        as.character(order[1]), ",", as.character(order[2]), ")", sep = "")
    fit$order = order
    fit$x = x
    fit$se.coef = fit$asy.se.coef
    # fit$sigma2 = var(na.remove(fit$residuals))
    fit$sigma2 = var(as.vector(na.omit(fit$residuals)))
    
    # Retrun Value:
    class(fit) = "fGARCH"
    fit
}


# ------------------------------------------------------------------------------


print.fGARCH = 
function(x, ...)
{   # A function implemented by D. Wuertz

    # Description:
    #   Print method for an object of class "fGARCH".
    
    # Notes:
    #   Required:
    #   object$call
    #   object$coef
    
    # FUNCTION:
    
    # Call:
    cat("\nCall:\n")
    cat(paste(deparse(x$call), sep = "\n", collapse = "\n"), 
        "\n\n", sep = "")
    
    # Coefficients:
    cat("Coefficient(s):\n")
    digits = max(4, getOption("digits") - 4) 
    print.default(format(x$coef, digits = digits), print.gap = 2, 
        quote = FALSE)
    cat("\n")
        
    # Return Value
    invisible(x)
}


# ------------------------------------------------------------------------------


summary.fGARCH = 
function (object, ...) 
{   # A function implemented by D. Wuertz

    # Description:
    #   Summary method for an object of class "fGARCH".
    
    # FUNCTION:

    # Initialize:
    ans = NULL
    
    # Fit Call and Model:
    ans$call = object$call
    ans$model = object$model
    
    # Calculate Residuals and Variance:
    # ans$residuals = na.remove(object$residuals)
    ans$residuals = as.vector(na.omit(object$residuals))
    ans$var = var(ans$residuals)
    
    # Generate Coef Matrix:
    tval = object$coef/object$se.coef
    prob = 2 * (1 - pnorm(abs(tval)))
    ans$coef = cbind(object$coef, object$se.coef, tval, prob)
    dimnames(ans$coef) = list(names(object$coef), 
        c(" Estimate", " Std. Error", " t value", "Pr(>|t|)"))
   
    # Fit:
    ans$aic = (object$n.used * (1 + log(2 * pi)) + object$n.used * 
        log(ans$var) + 2 * length(object$coef))
        
    # Return Value:   
    class(ans) = "summary.fGARCH"
    ans
}


# ------------------------------------------------------------------------------


print.summary.fGARCH = 
function(x, ...)
{   # A function implemented by D. Wuertz

    # Description:
    #   Print summary method for an object of class "fGARCH".
    
    # FUNCTION:

    # Call and Model:
    object = x
    cat("\nCall:\n", deparse(object$call), "\n", sep = "")
    cat("\nModel:\n", object$model, "\n", sep = "")
    
    # Residuals:
    cat("\nResiduals:\n")
    rq = structure(quantile(object$residuals), 
        names = c("Min", "1Q", "Median", "3Q", "Max"))
    digits = max(4, getOption("digits") - 4)
    print(rq, digits=digits, ...)
    
    # Coef Matrix:
    cat("\nCoefficient(s):\n")
    signif.stars = getOption("show.signif.stars")
    printCoefmat(object$coef, digits = digits, signif.stars = signif.stars, ...)
    
    # Fit:
    cat("\nFit:\n")
    cat("sigma^2 estimated as:       ", 
        format(object$var, digits = digits), "\n")
    cat("AIC Criterion:              ", 
        format(round(object$aic, 2)), "\n")
    
    # Return Value:
    cat("\n")
    invisible()
}


# ------------------------------------------------------------------------------ 


plot.fGARCH = 
function(x, gof.lag = 10, ...)
{   # A function implemented by D. Wuertz

    # Description:
    #   Plot Method for an object of class "fGARCH".
    
    # FUNCTION:

    # Standardized residuals:
    object = x
    rs = object$residuals
    stdres = rs/sqrt(object$sigma2)
    plot(stdres, type = "h", main = "Standardized Residuals", 
        ylab = "Residuals")
    abline(h = 0)
    
    # Plot ACF:
    acf(abs(object$residuals), plot = TRUE, main = "ACF of |Residuals|", 
        na.action=na.pass)
    
    # QQ Plot of Residuals:
    qqnorm(stdres, xlab = "Normal Quantiles", ylab = "Residual Quantiles", 
        main = "QQ-Plot of Residuals")
    qqline(stdres)
    
    # LB Statistic:
    nlag = gof.lag
    pval = numeric(nlag)
    for (i in 1:nlag) pval[i] = Box.test(abs(rs), i, type = "Ljung-Box")$p.value
    plot(1:nlag, pval, xlab = "lag", ylab = "p value", ylim = c(0, 1), 
        main = "Ljung-Box p-values")
    title(sub = "Test on absolute Values of Residuals")
    abline(h = 0.05, lty = 2, col = "blue")

    # Return Value:
    invisible()
}


# ------------------------------------------------------------------------------


fitted.values.fGARCH = 
function(object, ...)
{   # A function implemented by D. Wuertz

    # Description:
    #   Fitted values method for an object of class "fGARCH".
    
    # FUNCTION:
    
    # Fitted Values:
    ans = as.vector(object$fitted.values)
    
    # Return Value:
    ans
}


# ------------------------------------------------------------------------------


residuals.fGARCH = 
function(object, ...)
{   # A function implemented by D. Wuertz

    # Description:
    #   Residuals method for an object of class "fGARCH".
    
    # FUNCTION:
    
    # Residuals:
    ans = as.vector(object$residuals)
    
    # Return Value:
    ans
}


# ******************************************************************************


fAPARCH = 
function(x, ...)
{   # A function implemented by D. Wuertz

    UseMethod("fAPARCH")
}


# ------------------------------------------------------------------------------


aparchSim = 
function(model = list(omega=1e-6, alpha=0.1, gamma=0, alpha.lags=1,
beta = 0.8, beta.lags = 1, delta = 1), n = 100, innov = rand.gen(n, ...), 
n.start = 100, start.innov = NULL, rand.gen = rnorm, ...)
{   # A function implemented by D. Wuertz   

    # Description:
    #   Allows to simulate time series from the asymmetric
    #   power ARCH family, e.g. ARCH, GARCH, GJR, ...
    #   sigma(t)**delta = omega
    #         + alpha(i) * ||eps(t-i)|-gamma(i)*eps(t-i)|**d
    #         + beta(j) * sigma(t-j)**delta
    
    # Note:
    #   The default model is a symmetric Taylor-Schwert(1,1) 
    #   GARCH model with Gaussian innovations!

    # Arguments:
    #   x
    #   model
    #     omega            Variance coefficient
    #     alpha            Autoregressive coefficients
    #       alpha.lags       integer vector denoting to which lags the 
    #                          alpha coefficients belong
    #       gamma            Asymmetry coefficients, same length as alpha
    #     beta             Moving average coefficients
    #       beta.lags        integer vector denoting to which lags the 
    #                          beta coefficients belong
    #     delta            delta exponent
    #   innov              innovations from a random number generator
    #   start.innov        innovations for starting values
    
    # FUNCTION: 
    
    # Innovations:      
    if (is.null(innov)) innov = rand.gen(n, ...)
    if (is.null(start.innov)) start.innov = rand.gen(n.start, ...)   
    
    # subroutine aparchsim( z, x, h, nt,
    #  &   omega, alpha, gamma, lagsa, na, beta, lagsb, nb, delta)
    #  z   innovations
    #  x   time series,
    #  h   std's
    #  nt  nx vector length
    #  omega, alpha, gamma, alpha.lags, na=length(alpha),
    #         beta, nbeta=beta.lags, nb=length(beta), delta
    # Call Fortran Routine:
    result = .Fortran( "aparchsim",
        as.double(c(start.innov, innov)),
        as.double(rep(0, length(start.innov)+length(innov))),
        as.double(rep(0, length(start.innov)+length(innov))),
        as.integer(length(start.innov)+length(innov)),
        as.double(model$omega),
        as.double(model$alpha),
        as.double(model$gamma),
        as.integer(model$alpha.lags),
        as.integer(length(model$alpha.lags)),
        as.double(model$beta),
        as.integer(model$beta.lags),
        as.integer(length(model$beta.lags)),
        as.double(model$delta),
        PACKAGE = "fSeries")
    
    # Time Series:
    start = length(start.innov)+1
    end = start - 1 + length(innov)
    x = result[[2]][start:end]
    
    # Return Value:
    as.ts(x)
}


# ------------------------------------------------------------------------------


aparchFit = 
function(x, 
order = list(alpha.lags = 1, beta.lags = 1, delta = 2), 
opt = list(gamma = FALSE, delta = FALSE, disparm = FALSE),
distribution = c("norm", "t", "symstb"), disparm = c(1, 4, 1.9), 
n.cond = NULL, doprint = TRUE, method = "Nelder-Mead", ...)
{   # A function written by Diethelm Wuertz

    # Description:
    #   Fits the parameter of an APARCH Model.
    
    # Note:
    #   This program is still under construction !!
    
    # FUNCTION:
        
    # Settings:         
    nt = length(x)             # time series length
    h = 0 * x                  # variances  
    laga = order$alpha.lags    # lag vector alpha
    lagb = order$beta.lags     # lag vector beta
    na = length(laga)          # number of coefficients alpha
    nb = length(lagb)          # number of coefficients beta    
    np = laga[length(laga)]    # highest order      
    nq = lagb[length(lagb)]    # highest order  
    if (is.null(n.cond)) n.cond = np + nq
    delta = order$delta
    
    # Select always the first!
    distribution = distribution[1]
    disparm = disparm[1]
    
    # Starting Variables:
    tlog = function(x) {-log((1-x)/x) }
    par = -log( mean(abs(x))^delta )               # omega            
    par = c(par, tlog(rep(0.1000/na, times=na)))   # alphas
    par = c(par, tlog(rep(0.8000/nb, times=nb)))   # betas                      
    if (opt$gamma) par = c(par, rep(0, times=na))  # gamma  
    if (opt$delta) par = c(par, -log(delta))       # delta
    if (opt$disparm) par = c(par, disparm)         # disparm
        
    # Global Save:
    h <<- h; nt <<- nt
    laga <<- laga; lagb <<- lagb
    na <<- na; nb <<- nb
    delta <<- delta
    disparm <<- disparm
    distribution <<- distribution
    n.cond <<- n.cond
    opt <<- opt
    doprint <<- doprint
    
    # Likelihood Function:
    "fn" = function(par) {      
        # Needs: x, nt, n.cond, laga, lagb, na, nb, gamma, delta
                
        # Parameters:
        omega <<- exp(-par[1])
        alpha <<- 1/(1 + exp(-par[2:(na+1)]))
        beta  <<- 1/(1 + exp(-par[(na+2):(1+na+nb)]))
        if (opt$gamma) { 
            ng = na
            gamma <<- par[(1+na+nb+1):(1+na+nb+ng)]  }
        else { 
            ng = 0
            gamma <<- rep(0, times = na) }
        if (opt$delta) { 
            nd = 1
            delta <<- exp(-par[1+na+nb+ng+nd]) }
        else { 
            nd = 0
            delta <<- delta }
        if (opt$disparm) { 
            disparm <<- par[1+na+nb+ng+nd+1] }
        else { 
            disparm <<- disparm }       
            
        if (TRUE) {
        # Iteration Loop:   
        # EVALUATION OF LOG-LIKELIHOOD FUNCTION FOR APARCH(p,q) MODELS
        #   x(t) = z(t) * h(t)**(1/delta)
        #   h(t) = omega
        #        + A(L) * [(|eps(t)|-gamma(t)*eps(t))**delta]   
        #        + B(L) * h(t)          
        result = .Fortran("sumllh",
            as.double(x),
            as.double(h),
            as.integer(nt),
            as.double(omega),
            as.double(alpha),
            as.double(gamma),
            as.integer(laga),
            as.integer(na),
            as.double(beta),
            as.integer(lagb),
            as.integer(nb),
            as.double(delta),
            as.integer(n.cond),
            PACKAGE = "fSeries")            
        h <<- result[[2]]   
        
        # Log Likelihood:
        hh = abs(h[(n.cond+1):nt])^(1/delta)
        if (distribution == "norm") 
            llh <<- -sum(log(dnorm(x[(n.cond+1):nt]/hh)/hh))
        if (distribution == "t") 
            llh <<- -sum(log(dt(x[(n.cond+1):nt]/hh, df=disparm)/hh))
        if (distribution == "symstb") 
            llh <<- -sum(log(dsymstb(x[(n.cond+1):nt]/hh, alpha=disparm)/hh))}
    
        # Trace:
        if (doprint) print(c(llh, omega, alpha, gamma, beta, delta, disparm))
    
        # Return Result:
        llh }
        
    # Optimization:
    fit = optim(par = par, fn = fn, method = method, ...)
    # print(c(llh, omega, alpha, gamma, beta, delta, disparm))
    
    # Add Results:
    fit$omega = omega
    fit$alpha = alpha
    fit$gamma = gamma
    fit$beta = beta
    fit$delta = delta
    fit$disparm = disparm
    
    fit$call = match.call()
    fit$coef = c(omega, alpha, gamma, beta, delta, disparm)
    
    # Return Value
    class(fit) = "fAPARCH"
    fit    
}


# ------------------------------------------------------------------------------


print.fAPARCH = 
function(x, ...)
{   # A function imlemented by Diethelm Wuertz

    # Description:
    #   Print method for an object of class "fAPARCH".
    
    # Note:
    #   Required:
    #   object$call
    #   object$coef
    
    # FUNCTION:
    
    # Call:
    object = x
    cat("\nCall:\n")
    cat(paste(deparse(object$call), sep = "\n", collapse = "\n"), 
        "\n\n", sep = "")
    
    # Coefficients:
    cat("Coefficient(s):\n")
    digits = max(4, getOption("digits") - 4) 
    print.default(format(object$coef, digits = digits), print.gap = 2, 
        quote = FALSE)
        
    # Return Value
    cat("\n")
    invisible(object)
}


# ------------------------------------------------------------------------------


summary.fAPARCH = 
function(object, ...)
{   # A function imlemented by Diethelm Wuertz

    # Description:
    #   Summary method for an object of class "fAPARCH".
    
    # Note:
    #   Sorry, not yet implemented.

    # FUNCTION:
    
    # Return value:
    print(object)
}


# ------------------------------------------------------------------------------


print.summary.fAPARCH = 
function(x, ...)
{   # A function imlemented by Diethelm Wuertz

    # Description:
    #   Print summary method for an object of class "fAPARCH".
    
    # Note:
    #   Sorry, not yet implemented.

    # FUNCTION:
    
    # Return Value:
    print(x)
}


# ------------------------------------------------------------------------------


# fitted.values.fAPARCH = 
# function(object, ...)
# {   # A function imlemented by Diethelm Wuertz

    # Description:
    #   Fitted values method for an object of class "fAPARCH".
    
    # FUNCTION:

    # Return Value:
#   as.vector(object$fitted.values)
#}


# ------------------------------------------------------------------------------


# residuals.fAPARCH = 
# function(object, ...)
# {   # A function imlemented by Diethelm Wuertz

    # Description:
    #   Residuals method for an object of class "fAPARCH".
    
    # FUNCTION:
    
    # Return Value:
#    as.vector(object$residuals)
#}


################################################################################


