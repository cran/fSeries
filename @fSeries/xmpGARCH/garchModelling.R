
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

# Copyright (C) 1998-2003 by Diethelm Wuertz


################################################################################
# FUNCTION:              CLASS REPRESENTATION:  
#  setClass               S4: Class representation for 'garchSpec'
#  garchSpec              S4: Creates a 'garchSpec' object from scratch
#  print.garchSpec        S3: Print method for an object of class 'garchSpec'
# FUNCTION:              DESCRIPTION - SIMULATION:  
#  garchSim               Simulates a GARCH process
#  garch.sim              Simulates a GARCH process
#  aparch.sim             Simulates an APARCH process
# FUNCTION:              DESCRIPTION - PARAMETER ESTIMATION:    
#  setClass               S4: Class representation for 'fGARCH'
#  garchFit               Fits GARCH and APARCH processes
#  .garchGlobal           Internal undocumented function
#  .tgarchCondVariances
#  .aparchCalcCondVariances
# FUNCTION:              DESCRIPTION - METHODS: 
#  print.fGARCH           Print method for an object of class fGARCH
#  plot.fGARCH            Plot method for an object of class fGARCH
#  summary.fGARCH         Summary method for an object of class fGARCH
#  print.summary.fGARCH   Print summary method for an object of class fGARCH
# PROGRAM STRUCTURE:     FOR FUNCTION - garchFit:
#  - CheckSyntax          Checks syntax of input model arguments
#  - SaveGlobally         Saves arguments globally
#  - FitModel             Fits the parameters of the model
#  - ReturnValue          Returns the parameter estimate
# INTERNAL PROGRAMS:     FOR FUNCTION - .garchGlobal:
#  .garchCheckModel       Checks for Model Consistency
#  .garchGetOrder         Gets the Order of te Model
#  .garchInitTimeSeries   Initializes Time Series
#  .garchSetCondDist      Selects Conditional Density Function
#  .garchInitParms        Initializes all Parameters to be Optimized
#  .garchLLH              Computes Log-Likelihood Function
#  .garchOptimizeLLH      Opimizes Log-Likelihood Function
#  .garchCalcHessian      Calculates LLH Hessian
#  .GarchPrintCoef        Prints estimated Coefficients
#  .garchOutSeries        Outputs Series
################################################################################


require(methods)


# ------------------------------------------------------------------------------


setClass("garchSpec", 
    representation(
        call = "call",
        formula = "formula",        
        model = "list",
        distribution = "character",
        presample = "matrix",
        title = "character",
        description = "character")  
)
        
        
# ------------------------------------------------------------------------------


garchSpec =
function (model = list(omega = 1e-6, alpha = 0.1, beta = 0.8), 
distribution = c("norm", "snorm", "ged", "sged", "std", "sstd"), 
presample = NA, title = NULL, description = NULL)
{   # A function implemented by Diethelm Wuertz

    # Description:
    #   Creates a "garchSpec" object from scratch.
    
    # Arguments:
    #   model - a list with the model parameters as entries
    #     ar - a vector of autoregressive coefficients of length m for
    #       the ARMA specification,
    #     ma - a vector of moving average coefficients of length n for
    #       the ARMA specification,
    #     omega - the variance value for GARCH/APARCH specification,
    #     alpha - a vector of autoregressive coefficients of length p
    #       for the GARCH/APARCH specification,
    #     gamma - a vector of leverage coefficients of length p for 
    #       the APARCH specification,
    #     beta - a vector of moving average coefficients of length q 
    #       for the GARCH/APARCH specification,
    #     mu - the mean value for ARMA specification,
    #     delta - the exponent value used in the variance equation.
    #     parm - a vector with the distributional parameters.
    #   distribution - a character string naming the distribution 
    #     function.
    #   presample - either a multivariate "timeSeries", a multivariate 
    #     "ts", a "data.frame" object or a numeric "matrix" with 3 columns 
    #     and at least max(m,n,p,q) rows. The first culumn ...
    #   title - a title string.
    #   description - a description string.
    
    # Details:
        #   formula - a formula object describing the model, e.g. 
    #       ARMA(m,n) + GARCH(p,q). ARMA can be missing or 
    #       specified as AR(m) or MA(n) in the case of pure 
    #       autoregressive or moving average models. GARCH may 
    #       alternatively specified as ARCH(p) or APARCH(p,q).
    #       If formula is set to "NA", the formula is constructed
    #       from the "model" list.
   
    # FUNCTION:
        
    # Formula:  
    formula.mean = ""
    formula.var = "garch"
    
    # Model:
    # Autoregressive Coefficients:
    if (is.null(model$ar)) {
        model$ar = 0
        order.ar = 0 }
    else {
        order.ar = length(model$ar) }   
    
    # Moving Average Coefficients:
    if (is.null(model$ma)) {
        model$ma = 0
        order.ma = 0 }
    else {
        order.ma = length(model$ma) }
    
    # Mean value:
    if(is.null(model$mu)) 
        model$mu = 0  
    
    # Alpha Coefficients:
    if (is.null(model$alpha)) {
        model$alpha = 0
        order.alpha = 0 }
    else {
        order.alpha = length(model$alpha) }    
    
    # Gamma Coefficients:
    if (is.null(model$gamma)) {
        model$gamma = rep(0, times = length(model$alpha)) }
    else {
        formula.var = "aparch" }    
    
    # Beta Coefficients:
    if(is.null(model$beta)) {
        model$beta = 0
        order.beta = 0 }
    else {
        order.beta = length(model$beta) }    
    
    # Delta Exponent:
    if (is.null(model$delta)) {
        model$delta = 2 }
    else {
        formula.var = "aparch" }

    # Distributional Parameters:
    if(is.null(model$parm)) {                      #   nu xi 
        if (distribution[1] == "norm")  model$parm = c(NA, 1)
        if (distribution[1] == "ged")   model$parm = c(2,  1)
        if (distribution[1] == "std")   model$parm = c(4,  1)
        if (distribution[1] == "snorm") model$parm = c(NA, 1)
        if (distribution[1] == "sged")  model$parm = c(2,  1)
        if (distribution[1] == "ssdt")  model$parm = c(4,  1) }
    else { 
        model$parm = model$parm }
    
    # Distribution:
    distribution = distribution[1]
    
    # Presample:
    # First Column:  pre Innovations - z
    # Second Column: pre Sigmas      - h
    # Third Column:  pre Series      - y
    order.max = max(order.ar, order.ma, order.alpha, order.beta)
    iterate = TRUE
    if (!is.matrix(presample)) {
        if (is.na(presample)) {
            iterate = FALSE
            n.start = order.max }
        else {
            n.start = presample }
        z = rnorm(n = n.start)
        # GARCH(p,q)
        h = rep(model$omega/(1-sum(model$alpha)-sum(model$beta)), 
            times = n.start)
        y = rep(model$mu/(1-sum(model$ar)), times = n.start) }
    presample = cbind(z, h, y)
    # Presample Iteration:
    if (iterate) {
        n.iterate = length(z) - order.max
        deltainv = 1/model$delta
        for (i in n.iterate:1) {
            h[i] = model$omega +    
                sum(model$alpha*(abs(abs(y[i+(1:order.alpha)]) - 
                    model$gamma*y[i+(1:order.alpha)])^model$delta)) +
                sum(model$beta*h[i+(1:order.beta)]) 
            y[i] = model$mu  +  
                sum(model$ar*y[i+(1:order.beta)]) +
                sum(model$ma*(h[i+(1:order.beta)]**deltainv)) +
                h[i]^deltainv * z[i] }
        presample = cbind(z, h, y) }
    
    # Complete Formula:
    if (order.ar > 0 && order.ma == 0) 
        formula.mean = paste ("~ ar(", as.character(order.ar), ")", 
            sep = "")
    if (order.ar == 0 && order.ma > 0) 
        formula.mean = paste ("~ ma(", as.character(order.ma), ")", 
            sep = "")
    if (order.ar > 0 && order.ma > 0) 
        formula.mean = paste ("~ arma(", as.character(order.ma), ", ",
        as.character(order.ma), ")", sep = "")
    if (formula.mean == "") {
        formula.mean = "~ " }
    else {
        formula.mean = paste(formula.mean, " + ") }         
    # Variance Formula:
    formula.var = paste(formula.var, "(", as.character(order.alpha), ", ",
        as.character(order.beta), ")", sep = "")    
    # Common Formula:
    formula = paste(formula.mean, formula.var)
    
    # Title:
    if (is.null(title)) 
        title = "ARMA - GARCH/APARCH Specification"
        
    # Description:
    if (is.null(description)) 
        description = paste("Specification as of ", date())
    
    # Return Value:
    new("garchSpec",     
        call = as.call(match.call()),
        formula = as.formula(formula), 
        model = list(ar = model$ar, ma = model$ma, mu = model$mu, 
            omega = model$omega, alpha = model$alpha, gamma = model$gamma, 
            beta = model$beta, delta = model$delta, parm = model$parm), 
        distribution = as.character(distribution), 
        presample = as.matrix(presample),
        title = as.character(title), 
        description = as.character(description) )      
}


# ------------------------------------------------------------------------------


print.garchSpec =
function(x)
{   # A function implemented by Diethelm Wuertz

    # Description:
    #   S3 Print method for an object of class "garchSpec"
    
    # Arguments:
    #   x - an object of class "garchSpec"
    
    # FUNCTION:
    
    # Call:
    cat("\nCall:\n", deparse(x@call), "\n", sep = "")
    
    # Title:
    cat("\n          Title: ", as.character(x@title), "\n")
    
    # Formula:
    cat("\n        Formula: ", as.character(x@formula), "\n")
    
    # Model:
    cat("\n          Model: ")
    if (sum(abs(x@model$ar)) != 0) 
        cat("\n             ar: ", x@model$ar)
    if (sum(abs(x@model$ma)) != 0)    
        cat("\n             ma: ", x@model$ar)
    if (x@model$mu != 0)              
        cat("\n             mu: ", x@model$mu)
    if (x@model$omega != 0)           
        cat("\n          omega: ", x@model$omega)
    if (sum(abs(x@model$alpha)) != 0) 
        cat("\n          alpha: ", x@model$alpha)
    if (sum(abs(x@model$gamma)) != 0) 
        cat("\n          gamma: ", x@model$gamma)
    if (sum(abs(x@model$beta)) != 0)  
        cat("\n           beta: ", x@model$beta)
        cat("\n          delta: ", x@model$delta, "\n")
    
    # Distribution: 
    cat("\n   Distribution: ", x@distribution)  
    cat("\n      Parameter: ", x@model$parm, "\n")
    
    # Presample:
    cat("\n      Presample: \n")
    n = -(length(x@presample[,1])-1)
    time = 0:n
    print(data.frame(cbind(time, x@presample)))
     
    # Description:
    cat("\n    Description:")
    cat("\n   ", as.character(x@description))
    cat("\n\n")
}


# ******************************************************************************


garchSim =
function(n, spec = NULL, model = list(omega = 1e-6, alpha = 0.1, beta = 0.8), 
distribution = c("norm", "snorm", "ged", "sged", "std", "sstd"),
presample = 10, title = NULL, description = NULL))
{   # A function implemented by Diethelm Wuertz

    # Description:
    #   Simulates a garch time series Process
    
    # Arguments:
    #   n - the length of the series to be simulated
    #   spec - the Garch specification, an object of class "garchSpec"
    #   model - a list with the model parameters as entries
    #     ar - a vector of autoregressive coefficients of length m for
    #       the ARMA specification,
    #     ma - a vector of moving average coefficients of length n for
    #       the ARMA specification,
    #     omega - the variance value for GARCH/APARCH specification,
    #     alpha - a vector of autoregressive coefficients of length p
    #       for the GARCH/APARCH specification,
    #     gamma - a vector of leverage coefficients of length p for 
    #       the APARCH specification,
    #     beta - a vector of moving average coefficients of length q 
    #       for the GARCH/APARCH specification,
    #     mu - the mean value for ARMA specification,
    #     delta - the exponent value used in the variance equation.
    #     parm - a vector with the distributional parameters.
    #   distribution - a character string naming the distribution 
    #     function.
    #   presample - either a multivariate "timeSeries", a multivariate 
    #     "ts", a "data.frame" object or a numeric "matrix" with 3 columns 
    #     and at least max(m,n,p,q) rows. The first culumn ...
    #   title - a title string.
    #   description - a description string.
    
    # Details:
    #   If the specification structure is NULL, then the arguments
    #   model, distribution, presample, title, and description are used
    #   to create internally a GARCH specification structure. Otherwise,
    #   these arguments will be ignored and overwritten form the values
    #   defined through the GARCH specification structure. By default a 
    #   GARCH(1,1) model is specified. 
    
    # FUNCTION:
    
    # Create Garch Specification, if required:
    if (is.null(spec)) spec = garchSpec(model = model, distribution = 
        distribution, presample = presample, title = title, description = 
        description)
    
    # Retrive Order:
    order.ar = order.ma = order.alpha = order.gamma = order.beta = 1    
    if (sum(abs(spec@model$ar)) != 0)    
        order.ar = length(spec@model$ar) 
    if (sum(abs(spec@model$ma)) != 0) 
        order.ma = length(spec@model$ma)
    if (sum(abs(spec@model$alpha)) != 0)
        order.alpha = length(spec@model$alpha)
    if (sum(abs(spec@model$gamma)) != 0)
        order.gamma = length(spec@model$gamma)
    if (sum(abs(spec@model$beta)) != 0)
        order.beta = length(spec@model$beta)
    order.ar; order.ma; order.alpha; order.gamma; order.beta    
    
    # Create Innovations:
    if (spec@distribution == "norm")     
        z = rnorm(n)
    if (spec@distribution == "ged")      
        z = rged(n, nu=spec@model$parm[1])
    if (spec@distribution == "st")       
        z = rst(n, nu=spec@model$parm[1])
    if (spec@distribution == "skewnorm") 
        z = rskewnorm(n, xi = spec@model$parm[2])
    if (spec@distribution == "skewged")  
        z = rskewged(n, nu = spec@model$parm[1], xi = spec@model$parm[2])
    if (spec@distribution == "skewst")   
        z = rskewst(n, nu = spec@model$parm[1], xi = spec@model$parm[2])
    
    # Expand to Whole Sample:
    delta = spec@model$delta
    z = c(rev(spec@presample[,1]), z)
    h = c(rev(spec@presample[,2])^delta, rep(NA, times=n))
    y = c(rev(spec@presample[,3]), rep(NA, times=n))
    m = length(spec@presample[,1])
        
    # Iterate Model:
    deltainv = 1/delta
    for (i in m:(n+m)) {
        h[i] =  spec@model$omega +  
            sum(spec@model$alpha*(abs(abs(y[i-(1:order.alpha)]) - 
                spec@model$gamma*y[i-(1:order.alpha)])^spec@model$delta)) +
            sum(spec@model$beta*h[i-(1:order.beta)]) 
        y[i] =  spec@model$mu  +    
            sum(spec@model$ar*y[i-(1:order.beta)]) +
            sum(spec@model$ma*(h[i-(1:order.beta)]**deltainv)) +
            h[i]^deltainv * z[i] }
    
    # Bind Sample:       
    data = cbind(
        z = z[(m):(n+m)], 
        h = h[(m):(n+m)]^deltainv, 
        x = y[(m):(n+m)])
        
    # Return Value: 
    list(data = data, spec = spec)
}


# ******************************************************************************


garch.sim =
function(model = list(omega = 1e-6, alpha = 0.1, beta = 0.8, mu = 0), n = 100, 
innov = NULL, n.start = 100, start.innov = NULL, rand.gen = rnorm, ...)
{   # A function implemented by Diethelm Wuertz   
    
    # Description:
    #   Generates a GARCH(p,q) process. It is possible to generate
    #   ARCH(p) processes as GARCH(p,0) processes.
    
    # Arguments:
    #   model - a list to specify the GARCH(p,q) model with the
    #       following elements:
    #    omega - Constant variance parameter, a numeric value.
    #    alpha - Alpha parameter which iterates the squared time series
    #       values "x^2", a numeric value or a numeric vector of length 
    #       "p".
    #    beta - Beta parameter which iterates the conditional variances 
    #       "h", a numeric value or a numeric vector of length "q".
    #   mu - Constant mean value.
    #   n - Number of time series points to be generated.
    #   innov - The innovations for the generation of the time series 
    #       "x", a vector of length n, If "innov" is set to "NULL", 
    #       the innovations will be generated with the help of the  
    #       random number generator "rand.gen=rnorm".
    #   n.start - Number of start values of innovations for the 
    #       generation of the time series "x".
    #   start.innov - The innovations for the generation of the time
    #       series "x", a vector of length n. If "start.innov" is 
    #       set to "NULL", the innovations will be generated with
    #       the help of the random number generator "rand.gen=rnorm".
    #   rand.gen - a function which is called to generate the innovations. 
    #       Usually, "rand.gen" will be a random number generator. 
    #   ... - Optional Arguments passed to the function "rand.gen".
    
    # Details:
    #   This is a S-Plus like function call.
    
    # Note:
    #   To generate the subclass of ARCH models, use GARCH(p,0).
    
    # FUNCTION:
    
    # If Missing from Model List - Add to List:
    if (!exists("model$alpha")) model$alpha = 0
    if (!exists("model$beta")) model$beta = 0
    if (!exists("model$mu")) model$mu = 0
    
    # Innovations:
    max.order = max(length(model$alpha), length(model$beta))
    if (n.start < max.order)
        stopt("n.start must be greater or equal max(alpha,beta)")   
    if (is.null(start.innov)) start.innov = rand.gen(n.start, ...)  
    if (is.null(innov)) innov = rand.gen(n, ...)
    
    # Setting Start Values and Vectors:
    h = x = z = c(start.innov, innov)
    deltainv = 1/model$delta    
        
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


aparch.sim =
function(model = list(omega = 1e-6, alpha = c(0.1, 0.1), gamma = c(0, 0), 
beta = 0.7, mu = 0, delta = 2), n = 100, innov = NULL, n.start = 100, 
start.innov = NULL, rand.gen = rnorm, ...)
{   # A function implemented by Diethelm Wuertz   

    # Description:
    #   Generates a APARCH time series process. 
    
    # Arguments:
    #   model - a list to specify the GARCH(p,q) model with the
    #       following elements:
    #    omega - Constant variance parameter, a numeric value.
    #    alpha - Alpha parameter which iterates the squared time series 
    #       values "x^2", a numeric value or a numeric vector of length 
    #       "p".
    #    gamma - Assymmetry parameters, a numeric value or numeric
    #       vector of the same length as "alpha".
    #    beta - Beta parameter which iterates the conditional variances 
    #       "h", a numeric value or a numeric vector of length "q".
    #    mu - Constant mean value.
    #    delta - The exponent of the conditional variances, by default
    #       2 which means that variances are iterated, a numeric value.
    #   n - Number of time series points to be generated.
    #   innov - The innovations for the generation of the time series 
    #       "x", a vector of length n, If "innov" is set to "NULL", 
    #       the innovations will be generated with the help of the 
    #       random number generator "rand.gen=rnorm".
    #   n.start - Number of start values of innovations for the 
    #       generation of the time series "x".
    #   start.innov - The innovations for the generation of the time
    #       series "x", a vector of length n. If "start.innov" is set 
    #       to "NULL", the innovations will be generated with the help 
    #       of the random number generator "rand.gen=rnorm".
    #   rand.gen - a function which is called to generate the innovations. 
    #       Usually, "rand.gen" will be a random number generator. 
    #   ... - Optional Arguments passed to the function "rand.gen".
    
    # Value:
    
    # Details:
    #   This is a S-Plus like function call.
    
    # Note:
    #   Although it is possible to generate as sublcasses of this family
    #   ARCH and GARCH processes, we recommend to use in this case the 
    #   function 'garchSim'.
    #   To generate the subclass of ARCH models, use GARCH(p,0).
    
    # FUNCTION:
    
    # Innovations:
    max.order = max(length(model$alpha), length(model$beta))
    if (n.start < max.order)
        stopt("n.start must be greater or equal max(alpha,beta)")       
    if (is.null(start.innov)) start.innov = rand.gen(n.start, ...)  
    if (is.null(innov)) innov = rand.gen(n, ...)
    
    # Setting Start Values and Vectors:
    h = x = z = c(start.innov, innov)
    deltainv = 1/model$delta    
        
    # Initialization:
    for (i in 1:max.order) { 
        h[i] = model$omega/(1-sum(model$alpha)-sum(model$beta))
        x[i] = sqrt(h[i]) * z[i] + model$mu }
        
    # Iteration:
    n.alpha = length(model$alpha)
    n.beta = length(model$beta)
    for (i in (max.order+1):(n.start+n)) {        
        h[i] = model$omega +    
            sum(model$alpha*(abs(abs(x[i-(1:n.alpha)]) - 
                model$gamma*x[i-(1:n.alpha)])^model$delta)) +
            sum(model$beta*h[i-(1:n.beta)]) 
        x[i] = h[i]**deltainv * z[i] + model$mu }
        
    # Return Value:
    as.ts(x[-(1:n.start)])
}


# ******************************************************************************


setClass("fGARCH", 
    representation(
        call = "call",
        formula = "formula",        
        model = "list",
        distribution = "character",
        presample = "matrix",
        title = "character",
        description = "character")  
)


# ------------------------------------------------------------------------------


garchFit =
function(formula.mean = ~arma(0, 0), formula.var = ~garch(1, 1), series = x,
fixed = NA, pre.mean = NA, pre.var = NA, h.start = NA, llh.start = NA, 
init = NA, delta.est = FALSE, delta.par = NA, cond.dist = c("dnorm", "dt"), 
dist.est = FALSE, dist.par = NA, trace = FALSE, algorithm = c("nlm", "optim"), 
title = NULL, description = NULL, ...)
{   # A function implemented by Diethelm Wuertz

    # Description
    #   Fit parameters to a GARCH model.
    
    # Arguments:
    #   formula.mean - ARMA(u,v) specification, not yet implemented
    #   formula.var - GARCH(p,q) specification
    #   series - names time series "x" of length nt
    #   fixed - which of the coefficients of the ARMA, GARCH, 
    #       "delta", and "dist.par" parameters should be fixed?
    #   pre.mean - a vector of length mx (at lest max(p,q, u,v ), 
    #        added in front of x
    #   pre.var - a vector of length mh=mx+h.start the first mx 
    #       elements are in front of h pre.var[mx+1:mx+h.start] 
    #       replaces the first h's
    #   h.start - where the iteration of h starts
    #   llh.start - where the contributions to the log likelihood
    #       function starts
    #   init - optional starting values for the optimization
    #   delta - power exponent of conditional variance [h]
    #   cond.dist - conditional distribution for log-likelihood
    #       a one-parametric distribution
    #   dist.est - should the distributional parameter be estimated?
    #   dist.par - fixed or init value of the distribution
    #   trace - should the optimization be traced?
    #   algorithm - optimization algorithm
    #   ... - attitional arguments passed to the optimization
    #       algorithm
    
    # Value:
        
    # FUNCTION:
    
    # Call:
    # DEBUG <<- TRUE    
    call = match.call()
    
    # Check Formula Syntax:
    m = length(formula.var)
    if (m != 2) stop("Formula misspecified")
    end = regexpr("\\(", as.character(formula.var[m])) - 1
    model <<- substr(as.character(formula.var[m]), 1, end)
    
    # Save Arguments Globally:
    formula.mean <<- formula.mean
    formula <<- formula.var
    series <<- series
    fixed <<- fixed
    pre.mean <<- pre.mean
    pre.var <<- pre.var 
    h.start <<- h.start
    llh.start <<- llh.start   
    init <<- init
    delta.est <<- delta.est
    delta.par <<- delta.par
    cond.dist <<- cond.dist[1]
    dist.est <<- dist.est
    dist.par <<- dist.par
    trace <<- trace
    algorithm <<- algorithm[1]
    asym <<- FALSE
    
    # Fit the Model:
    .garchGlobal(...)
    
    # Create GARCH Specification Structure:
    #   call = "call"
    #   formula = "formula"       
    #   model = "list"
    #   distribution = "character"
    #   presample = "matrix",
    #   title = "character"
    #   description = "character"
    
    # Return Value:
    result = NULL
    result$call = call
    result$algorithm = algorithm
    result$x = x
    result$h = h
    result$res = res
    result$model = tsmodel
    result$cond.dist = cond.dist[1]
    result$arma.order = c(0, 0)
    result$garch.order = c(p, q)
    result$coef = coef
    result$hessian = hessian
    class(result) = "fGARCH"
    result
}


# ------------------------------------------------------------------------------


.garchGlobal =
function(...)
{   # A function written by Diethelm Wuertz

    # Description:
    #   The Global Garch Optimizer
    
    # Arguments:
    #   ... - passed to the optimization algorithms
    
    # Ddetails:
    #   The following steps will be performed:
    #   .garchCheckModel()      Chceck for Model Consistency
    #   .garchGetOrder()        Get the Order of te Model
    #   .garchInitTimeSeries()  Initialize Time Series
    #   .garchSetCondDist()     Select Conditional Density Function
    #   .garchInitParms()       Initialize all Parameters to be Optimized
    #   .garchLLH()             Compute Log-Likelihood Function
    #   .garchOptimizeLLH()     Opimize Log-Likelihood Function
    #   .garchCalcHessian()     Calculate LLH Hessian
    #   .GarchPrintCoef()       Print estimated Coefficients
    #   .garchOutSeries()       Output Series
    
    # FUNCTION:
    
    # Check Variance Model:
    .garchCheckModel <<- function () {
        # ARMA - Check Mean Formula:
          mm = length(formula.mean)
          if (mm != 2) stop("Mean Formula misspecified")
          end = regexpr("\\(", as.character(formula.mean[mm])) - 1
          model.arma <<- substr(as.character(formula.mean[mm]), 1, end)
          if (trace) {cat("\nMean Model:\n"); print(formula.mean) } 
        # GARCH - Check Variance Formula:
          mv = length(formula)
          if (mv != 2) stop("Variance Formula misspecified")
          end = regexpr("\\(", as.character(formula[mv])) - 1
          tsmodel <<- substr(as.character(formula[mv]), 1, end)
          if (tsmodel != model) stop("Invalid Model Specification")
          if (trace) {cat("\nVariance Model:\n"); print(formula) } 
        # Symmetry:
          if(model == "garch") asym <<- FALSE 
          if(model == "tgarch") asym <<- TRUE
          if(model == "aparch") asym <<- TRUE
          if (DEBUG) cat("\n .garchCheckModel - symmetry: ", !asym, "\n")
        # Return Value:
          invisible() }
  
                
    # Get the Model Order:
    .garchGetOrder <<- function () {
        # ARMA - Determine Mean Order:
          order = as.numeric(strsplit(strsplit(strsplit(as.character(
            formula.mean), "\\(")[[2]][2], "\\)")[[1]], ",")[[1]])
          u <<- order[1]
          v <<- 0; if (length(order) == 2) v <<- order[2]
          maxuv <<- max(u, v)
          if (DEBUG) cat("\n .garchGetOrder - ARMA Order : ", u, v, "\n")
          if (DEBUG) cat("\n .garchGetOrder - Max Order  : ", maxuv, "\n")
        # GARCH - Determine Variance Order:
          order = as.numeric(strsplit(strsplit(strsplit(as.character(
            formula), "\\(")[[2]][2], "\\)")[[1]], ",")[[1]])
          p <<- order[1]
          q <<- 0; if (length(order)==2) q <<- order[2]
          maxpq <<- max(p, q)
          if (DEBUG) cat("\n .garchGetOrder - GARCH Order: ", p, q, "\n")
          if (DEBUG) cat("\n .garchGetOrder - Max Order  : ", maxpq, "\n")
        # Return Value:
          invisible() }
    
          
    # Initialize Time Series "x" and "h" and Start Values:
    .garchInitTimeSeries <<- function() {
        # Time Series "x" Properties:
          x <<- series
          nt <<- length(x)
          x.mean <<- mean(x)
          x.var <<- var(x)
        # Put in front pre.mean:   
          if (is.na(pre.mean[1])) pre.mean <<- rep(x.mean, maxpq)
        # Initialize time series vector x:
          x <<- c(pre.mean, x)
          pre <<- length(x) - nt
          nt.max <<- length(x)
        # Initialize Variance Vector "h":
          if (is.na(h.start)) h.start <<- maxpq+1+length(pre.mean)
          else h.start <<- h.start + length(pre.mean)
          h <<- rep(x.var, nt.max)  
          if(!is.na(pre.var[1])) h[1:length(pre.var)] <<- pre.var
          if (DEBUG) cat("\n .garchGlobal - h.start: ", h.start, "\n")
        # llh Start Value, Previous LLH, and Counter:
          if (is.na(llh.start)) llh.start <<- h.start
          else llh.start <<- llh.start + length(pre.mean)
          llh.previous <<- Inf
          counter <<- -1 
        # Pre-Values:
          p.start = max(h.start, llh.start) 
          timeIndex = (1:p.start)-length(pre.mean)
          xSeries <<- x[1:p.start]
          hSeries <<- h[1:p.start]
          preValues <<- cbind(timeIndex, xSeries, hSeries)
          rownames(preValues) <<- rep("", p.start)
          if (trace) {
             cat("\nPre-Values:")
             cat("\n Order:     ", p, q)
             cat("\n preMean's: ", -length(pre.mean))
             cat("\n h.start:   ", h.start-length(pre.mean))
             cat("\n llh.start: ", llh.start-length(pre.mean))
             cat("\n mean(x):   ", x.mean)
             cat("\n var(x):    ", x.var, "\n\n")
             print(preValues) }
        # Return Value:
          invisible() }
    
          
    # Select Conditional Density Function:
    .garchSetCondDist <<- function() {
        # Set Conditional Distribution Parameter:
          if (cond.dist[1] == "dnorm") {
             dist.est <<- FALSE
             if (is.na(dist.par)) dist.par <<- 1 }
          if (cond.dist[1] == "dt") { 
          if (is.na(dist.par)) dist.par <<- 4 }
        # NORMAL DISTRIBUTION:
          if (cond.dist == "dnorm") {
             .garchDist <<- function(z, hh, parm=NA){
                .Internal(dnorm(x=z/hh, mean=0, sd=1, log=FALSE)) / hh } }
        # STUDENT-T DISTRIBUTION:
          if (cond.dist == "dt") {
             .garchDist <<- function(z, hh, parm=4){ 
                hh = hh*sqrt((parm-2)/parm)
                .Internal(dt(x=z/hh, df=parm, log=FALSE)) / hh } }
        # GED DISTRIBUTION:
          if (cond.dist == "ged") {
             .garchDist <<- function(z, hh, parm=NA){ 
                NA } }
        # DOUBLE-EXPONENTIAL DISTRIBUTION:
          if (cond.dist == "double.exp") {
             .garchDist = function(z, hh, parm=NA){ 
                sqrt2 = sqrt(2)
                exp(-sqrt2*abs*(z/hh)) / (hh*sqrt2) } }
        # Print:
          if (trace) cat("\nConditional Distribution:", cond.dist, "\n")
          if (DEBUG) cat("\n .garchSetCondDist - dist.par: ", dist.par, "\n") 
          if (DEBUG) cat("\n .garchSetCondDist - dist.est: ", dist.est, "\n")
        # return Value:
          invisible() }
    
              
    # Initialize all Parameters to be Optimized:
    .garchInitParms <<- function() {
        # The parameters are sorted in the following order:
        # omega - alpha - gamma - beta - mu - ar - ma - delta.par - dist.par
        # 1       q       asym*q  p      1    v    u    1           1
          positions <<- cumsum(c(1, q, asym*q, p, 1, v, u, 1, 1))
        # If the model is symmetric, then the "gamma"'s are set to zero
        #   and the parameters are fixed.
        # The distribution parameter "parm" is set to 1 for the normal
        #   distribution and to 4 or "dist.par" if specified in the
        #   argument list
        
        # Power Exponent:
          if (is.na(delta.par)) delta.par <<- 2
        # Transformations:
          tlog <<- function(z) { -log((1-z)/z) }
          texp <<- function(z) { 1/(1+exp(-z)) }
        # NOTE: "par.init" is the vector of all model parameters including the
        # power exponent "delta.par" and the distribution parameter "dist.par"
        # MODEL: omega - alpha - omega - beta - mu - ar - ma - parm -
          alpha = 0.1
          beta = 0.8
          if (is.na(init[1])) par.init <<- c(
             1-alpha-beta - var(x)*alpha*mean(x)^2, 
             if (q > 0) tlog(rep(alpha/q, q)),
             if (q > 0 && asym) rep(0, q),  
             if (p > 0) tlog(rep(beta/p, p)), 
             1,
             if (v > 0) rep(0, v),  
             if (u > 0) rep(0, u),
             delta.par,
             dist.par )                             
          else par.init <<- c(
             init[1]/x.var, 
             if (q > 0) tlog(init[2:(q+1)]),
             if (q > 0 && asym) init[(q+2):(2*q+1)],
             if (p > 0) tlog(init[(q+2):(q+p+1)]), 
             init[q+p+2]/x.mean,
             if (v > 0) init[(q+p+2+1):(q+p+2+v)],  
             if (u > 0) init[(q+p+2+v+1):(q+p+2+v+u)],  
             delta.par,
             dist.par )  
          names(par.init) <<- c(
             "omega", 
             if (q > 0) paste("alpha", 1:q, sep=""), 
             if (asym)  paste("gamma", 1:q, sep=""), 
             if (p > 0) paste("beta", 1:p, sep=""), 
             "mu",
             if (v > 0) paste("ar", 1:v, sep=""), 
             if (u > 0) paste("ma", 1:u, sep=""),
             "delta.par",
             "dist.par" ) 
        # NOTE: The vector "fix" determines which of the parameters in the
        # will be fixed (unchanged) during the LLH optimization. The length
        # of this vector is the same as the length of the vector "par.init".
        # Setting up the vector "fix" ...
          if (is.na(fixed[1])) {
             fix <<- c(rep(FALSE, 1+q+asym*q+p+1+v+u), !delta.est, !dist.est) }
          else {
             fix <<- c(fixed, delta.est, dist.est) }
        # if (trace) cat("\nFixed Coefficients/Parameters:", fix, "\n")   
        # NOTE: The vector "par" holds all parameters which have to be
        # optimized. This vector has the length of "FALSE" entries in the
        # vector "fix". In addition "par.index" indexes the position of the 
        # paramters to be optimized with respect to the vector of all
        # parameters. 
        # Excluding the fixed values  ...
          par <<- par.init[!fix]
          par.index <<- (1:length(par.init))[fix == FALSE]
          names(par) <<- names(par.init)[!fix]
        # Print:
          if(DEBUG) cat("\n .garchInitParms - par.init: ", par.init, "\n")
          if(DEBUG) cat("\n .garchInitParms - par: ", par, "\n")
          if(DEBUG) cat("\n .garchInitParms - par.index: ", par.index, "\n")
        # Return Value:
          invisible() }
    
              
    # Compute Log-Likelihood Function used by ".garchLLH" - Par:
    .garchLLH <<- function(par) {   
        # Save Globally:
          par <<- par
        # To evaluate the conditional variances we first extract the
        # parameters by name and back-transform the "alpha" and "beta"
        # coefficients. "par.init" is overwritten by the current values.
        # Extracting the parameters by name ...
          par.init[par.index] <<- par
          omega <<- par.init["omega"]*x.var
          alpha <<- texp(par.init[substr(names(par.init), 1, 5) == "alpha"])   
          beta <<-  texp(par.init[substr(names(par.init), 1, 4) == "beta"])
          mu <<- par.init["mu"]*x.mean
          if(v > 0) ar <<- par.init["ar"]
          if(u > 0) ma <<- par.init["ma"]
        # Iterate Conditional Variances h:  
          ARMA <<- mu
          for (i in (h.start):(nt.max)) h[i] <<- omega + 
            sum(alpha*(x[i-(1:q)]-ARMA)^2) + sum(beta*h[i-(1:p)]) 
        # Calculate Log Likelihood:
          dist.par <<- par.init["dist.par"]
          hh = sqrt(abs(h[llh.start:nt.max]))
          llh <<- -sum(log(.garchDist((x[llh.start:nt.max]-ARMA), 
            hh, dist.par))) 
        # Trace:
          counter <<- counter + 1
          if (trace && counter%%10 == 0) .garchPrintCoef()
        # Save previous llh:
          llh.previous <<- llh  
        # Return Value:
          llh }
    
                
    # Opimize Log-Likelihood Function:
    .garchOptimizeLLH <<- function(...) {
        # Initialization:
          cat("\nInitialization:\n")
          coef = par.init
          names(coef) = substr(names(par.init), 1, 2)
          coef["om"] = coef["om"]*x.var
          coef["al"] = texp(coef["al"])
          coef["be"] = texp(coef["be"])
          coef["mu"] = coef["mu"]*x.mean
          names(coef) = names(par.init)
          print(data.frame(coef, fix, par.init))
        # Optimize:
          cat("\nIteration Path:\n") 
          if (algorithm == "nlm") {
             fit <<- nlm(.garchLLH, par, hessian = FALSE, ...)
             par <<- fit$estimate }
          if (algorithm == "optim") {
             fit <<- optim(par, .garchLLH, hessian = FALSE, ...)
             par <<- fit$par }
        # Final Estimate:
          if (trace) {
            cat("\nFinal Estimate:\n")
            .garchPrintCoef()}
        # Print:
          if (DEBUG) cat("\n .garchOptimizeLLH: ", algorithm, "\n")
        # Save Coefficients:
          coef.init <<- par.init
          coef.init[par.index] <<- par
          names(coef.init) <<- substr(names(par.init), 1, 2)
          coef.init["om"] <<- coef.init["om"]*x.var
          coef.init["al"] <<- texp(coef.init["al"])
          coef.init["be"] <<- texp(coef.init["be"])
          coef.init["mu"] <<- coef.init["mu"]*x.mean  
          coef.init["di"] <<- coef.init["di"]
          coef <<- coef.init[fix == FALSE]
        # Return value:
          invisible() }
    
          
    # Compute LLH Hessian used by ".garchCalcHessian" - Coef:
    .garchHessian <<- function(coef) {  
        # Start with fixed Parameters:
          coef.init[par.index] <<- coef
          names(coef.init) <<- substr(names(par.init), 1, 2)
          omega <<- coef.init["om"]
          alpha <<- coef.init["al"]
          beta  <<- coef.init["be"]
          mu    <<- coef.init["mu"]
          dist.par <<- coef.init["di"]
          names(coef) <<- names(par.init)[fix == FALSE]
        # Iterate Conditional Variances h:      
          ARMA <<- mu
          for (i in (h.start):(nt.max)) h[i] <<- omega + 
            sum(alpha*(x[i-(1:q)]-ARMA)^2) + sum(beta*h[i-(1:p)])                       
        # Calculate Log Likelihood:
          hh = sqrt(abs(h[llh.start:nt.max]))
          llh <<- -sum(log(.garchDist((x[llh.start:nt.max]-ARMA), hh, dist.par)))   
        # Trace:
          if (DEBUG) { print(c(llh, coef)) }
        # Return Value:
          llh }

                
    # Calculate Hessian:
    .garchCalcHessian <<- function (...) {
        # Print Header:
          if (trace) cat("\nHessian Matrix:\n")
        # Calculate Hessian:
          if (algorithm == "nlm") {
             fit.final = nlm(.garchHessian, coef, hessian = TRUE, ...)
             hessian <<- fit.final$hessian 
             estimate2 <<- fit.final$estimate }
          if (algorithm == "optim") {
             fit.final = optim(coef, .garchHessian, hessian = TRUE, ...)
             hessian <<- fit.final$hessian 
             estimate2 <<- fit.final$par }
        # Give Names to Hessian and Coefficients:
          names(estimate2) <<- rownames(hessian) <<- colnames(hessian) <<- 
             names(par.init)[fix == FALSE]
          names(coef) <<- names(par.init)[fix == FALSE]
        # Print:
          if (trace) print(hessian)
          if (DEBUG) {cat("\nHessian Estimate:\n"); print(estimate2); cat("\n")} 
        # return Value:
          invisible() }
    
          
    # Some Printing:
    .garchPrintCoef <<- function() {
        # Printing for Iteration Paths:
          coef <<- par.init
          coef[par.index] <<- par
          names(coef) <<- substr(names(par.init), 1, 2)
          coef["om"] <<- coef["om"]*x.var
          coef["al"] <<- texp(coef["al"])
          coef["be"] <<- texp(coef["be"])
          coef["mu"] <<- coef["mu"]*x.mean
          names(coef) <<- names(par.init)     
          estimate  <<- c(llh, coef[fix == FALSE])
          names(estimate) <<- c("llh", names(par.init)[fix == FALSE])  
        # Print:
          print(estimate)
        # Return Value:
          invisible() }     
        
                
    # Output Series:
    .garchOutSeries <<- function() {
        # Prepare Series for Output:
          x <<- x[(length(pre.mean)+1):nt.max]
          h <<- h[(length(pre.mean)+1):nt.max]
          res <<- (x - coef["mu"])/sqrt(h) 
        # Return Value:
          invisible() } 
    
    # Run:
    .garchCheckModel()
    .garchGetOrder()
    .garchInitTimeSeries()
    .garchSetCondDist()
    
    # Optimize:
    .garchInitParms()
    .garchOptimizeLLH()
    .garchCalcHessian()
    
    # Output:
    .garchOutSeries()
        
    # Return Value:
    cat("\n")
    invisible()
}


# ******************************************************************************


.tgarchCondVariances <<- 
function() 
{
    # FUNCTION:
    
    S <<- 0.5*(1 - sign(x-mu))
    if(DEBUG) cat("\n .garchCondVariances - h[i]: summing up ...\n")
    for (i in (maxpq+1):(nt+maxpq)) h[i] <<- omega + 
        sum( (alpha + gamma*S[i-(1:q)]) * (x[i-(1:q)]-mu)^2 ) + 
        sum( beta * h[i-(1:p)] )  
}       

            
# ------------------------------------------------------------------------------


.aparchCalcCondVariances <<-
function() 
{
    # FUNCTION:
    
    hd = h^(delta/2)  
    # Note, h is sigma^2, but we have to iterate sigma^delta!
    if(DEBUG) cat("\n .garchCondVariances - h[i]: summing up ...\n")
    for (i in (maxpq+1):(nt+maxpq)) hd[i] = omega + 
        sum(alpha*(abs(x[i-(1:q)]-mu)-gamma*(x[i-(1:q)]-mu) )^delta) + 
        sum(beta*hd[i-(1:p)]) 
    h <<- hd^(2/delta) 
}


# ******************************************************************************
                

print.fGARCH = 
function(object) 
{   # A function implemented by Diethelm Wuertz
    
    # Description:
    #   Print method for an object of class "fGARCH"
    
    # FUNCTION:
    
    # Check object:
    # if (!inherits(object, "fGARCH")) 
    #    stop("method is only for fGARCH objects")
    
    # Call:
    cat("\nCall:\n")
    cat(paste(deparse(object$call), sep = "\n", collapse = "\n"), 
        "\n\n", sep = "")
    
    # Mean Equation:
    cat("\nMean Equation: ~ arma(", object$arma.order[1], ", ",
        object$arma.order[2], ")\n", sep = "")
    
    # Conditional Variance Equation: 
    if (object$garch.order[1] == 0)
        cat("\nConditional Variance Equation: ~ ", object$model, "(", 
            object$garch.order[2], ")\n", sep = "")
    else 
        cat("\nConditional Variance Equation: ~ ", object$model, "(", 
            object$garch.order[1], ", ", object$garch.order[2], ")\n", 
            sep = "")
        
    # Conditional Distribution:
    cat("\nConditional Distribution: ", object$cond.dist, "\n", sep = "")
  
    # Coefficients:
    cat("\nCoefficient(s):\n")
    digits = max(4, getOption("digits") - 4)
    print.default(format(object$coef, digits = digits), print.gap = 2, 
        quote = FALSE)    
   
    # Return Value:
    cat("\n")
    invisible(object)
}


# ------------------------------------------------------------------------------


plot.fGARCH = 
function(object) 
{   # A function implemented by Diethelm Wuertz
    
    # Description:
    #   Plot method for an object of class "fGARCH"

    # FUNCTION:
    
    # Check Object:
    if (!inherits(object, "fGARCH")) 
        stop("method is only for fGARCH objects")
        
    # Plot Time Series"
    plot(object$x, type="l", main="Time Series")
    
    # Conditional Variances:
    plot(object$h, type="l", main="Conditional Variances")
    
    # Autocorrelation Functions: 
    acf(object$res)
    acf(object$res^2)  
}


# ------------------------------------------------------------------------------


summary.fGARCH = 
function(object)
{   # A function implemented by Diethelm Wuertz
    
    # Description:
    #   Summary method for an object of class "fGARCH"

    # Requirements:
    #   Requires R's tseries Package!
    
    # FUNCTION:
    
    # Check object:
    # if (!inherits(object, "fGARCH")) 
    #    stop("method is only for fGARCH objects")
    ans = object
    
    # Residuals: 
    ans$cvar = solve(ans$hessian)
    ans$se.coef = sqrt(diag(ans$cvar))
    ans$tval = ans$coef/ans$se.coef
    ans$matcoef = cbind(ans$coef, ans$se.coef, 
        ans$tval, 2*(1-pnorm(abs(ans$tval))))
    dimnames(ans$matcoef) = list(names(ans$tval), c(" Estimate", 
        " Std. Error", " t value", "Pr(>|t|)"))
    
    # Tests:
    ans$jbtest = jarque.bera.test(ans$res)
    ans$lbtest = Box.test(ans$res^2, type = "Ljung-Box")
    
    # Return Value:
    class(ans) = "summary.fGARCH"
    return(ans)
}


# ------------------------------------------------------------------------------


print.summary.fGARCH = 
function (object) 
{   # A function implemented by Diethelm Wuertz
    
    # Description:
    #   Print summary method for an object of class "fGARCH"
    
    # FUNCTION:
    
    # Check Object:
    if (!inherits(object, "summary.fGARCH")) 
        stop("method is only for summary.garch objects")
    
    # Header:
    digits = max(4, getOption("digits") - 4)
    print.fGARCH(object)
    
    # Residuals:
    cat("Residuals:\n")
    rq = structure(quantile(object$res), names = c("Min", "1Q", 
        "Median", "3Q", "Max"))
    print(rq, digits = digits)
    
    # Coefficients:
    signif.stars = getOption("show.signif.stars")
    cat("\nCoefficient(s):\n")
    print.coefmat(object$matcoef, digits = digits, 
        signif.stars = signif.stars)
    
    # 1. Diagnostic Test:
    cat("\nJarque Bera Test of Residuals:\n")
    out1 = paste(
        names(object$jbtest$statistic), " = ", 
        format(round(object$jbtest$statistic, 4)),
        ", ", sep="")
    out2 = paste(
        names(object$jbtest$parameter), " = ", 
        format(round(object$jbtest$parameter, 3)),
        ", ", sep="")
    out3 = paste(
        "p-value =", 
        format.pval(object$jbtest$p.value, digits = 4))
    cat(out1, out2, out3, "\n")
    
    # 2. Diagnostic Test:
    cat("\nLjung-Box Test of Squared Residuals:\n")
    out1 = paste(names(object$lbtest$statistic), " = ", 
        format(round(object$lbtest$statistic, 4)), ", ", sep="")
    out2 = paste(names(object$lbtest$parameter), " = ", 
        format(round(object$lbtest$parameter, 3)), ", ", sep="")
    out3 = paste("p-value =", format.pval(object$lbtest$p.value, digits = 4))
    cat(out1, out2, out3, "\n")
        
    # Return Value:
    invisible(x)
}


# ------------------------------------------------------------------------------