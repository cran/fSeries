
#
# Example: 
#   Write functions for the modelling of GARCH(p,q) processes
#   
# Details:
#   The orders p and q must both be greater or equal 1. 
#
# Description: 
#   This is an exercise which shows how to write functions for the
#   modelling of GARCH(p,q) processes. The functions can be devided
#   into the following classes including
#   PART I:     Model Specification
#   PART II:    Time Series Simulation
#   PART III:   Parameter Estimation
#   PATY IV:    Forecasting Volatility  
# 
# Examples:
#   DEMGBP Benchmark Data
#   Garch(p,q) Specification
#   Garch(p,q) Simulation
#   Garch(1,1) Parameter Estimation - With Mean
#   Garch(1,1) Parameter Estimation - Zero Mean
#   GARCH(1,2) AND GARCH(2,1) Processes
#   Garch(1,1) Forecasting
#
# Author:
#   (C) 2002, Diethelm Wuertz, GPL
#
    
    
################################################################################
# Part I: Model Specification

    
    # Class Representation:
    setClass("garchpqSpec", 
        representation(
            formula = "formula", 
            model = "list", 
            presample = "data.frame",
            distribution = "character",
            h = "numeric" ) 
         ) 
        
    
    # Specification Function:
    garchpqSpec = function(
    model = list(omega = 1.0e-6, alpha = 0.1, beta = 0.8), presample = NULL)
    {
        # Add Missing Coefficients:
        if (!any(names(model) == "omega")) model$omega = 1.0e-6
        if (!any(names(model) == "alpha")) model$alpha = 0.1
        if (!any(names(model) == "beta")) model$beta = 0.8
        
        # Model Order:
        order.alpha = length(model$alpha) 
        order.beta = length(model$beta) 
        
        # Compose Formula Object from Model Argument:
        formula = paste("~ garch(", as.character(order.alpha), ", ",
            as.character(order.beta), ")", sep = "")
    
        # Define Missing Presample:
        order.max = max(order.alpha, order.beta)
        if (is.null(presample)) {
            h = model$omega/(1-sum(model$alpha)-sum(model$beta))            
            presample = data.frame(
                z = rep(0, times = order.max), 
                h = rep(h, times = order.max),
                x = rep(0, times = order.max))
        }
        
        # Return Value:
        new("garchpqSpec",
            formula = as.formula(formula),
            model = model, 
            presample = presample,
            distribution = "dnorm",
            h = numeric() )
    }

    
    # Print Method:
    print.garchpqSpec = function(x, ...)
    {
        # Print Title:
        cat("\nTitle: \n ")
        cat("GARCH(p,q) Modelling\n")
        
        # Print Formula:
        cat("\nFormula: \n ")
        cat(as.character(x@formula))

        # Print Model Parameters:
        cat("\n\nModel Parameters:")
        cat("\n omega:", x@model$omega)
        cat("\n alpha:", x@model$alpha)
        cat("\n beta: ", x@model$beta)

        # Print Presample:
        cat("\n\nPresample: \n")
        print(as.data.frame(x@presample))

        # Print Distribution:
        cat("\nDistribution: \n")
        cat("", x@distribution, "\n\n")

        # Return Value:
        invisible()
    }
    
    
################################################################################
# Part II: Time Series Simulation


    # Compare with:
    
    #   armaSim = function (
    #       model = list(ar = c(0.5, -0.5), d = 0, ma = 0.1), 
    #       n = 100, innov = NULL, n.start = 100, start.innov = NULL, 
    #       rand.gen = rnorm, rseed = NULL, ...) 
    
    #   garchpqSim = function(
    #       model = list(omega = 1.0e-6, alpha = 0.1, beta = 0.8), 
    #       n = 100, innov = NULL, n.start = 100, start.innov = NULL,
    #       presample = NULL, rand.gen = "rnorm", rseed = NULL, ...)  


    # Simulate Artificial GARCH Process:
    garchpqSim =
    function(model = list(omega = 1.0e-6, alpha = 0.1, beta = 0.8), 
    n = 100, innov = NULL, n.start = 100, start.innov = NULL,
    presample = NULL, rand.gen = "rnorm", rseed = NULL, ...) 
    {
        # Specification:
        spec = garchpqSpec(model= model, presample = presample)
       
        # Determine Orders:
        order.alpha = length(spec@model$alpha)
        order.beta = length(spec@model$beta)

        # Create Innovations:
        if (is.integer(rseed)) set.seed(rseed)
        FUN = match.fun(rand.gen)
        if (is.null(start.innov) & is.null(innov))
            z = FUN(n.start + n, ...)
        if (is.null(start.innov) & is.numeric(innov)) 
            z = c(FUN(n.start, ...), innov)
        if (is.numeric(start.innov) & is.null(innov)) 
            z = c(start.innov, FUN(n, ...))
        if (is.numeric(start.innov) & is.numeric(innov)) 
            z = c(start.innov, innov)
            
        # Expand to Whole Sample:
        z = c(spec@presample[, 1], z)
        h = c(spec@presample[, 2]^2, rep(NA, times = n.start + n))
        y = c(spec@presample[, 3], rep(NA, times = n.start + n))

        # Iterate GARCH(p,q) Model:
        m = length(spec@presample[, 1])
        for (i in (m+1):(m+n.start+n)) {
            h[i] = spec@model$omega +
                sum(spec@model$alpha*(y[i-(1:order.alpha)]^2)) +
                sum(spec@model$beta*h[i-(1:order.beta)] )
            y[i] = sqrt(h[i]) * z[i]
        }
        
        # Update Specification Structure:
        spec@presample = 
            data.frame(z = z, h = sqrt(h), x = y)[(n.start+1):(m+n.start), ]
        spec@h = sqrt(h[-(1:(m+n.start))]) 
        spec@distribution = rand.gen

        # Cut Result:
        ans = y[-(1:(m+n.start))]
        attr(ans, "spec") = spec
       
        # Return Value:
        ans      
    }
 

################################################################################
# Part III: Parameter Estimation
    
    # Compare with:
    # setClass("fARMA", 
    #    representation(
    #        call = "call",
    #        formula = "formula",
    #        method = "character",
    #        parameter = "list",
    #        data = "list",
    #        fit = "list",
    #        residuals = "numeric",
    #        fitted.values = "numeric",
    #        title = "character",
    #        description = "character")  
    #)

    
    # Class Representation:
    setClass("fGARCHPQ", 
        representation(
            call = "call",
            formula = "formula",
            method = "character",
            parameter = "list",
            data = "list",
            fit = "list",
            residuals = "numeric",
            fitted.values = "numeric",
            variances = "numeric",
            title = "character",
            description = "character")  
    )
        
    
    # The function garchpqFit consists of three parts. First we 
    # initialize the time series information using the internal
    # function ".garchpqInitSeries" and store the result in 
    # a global variable named ".series"; second we intitialize
    # the model parameters using the internal function 
    # ".garchpqInitParams" and store the result in a global 
    # variable named :\".params", and third we optimize the 
    # Log-Likelihood function using the function ".garchpqOptimizeLLH".

    # As additional functions we need internal functions to perform
    # the transformation of parameters. These are the functions
    # ".transLog" and ".transExp". We also need internal functions
    # to compute the conditional distribution function and the 
    # Log-Likelihood function. These are the functions "garchpqDist"
    # and "garchpqLLH".


    # Fit GARCH(p,q) Model:
    garchpqFit = 
    function(formula.var = ~ garch(1, 1), series = x, presample = NULL, 
        init = NULL, fixed = NULL, trace = TRUE, title = NULL, 
        description = NULL, ...)
    { 
        # Initialize Time Series Information:            
        .series <<- .garchpqInitSeries(formula.var, series,
            presample, trace = trace)
            
        # Initialize Model Parameters:
        .transLog <<- function(z, U, V)  { -log( (V-z)/(z-U) ) } 
        .transExp <<- function(z, U, V)  {  U + (V-U)/(1+exp(-z)) }
        .params <<- .garchpqInitParams(formula.var, init, fixed, trace = trace) 
        
        # Optimize:         
        fit = .garchpqOptimizeLLH(trace = trace, ...)  

        # Add Title and Description:
        if (is.null(title)) title = "GARCH(p,q) Modelling"
        if (is.null(description)) description = as.character(date())
              
        # Return Value:
        new("fGARCHPQ",     
            call = as.call(match.call()),
            formula = as.formula(formula.var), 
            method = "Max Log-Likelihood Estimation", 
            parameter = list(presample = presample),
            data = list(x = series),
            fit = fit,        
            residuals = series - fit$par["mu"],
            fitted.values = fit$par["mu"],
            variances = .series$h[-(1:length(.series$presample[, 1]))],
            title = as.character(title),
            description = as.character(description) 
        )
    }
    
    
    # Initialize Time Series Information:
    .garchpqInitSeries = 
    function(formula.var, series, presample, trace = TRUE)
    {
        # Extract Order from Formula Object:
        order = as.numeric(strsplit(strsplit(strsplit(as.character(
            formula.var), "\\(")[[2]][2], "\\)")[[1]], ",")[[1]])

        # Generate Presample:
        max.order = max(order)
        if (is.null(presample)) {
            pre.z = rep(0, max.order)
            pre.h = rep(var(series), max.order)
            pre.x = rep(mean(series), max.order)
            presample = data.frame(z = pre.z, h = pre.h, x = pre.x)
        }

        # Save Augmented Series:
        z = c(presample[, 1], rep(0, times = length(series)))
        h = c(presample[, 2], rep(var(series), times = length(series)))
        x = c(presample[, 3], series)
      
        # Trace:
        if (trace) {
            cat("\nExtracted Order:           ", order)
            cat("\nLength of Series:          ", length(z))
            cat("\n\nPresample Series: \n")
            print(presample)
            cat("\n")
        }

        # Return Value:
        list(order = c(p = order[1], q = order[2]), z = z, h = h, x = x, 
            presample = presample)
    }
    
    
    # Initialize Model Parameters:
    .garchpqInitParams = 
    function(formula.var, init, fixed, trace = TRUE)
    {
        # Extract Order from Formula Object:
        order = as.numeric(strsplit(strsplit(strsplit(as.character(
            formula.var), "\\(")[[2]][2], "\\)")[[1]], ",")[[1]])
        p = order[1]
        q = order[2]
        
        # Names for Initial Parameters:
        Names = c("mu", "omega",
            if (p > 0) paste("alpha", 1:p, sep = ""),
            if (q > 0) paste("beta", 1:q, sep = ""))      
        
        # Set Limits for Transformations:       
        U = c(-10*abs(mean(.series$x)), 0, 
            if (p > 0) rep(0, times = p),
            if (q > 0) rep(0, times = q))
        V = c(10*abs(mean(.series$x)), 10*var(x), 
            if (p > 0) rep(1, times = p),
            if (q > 0) rep(1, times = q))
        names(U) = names(V) = Names
        
        # Initialize Parameters:
        if (is.null(init)) {
            mu = mean(x)
            alpha.start = 0.1
            beta.start = 0.8
            omega = var(x) * (1 - alpha.start - beta.start)
            params = c(mu, omega,
                if (p > 0) rep(alpha.start/p, p),
                if (q > 0) rep(beta.start/q, q) )
        } else {
            params = init
        }
        names(params) = Names
            
        # Which of the Parameters Should be Fixed?
        if (is.null(fixed)) fixed = rep(FALSE, 1+1+p+q)
        names(fixed) = Names 
        index = (1:length(params))[fixed == FALSE]
        names(index) = names(params)[fixed == FALSE]
      
        # Trace:
        if (trace) {
            cat("\nInitial Parameters:\n")
            print(params)
            cat("\nLimits of Transformations:\n")
            print(rbind(U, V))
            cat("\nWhich Parameters are Fixed?\n")
            print(fixed)
            cat("\nIndex List of Parameters to be Optimized:\n")
            print(index)
            cat("\n")
        }

        # Return Value:
        list(params = params, U = U, V = V, fixed = fixed, index = index)
    }
      
     
    # Define Conditional Density Function:
    .garchpqDist = function(z, hh) 
    {
        # Use the Normal Distribution
        ans = dnorm(x = z/hh, mean = 0, sd = 1) / hh 
        
        # Return Value:
        ans
    }
            
    
    # Define log-Likelihood Function:
    .garchpqLLH = function(par, trace = FALSE, trans = FALSE)
    {       
        # Backtransform Parameters:
        if (trans) {
            .params$params[.params$index] = .transExp(par, 
                .params$U[.params$index], .params$V[.params$index])
        } else {
            .params$params[.params$index] = par
        }  
            
        # Extract Parameters by Name:
        .par = .params$params
        names(.par) = names(.params$params)     
        mu    = .par[substr(names(.par), 1, 2) == "mu"]
        omega = .par[substr(names(.par), 1, 2) == "om"]
        alpha = .par[substr(names(.par), 1, 2) == "al"]
        beta  = .par[substr(names(.par), 1, 2) == "be"]

        # Iterate and Save Conditional Variances:
        p = .series$order[1]
        q = .series$order[2]
        x = .series$x
        h = .series$h
        N = length(x)
        z = x - mu
        M = length(.series$presample[,1]) + max(p,q) + 1
        for ( i in M:N )
            h[i] = omega + sum(alpha*z[i-(1:p)]^2) + sum(beta*h[i-(1:q)])
        .series$h <<- h
        
        # Calculate log-Likelihood:
        hh = sqrt(abs(h[M:N]))
        llh = -sum(log(.garchpqDist(z[M:N], hh)))
        if (trace) {
            cat("\n LLH: ", llh, "\n")
            # print(.params$params[.params$index])
            print(.par)
        }

        # Return Value:
        llh
    }
   
     
    # Optimize Log Likelihood Function:
    .garchpqOptimizeLLH = 
    function(trace = TRUE, ...)
    {
        # Transform Parameters:
        par = .transLog(.params$params[.params$index], 
            .params$U[.params$index], .params$V[.params$index])
            
        # Estimate Parameters:
        if (trace) cat("\nIteration Path:\n")
         fit = optim(par = par, fn = .garchpqLLH, trace = trace, 
            trans = TRUE, ...) 
        # fit = nlm(f = .garchpqLLH, p = par, trace = trace, trans = TRUE, ...) 
        # fit$par = fit$estimate
        # fit$value = fit$minimum
        
        # Backtransform Parameters:
        .params$params[.params$index] <<- .transExp(fit$par, 
            .params$U[.params$index], .params$V[.params$index]) 
        fit$par = .params$params[.params$index]
        names(fit$par) = names(.params$params[.params$index])
        
        # Compute Simple Hessian:
        eps = 0.0001 * fit$par
        n = length(fit$par)
        m = matrix(0, ncol = n, nrow = n)
        f = .garchpqLLH
        for (i in 1:n) {
            for (j in 1:n) {
                x1 = x2 = x3 = x4 = fit$par
                x1[i] = x1[i] + eps[i]
                x1[j] = x1[j] + eps[j]
                x2[i] = x2[i] + eps[i]
                x2[j] = x2[j] - eps[j]
                x3[i] = x3[i] - eps[i]
                x3[j] = x3[j] + eps[j]
                x4[i] = x4[i] - eps[i]
                x4[j] = x4[j] - eps[j]
                m[i, j] = (f(x1) - f(x2) - f(x3) + f(x4))/(4 * eps[i] * eps[j])
            }
        }
        fit$hessian = m 
            
        
        # Final Print Output:
        if (trace) {
            cat("\nHessian:\n")
            print(fit$hessian)
            cat("\nLog-Likelihood: \n", fit$value, "\n") 
            cat("\nFinal Estimate:\n")
            print(fit$par)
        }   
        
        # Return Value:
        fit
    }
    
    
    # S3 Print Method - Copy and adapt print.fARMA:
    print.fGARCHPQ = 
    function(x, ...)
    {           
        # Title:
        cat("\nTitle:\n")
        cat(x@title, "\n")
        
        # Call:
        cat("\nCall:\n")
        cat(paste(deparse(x@call), sep = "\n", collapse = "\n"), 
            "\n", sep = "")
          
        # Model: 
        cat("\nModel:\n", as.character(x@formula), "\n", sep = "")
        
        # Coefficients:
        coef = x@fit$par
        se.coef = sqrt(diag(solve(x@fit$hessian)))
        tval = coef/se.coef
        matcoef = cbind(coef, se.coef, tval, 2*(1-pnorm(abs(tval))))
        dimnames(matcoef) = list(names(tval), c(" Estimate", 
            " Std. Error", " t value", "Pr(>|t|)"))
        signif.stars = getOption("show.signif.stars")
        # digits = max(4, getOption("digits") - 4) 
        digits = 6
        cat("\nCoefficient(s):\n")
        printCoefmat(matcoef, digits = digits, signif.stars = signif.stars)
       
        # Log Likelihood: 
        cat("\nLog-Likelihood:\n", x@fit$value, "\n", sep = "")
        
        # Persistence: 
        persistence = 
            sum(x@fit$par[substr(names(x@fit$par), 1, 5) == "alpha"]) +
            sum(x@fit$par[substr(names(x@fit$par), 1, 4) == "beta"])
        cat("\nPersistence:\n", persistence, "\n", sep = "")
                    
        # Description:
        cat("\nDescription:\n")
        cat(x@description, "\n\n")
            
        # Return Value:
        invisible()
    }
    

################################################################################
# PART IV: Forecasting


    predict.fGARCHPQ = 
    function(object, n.ahead = 10, ...)
    {  
        
        # Extract Parameters by Name:
        .par = object@fit$par   
        omega = .par[substr(names(.par), 1, 2) == "om"]
        alpha = .par[substr(names(.par), 1, 2) == "al"]
        beta  = .par[substr(names(.par), 1, 2) == "be"]
        mu    = .par[substr(names(.par), 1, 2) == "mu"]
 
        # Get Series:
        M = n.ahead
        N = length(object@data$x)
        xx = c(object@data$x, rep(mu, M))
        hh = c(object@variances, rep(0, M))
        zz = xx - mu
        
        # Forecast:
        p = length(alpha)
        q = length(beta)
        for (i in 1:M ) {
            hh[N+i] = omega  + sum(beta*hh[N+i-(1:q)])
            for ( j in 1:p ) {
                if (i-j > 0) {
                    s = hh[N + i - j]
                } else { 
                    s = zz[N + i - j]^2
                }
                hh[N+i] = hh[N+i] + alpha * s
            }
        }
            
        # Result:
        forecast = data.frame(Mean = xx, Variance = hh)[-(1:N), ]
        
        # Return Value:
        forecast
    }
   
   
################################################################################


    # tseries - Wrapper:
    garchTsFit = function(x, order) 
    {
        require(tseries)
        fit = garch(x, order)
        llh = -(fit$n.likeli + length(x)*log(2*pi)/2)
        ans = cbind(Estimate = fit$coef, "Std. Error" = fit$asy.se.coef, 
            "t value" = fit$coef/fit$asy.se.coef)
        rownames(ans) = names(fit$coef)
        round(ans, digits = 6)
    }
    
    
    # Splus - Wrapper:
    garchSplusFit = function(...)
    {
        fit = garch(..., control = bhhh.control(tol = 1.0e-12, n.iter = 10000))
        ans = cbind(fit$coef, sqrt(diag(solve(-fit$cov$A))), 
        fit$coef/sqrt(diag(solve(-fit$cov$A))))
        colIds(ans) = c("Estimate", "Std.Error", "t.value")
        round(ans, digits = 6)
    }
    
    
################################################################################
    

# DEMGBP Benchmark Data:

    # Bollerslev and Ghysels DEMGBP Rates:
    require(fSeries)
    data(dem2gbp)
    x = dem2gbp[, 1]
  

# ------------------------------------------------------------------------------
# Garch(p,q) Specification:
    

    # Default GARCH(1,1) Model:
    garchpqSpec()
    
    # GARCH(2,1) Model - Omega Missing:
    garchpqSpec(model = list(alpha = c(0.15, 0.05), beta = 0.75))
    
    # GARCH(1,1) with User Defined Presample:
    garchpqSpec(model = list(omega = 2e-6, alpha = 0.12, beta = 0.85),
        presample = data.frame(cbind(z = 0.12, h = 0.05, y = -0.02)))
        

# ------------------------------------------------------------------------------
# Garch(p,q) Simulation:


    # Default GARCH(1,1) Model:
    garchpqSim()
    
    # GARCH(1,1) with User Defined Presample:
    garchpqSim(model = list(omega = 2e-6, alpha = 0.12, beta = 0.85),
        presample = data.frame(cbind(z = 0.12, h = 0.05, y = -0.02)))
        
    # GARCH(1,1) with User Defined Presample - no warm up:
    garchpqSim(presample = data.frame(cbind(z = 0.12, h = 0.05, y = -0.02)),
        n.start = 0)


# ------------------------------------------------------------------------------
# Garch(1,1) Parameter Estimation - With Mean:
    

    # Data:
    require(fSeries)
    data(dem2gbp)
    x = dem2gbp[, 1]
   
    # -0.006190  0.01076   0.1531   0.8060
    #  0.008462  0.002852  0.02652  0.03355
        
    # Rmetrics:
    fit1 = garchpqFit(formula.var = ~ garch(1,1))
    fit1
    #              Estimate  Std. Error   t value    
    #   mu        -0.006209    0.008469    -0.733   
    #   omega      0.010760    0.002853     3.772    
    #   alpha1     0.153403    0.026579     5.771    
    #   beta1      0.805885    0.033564    24.010     

    # Ox - R/Interface:
    fit2 = garchOxFit(formula.var = ~ garch(1, 1))
    fit2
    #                 Value  Std. Error   t value
    #   Cst(M)    -0.006183   0.0084616    -0.731
    #   Cst(V)     0.010761   0.0028506     3.775
    #   ARCH(1)    0.153410   0.0265680     5.774
    #   GARCH(1)   0.805880   0.0335420    24.026
    
    # S-Plus:
    x = scan("dem2gbp.csv", skip = 1)
    garchSplusFit(formula.var =  ~ garch(1, 1), series = x)
    #              Estimate  Std. Error   t value 
    #        C    -0.006091   0.0084693    -0.719
    #        A     0.010848   0.0028900     3.754
    #  ARCH(1)     0.153982   0.0267438     5.758
    # GARCH(1)     0.804816   0.0338972    23.743
    

# ------------------------------------------------------------------------------
# Parameter Estimation - Zero Mean:

    
    # Rmetrics:
    fit1 = garchpqFit(formula.var = ~ garch(1,1),  
        init = c(0, 0.01, 0.10, 0.80), fixed = c(T, F, F, F))
    fit1  
    #              Estimate  Std. Error   t value
    #   omega      0.010866    0.002888     3.763   
    #   alpha1     0.154596    0.026784     5.772    
    #   beta1      0.804432    0.033857    23.760   
    
    # Ox - R/Interface:
    fit2 = garchOxFit(formula.var = ~ garch(1, 1), include.mean = FALSE)
    fit2
    #              Estimate  Std. Error   t value
    #   Cst(V)     0.010873    0.002887     3.766
    #   ARCH(1)    0.154640    0.026778     5.775
    #   GARCH(1)   0.804350    0.033847    23.765
    
    # R tseries:
    fit3 = garchTsFit = function(x, order = c(1,1)
    fit3
    #              Estimate  Std. Error   t value
    #   omega      0.010784    0.001288     8.372
    #   alpha      0.154074    0.013823    11.147
    #   beta       0.805295    0.015966    50.439
    
 
    
# ------------------------------------------------------------------------------
# GARCH(1,2) AND GARCH(2,1) Processes:


    # Note: Garch(p,q) 
    #   in Rmetrics / S-Plus  p is ARCH-order
    #   in Ox / tseries       q is ARCH-order !

    # Rmetrics:
    garchpqFit(formula.var = ~ garch(2, 1))
    #              Estimate  Std. Error   t value
    # mu          -0.006279    0.008477    -0.741  
    # omega        0.010786    0.002857     3.775  
    # alpha1       0.153372    0.002655     5.776  
    # alpha2       3.17e-10          NA        NA    
    # beta1        0.805781    0.003356    24.012 
    garchpqFit(formula.var = ~ garch(1, 2))  
    #              Estimate  Std. Error   t value
    # mu          -0.005062    0.008524    -0.594  
    # omega        0.011249    0.002983     3.771  
    # alpha1       0.168611    0.027673     6.093  
    # beta1        0.489901    0.130747     3.747  
    # beta2        0.297293    0.12587      2.362  

    # Ox:
    garchOxFit(formula.var = ~ garch(1, 2), series = dem2gbp[, 1])
    #              Estimate  Std. Error   t value
    # Cst(M)      -0.006218    0.008484    -0.733  
    # Cst(V)       0.010759    0.007692     1.399   
    # ARCH(Alpha1) 0.153439    0.026693     5.748  
    # ARCH(Alpha2) 0.000000    0.083615     0.000   
    # GARCH(Beta1) 0.805868    0.11124      7.245   
    garchOxFit(formula.var = ~ garch(2, 1), series = dem2gbp[, 1])
    #              Estimate  Std. Error   t value
    # Cst(M)      -0.005035    0.008510    -0.592   
    # Cst(V)       0.011247    0.002980     3.774  
    # ARCH(Alpha1) 0.168616    0.027663     6.095  
    # GARCH(Beta1) 0.489921    0.13074      3.747   
    # GARCH(Beta2) 0.297278    0.12586      2.362   
            
    # S-Plus:
    x = scan("dem2gbp.csv", skip = 1)
    garchSplusFit(formula.var = ~ garch(2, 1), series = x)
	#              Estimate  Std. Error   t value
	#        C    -0.002247    0.008514    -0.264
	#        A     0.001875    0.000840     2.231
	#  ARCH(1)     0.221632    0.033738     6.569
	#  ARCH(2)    -0.176475    0.035500    -4.971
	# GARCH(1)     0.947134    0.014626    64.757
    garchSplusFit(formula.var = ~ garch(1, 2), series = x)
	#              Estimate  Std. Error   t value
	#        C    -0.004895    0.008515    -0.575
	#        A     0.011198    0.002976     3.763
	#  ARCH(1)     0.168278    0.027592     6.099
	# GARCH(1)     0.486852    0.130312     3.736
	# GARCH(2)     0.300699    0.125557     2.395
    
    # R tseries:
    garchTsFit(x, order = c(1, 2))
    #              Estimate  Std. Error   t value
    # a0           0.011361    0.001568     7.245
    # a1           0.163420    0.018222     8.968
    # a2           0.000000    0.025379     0.000
    # b1           0.795568    0.022755    34.963
    garchTsFit(x, order = c(2, 1))      
    #              Estimate  Std. Error   t value
    # a0           0.011237    0.001497     7.508
    # a1           0.169283    0.016457    10.286
    # b1           0.484615    0.110617     4.381
    # b2           0.302105    0.100985     2.992
    
    
# ------------------------------------------------------------------------------       
# Garch(1,1) Forecasting:

    
    # Data:
    require(fSeries)
    data(dem2gbp)
    x = dem2gbp[, 1]
    
    # Rmetrics:
    fit = garchpqFit(formula.var = ~ garch(1, 1))   
    predict(fit, 10)
    #            Mean   Variance
    # 1975  -0.006216     0.1470
    # 1976  -0.006216     0.1518
    # 1977  -0.006216     0.1564
    # 1978  -0.006216     0.1607
    # 1979  -0.006216     0.1650
    # 1980  -0.006216     0.1690
    # 1981  -0.006216     0.1729
    # 1982  -0.006216     0.1766
    # 1983  -0.006216     0.1802
    # 1984  -0.006216     0.1836
        
    # Ox - R/Interface:
    fit = garchOxFit(formula.var = ~ garch(1, 1))  
    # Horizon   Mean    Variance
    #   1   -0.006185     0.1470
    #   2   -0.006185     0.1517
    #   3   -0.006185     0.1563
    #   4   -0.006185     0.1607
    #   5   -0.006185     0.1648
    #   6   -0.006185     0.1689
    #   7   -0.006185     0.1727
    #   8   -0.006185     0.1764
    #   9   -0.006185     0.1800
    #  10   -0.006185     0.1834
    
    # Splus:
    fit = garch(formula.var = ~ garch(1, 1), series = x,
    	control = bhhh.control(tol = 1.0e-12, n.iter = 10000))
    predict(fit. 10)
    $series.pred:
	 [1] -0.00609 -0.00609  -0.00609  -0.00609  -0.00609
 	 [6] -0.00609 -0.00609  -0.00609  -0.00609  -0.00609
	$sigma.pred:
	 [1] 0.3836317 0.3898176 0.3956579 0.4011777 0.4063997 
	 [6] 0.4113442 0.4160299 0.4204735 0.4246903 0.4286945
	   
    
# ------------------------------------------------------------------------------   
        
 