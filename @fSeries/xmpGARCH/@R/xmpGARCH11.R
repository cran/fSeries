
## SETTINGS:
    
   # Load Libraries:
   library(fSeries)
   require(tseries)
   DEBUG = FALSE

   # Load Benchmark Data:
   data(dem2gbp)
   x <<- dem2gbp[,1]

   # Fit Benchmark Model:
   f1 = garchFit(formula.var = ~ garch(1, 1), series = x, cond.dist = "dt",
        dist.est = TRUE, dist.par = 4, algorithm = "optim", trace = TRUE)
    
   # % MATLAB CODE
   # load garchdata;
   # [coeff, errors] = garchfit(garchset('P',1,'Q',1), price2ret(DEM2GBP));
   # garchdisp(coeff, errors);
    

## BENCHMARK RESULTS:

   # COEFFICIENTS
   # FROM:
   # ------------------ mu ------- omega ----- alpha ------- beta ------
   # FCP          -0.00619       0.01076      0.1531       0.8060
   # G@RCH-OX     -0.00618       0.01076      0.1534       0.8059
   # Eviews       -0.00541       0.00958      0.1423       0.8213
   # PcGive       -0.00625       0.01076      0.1534       0.8059
   # TSP          -0.00619       0.01076      0.1531       0.8060
   # S-Plus       -0.00919       0.01170      0.1543       0.8003  
   # -------------------------------------------------------------------   
    
   # COEFFICIENTS AND T-VALUES
   # FROM:
   # ------------------ mu ------- omega ----- alpha ------- beta ------
   # E-VIEWS      -0.00540  -0.64 0.0096  8.01 0.143  11.09 0.821  53.83
   # GAUSS-FANPAC -0.00600  -0.75 0.0110  3.67 0.153   5.67 0.806  23.71
   # LIMDEP       -0.00619  -0.71 0.0108  3.44 0.153   5.61 0.806  26.73
   # MATLAB       -0.00619  -0.73 0.0108  8.13 0.153  10.96 0.806  48.67
   # MICROFIT     -0.00621  -0.73 0.0108  3.78 0.153   5.78 0.806  24.02
   # SAS          -0.00619  -0.74 0.0108  8.15 0.153  10.97 0.806  48.60
   # SHAZAM       -0.00613  -0.73 0.0107  5.58 0.154   7.91 0.806  36.96
   # RATS         -0.00625  -0.71 0.0108  3.76 0.153   5.79 0.806  23.93
   # TSP          -0.00619  -0.67 0.0108  1.66 0.153   2.86 0.806  11.11
   # -------------------------------------------------------------------
    

## MODIFY INITIAL CONDITIONS:
    
   f1 = garchFit(formula.var = ~ garch(1, 1), h.start = 2, llh.start = 2, 
        pre.mean = 0, pre.var = rep(var(x), 2), trace = TRUE) # DEFAULT   
   f2 = garchFit(formula.var = ~ garch(1, 1), h.start = 2, llh.start = 1, 
        pre.mean = 0, pre.var = rep(var(x), 2), trace = TRUE)
   f3 = garchFit(formula.var = ~ garch(1, 1), h.start = 1, llh.start = 1, 
        pre.mean = mean(x), pre.var = var(x), trace = TRUE)    
   f4 = garchFit(formula.var = ~ garch(1, 1), h.start = 1, llh.start = 1, 
        pre.mean = 0, pre.var = var(x), trace = TRUE)
         
   f11 = garchFit(formula.var = ~ garch(1, 1), series = x, h.start = 2, 
         llh.start = 2, pre.mean = 0, pre.var = rep(var(x), 2), trace = TRUE) 
         
   # Pre-Value Scheme:
   #   m0 x1 x2 x3 ...
   #   v0 v1 h2 h3 ...
   #         h.start
   #         llh.start
     
   #                  omega    alpha1      beta1           mu
   # --------------------------------------------------------
   # BENCHMARK:    0.010761   0.15313   0.805974   -0.0061904
   # f1:           0.010759   0.15340   0.805896   -0.0062092
   # f2:           0.010760   0.15340   0.805886   -0.0061665
   # f3:           0.010741   0.15356   0.805857   -0.0061326
   # f4:           0.010741   0.15356   0.805855   -0.0061278
   # --------------------------------------------------------

    
## USE DIFFERENT OPTIMIZATION ALGORITHMS:

   g1 = f1 # "nlm"
   g2 = garchFit(x ~ garch(1, 1), algorithm = "optim", trace = TRUE)
   g3 = garchFit(x ~ garch(1, 1), algorithm = "optim", method = "BFGS", 
        trace = TRUE)
    
   #                  omega    alpha1      beta1           mu
   # --------------------------------------------------------
   # BENCHMARK:    0.010761   0.15313   0.805974   -0.0061904
   # g1:           0.010759   0.15340   0.805896   -0.0062092
   # g2:           0.010756   0.15353   0.805813   -0.0062256
   # g3:           0.010778   0.15353   0.805687   -0.0062053
   # --------------------------------------------------------
    

## INCREASE TOLERANCE:

   h1 = f1 # "nlm"
   h2 = garchFit(x ~ garch(1, 1), trace = TRUE, ndigit = 12, 
        gradtol = 1e-9, steptol = 1e-9, iterlim = 1000)
   h3 = garchFit(x ~ garch(1, 1), trace=TRUE, ndigit=18,
        gradtol = 1e-12, steptol = 1e-12, iterlim = 1000)

   #                  omega    alpha1      beta1           mu  
   # --------------------------------------------------------
   # BENCHMARK:    0.010761   0.15313   0.805974   -0.0061904  
   # h1:           0.010759   0.15340   0.805896   -0.0062092    
   # h2:           0.010760   0.15340   0.805884   -0.0062090
   # h3:           0.010760   0.15340   0.805884   -0.0062090
   # --------------------------------------------------------
    
   # 0.01076003  0.15340619  0.80588173 -0.00616632 

        
## COMPARE WITH:

   k1 = f1 # "nlm" 
   k2 = garchFit(x ~ garch(1, 1), trace = TRUE, 
        init=c(0.1*var(x), 0.1, 0.8, -0.0061904), 
        fixed = c(F, F, F, T))  # fixed mean to bechmark mean 
   k3 = garch(x + 0.0061904, c(1, 1)) # from "tseries"
   k4 = garchOxFit(formula.var = ~ garch(1, 1), series = x) # from "Ox"

   #                  omega    alpha1      beta1           mu
   # --------------------------------------------------------  
   # BENCHMARK:    0.010761   0.15313   0.805974   -0.0061904
   # k1:           0.010759   0.15340   0.805896   -0.0062092
   # k2:           0.010760   0.15341   0.805880   -0.0061904
   # k3:           0.010720   0.15320   0.806256        FIXED
   # k4:           0.010770   0.15338   0.805848   -0.006144
   # --------------------------------------------------------

    
## ROUND TIME SERIES X TO 6 DIGITS:

   m1 = f1 # "nlm" 
   x  = round(x, digits=6)
   m2 = garchFit(x ~ garch(1, 1), trace=TRUE)
    
   #                  omega    alpha1      beta1           mu  
   # --------------------------------------------------------
   # BENCHMARK:    0.010761   0.15313   0.805974   -0.0061904
   # m1:           0.010759   0.15340   0.805896   -0.0062092
   # m2:           0.010739   0.15355   0.805877   -0.0061328
   # --------------------------------------------------------
    
    
## SUMMARY:

   # f1: Coefficient(s):
   #             Estimate    Std.Error   t.value  Pr(>|t|)    
   # ----------------------------------------------------------
   # omega       0.010759    0.003081      3.492   0.000479 ***
   # alpha1      0.153396    0.027775      5.523   3.34e-08 ***
   # beta1       0.805896    0.035767     22.532    < 2e-16 ***
   # mu         -0.006209    0.008472     -0.733   0.463612   
   # ----------------------------------------------------------
        
   # k3: Coefficient(s):
   #             Estimate   Std.Error    t.value   Pr(>|t|)    
   # ----------------------------------------------------------
   # a0          0.010720    0.001283      8.355    <2e-16 ***
   # a1          0.153196    0.013787     11.111    <2e-16 ***
   # b1          0.806256    0.015961     50.513    <2e-16 ***
   # ----------------------------------------------------------
        
   # k4: Coefficient(s):
   #                Value   Std.Error    t.value
   # -------------------------------------------
   # Cst(M)     -0.006144   0.0084372    -0.7282
   # Cst(V)      0.010770   0.0013244     8.1318
   # ARCH(1)     0.153380   0.0140020    10.9540
   # GARCH(1)    0.805850   0.0165750    48.6190
   # -------------------------------------------
    
   # T-VALUES         mu   omega  alpha1   beta1
   # -------------------------------------------
   # f1:           -0.73    3.49    5.52    22.5
   # k4:           -0.73    8.13   10.95    48.6
   # k3:           FIXED    8.36   11.11    50.5
   # -------------------------------------------
   # E-VIEWS       -0.64    8.01   11.09    53.8
   # GAUSS-FANPAC  -0.75    3.67    5.67    23.7
   # LIMDEP        -0.71    3.44    5.61    26.7
   # MATLAB        -0.73    8.13   10.96    48.7
   # MICROFIT      -0.73    3.78    5.78    24.0
   # SAS           -0.74    8.15   10.97    48.6
   # SHAZAM        -0.73    5.58    7.91    37.0
   # RATS          -0.71    3.76    5.79    23.9
   # TSP           -0.67    1.66    2.86    11.1
    

# ******************************************************************************    
        
##   ARCH(1):
    fit11 <- garchFit(x ~ arch(1))
    fit12 <- garch(x-fit1$coef[3], order=c(0, 1), grad="analytical")
    fit13 <- garch(x-fit1$coef[3], order=c(0, 1), grad="numerical")
    fit14 <- garchOxFit(formula.var = ~garch(1, 1), series=x)

    
# ------------------------------------------------------------------------------


#   GARCH(1,1):
    
    
    
#   MORE GARCH(1,1) ...
    # GARCH(1,1) - with initial values:
    garchFit(x ~ garch(1, 1), init=c(0.0108, 0.153, 0.806, -0.00621))
    # GARCH(1,1) - mu fixed to mean(x):
    garchFit(x ~ garch(1, 1), fixed=c(F, F, F, T))
    # GARCH(1,1) - mu fixed to benchmark value -0.00619041:
    garchFit(x ~ garch(1, 1),
        init=c(var(x)/10, 0.1, 0.8, -0.00619041), fixed=c(F, F, F, T))

        
#   GARCH(2,1):
    garchFit(x ~ garch(2, 1))

    
# ------------------------------------------------------------------------------

    
#   GJR-ARCH(1,1)
    
    garchFit(x ~ gjrarch(1,1), doprint=TRUE, trace=TRUE,
        init=c(0.011231, 0.140777, 0.028321, 0.801357, -0.007906) )
    #    omega     alpha1     gamma1      beta1         mu  
    # 0.011232   0.140773   0.028358   0.801348  -0.007936  
    
    garchOxFit(formula.mean = ~ arma(0, 0), formula.var = ~ gjr(1, 1), 
        series = x, trace = TRUE)
    #    omega     alpha1     gamma1      beta1         mu  
    # 0.011231   0.140777   0.028321   0.801357  -0.007906


# ------------------------------------------------------------------------------


#   APARCH(1,1):
    garchFit(x ~ aparch(1,1), doprint=TRUE, trace=TRUE)
    #   omega    alpha1    gamma1     beta1        mu     delta  
    # 0.02318   0.17483   0.09533   0.79700  -0.00941   1.35396  
    #     llh -1102.616
    
    garchFit(x ~ aparch(1,1), doprint=TRUE, trace=TRUE,
        ndigit=16, gradtol=1e-12, steptol=1e-12, iterlim=1000)
    #   omega    alpha1    gamma1     beta1        mu     delta  
    # 0.02318   0.17484   0.09534   0.79699  -0.00941   1.35363  
    #     llh -1102.615
    
    garchFit(x ~ aparch(1,1), doprint=TRUE, trace=TRUE, 
        init = c(0.0243, 0.173, 0.0953, 0.797, -0.00961, 1.29) )        
    #   omega    alpha1    gamma1     beta1        mu     delta  
    # 0.02318   0.17484   0.09533   0.79700  -0.00941   1.35396  
    #     llh -1102.616

        
    garchOxFit(formula.mean = ~ arma(0, 0), formula.var = ~ aparch(1, 1), 
        series = x, trace = TRUE)
    
    Cst(M)              -0.00961     -0.00941
    Cst(V)               0.0243       0.0232 
    ARCH(Alpha1)         0.173        0.175  
    GARCH(Beta1)         0.800        0.797
    APARCH(Gamma1)       0.100        0.0953
    APARCH(Delta)        1.29         1.35
                        -1101.849    -1102.616
                        
                        
                        
           llh      omega     alpha1     gamma1      beta1         mu      delta 
   -1102.98474    0.02430    0.17300    0.09530    0.79700   -0.00961    1.29000

   
# ------------------------------------------------------------------------------

"garchOxFit" <-
#     function(formula.mean = ~ arma(0, 0), formula.var = ~ garch(1, 1), 
#     series = x, cond.dist = c("gaussian", "t", "ged", "skewed-t"), 
#     arfima=FALSE, arch.in.mean=0, truncation=100, trace = TRUE)




    f1 <- garchFit(formula.var = ~ garch(1, 1), h.start=2, llh.start=2, 
        series=dem2gbp[,1], cond.dist="dnorm", dist.est=TRUE, trace=TRUE)
        
        
    f2 <- garchFit(formula.var = ~ garch(1, 1), h.start=2, llh.start=2, 
        series=dem2gbp[,1], cond.dist="dnorm", dist.est=TRUE, trace=TRUE)
        
        
        
optim <<-
function (par, fn, gr = NULL, method = c("Nelder-Mead", "BFGS", 
    "CG", "L-BFGS-B", "SANN"), lower = -Inf, upper = Inf, control = list(), 
    hessian = FALSE, ...) 
{
    fn1 <- function(par) fn(par, ...)
    gr1 <- if (!is.null(gr)) 
        function(par) gr(par, ...)
    method <- match.arg(method)
    if ((length(lower) > 1 || length(upper) > 1 || lower[1] != 
        -Inf || upper[1] != Inf) && method != "L-BFGS-B") {
        warning("bounds can only be used with method L-BFGS-B")
        method <- "L-BFGS-B"
    }
    con <- list(trace = 0, fnscale = 1, parscale = rep(1, length(par)), 
        ndeps = rep(0.001, length(par)), maxit = 100, abstol = -Inf, 
        reltol = sqrt(.Machine$double.eps), alpha = 1, beta = 0.5, 
        gamma = 2, REPORT = 10, type = 1, lmm = 5, factr = 1e+07, 
        pgtol = 0, tmax = 10, temp = 10)
    if (method == "Nelder-Mead") 
        con$maxit <- 500
    if (method == "SANN") 
        con$maxit <- 10000
    con[(namc <- names(control))] <- control
    if (con$trace < 0) 
        warning("read the documentation for `trace' more carefully")
    if (method == "L-BFGS-B" && any(!is.na(match(c("reltol", 
        "abstol"), namc)))) 
        warning("Method L-BFGS-B uses `factr' (& `pgtol') instead of `reltol' and `abstol'")
    npar <- length(par)
    if (npar == 1 && method == "Nelder-Mead") 
        warning("one-diml optimization by Nelder-Mead is unreliable: use optimize")
    lower <- as.double(rep(lower, , npar))
    upper <- as.double(rep(upper, , npar))
    res <- .Internal(optim(par, fn1, gr1, method, con, lower, 
        upper))
    names(res) <- c("par", "value", "counts", "convergence", 
        "message")
    nm <- names(par)
    if (!is.null(nm)) 
        names(res$par) <- nm
    names(res$counts) <- c("function", "gradient")
  
    
    if (hessian) {
    print(res$par)
    print(fn1)
    print(gr1)
    print(con)
    
        hess <- .Internal(optimhess(res$par, fn1, gr1, con))
        hess <- 0.5 * (hess + t(hess))
        if (!is.null(nm)) 
            dimnames(hess) <- list(nm, nm)
        res$hessian <- hess
    }
    res
}


"hess"<-
function(f, x)
{   # A function by Stuart Coles
    ep <- 0.0001
    eps <- ep * x
    n <- length(x)
    m <- matrix(0, ncol = n, nrow = n)
    for(i in 1:n) {
        for(j in 1:n) {
            x1 <- x
            x1[i] <- x1[i] + eps[i]
            x1[j] <- x1[j] + eps[j]
            x2 <- x
            x2[i] <- x2[i] + eps[i]
            x2[j] <- x2[j] - eps[j]
            x3 <- x
            x3[i] <- x3[i] - eps[i]
            x3[j] <- x3[j] + eps[j]
            x4 <- x
            x4[i] <- x4[i] - eps[i]
            x4[j] <- x4[j] - eps[j]
            m[i, j] <- (f(x1) - f(x2) - f(x3) + f(x4))/(4 * eps[i] * 
                eps[j])
        }
    }
    solve(m)
}