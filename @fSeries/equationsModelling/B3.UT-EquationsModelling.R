
## Load Data:

   kmenta = read.table("kmenta.csv", header = TRUE, sep = ";")


## OLS Estimation with 2 Restrictions:  
 
   formulas = list( demand = q ~ p + d, supply = q ~ p + f + a )
   FITOLS = eqnsFit(formulas, data = kmenta)
   FITOLS
   
   print(FITOLS)
   plot(FITOLS)
   summary(FITOLS)
   
   coef(FITOLS)
   fitted(FITOLS)
   residuals(FITOLS)
   vcov(FITOLS)
   
   predict(FITOLS)
   

## SUR Estimation:

   formulas = list( demand = q ~ p + d, supply = q ~ p + f + a )
   FITSUR = eqnsFit(formulas, data = kmenta, method = "SUR")
   FITSUR
   
   coef(FITSUR)
   fitted(FITSUR)
   residuals(FITSUR)
   vcov(FITSUR)
   
   FITSUR@fit$coef
   FITSUR@fit$fitted
   FITSUR@fit$residuals
   FITSUR@fit$vcov
     
   
## Iterated SUR Estimation:

   formulas = list( demand = q ~ p + d, supply = q ~ p + f + a )
   FITITSUR = eqnsFit(formulas, data = kmenta, method = "SUR", maxit = 100)
   FITITSUR
   
   
## 2SLS Estimation:

   formulas = list( demand = q ~ p + d, supply = q ~ p + f + a )
   inst = ~ d + f + a
   FIT2SLS = eqnsFit(formulas, data = kmenta, method = "2SLS", inst = inst)
   FIT2SLS
   
   FIT2SLS@fit$coef
   FIT2SLS@fit$fitted
   FIT2SLS@fit$residuals
   FIT2SLS@fit$vcov
   
   
## 3SLS Estimation:

   formulas = list( demand = q ~ p + d, supply = q ~ p + f + a )
   inst = ~ d + f + a
   FIT3SLS = eqnsFit(formulas, data = kmenta, method = "3SLS", inst = inst)
   FIT3SLS
   
   FIT3SLS@fit$coef
   FIT3SLS@fit$fitted
   FIT3SLS@fit$residuals
   FIT3SLS@fit$vcov
   
   cov2cor(FIT3SLS@fit$vcov)
   
## The formula for calculating the 3SLS estimator: 
   
   # Use "formula3sls", one of "GLS", "IV", "GMM", "Schmidt" or "EViews".
   # Note, the formulas to calculate the 3SLS estimator lead to identical 
   # results if the same instruments are used in all equations. If 
   # different instruments are used in the different equations, only 
   # the GMM-3SLS estimator ("GMM") and the 3SLS estimator proposed by 
   # Schmidt (1990) ("Schmidt") are consistent, whereas "GMM" is efficient 
   # relative to "Schmidt" (see Schmidt, 1990). 

 
   FIT3SLS = eqnsFit(formulas, data = kmenta, method = "3SLS", 
   	 formula3sls = "IV", inst = inst)
   FIT3SLS = eqnsFit(formulas, data = kmenta, method = "3SLS", 
   	 formula3sls = "GMM", inst = inst)
   	 
   	 