
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
# for the Rmetrics-port: 
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
# FUNCTION:
#  BIpolymars
#  predictBIpolymars
#  summaryBIpolymars
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
 

BIpolymars = 
function(responses, predictors, maxsize, gcv = 4., additive = FALSE, 
startmodel, weights, no.interact, knots, knot.space = 3, ts.resp, ts.pred, 
ts.weights, classify, factors, tolerance = 1e-05, verbose = FALSE)
{	# A slightly modified copy from R's contributed package polspline

	# Arguments:
	#  responses   - a vector (or matrix) of responses. (Can be a a vector of 
	#                characters for classification)
	#  predictors  - a matrix of predictors with same number of cases as 
	#                response. Columns are predictors.
	
	# Optional Arguments:
	#  maxsize     - maximum number of basis function the model can contain 
	#  gcv         - parameter for overall best model seletion
	#  additive    - boolean, is the model to be additive
	#  startmodel  - either a matrix (m*4 or m*5) or a polymars object from 
	#                a previous call to polymars 
	#                an initial model the procedure should start with in model 
	#                selection
	#  weights     - a vector of length equal to the number of cases
	#  no.interact - a 2*l matrix of columns numbers of the predictor 
	#                matrix (each row pair cannot have interaction terms)
	#  knots       - a vector specifying many knots per predictor are 
	#                wanted (with -1 for categorical variables) 
	#                ncol(predictors)==length(knots), or a matrix with 
	#                ncol(predictors) == ncol(knots) with actual knot 
	#                specified and filled out with NA's.
	#                Can also be a single number - "knots" number of knots 
	#                per predictor                    
	#  knot.space  - minimum number of order statistics between knots
	#  ts.resp     - testset reponses, same format as responses
	#  ts.pred     - testset predictors, same format as predictors
	#  ts.weights  - testset weights, same format as weights
	#  classify    - whether classification is to be done, set = TRUE if the 
	#                response vector is integer, if 
	#                if character classify is automatically true
	#  factors     - a vector of column numbers of the predictor matrix of 
	#                categorical variables
	#  tolerance   - a numerical parameter which may need to be made smaller 
	#                if the program crashes store the call to the polymars 
	#                function
   
   
	call = match.call()
	ism0 = missing(classify)
	ism1 = missing(ts.resp)
	ism2 = missing(maxsize)
	ism3 = missing(ts.pred)
	ism4 = missing(ts.weights)
	ism5 = missing(knots)
	ism6 = missing(factors)
	ism7 = missing(startmodel)
	ism8 = missing(weights)
	ism9 = missing(no.interact)
	
	# Internal added by DW:
	unstrip = function(x) {
       dd = dim(x)
       y = x
       if (length(dd)==2){
          dd2 = dd[2]
          if (dd2==1) y<- c(x[,1])
          if (dd2==2) y<- cbind(c(x[,1]),c(x[,2]))
          if (dd2>2) y<- cbind(c(x[,1]),c(x[,2]),c(x[,3]))
          if (dd2>3)for (i in 4:dd2) y = cbind(y,c(x[,i]))
          y }
       if (length(dd)==1 || length(dd)==0){
          y = c(unlist(c(unlist(x))))
          names(y) = NULL }
       y}

	if (!missing(responses))
	  responses = unstrip(responses)
	if (!missing(predictors))
	  predictors = unstrip(predictors)
	if (!missing(weights))
	  weights = unstrip(weights)
	if (!missing(no.interact))
	  no.interact = unstrip(no.interact)
	if (!missing(knots))
	  knots = unstrip(knots)
	if (!missing(ts.resp))
	  ts.resp = unstrip(ts.resp)
	if (!missing(ts.pred))
	  ts.pred = unstrip(ts.pred)
	if (!missing(ts.weights))
	  ts.weights = unstrip(ts.weights)
	if (!missing(factors))
	  factors = unstrip(factors)
	
	responses = as.matrix(responses)
	predictors = data.matrix(predictors)
	nresponses = ncol(responses)
	npredictors = ncol(predictors)
	ncases = nrow(predictors)
	if (ism0) classify = FALSE
	if (mode(responses) == "character" || classify == TRUE) {
      if (ncol(responses) > 1) {
        stop(paste(
            "When using character responses or classify = TRUE", 
            "only 1 response per case is allowed\n") ) }
      char.responses = responses
      int.responses = as.integer(as.factor(responses))
      nresponses = length(unique(responses))
      responses = matrix(ncol = nresponses, nrow = ncases, data = 
         int.responses)
      for (i in 1:nresponses) {
         responses[, i] = (responses[, i] == (unique(int.responses)[i])) }
      conversion = matrix(ncol = 2, nrow = nresponses, c(unique(
         char.responses), unique(int.responses)))
      classify = TRUE
      if (!ism1) {
         char.responses.test = ts.resp
         ts.resp = matrix(ncol = nresponses, nrow = length(
            char.responses.test), data = 0)
         for (i in 1:nresponses) {
            ts.resp[, i] = as.integer(char.responses.test ==
               conversion[i, 1]) } } }
   else {
      conversion = FALSE
      classify = FALSE }
      
	# Maxsize that the model can grow to	
	if (ism2) maxsize = ceiling(min(6 * (ncases^(1/3)), ncases/4, 100))
	
	# If a testset is to be used in model selection
	if (!ism1 || !ism3) {
      if (ism1 || ism3) {
      stop(paste(
        "Both ts.resp (testsets responses) and", 
        "ts.pred (testset predictors) should be specified\n") )}
      if (!is.matrix(ts.resp))
         ts.resp = as.matrix(ts.resp)
      if (!is.matrix(ts.pred))
         ts.pred = as.matrix(ts.pred)
      if (ncol(ts.resp) != nresponses) {
         stop(paste(
            "Testset should have the same number of responses",
            " as the training set\n ")) }
      if (ncol(ts.pred) != npredictors) {
         stop(paste(
            "Testset should have the same number of predictors",
            "  as the training set\n " )) }
      if (nrow(ts.resp) != nrow(ts.pred)) {
         stop(paste(
            "Testset ts.pred and ts.resp should have the same",
            " number of cases (rows)") ) }
      testsetmatrix = cbind(ts.resp, ts.pred)
      testsetcases = nrow(testsetmatrix)
      testset = TRUE
      if (!ism4) {
         if (length(ts.weights) != testsetcases) {
            stop("length of testset weights misspecified\n") }
         testset.weighted = TRUE
         testsetmatrix = cbind(ts.resp * ts.weights, ts.pred) }
      else {
         testset.weighted = FALSE
         ts.weights = 0 } }
   else {
      testsetmatrix = 0
      testsetcases = 0
      testset = FALSE
      testset.weighted = FALSE
      ts.weights = 0 }
      
    # If the mesh is specified by the knots arguement this will be 
    # changed to true later
    mesh.specified = FALSE
    mesh.vector = 0
    if (nrow(responses) != nrow(predictors)) {
      stop(paste(
        "The number of rows (cases) of the response and",
        " predictor matricies should be the same" )) }
    if (!ism5 && !is.matrix(knots) && length(knots) != npredictors && length(
      knots) != 1) {
      stop(paste(
        "Length of vector of `knots per predictor' should",
        " be equal to number of predictors or 1\n" )) }
    if (!ism5) {
      if (!is.matrix(knots)) {
         # if knots is specified as a single number it is expanded  
         # to a vector length npredictors
         if (length(knots) == 1) {
            knots = rep(knots, npredictors)
            if (!ism6) {
               for (i in 1:length(factors)) {
                  if (!is.vector(factors)) {
                     stop(paste(
                     "`factors' should be a vector whose elements",
                     " are indicies of predictors that are factors\n" )) }
                     
                  # in knots the number of knots(per predictor) is specified
                  # or -1 if the predictor is a factor and all it values are 
                  # levels  
                  
                  knots[factors[i]] = -1 } } } }
    else {
         mesh = knots
         mesh.vector = 
            vector(length = ncol(mesh) * nrow(mesh), mode = "double")
         knots = vector(length = npredictors, mode = "integer")
         k = 0
         for (i in 1:npredictors) {
            knots[i] = length(unique(mesh[is.na(mesh[
               , i]) == FALSE, i]))
            for (j in 1:knots[i]) {
               k = k + 1
               mesh.vector[k] = unique(mesh[!is.na(
                  mesh[, i]), i])[j] } }
         if (!ism6) {
            for (i in 1:length(factors)) {
               if (!is.vector(factors)) {
                  stop(paste(
                    "`factors' should be a vector whose elements are",
                    " indicies of predictors that are factors\n" )) }
                    
               # in knots the number of knots(per predictor) is specified
               # or -1 if the predictor is a factor and all it values are 
               # levels
                 
               knots[factors[i]] = -1 } }
         mesh.specified = TRUE } }
   
    if (ism5) {
      knots = rep(min(20, round(ncases/4)), npredictors)
      if (!ism6) {
         for (i in 1:length(factors)) {
            if (!is.vector(factors)) {
               stop(paste(
                "`factors' should be a vector whose elements",
                " are indicies of predictors that are factors\n" )) }
                
            # in knots the number of knots(per predictor) is specified
            # or -1 if the predictor is a factor and all it values are 
            # levels 
            
            knots[factors[i]] = -1 } } }
            
    startmodelsize = 1
    
    # A starting model must be specified as a object of class 
    # polymars or a matrix with 4 or 5 columns
    
    no.remove = 0
    no.remove.size = 0
    if (!ism7) {
      if (is.vector(startmodel))
         startmodel = t(as.matrix(startmodel))
      v1 = (class(startmodel) == "polymars")
      if (length(v1) == 0)
         v1 = FALSE
      if (!(is.matrix(startmodel) || v1) || (is.matrix(startmodel) &&
         (ncol(startmodel) != 4 && (ncol(startmodel) != 5)))) {
         stop(paste(
            "startmodel should be a matrix with each row corresponding to",
            "a function with number of columns = 4 (or 5 for extra boolean\n",
            "column specifying predictors which cannot be removed)",
            "or startmodel should be a polymars object\n")) }
            
      if (is.matrix(startmodel)) {
         # Fifth column denotes which basis functions must remain in 
         # the model at all times
         if (ncol(startmodel) == 5) {
            no.remove = vector(length = (nrow(startmodel)) )
            j = 0
            for (i in 1:nrow(startmodel)) {
               if (startmodel[i, 5] == TRUE) {
                  j = j + 1
                  no.remove[j] = i } }
            no.remove.size = j }
            
         # The startknots are taken from the startmodel and put into a vector
         # The startmodel becomes a 4*n matrix with a "1" in the 2nd and 4th 
         # columns where knots appear
         
         startknots = as.vector(t(cbind(startmodel[, 2], startmodel[, 4])))
         startknots[is.na(startknots)] = 0.
         startmodel = matrix(startmodel[, 1:4], ncol = 4)
         startmodel[!is.na(startmodel[, 2]), 2] = 1
         startmodel[is.na(startmodel[, 2]), 2] = 0
         startmodel[is.na(startmodel[, 3]), 3] = 0
         startmodel[startmodel[, 3] == 0, 4] = 0
         for (i in 1:nrow(startmodel)) {
            if ((!is.na(startmodel[i, 4])) && startmodel[
               i, 3] != 0)
               startmodel[i, 4] = 1 }
         startmodel[is.na(startmodel[, 4]), 4] = 0
         startmodelsize = nrow(startmodel) + 1 }
      else {
         startmodelsize = startmodel$model.size
         startmodel = startmodel$model[-1,  ]
         startknots1 = startmodel$knot1
         startknots2 = startmodel$knot2
         L1 = FALSE
         if (!is.null(startmodel$level1)) {
            L1 = TRUE
            level1 = startmodel$level1 }
         if (L1) {
            startmodel$knot1[!is.na(level1)] = 1
            startknots1[!is.na(level1)] = level1[!is.na(
               level1)] }
         startknots = cbind(startknots1, startknots2)
         startknots = as.vector(t(startknots))
         startknots[is.na(startknots)] = 0.
         startmodel = cbind(startmodel[, "pred1"], startmodel[
            , "knot1"], startmodel[, "pred2"], startmodel[
            , "knot2"])
         startmodel[, 2] = !is.na(startmodel[, 2])
         startmodel[, 4] = !is.na(startmodel[, 4]) } }
    else {
        startmodel = 0
        startknots = 0 }
    if (!ism8) {
        if (length(weights) != ncases) {
            stop("Number of weights not equal to the numnber of cases\n" ) }
        weighted = TRUE
        responses = responses * weights }
    else {
        weighted = FALSE
        weights = 0 }
        
    datamatrix = cbind(responses, predictors)
   
    # Predictors which cannot interact together in the model are  
    # specified by a 2*n matrix of predictor indicies
   
    if (!ism9) {
        if (!is.matrix(no.interact) || ncol(no.interact) != 2) {
            stop(paste(
                "list of interactions disallowed has been misspecified",
                " must be a 2*n matrix") ) }
            no.interact = t(no.interact)
            no.interact.size = ncol(no.interact) }
    else {
        no.interact.size = 0
        no.interact = 0 }
    if (startmodelsize > maxsize) {
        stop(paste(
            "start model should not be of greater size than",
            " the max model size\n" )) }
   
    # Some error checking on the startmodel
    
    if (startmodelsize != 1) {
      for (i in 1:(startmodelsize - 1)) {
         if (startmodel[i, 1] == 0) {
            stop("first column of startmodel cannot be zero\n" ) }
         if (startmodel[i, 2] == 1) {
            if (startknots[(i * 2) - 1] < min(predictors[
               , startmodel[i, 1]]) || startknots[
               (i * 2) - 1] > max(predictors[, 
               startmodel[i, 1]])) {
               stop("Knot out of range of its predictor \n" ) } }
         if (startmodel[i, 4] == 1) {
            if (startknots[(i * 2)] <= min(predictors[,
               startmodel[i, 3]]) || startknots[(
               i * 2)] >= max(predictors[, startmodel[
               i, 3]])) {
               stop("Knot out of range of its predictor\n" ) } } }
      if (max(startmodel[, c(1, 3)] > npredictors)) {
         stop("Initial model misspecified on input\n")  } }
   
    startmodel = t(startmodel)
    resultmodelsize = 0
    end.state = 0
    step.count = 0
   
    z = .C("polymars", as.integer(npredictors),
        as.integer(nresponses), as.integer(ncases),
        as.double(datamatrix), as.integer(knots),
        as.double(mesh.vector), as.integer(mesh.specified),
        as.integer(maxsize), as.double(gcv),
        as.integer(additive), as.integer(startmodelsize),
        start.model = as.integer(startmodel),
        start.knots = as.double(startknots),
        as.integer(weighted), as.double(weights),
        as.integer(no.interact.size), as.integer(no.interact),
        as.integer(no.remove.size), as.integer(no.remove),
        as.integer(knot.space), as.integer(testset),
        as.double(testsetmatrix), as.integer(testsetcases),
        as.integer(testset.weighted), as.double(ts.weights),
        as.integer(classify), as.double(tolerance),
        as.integer(verbose),
        best.model = as.integer(matrix(nrow = maxsize, ncol = 4, data =
            rep(0, maxsize * 4))),
        coefficients = as.double(matrix(nrow = maxsize, ncol = 
            nresponses, data = rep(0., maxsize * nresponses))),
        steps = as.integer(matrix(nrow = maxsize * 2, ncol = 2,
        data = rep(0, maxsize * 4))),
        rss.gcv = as.double(matrix(nrow = maxsize * 2, ncol = 
            nresponses + 1, data = rep(0., maxsize * 2 * (
            nresponses + 1)))),
        modelsize = as.integer(resultmodelsize),
        modelknots = as.double(matrix(nrow = maxsize, ncol = 2, data = 
            rep(0., maxsize * 2))),
        coefficient.se.term = as.double(rep(0., maxsize)),
        end.state = as.integer(end.state),
        step.count = as.integer(step.count),
        PACKAGE = "fSeries")
   
    # The C function returns information about how it ended
   
    if (z$end.state != 0 && z$end.state != 5) {
      switch(z$end.state,
         stop("Mis-specification of initial model\n"),
         stop(paste(
            "Initial model with non-linear function must contain",
            " the corresponding linear function\n" )),
         stop(paste(
            "Initial model contains two-predictor functions that", 
            " require prerequisite functions\n" ))) }
    else {
      model = matrix(z$best.model[1:((z$modelsize - 1) * 4)], ncol
          = 4, byrow = TRUE)
      knot.values = matrix(z$modelknots[1:((z$modelsize - 1) * 2)],
         ncol = 2, byrow = TRUE)
      for (i in 1:nrow(model)) {
         if (model[i, 2] != 0) {
            model[i, 2] = knot.values[i, 1] }
         else {
            model[i, 2] = NA }
         if (model[i, 4] != 0) {
            model[i, 4] = knot.values[i, 2] }
         else {
            model[i, 4] = NA } }
      if (length(knots[model[, 1]]) != 0 && min(knots[model[, 1]]) <  0) {
         factor1 = TRUE
         levels1 = rep(NA, z$modelsize - 1)
         factor.variables = unique(model[knots[model[, 1]] < 0, 1])
         for (i in 1:length(factor.variables)) {
            for (j in 1:length(model[, 1])) {
               if (model[j, 1] == factor.variables[i]) {
                  levels1[j] = model[j, 2] } }
            model[model[, 1] == factor.variables[i], 2] = NA }
         levels1 = c(NA, levels1) }
      else {
         factor1 = FALSE }
      coefs = matrix(z$coefficients[1:(z$modelsize * nresponses)],
         ncol = nresponses)
      
	  # The model that the C-function returns does not explicitly  
	  # contain an intercept so in formatting the output one is added
        
      if (z$modelsize > 1) {
         if (factor1 == FALSE) {
            model = rbind(c(0, NA, 0, NA), model)
            model = data.frame(model, coefs)
            if (nresponses == 1) {
               dimnames(model) = list(1:z$modelsize,
                  c("pred1", "knot1", "pred2",
                  "knot2", "coefs")) }
            else {
               dimnames(model) = list(1:z$modelsize,
                  c("pred1", "knot1", "pred2",
                  "knot2", paste("Coefs", 1:
                  nresponses))) } }
         if (factor1 == TRUE) {
            model[(knots[model[, 1]] < 0), 2] = NA
            model = rbind(c(0, NA, 0, NA), model)
            model = data.frame(model[, 1:2], levels1,
               model[, 3:4], coefs)
            if (nresponses == 1) {
               dimnames(model) = list(1:z$modelsize,
                  c("pred1", "knot1", "level1",
                  "pred2", "knot2", "coefs")) }
            else {
               dimnames(model) = list(1:z$modelsize,
                  c("pred1", "knot1", "level1",
                  "pred2", "knot2", paste("Coefs",
                  1:nresponses))) } } }
      else {
         model = data.frame(0, NA, 0, NA, coefs)
         if (nresponses == 1) {
            dimnames(model) = list(1:z$modelsize, c(
               "pred1", "knot1", "pred2", "knot2",
               "coefs")) }
         else {
            dimnames(model) = list(1:z$modelsize, c(
               "pred1", "knot1", "pred2", "knot2",
               paste("Coefs", 1:nresponses))) } }
               
      # For later plotting the ranges and medians of the 
      # predictors are stored
      
      ranges.and.medians = matrix(ncol = npredictors, nrow = 3, data = 0)
      for (i in 1:npredictors) {
         ranges.and.medians[1, i] = min(predictors[, i]) }
      for (i in 1:npredictors) {
         ranges.and.medians[2, i] = max(predictors[, i]) }
      for (i in 1:npredictors) {
         ranges.and.medians[3, i] = median(predictors[, i]) }
         
      # A table with information from the fitting is formatted here
      
      steps = matrix(z$steps[1:(2 * (z$step.count + 1))], ncol = 2,
         byrow = TRUE)
      rss.gcv = matrix(z$rss.gcv[1:((nresponses + 1) * (z$step.count +
         1))], ncol = nresponses + 1, byrow = TRUE)
      fitting = data.frame(steps, rss.gcv)
      if (testset == FALSE) {
         if (nresponses == 1) {
            dimnames(fitting) = list(1:(nrow(fitting)),
               c("0/1", "size", "RSS", "GCV")) }
         else {
            dimnames(fitting) = list(1:nrow(fitting),
               c("0/1", "size", paste("RSS", 1:nresponses), "GCV")) } }
      else {
         if (classify == FALSE) {
            if (nresponses == 1) {
               dimnames(fitting) = list(1:(nrow(
                  fitting)), c("0/1", "size", "RSS", "T.S. RSS")) }
            else {
               dimnames(fitting) = list(1:nrow(
                  fitting), c("0/1", "size",
                  paste("RSS", 1:nresponses), "T.S. RSS")) } }
         else {
            if (nresponses == 1) {
               dimnames(fitting) = list(1:(nrow(
                  fitting)), c("0/1", "size", "RSS", "T.S.M.C.")) }
            else {
               dimnames(fitting) = list(1:nrow(
                  fitting), c("0/1", "size",
                  paste("RSS", 1:nresponses),
                  "T.S.M.C.")) } } }
                  
      # Calculates fitted values and residual of the data according 
      # to the model returned 
      
      if (z$modelsize > 1) {
         temp = list(model = model, model.size = z$modelsize,
            ranges.and.medians = ranges.and.medians, 
            responses = nresponses)
         class(temp) = "polymars"
         fitted = matrix(ncol = nresponses, nrow = ncases,
            data = rep(0, nresponses * ncases))
         residuals = matrix(ncol = nresponses, nrow = ncases,
            data = rep(0, nresponses * ncases))
         fitted = predictBIpolymars(temp, x = predictors)
         residuals = responses - fitted }
      else {
         fitted = matrix(ncol = nresponses, nrow = ncases,
            data = coefs[1, 1])
         residuals = matrix(ncol = nresponses, nrow = ncases,
            data = responses - coefs[1, 1]) }
            
      # If their are factors present in the model the factors  
      # must be stored for use during plotting
      
      if (factor1 == TRUE) {
         model2 = model[-1,  ]
         factors.in.model = unique(model2[knots[model2[, 1]] <
            0, 1])
         maxfactors = 0
         for (i in 1:length(factors.in.model)) {
            maxfactors = max(maxfactors, length(unique(
               predictors[, factors.in.model[i]]))) }
         factor.matrix = matrix(ncol = length(factors.in.model),
            nrow = maxfactors + 2, data = NA)
         for (i in 1:length(factors.in.model)) {
            factor.matrix[1, i] = factors.in.model[i]
            factor.matrix[2, i] = length(unique(predictors[
               , factors.in.model[i]]))
            for (j in 3:(length(unique(predictors[, 
               factors.in.model[i]])) + 2)) {
               factor.matrix[j, i] = unique(
                  predictors[, factors.in.model[
                  i]])[j - 2] } } }
      else {
         factor.matrix = 0 }
      if (nresponses == 1) {
         SE = round(sqrt((sum(residuals^2)/(ncases - z$
            modelsize)) * z$coefficient.se.term[1:z$
            modelsize]), 4)
         model = cbind(model, SE)
         dimnames(model)[[2]][length(dimnames(model)[[2]])] <-
            "SE" }
      else {
         for (i in 1:nresponses) {
            SE = round(sqrt((sum(residuals[, i]^2)/(ncases -
               z$modelsize)) * z$coefficient.se.term[
               1:z$modelsize]), 4)
            model = cbind(model, SE)
            dimnames(model)[[2]][length(dimnames(model)[[
               2]])] = paste("SE", i) } }
      if (nresponses == 1) {
         Rsquared = 1 - (sum(residuals^2)/sum((responses - mean(
            responses))^2)) }
      else {
         Rsquared = NULL }
      result = list(model = model, fitting = fitting, model.size = z$
         modelsize, fitted = fitted, responses = nresponses,
         residuals = residuals, ranges.and.medians = 
         ranges.and.medians, call = call, conversion = 
         conversion, factor.matrix = factor.matrix, Rsquared = 
         Rsquared)
      class(result) = "polyMARS"
      return(result)
   }
}


# ------------------------------------------------------------------------------


predictBIpolymars = 
function(object,x,classify=FALSE,intercept,...)
{	# A slightly modified copy from R's contributed package polspline

    # Description:
    #   Produces predicted values for a polymars object
    
    # Arguments:
    #   pmars.model  
    #           an object returned from a call to polymars
    #   x       
    #           a matrix with number of columns equal to number of columns 
    #           of predictor matrix in 
    #           original call to polymars and predictor values in the 
    #           corresponding columns. Can 
    #           also be a matrix with number of column equal to the number 
    #           of predictors in the 
    #           model, in the order of the original dataset.
    #   classify     
    #           If the original call to polymars was for classification 
    #           setting  classify=TRUE will 
    #           the new data otherwise it will return the multi-response 
    #           fitted values.
    #   intercept    
    #           By default TRUE. The full intercept is included; or when 
    #           FALSE the intercept is left out.
    #           Can also be givebn a numerical value
    
    # Internal Function added by DW:
    
    unstrip = function(x) {
		dd = dim(x)
        y = x
        if (length(dd)==2){
            dd2 = dd[2]
            if (dd2==1) y<- c(x[,1])
            if (dd2==2) y<- cbind(c(x[,1]),c(x[,2]))
            if (dd2>2) y<- cbind(c(x[,1]),c(x[,2]),c(x[,3]))
            if (dd2>3)for (i in 4:dd2) y = cbind(y,c(x[,i]))
            y }
        if (length(dd)==1 || length(dd)==0){
            y = c(unlist(c(unlist(x))))
            names(y) = NULL }
        y }
 
    if (missing(intercept)) { intercept<-TRUE }
    if (!missing(x))x = unstrip(x)
    
    # Some error checking:
    
    if (class(object)!="polymars") stop("object is not a polymars object")
    pmars.model = object
    
    # The x matrix number of columns can be of length equal to the 
    # number of predictors in the original model or shorten to the 
    # number of predictors in the model in `pmars.model'
    
    if (!(is.matrix(x)))  {
    if (length(unique(pmars.model$model[, "pred1"]))== 1 ||  
        ncol(pmars.model$ranges.and.medians)== 1  ) {
        x<-matrix(data=x,ncol=1) } }
    if ((is.matrix(x) && ncol(x) 
        != length(unique(pmars.model$model[,"pred1"])))) {
    if (ncol(x) != ncol(pmars.model$ranges.and.medians)) {    
        stop(paste(
            "Input should be a matrix with number of columns",
            " equal to either number of original predictors or",
            " number of predictors in model\n")) } }
    
    # If the number of columns of the matrix is not length equal to 
    # number of predictors it is expanded to that size.
    
    if (is.matrix(x) && ncol(x) == 
        length(unique(pmars.model$model[, "pred1"])) && 
        ncol(x) != ncol(pmars.model$ranges.and.medians)) {
            tempmatrix<-x
    
    x<-matrix(nrow=nrow(tempmatrix), 
        ncol=ncol(pmars.model$ranges.and.medians),data = 0)
    for (i in 1:length(unique(pmars.model$model[, "pred1"])))  {
     for (j in 1:nrow(tempmatrix)) {
       x[j,sort(unique(pmars.model$model[,"pred1"]))[i]]<-x[j] } } }
       
    # If x is a vector put it into matrix form expanding it if it is  
    # of length equal to only the number of predictors in the model 
    # in `pmars.model'
    
    if (!(is.matrix(x))) {
    if (!(length(x) == ncol(pmars.model$ranges.and.medians) || 
        length(x) == unique(pmars.model$model[, "pred1"]))) {
        stop(paste(
            "The vector of values must be equal in length to",
            "  either the number of original predictors or",
            " predictors in the model\n") ) }
        if (length(x) == unique(pmars.model$model[, "pred1"]) && 
            length(x) != ncol(pmars.model$ranges.and.medians)) {
            x = rep(0, ncol(pmars.model$ranges.and.medians))
            for (i in 1:length(unique(pmars.model$model[, "pred1"]))) {
                x[sort(unique(pmars.model$model[, "pred1"]))[i]]<-x[i] } }
                x = t(as.matrix(x)) }
 
    # Checking to see if there are factor variables in the model
    
    if (dimnames(pmars.model$model)[[2]][3] == "level1") {
        level1<-TRUE
        pmars.model$model<-pmars.model$model[,c(1:(5+pmars.model$responses))]

    # if(dimnames(pmars.model$model)[[2]][6] == 
    # "level2"){level2<-TRUE}else{level2<-FALSE}
    
    } else {
       level1 <- FALSE
       pmars.model$model = pmars.model$model[, c(1:(4+pmars.model$responses))]
   
    # if(dimnames(pmars.model$model)[[2]][5] == "level2")
    # {level2<-TRUE}else{level2<-FALSE}
   
    }
    
    # Setting up the fitted responses matrix
    
    responses<-pmars.model$responses
    Y = matrix(ncol = responses, nrow = nrow(x), data = rep(0, nrow(x)))
    Y1 = matrix(ncol = 1, nrow = nrow(x), data = rep(0, nrow(x)))
    Y2 = matrix(ncol = 1, nrow = nrow(x), data = rep(0, nrow(x)))
    if (is.logical(intercept)){
        if (intercept == TRUE) {
            for (i in 1:responses)Y[,i] = 
            pmars.model$model[1,ncol(pmars.model$model)-responses+i] }
    else {
        if (intercept == FALSE) {
            for (i in 1:responses)Y[,i] = 0.0 } } }
	else {
	   if (is.numeric(intercept)) {
	     if (length(intercept)==responses) {
	       for (i in 1:responses)Y[,i] = intercept[i] }
	     else {
	       if (length(intercept) != 1) {
	                stop("Intercept arguement mispecified \n") }
	       for (i in 1:responses)Y[,i] = intercept } } }
       
    # Computing fitted values
    
	if (pmars.model$model.size>1) {
	   for (i in 2:pmars.model$model.size)  { 
	     Y2[] = 1   
	     Y1[] = x[,pmars.model$model[i, "pred1"]]
	     if (!is.na(pmars.model$model[i, "knot1"])) {
	       Y1 = Y1 - pmars.model$model[i,"knot1"]
	       Y1[Y1 < 0,] = 0 }
	     if (level1) {
	     if (!is.na(pmars.model$model[i, "level1"])){
	       Y1<- (Y1 == pmars.model$model[i, "level1"]) } }
	   if (!is.na(pmars.model$model[i, "pred2"]) & 
	    pmars.model$model[i, "pred2"] != 0)  {
	     Y2[] = x[,pmars.model$model[i,"pred2"]]
	     if (!is.na(pmars.model$model[i,"knot2" ]))  {
	       Y2 = Y2 - pmars.model$model[i,"knot2"]
	       Y2[Y2 < 0,] = 0 } }
	   
    for (j in 1:responses){Y[,j]<-Y[,j]+(Y1 * Y2 * 
        pmars.model$model[i,ncol(pmars.model$model)-responses+j])} } }
        
	# If classification is to be used the original polymars fitting 
	# expanded the response into a vector of indicator variables. The 
	# largest of the responses correspondes to the fitted class for 
	# each case.
	
	if (classify == TRUE) {
   	for (i in 1:nrow(Y)) {
     	Y[i,]<-Y[i,]==max(Y[i,]) }
   	if (is.matrix(pmars.model$conversion))
   	Z<-Y
   	Y<-matrix(ncol=1,nrow=nrow(Z))
   	for (i in 1:nrow(Y)) {
     	for (j in 1:ncol(Z)) {  
       if (Z[i,j] == 1) Y[i,] = pmars.model$conversion[j] } }}
        
    # Return Value:
    return(Y)
}


# ------------------------------------------------------------------------------


summaryBIpolymars = 
function(object,...)
{	# A slightly modified copy from R's contributed package polspline

    ## if (class(object)!="polymars")
    ## stop("object is not a polymars object")
    
    pmars.model = object
    cat("Call:\n")
    print(pmars.model$call)
    cat("\nModel fitting\n\n")
    print(pmars.model$fitting)
    cat("\n\nModel produced\n\n")
    print(pmars.model$model)
    if (pmars.model$responses != 1)
            cat("\nRESPONSES :", pmars.model$responses, "\n")
    if (!is.null(pmars.model$Rsquared))
            cat("\nRsquared :",round(pmars.model$Rsquared,3),"\n")
            
    # Return Value:
    invisible()
}


# ------------------------------------------------------------------------------


