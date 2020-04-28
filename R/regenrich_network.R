#' @rdname regenrich_network
#' @export
setGeneric("regenrich_network", function(object, ...) standardGeneric("regenrich_network"))
.regenrich_network = function(object, ...) {
  argsIn = list(...)
  
  mustInArgs = c("networkConstruction", "reg", "rowSample",
                 "softPower", "networkType", "TOMDenom", "RsquaredCut",
                 "edgeThreshold", "K", "nbTrees", "importanceMeasure",
                 "trace", "BPPARAM", "minR", "topNetPercent", "directed")
  
  object = checkParams(object, argsIn, mustInArgs)
  
  topNetPercent = getParamsIn(object, "topNetPercent")
  stopifnot(is.numeric(topNetPercent))
  
  networkConstruction = match.arg(getParamsIn(object, "networkConstruction"),
                                  c("COEN", "GRN", "new"))  # COEN/GRN/new
  if (networkConstruction == "new") {
    if (!(topNetPercent <= 100 && topNetPercent > 0)) {
      stop("'topNetPercent' must be between 0 (excluding) and 100 ",
           "(including), when 'networkConstruction' is 'new'.")
    }
  } else {
    if (!(topNetPercent <= 10 && topNetPercent > 0)) {
      stop("'topNetPercent' must be between 0 (excluding) and 10 ",
           "(including), when 'networkConstruction' is not 'new'.")
    }
  }
  reg = getParamsIn(object)[["reg"]]
  if (networkConstruction == "COEN") {
    params = getParamsIn(object)[c("rowSample", "softPower",
                                   "networkType", "TOMDenom", "RsquaredCut", 
                                   "edgeThreshold", "trace")]
    resCall = as.call(c(list(quote(COEN), expr = quote(assay(object))),
                        reg = quote(reg), params))
    directed = FALSE
  } else if (networkConstruction == "GRN") {
    params = getParamsIn(object)[c("rowSample", "K",
                                   "nbTrees", "importanceMeasure", "trace",
                                   "BPPARAM", "minR")]
    resCall = as.call(c(list(quote(GRN), expr = quote(assay(object))), 
                        reg = quote(reg), params))
    directed = TRUE
  } else if (networkConstruction == "new") {
    stop("When networkConstruction = 'new', use 'regenrich_network<-' to ",
         "assign the network.")
  } else {
    stop("Unknown networkConstruction = \n", networkConstruction)
  }
  net = eval(resCall)
  
  object@network = newTopNetwork(networkEdgeTable = net$weightHi,
                                 networkConstruction = networkConstruction, 
                                 directed = directed,
                                 reg = reg, percent = 100)
  
  topNetP = topNet(object@network@elementset, percent = topNetPercent,
                   directed = object@network@directed, reg = reg)
  object@topNetwork = newTopNetwork(as.data.frame(topNetP$netTable),
                                    networkConstruction = networkConstruction, 
                                    directed = TRUE,
                                    reg = reg, percent = topNetPercent)
  
  # Since network changes, the enrichment and regulator ranking
  # must change
  object@resEnrich = newEnrich()
  object@resScore = newScore()
  object@paramsOut = list(DeaMethod = object@paramsOut$DeaMethod,
                          networkType = networkConstruction, 
                          enrichTest = NULL,
                          percent = topNetPercent)
  return(object)
}
#' Regulator-target network inference step
#'
#' As the second step of RegEnrich analysis, network inference
#' is followed by differential expression analysis (regenrich_diffExpr).
#'
#' @param object a `RegenrichSet` object, to which
#' \code{\link{regenrich_diffExpr}} function has been already applied.
#'
#' @param ... arguments for network inference.
#' After constructing a `RegenrichSet` object using \code{\link{RegenrichSet}}
#' function, all arguments for RegEnrich analysis have been initialized and
#' stored in `paramsIn`` slot. The arguments for network inference can be
#' re-specified here.\cr\cr
#' These arguments include 'networkConstruction', 'reg', 'rowSample',
#' 'softPower', 'networkType', 'TOMDenom', 'RsquaredCut', 'edgeThreshold',
#' 'K', 'nbTrees', 'importanceMeasure', 'trace',
#' 'BPPARAM', 'minR', 'topNetPercent', and 'directed'.\cr\cr
#' See \code{\link{RegenrichSet}} function for more details about these
#' arguments.
#'
#' @return This function returns a `RegenrichSet` object with an updated
#' `network` and `topNetP` slots, which are `TopNetwork` objects, and
#' an updated `paramsIn` slot.
#' See \code{\link{TopNetwork-class}} class for more details.
#'
#' @examples
#' # library(RegEnrich)
#' data("Lyme_GSE63085")
#' data("TFs")
#' 
#' data = log2(Lyme_GSE63085$FPKM + 1)
#' colData = Lyme_GSE63085$sampleInfo
#' 
#' # Take first 2000 rows for example
#' data1 = data[seq(2000), ]
#'
#' design = model.matrix(~0 + patientID + week, data = colData)
#' 
#' # Initializing a 'RegenrichSet' object
#' object = RegenrichSet(expr = data1,
#'                       colData = colData,
#'                       method = 'limma', minMeanExpr = 0,
#'                       design = design,
#'                       contrast = c(rep(0, ncol(design) - 1), 1),
#'                       networkConstruction = 'COEN',
#'                       enrichTest = 'FET')
#'
# \donttest{
#' # Differential expression analysis
#' (object = regenrich_diffExpr(object))
#'
#' # Network inference using 'COEN' method
#' (object = regenrich_network(object))
# }
#' @rdname regenrich_network
#' @include regenrichClasses.R
#'
#' @seealso Previous step \code{\link{regenrich_diffExpr}},
#' and next step \code{\link{regenrich_enrich}}. User defined
#' network \code{\link{regenrich_network<-}}
#'
#' @export
#'
setMethod("regenrich_network", signature = "RegenrichSet", .regenrich_network)



#' @rdname regenrich_network
#' @export
setGeneric("regenrich_network<-", function(object, value) standardGeneric("regenrich_network<-"))

.regenrich_network_Top = function(object, value) {
  # if (!is(results_DEA(object), "DeaSet")) {
  if (isEmptyPFC(results_DEA(object))) {
    stop("Differential expression analysis needs to ",
         "be performed first.")
  }
  networkConstruction = "new"
  
  # update network
  object@network = value
  
  # update topNetwork
  object@topNetwork = newTopNetwork(object@network@elementset,
                                    networkConstruction = networkConstruction,
                                    directed = TRUE, reg = object@paramsIn$reg, 
                                    percent = 100)
  
  object@paramsIn$topNetPercent = 100
  object@paramsIn$networkConstruction = networkConstruction
  object@paramsIn$directed = TRUE
  
  # Since network changes, the enrichment and regulator ranking
  # must change
  object@resEnrich = newEnrich()
  object@resScore = newScore()
  object@paramsOut = list(DeaMethod = object@paramsOut$DeaMethod,
                          networkType = networkConstruction, 
                          enrichTest = NULL, percent = 100)
  
  return(object)
}
#' Network provided by users.
#'
#' Provide a network to `RegenrichSet` object.
#'
#' @param object a `RegenrichSet` object, to which
#' \code{\link{regenrich_diffExpr}} function has been already applied.
#' @param value either a `TopNetwork` object or `data.frame` object.
#' If value is a `data.frame` object, then the number of columns of
#' @return This function returns a `RegenrichSet` object with an updated
#' `network` and `topNetP` slots, which are `TopNetwork` objects, and an
#' updated `paramsIn` slot.
#' See \code{\link{TopNetwork-class}} class for more details.
#'
#' @rdname regenrich_network
#' @export
setMethod(f = "regenrich_network<-", signature("RegenrichSet", "TopNetwork"), .regenrich_network_Top)


.regenrich_network_df = function(object, value) {
  if (isEmptyPFC(results_DEA(object))) {
    stop("Differential expression analysis needs to be performed first.")
  }
  value = as.data.frame(value)
  if (ncol(value) == 3) {
    if (!is.numeric(value[, 3])) {
      stop("The third column (weight) of 'value' is not numeric.")
    } else {
      colnames(value) = c("from.gene", "to.gene", "weight")
    }
  } else {
    stop("The number of columns of 'value' is equal to 3.")
  }
  
  networkConstruction = "new"
  
  # update network
  object@network = newTopNetwork(value,
                              networkConstruction = networkConstruction,
                              directed = TRUE, reg = object@paramsIn$reg, 
                              percent = 100)
  
  # update topNetwork
  object@topNetwork = newTopNetwork(object@network@elementset,
                                    networkConstruction = networkConstruction,
                                    directed = TRUE, reg = object@paramsIn$reg, 
                                    percent = 100)
  
  object@paramsIn$topNetPercent = 100
  object@paramsIn$networkConstruction = networkConstruction
  object@paramsIn$directed = TRUE
  
  # Since network changes, the enrichment and regulator ranking
  # must change
  object@resEnrich = newEnrich()
  object@resScore = newScore()
  object@paramsOut = list(DeaMethod = object@paramsOut$DeaMethod,
                          networkType = networkConstruction, 
                          enrichTest = NULL, percent = 100)
  
  return(object)
}

#' @rdname regenrich_network
#' @export
setMethod("regenrich_network<-", signature("RegenrichSet", "data.frame"), .regenrich_network_df)

