#' @rdname regenrich_diffExpr
#' @export
setGeneric("regenrich_diffExpr",
    function(object, ...) standardGeneric("regenrich_diffExpr"))

.regenrich_diffExpr = function(object, ...) {
    argsIn = list(...)

    deaParams = c("method", "minMeanExpr", "design", "reduced",
        "contrast", "coef", "name", "fitType", "sfType", "betaPrior",
        "minReplicatesForReplace", "useT", "minmu", "parallel",
        "BPPARAM", "altHypothesis", "listValues", "cooksCutoff",
        "independentFiltering", "alpha", "filter", "theta", "filterFun",
        "addMLE", "blind", "ndups", "spacing", "block", "correlation",
        "weights", "proportion", "stdev.coef.lim", "trend", "robust",
        "winsor.tail.p")
    object = checkParams(object, argsIn, deaParams)
    ParamsIn = getParamsIn(object)
    deaParams = ParamsIn[deaParams]

    deaParams$method = match.arg(deaParams$method, c("Wald_DESeq2",
        "LRT_DESeq2", "limma", "LRT_LM"))

    resCall = as.call(c(list(quote(DEA), expr = quote(object@assayData),
        pData = quote(object@phenoData)), deaParams))
    resDEA = eval(resCall)

    # expression after filtering (and normalization)
    object@assayData = resDEA$normExprHigh
    # differential expression results (DeaSet object)
    object@resDEA = DeaSet(resDEA$pFC, resDEA$fullRes, resDEA$fit)
    # update paramsOut
    object@paramsOut$method = deaParams$method

    # Since differential expression analysis changes, the
    # network, enrichment and regulator ranking must change
    object@network = NULL
    object@topNetP = NULL
    object@resEnrich = NULL
    object@resScore = NULL
    object@paramsOut = list(method = deaParams$method, network = NULL,
        enrichTest = NULL, percent = NULL)

    return(object)
}
#' Differential expression analysis step
#'
#' This is the first step of RegEnrich analysis.
#' differential expression analysis by this function needs to
#' be performed on a `RegenrichSet` object.
#'
#' @param object a `RegenrichSet` object, which is initialized by
#' \code{\link{RegenrichSet}} function.
#'
#' @param ... arguments for differential analysis.
#' After constructing a `RegenrichSet` object,
#' all arguments for RegEnrich analysis have been initialized and
#' stored in `paramsIn`` slot. while the arguments for differential analysis
#' can be re-specified here.\cr\cr
#' These arguments include 'method', 'minMeanExpr', 'design', 'reduced',
#' 'contrast',
#' 'coef', 'name', 'fitType', 'sfType', 'betaPrior', 'minReplicatesForReplace',
#' 'useT', 'minmu', 'parallel', 'BPPARAM', 'altHypothesis',
#' 'listValues', 'cooksCutoff', 'independentFiltering', 'alpha',
#' 'filter', 'theta', 'filterFun', 'addMLE', 'blind', 'ndups',
#' 'spacing', 'block', 'correlation', 'weights', 'proportion',
#' 'stdev.coef.lim', 'trend', 'robust', and 'winsor.tail.p'.\cr
#' See \code{\link{RegenrichSet}} function for more details about these
#' arguments.
#'
#' @return This function returns a `RegenrichSet` object with an updated
#' `resDEA` slot, which is a `DeaSet` object, and an updated `paramsIn` slot.
#' See \code{\link{DeaSet}} function for more details about `DeaSet` class.
#' If an argument not in the above list is specified in the regenrich_diffExpr
#' function, a warning or error will be raised.
#'
#' @examples
#' # library(RegEnrich)
#' # Initializing a 'RegenrichSet' object
#' data = log2(Lyme_GSE63085$FPKM + 1)
#' pData = Lyme_GSE63085$sampleInfo
#' x = apply(data, 1, sd)
#' data1 = data[seq_len(2000), ]
#'
#' pData$week = as.factor(pData$week)
#' pData$patientID = as.factor(sub('(\\d+)-(\\d+)', '\\1_\\2',
#'                             pData$patientID))
#'
#' design = model.matrix(~0 + patientID + week,
#'                       data = pData)
#' object = RegenrichSet(expr = data1,
#'                       pData = pData,
#'                       method = 'limma', minMeanExpr = 0,
#'                       design = design,
#'                       contrast = c(rep(0, ncol(design) - 1), 1),
#'                       networkConstruction = 'COEN',
#'                       enrichTest = 'FET')
#'
#' # Using the predifined parameters in the previous step
#' (object = regenrich_diffExpr(object))
#'
#' # re-specifying parameter 'minMeanExpr'
#' print(slot(object, 'paramsIn')$minMeanExpr)
#' (object = regenrich_diffExpr(object, minMeanExpr = 1))
#' print(slot(object, 'paramsIn')$minMeanExpr)
#'
#' # Unrecognized argument 'unrecognizedArg' (Error)
#' # (object = regenrich_diffExpr(object, minMeanExpr = 1,
#' #                              unrecognizedArg = 23))
#'
#' # Argument not for differential expression analysis (Warning)
#' print(slot(object, 'paramsIn')$networkConstruction)
#' (object = regenrich_diffExpr(object, minMeanExpr = 1,
#'                              networkConstruction = 'GRN'))
#' print(slot(object, 'paramsIn')$networkConstruction)
#' @rdname regenrich_diffExpr
#' @include regenrichClasses.R
#' @seealso Initialization of a `RegenrichSet` object
#' \code{\link{RegenrichSet}},and next step
#' \code{\link{regenrich_network}}.
#'
#' @export
#'
setMethod(f = "regenrich_diffExpr", signature = "RegenrichSet",
    definition = .regenrich_diffExpr)
