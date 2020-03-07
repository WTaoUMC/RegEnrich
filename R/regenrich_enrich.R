#' @rdname regenrich_enrich
#' @export
setGeneric("regenrich_enrich",
    function(object, ...) standardGeneric("regenrich_enrich"))
.regenrich_enrich = function(object, ...) {
    argsIn = list(...)
    mustInArgs = c("enrichTest", "namedScoresCutoffs", "minSize",
        "maxSize", "pvalueCutoff", "qvalueCutoff", "regAltName",
        "universe", "minSize", "maxSize", "pvalueCutoff", "nperm")
    object = checkParams(object, argsIn, mustInArgs)

    enrichTest = object@paramsIn$enrichTest
    enrichTest = match.arg(enrichTest, enrichTest)

    namedScores = stats::setNames(object@resDEA@pFC[, "p"],
        rownames(object@resDEA@pFC))
    if (enrichTest == "FET") {
        params = getParamsIn(object)[c("namedScoresCutoffs",
            "minSize", "maxSize", "pvalueCutoff", "qvalueCutoff",
            "regAltName", "universe")]
        resCall = as.call(c(list(quote(regFET), object = quote(object@topNetP),
            namedScores = quote(namedScores)), params))
    } else if (enrichTest == "GSEA") {
        params = getParamsIn(object)[c("minSize", "maxSize",
            "pvalueCutoff", "nperm")]
        resCall = as.call(c(list(quote(regSEA), object = quote(object@topNetP),
            namedScores = quote(namedScores)), params))
    } else {
        stop("'enrichTest' must be 'FET' or 'GSEA'.")
    }
    res = eval(resCall, envir = sys.frames())

    object@resEnrich = res
    object@paramsOut$enrichTest = enrichTest

    # Since enrichment changes, the regulator ranking must change
    object@resScore = NULL
    object@paramsOut = list(method = object@paramsOut$method,
        network = object@paramsOut$network, enrichTest = enrichTest,
        percent = object@paramsOut$percent)

    return(object)
}

#' Enrichment analysis step
#'
#' As the thrid step of RegEnrich analysis, enrichment analysis
#' is followed by differential expression analysis (regenrich_diffExpr),
#' and regulator-target network inference (regenrich_network).
#'
#' @param object a `RegenrichSet` object, to which
#' \code{\link{regenrich_diffExpr}}, and \code{\link{regenrich_network}},
#' functions have been already applied.
#'
#' @param ... arguments for enrichment analysis.
#' After constructing a `RegenrichSet` object using \code{\link{RegenrichSet}}
#' function, all arguments for RegEnrich analysis have been initialized and
#' stored in `paramsIn`` slot. The arguments for enrichment analysis can be
#' re-specified here.\cr\cr
#' These arguments include 'enrichTest', 'namedScoresCutoffs', 'minSize',
#' 'maxSize',
#' 'pvalueCutoff','qvalueCutoff', 'regAltName', 'universe',
#' 'minSize', 'maxSize', 'pvalueCutoff', and 'nperm'.\cr\cr
#' See \code{\link{RegenrichSet}} function for more details about these
#' arguments.
#'
#' @return This function returns a `RegenrichSet` object with an updated
#' `resEnrich` slots, which is `Enrich` objects, and an updated `paramsIn`
#' slot.
#' See \code{\link{Enrich-class}} function for more details about `Enrich`
#' class.
#' @rdname regenrich_enrich
#' @include regenrichClasses.R
#' @seealso Previous step \code{\link{regenrich_network}},
#' and next step \code{\link{regenrich_rankScore}}.
#' @export
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
#' # Differential expression analysis
#' object = regenrich_diffExpr(object)
#'
#' # Network inference using 'COEN' method
#' object = regenrich_network(object)
#'
#' # Enrichment analysis by Fisher's exact test (FET)
#' (object = regenrich_enrich(object))

setMethod(f = "regenrich_enrich", signature = "RegenrichSet",
    definition = .regenrich_enrich)


