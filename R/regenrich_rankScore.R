#' @rdname regenrich_rankScore
#' @export
setGeneric("regenrich_rankScore",
    function(object) standardGeneric("regenrich_rankScore"))

#' Regulator scoring and ranking
#'
#' As the fourth step of RegEnrich analysis, regulator ranking
#' is followed by differential expression analysis (regenrich_diffExpr),
#' regulator-target network inference (regenrich_network), and
#' enrichment analysis (regenrich_enrich).
#'
#' @param object a `RegenrichSet` object, to which
#' \code{\link{regenrich_diffExpr}}, \code{\link{regenrich_network}},
#' and \code{\link{regenrich_enrich}} functions all have been already applied.
#'
#' @return This function returns a `RegenrichSet` object
#' with an updated `resScore` slots, which is a `regEnrichScore` (also
#' `data.frame`) object, and an updated
#' `paramsIn` slot. In the `regEnrichScore` object there are five columns,
#' which are 'reg' (regulator), 'negLogPDEA' (-log10(p values of differential
#' expression analysis)), 'negLogPEnrich' (-log10(p values of enrichment
#' analysis), 'logFC' (log2 fold changes), and 'score' (RegEnrich ranking
#' score).
#' @rdname regenrich_rankScore
#' @include regenrichClasses.R
#' @seealso Previous step \code{\link{regenrich_enrich}}.
#' @export
#' @examples
#' \dontrun{
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
#' object = regenrich_enrich(object)
#'
#' # Regulators ranking
#' (object = regenrich_rankScore(object))
#' }
setMethod("regenrich_rankScore", signature = "RegenrichSet",
    definition = function(object) {
        enrichTest = object@paramsIn$enrichTest
        enrichTest = match.arg(enrichTest, enrichTest)

        resEnrich = object@resEnrich
        pFC = object@resDEA@pFC

        if (enrichTest == "FET") {
            resFET = resEnrich@allResult
            stopifnot(rownames(resFET) > 0)
            res = .rankScore(resFET = resFET, pDEA = pFC[, seq_len(2)],
                fcDEA = pFC[, c(1, 3)])
        } else if (enrichTest == "GSEA") {
            resSEA = resEnrich@allResult
            stopifnot(rownames(resSEA) > 0)
            res = .rankScore(resSEA = resSEA, pDEA = pFC[, seq_len(2)],
                fcDEA = pFC[, c(1, 3)])
        } else {
            stop("'enrichTest' must be 'FET' or 'GSEA'.")
        }

        object@resScore = res
        return(object)
    })



## The key function to calculate RegEnrich rank scores.
## @description The key function to generate a summary table
## for the results of differential analysis, enrichment
## analysis and RegEnrich rank scores.  @param resFET Result
## table by \code{regFET} @param resSEA Result table by
## \code{regSEA}.  If \code{resFET} is provided, this
## \code{resSEA} will not be used.  The default is
## \code{NULL}.  @param pDEA A data.frame of p-values for the
## regulators by differential expression analysis. The first
## column is the gene name/ID, the second column is the
## p-values.  @param fcDEA Same format as pDEA.  log2 fold
## change for the regulators by differential expression
## analysis (optional). The default is \code{NULL}.  @return
## A \code{regEnrichScore} object, including \code{reg}
## (regulators), \code{negLogPDEA} (-log10(pD)),
## \code{negLogPEnrich} (-log10(pE)), \code{logFC}
## (log2(fold change)), \code{score} (RegEnrich scores),
## @examples { \dontrun{ } }
.rankScore = function(resFET = NULL, resSEA = NULL, pDEA, fcDEA = NULL) {
    # stopifnot(!(is.null(resFET) | is.null(resSEA)))
    if (!is.null(resFET)) {
        enrichP = data.frame(reg = resFET$ID, pval = resFET$pvalue,
            padj = resFET$p.adjust, stringsAsFactors = FALSE)
    } else {
        if (!is.null(resSEA)) {
            enrichP = data.frame(reg = resSEA$regulator, pval = resSEA$pval,
                padj = resSEA$padj, stringsAsFactors = FALSE)
        } else {
            stop("Either 'resFET' or 'resSEA' should be provided!")
        }
    }

    # Calculate -log10(p)
    stopifnot(!any(duplicated(pDEA[, 1])))  # Gene duplicates are not allowed
    id = match(enrichP$reg, pDEA[, 1])
    negLogPDEA = -log10(pDEA[id, 2])
    negLogPDEA[is.na(negLogPDEA)] = 0

    negLogPEnrich = -log10(enrichP$pval)
    negLogPEnrich[is.infinite(negLogPEnrich)] =
        max(negLogPEnrich[!is.infinite(negLogPEnrich)])

    # Include fold change
    if (!is.null(fcDEA)) {
        id = match(enrichP$reg, fcDEA[, 1])
        logFC = fcDEA[id, 2]
    } else {
        logFC = Inf
    }

    # Calculate scores
    normFun = function(x) {
        mx = max(x[is.finite(x)])
        mnx = min(x[is.finite(x)])
        if ((mx - mnx) == 0) {
            return(rep(0, length(x)))
        } else {
            normx = (x - mnx)/(mx - mnx)
            normx[is.infinite(normx)] = 1
            return(normx)
        }
    }
    score = normFun(negLogPDEA) + normFun(negLogPEnrich)

    res = data.frame(reg = enrichP$reg, negLogPDEA = negLogPDEA,
        negLogPEnrich = negLogPEnrich, logFC = logFC, score = score,
        stringsAsFactors = FALSE)
    if (nrow(res) > 0) {
        res = sortDataframe(res, "score", decreasing = TRUE)
        rownames(res) = seq_len(nrow(res))
    }
    res = structure(res, class = c("regEnrichScore", "data.frame"))
    return(res)
}
