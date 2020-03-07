
#' @rdname regFET
# @export
setGeneric("regFET", def = function(object, ...) standardGeneric("regFET"))


.regFET = function(object, namedScores, namedScoresCutoffs = 0.05,
    minSize = 5, maxSize = 5000, pvalueCutoff = 0.05, pAdjustMethod = "BH",
    qvalueCutoff = 0.2, regAltName = NULL, universe = NULL) {
    netTable = object@netTable
    network = object@net
    tarReg = object@tarReg
    stopifnot(length(names(namedScores)) == length(namedScores))

    if (is.null(regAltName)) {
        regAltName = names(network)
    }

    stopifnot(is.numeric(namedScoresCutoffs))

    stopifnot(length(network) == length(regAltName))

    # The data format for enricher_internal()
    netENV = as.environment(list(PATHID2NAME = stats::setNames(names(network),
        regAltName), EXTID2PATHID = tarReg, PATHID2EXTID = network))

    topGene = names(namedScores[namedScores <= namedScoresCutoffs])

    # Fisher exact test
    enricher_internal = utils::getFromNamespace("enricher_internal", ns = "DOSE")
    y = enricher_internal(gene = topGene, pvalueCutoff = 1,
        pAdjustMethod = pAdjustMethod, universe = universe,
        minGSSize = minSize,
        maxGSSize = maxSize, qvalueCutoff = 1, USER_DATA = netENV)

    # The results to show
    new(Class = "Enrich", topResult = y@result[(y@result$pvalue <=
        pvalueCutoff) & (y@result$p.adjust <= pvalueCutoff) &
        (y@result$qvalue <= qvalueCutoff), ], allResult = y@result,
        gene = y@gene, namedScores = namedScores, type = "FET")
}

#' Enrichment analysi by Fisher's exact test.
#' @description Enrichment for regulators based on Fisher exact test.
#' @param object a topNetwork object, The result returned from
#' \code{\link{topNet}} function.
#' @param namedScores a named numeric vector of scores,
#' the names of the scores are the genes to perform enrichment analysis.
#' And the names should be the same as in the topNetwork object.
#' Here the scores are p-value of each gene.
#' @param namedScoresCutoffs the cutoff of \code{namedScores}.
#' @param minSize the minimum number (default 5) of target genes.
#' @param maxSize the maximum number (default 5000) of target genes.
#' @param pvalueCutoff numeric, the cutoff for adjusted enrichment p value.
#' This is used for obtaining the `topResult` slot in the final `Enrich`
#' object. Default is 0.05.
#' @param pAdjustMethod p adjust method, one of 'holm', 'hochberg',
#' 'hommel', 'bonferroni', 'BH', 'BY', 'fdr', 'none'.
#' @param qvalueCutoff cutoff of q-value.
#' @param regAltName alternative names of the regulators in the pathway.
#' @param universe a vector of charactors. Background target genes.
#' @param ... additional arguments.
#' @return An \code{enrichFET} object, consisting of 5 elements:\cr
#' \code{topResult}, the enrichment information of regulators that
#' pass minSize maxSize pvalueCutoff and qvalueCutoff cutoffs;\cr
#' \code{allResult}, the enrichment information of all regulators;\cr
#' \code{gene}, the genes pass the \code{namedScoresCutoffs}; \cr
#' \code{universe}, total target genes in the network;\cr
#' and \code{geneSets}, a list of regulators (the names of the list)
#' and their targets (the elements of the list).
#' @include regenrichClasses.R
#' @import DOSE
#' @importClassesFrom DOSE enrichResult
#' @rdname regFET
#' @seealso \code{\link{regSEA}}
# @export
setMethod("regFET", c("TopNetwork"), .regFET)



