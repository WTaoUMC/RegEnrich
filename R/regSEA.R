#' @rdname regSEA
# @export
setGeneric("regSEA", def = function(object, ...) standardGeneric("regSEA"))

.regSEA = function(object, namedScores, minSize = 5, maxSize = 5000,
    pvalueCutoff = 0.05, nperm = 10000, ...) {
    network = object@net
    stopifnot(length(names(namedScores)) == length(namedScores))

    ### scaled p values
    namedScores1 = stats::setNames(as.vector(scale(-log10(namedScores))),
        names(namedScores))

    if (utils::packageVersion("fgsea") >= as.package_version("1.13.2")) {
        fgseaFun = fgseaSimple
    } else {
        fgseaFun = fgsea
    }
    regseaRes = fgseaFun(pathways = network, stats = namedScores1,
        minSize = minSize, maxSize = maxSize, nperm = nperm,
        ...)

    # consider only one side
    regseaRes$pval = ifelse(sign(regseaRes$ES) >= 0, regseaRes$pval,
        1 - regseaRes$pval)
    regseaRes$nMoreExtreme = ifelse(sign(regseaRes$ES) >= 0,
        regseaRes$nMoreExtreme, nperm - regseaRes$nMoreExtreme)
    regseaRes$padj = stats::p.adjust(regseaRes$pval, method = "BH")
    regseaRes = regseaRes[order(regseaRes$pval), ]
    colnames(regseaRes)[1] = "regulator"

    # The results to show
    new(Class = "Enrich", topResult = as.data.frame(regseaRes[regseaRes$padj <
        pvalueCutoff & regseaRes$ES > 0]), allResult = as.data.frame(regseaRes),
        gene = names(namedScores), namedScores = namedScores,
        type = "GSEA")
}

#' Enrichment analysi by gene set enrichment analysis (GSEA).
#' @param object a \code{topNetwork} object, the result returned by
#' \code{topNet} funciton.
#' @param namedScores A named numeric vector of scores,
#' the names of the scores are the genes to perform enrichment analysis.
#' And the names should be the same as in the topNetwork object.
#' Here the scores are p-value of each gene.
#' @param minSize integer, the minimum number (default 5) of target genes.
#' @param maxSize integer, the maximum number (default 5000) of target genes.
#' @param pvalueCutoff numeric, the cutoff for adjusted enrichment p value.
#' This is used for obtaining the `topResult` slot in the final `Enrich`
#' object. Default is 0.05.
#' @param nperm integer, number of permutations. The minimial possible
#' nominal p-value is about 1/nperm. The default is 10000.
#' @param ... The rest parameters in \code{\link{fgsea}} function.
#'
#' @return a list of two elements: a table with GSEA results (see
#' \code{\link{fgsea}}) and a ggplot object.
#' @include regenrichClasses.R
#' @import fgsea
#' @rdname regSEA
#' @seealso \code{\link{regFET}}
# @export
#'
setMethod("regSEA", signature = "TopNetwork", definition = .regSEA)

