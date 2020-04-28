
#' Result accessor functions
#'
#' @description
#' \itemize{
#' \item results_expr accesses raw expression data.
#' \item results_DEA accesses results from differential expression analysis.
#' \item results_topNet accesses results from network inference.
#' \item retults_enrich accesses results from FET/GSEA enrichment analysis.
#' \item results_score accesses results from regulator scoring and ranking.
#' }
#' @param object RegenrichSet object.
#' @return results_expr retures an expression matrix.
#' @rdname results_expr
#' @export
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
#' object = regenrich_diffExpr(object)
#' results_expr(object)
#' results_DEA(object)
#'
#' # Network inference using 'COEN' method
#' object = regenrich_network(object)
#' results_topNet(object)
#'
#' # Enrichment analysis by Fisher's exact test (FET)
#' object = regenrich_enrich(object)
#' results_enrich(object)
#'
#' # Regulators ranking
#' object = regenrich_rankScore(object)
#' results_score(object)
# }
results_expr = function(object) {
    object@assayRaw
}

#' @return results_DEA returns a list result of differentila analysis.
#' @rdname results_expr
#' @export
results_DEA = function(object) {
    # object@resDEA
    mcols(object)
}

#' @return results_topNet returns a TopNetwork object.
#' @rdname results_expr
#' @export
results_topNet = function(object) {
    object@topNetwork
}

#' @return results_enrich returns an Enrich object by either FET or GSEA
#' method.
#' @rdname results_expr
#' @export
results_enrich = function(object) {
    object@resEnrich
}

#' @return results_score returns an data frame of summarized ranking scores
#' of regulators.
#' @rdname results_expr
#' @export
results_score = function(object) {
    object@resScore
}
