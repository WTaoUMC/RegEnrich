#' Co-expression network construction
#' @description Constructing co-expression network using Weighted gene
#' co-expression network analysis (WGCNA).
#' @param expr gene expression data, it is either a matrix or a data frame.
#' By default, each row represents a gene, each column represents a sample.
#' @param reg a vector of regulators. by default, these are transcription
#' (co-)factors defined by three literatures/databases, namely RegNet,
#' TRRUST, and Marbach2016.
#' @param rowSample logic, if TRUE, each row represents a sample.
#' Otherwise, each column represents a sample. The default is FALSE.
#' @param softPower numeric, a soft power to achieve scale free topology.
#' If not provided, the parameter will be picked automatically by
#' \code{\link{plotSoftPower}} function from the WGCNA package.
#' @param networkType network type. Allowed values are (unique abbreviations
#' of) 'unsigned' (default), 'signed', 'signed hybrid'.
#' See \code{\link{adjacency}}.
#' @param TOMDenom a character string specifying the TOM variant to be used.
#' Recognized values are 'min' giving the standard TOM described in Zhang
#' and Horvath (2005), and 'mean' in which the min function in the
#' denominator is replaced by mean. The 'mean' may produce better results
#' but at this time should be considered experimental.
#' @param RsquaredCut desired minimum scale free topology fitting index R^2.
#' The default is 0.85.
#' @param edgeThreshold numeric, the threshold to remove the low weighted
#' edges, the default is NULL, which means no edges will be removed.
#' @param trace logical. To show the progress or not (default).
#'
#' @return A list of two elements. weightHi: data frame of weighted edges
#' between regulators and targets; and networkType: 'COEN'.
#' @import WGCNA
#' @include globals.R
# @export
COEN = function(expr, reg = TFs$TF_name, rowSample = FALSE, softPower = NULL,
    networkType = "unsigned", TOMDenom = "min", RsquaredCut = 0.85,
    edgeThreshold = NULL, trace = FALSE) {
    # Enforce expression matrix to be the format that each row is
    # a sample.
    if (!rowSample) {
        expr = t(expr)
        rowSample = !rowSample
    }

    if (!any(reg %in% colnames(expr))) {
        stop("No expression data for the regulators.")
    }

    # Find soft power if not provided
    if (is.null(softPower)) {
        powerVector = c(seq_len(10), seq(12, 20, by = 2))
        tmp = utils::capture.output(sft <- plotSoftPower(expr,
            rowSample = rowSample, powerVector = powerVector,
            RsquaredCut = RsquaredCut, networkType = networkType,
            verbose = 0), file = NULL)  # suppress the output
        softPower = sft$powerEstimate
    }
    # stopifnot((!is.null(softPower)) & (!is.na(softPower)))
    if (is.null(softPower) || is.na(softPower)) {
        stop("RsquaredCut is too high to achieve.")
    }

    # Adjacency matrix and Topological Overlap Matrix (TOM)
    if (trace) {
        cat("Calculating Topological Overlap Matrix (TOM)\n")
    }
    tmp = utils::capture.output(TOM <- TOMsimilarity(adjMat = adjacency(expr,
        power = softPower, type = networkType), TOMType = networkType,
        TOMDenom = TOMDenom), file = NULL)
    dimnames(TOM) = list(colnames(expr), colnames(expr))
    if (trace) {
        cat("Collecting memory garbage\n")
    }
    collectGarbage()

    # Edges
    if (trace) {
        cat("Converting graph matrix to edges\n\n")
    }

    edge = mat2Edge(mat = TOM, mode = "upper")
    colnames(edge) = c("from.gene", "to.gene", "weight")

    nEdge = nrow(edge)
    if (nEdge == 0)
        stop("No edges were remained!")

    # Only remain the edges connecting regulators
    return(list(weightHi = edge[(edge$from.gene %in% reg) | (edge$to.gene %in%
        reg), ], networkType = "COEN"))
}
