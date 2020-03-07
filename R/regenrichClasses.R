
#' DeaSet class
#'
#' The `DeaSet` object is to store differnetial expression analysis
#' results, for further enrichment and regulator ranking.
#' @slot pFC a data frame of three columns, i.e. gene (gene names/IDs),
#' p (differential p values), and logFC (log2 expression fold changes).
#' @slot fullRes a full result table from differential analysis. The
#' contents in this table varies depending on the differential analysis
#' methods.
#' @slot fit either NULL or a linear model object from differential analysis.
#' @return an object of DeaSet class.
#' @import methods
#' @export
setClass("DeaSet", slots = c(pFC = "data.frame", fullRes = "ANY",
    fit = "ANY"), prototype = list(pFC = data.frame(gene = character(0),
    p = numeric(0), logFC = numeric(0)), fullRes = NULL, fit = NULL))

#' DeaSet object creator
#'
#' This is `DeaSet` object creator function.
#' @param pFC a data frame of three columns, i.e. gene (gene names/IDs),
#' p (differential p values), and logFC (log2 expression fold changes).
#' @param fullRes a full result table from differential analysis. The
#' contents in this table varies depending on the differential analysis
#' methods.
#' @param fit either NULL or a linear model object from differential analysis.
#' @return an object of DeaSet class.
#' @export
#' @examples
#' pFC = data.frame(gene = paste0('gene', seq_len(10)),
#'                  p = seq(0.01, 0.3, length.out = 10),
#'                  logFC = rnorm(10))
#' fullRes = NULL
#' fit = NULL
#' dSet = DeaSet(pFC, fullRes, fit)
DeaSet = function(pFC, fullRes, fit) {
    stopifnot(is.data.frame(pFC))
    new("DeaSet", pFC = pFC, fullRes = fullRes, fit = fit)
}

#' TopNetwork class
#'
#' The `TopNetwork` object is to store either a full network (the percentage
#' of top edges is 100%) or a sub-network (with percentage of top edges
#' between 0 to 10).
#' @slot net list, in which the names of elements are regulators,
#' and the elements are targets of the regulators indicated by the
#' element name.
#' @slot netTable data frame. A table showing regulator-target
#' relationship. It contains 3 columns, representing 'from.gene'
#' ('regulators'), 'to.gene' ('targets') and 'weight', respectively.
#' @slot validReg vector, a set of valid regulators in the network.
#' @slot tarReg list, in which the names of elements are targets,
#' and the elements are regulators of the targets indicated by the
#' element name.
#' @slot directed logical, whether the network is directed.
#' @slot networkConstruction character, by which method this network
#' is constructed. Either 'COEN' (coexpression network using WGCNA),
#' or 'GRN' (gene regulatory network using random forest), or 'new'
#' (a network provided by the user).
#' @slot percent numeric, what percentage of the top edges are remained.
#' The value must be between 0 (excluding) and 100 (including).
#' @export
setClass("TopNetwork", slots = c(net = "list", netTable = "data.frame",
    validReg = "vector", tarReg = "list", directed = "logical",
    networkConstruction = "character", percent = "numeric"),
    prototype = list(net = list(),
        netTable = data.frame(from.gene = character(),
        to.gene = character(), weight = numeric(), stringsAsFactors = FALSE),
        validReg = c(), tarReg = list(), directed = TRUE,
        networkConstruction = "new",
        percent = 5))


setValidity("TopNetwork", function(object) {
    if (object@percent > 100 || object@percent <= 0) {
        "@percent must be between 0 (excluding) and 100 (including)"
    } else {
        if (length(object@networkConstruction) != 1 ||
            !object@networkConstruction %in%
            c("COEN", "GRN", "new")) {
            "@networkConstruction must be one of 'COEN', 'GRN', and 'new'"
        } else {
            TRUE
        }
    }
})

#' TopNetwork object creator
#'
#' This function create `TopNetwork` object using 3-column edge table.
#'
#' @param networkEdgeTable a data frame of 3 columns, representing 'from.gene'
#' ('regulators'), 'to.gene' ('targets') and 'weight', respectively.
#' @param reg a vector of gene regulators.
#' @param directed logical, whether the network is directed. Default is TRUE.
#' @param networkConstruction the method to construct this network.
#' Possible can be:
#' 'COEN', coexpression network;
#' 'GRN', gene regulatory network by random forest;
#' 'new' (default), meaning a network provided by user, rather than
#' infered based
#' on the expression data.
#' @param percent the percentage of edges in the original whole network.
#' Default is 100, meaning 100\% edges in whole network.
#' @param extendObject logical, whether net and tarReg will be constructed.
#' Default is FALSE, which can save memory. This can be set FALSE for
#' construting 'network' slot
#' in RegenrichSet objects, but must be set TRUE when construting
#' 'topNetP' slot.
#' @return an object of topNetwork class.
#' @include globals.R
#' @export
#' @examples
#' edge = data.frame(from = rep(TFs$TF_name[seq_len(3)], seq_len(3)),
#'                   to = TFs$TF_name[11:16], weight = 0.1*(6:1))
#' object = TopNetwork(edge, networkConstruction = 'new', percent = 100)
TopNetwork = function(networkEdgeTable, reg = TFs$TF_name, directed = TRUE,
    networkConstruction = c("new", "COEN", "GRN"), percent = 100,
    extendObject = FALSE) {
    stopifnot(is.data.frame(networkEdgeTable) && ncol(networkEdgeTable) ==
        3)
    if (nrow(networkEdgeTable) == 0) {
        stop("The number of rows in networkEdgeTable must be not 0.")
    }
    if (!is.numeric(networkEdgeTable[, 3])) {
        stop("The third column ('weight') of networkEdgeTable must be numeric.")
    }
    colnames(networkEdgeTable) = c("from.gene", "to.gene", "weight")
    # retain edges which start from regulators
    networkEdgeTable = networkEdgeTable[networkEdgeTable$from.gene %in%
        reg, , drop = FALSE]
    reg = unique(as.character(networkEdgeTable$from.gene))

    # remove duplicated edges
    networkEdgeTable = networkEdgeTable[!duplicated(networkEdgeTable[,
        seq_len(2)]), , drop = FALSE]

    # convert potential factor to character
    networkEdgeTable$from.gene = as.character(networkEdgeTable$from.gene)
    networkEdgeTable$to.gene = as.character(networkEdgeTable$to.gene)

    networkConstruction = match.arg(networkConstruction)


    if (extendObject) {
        net = lapply(stats::setNames(reg, reg), function(x) {
            networkEdgeTable[networkEdgeTable$from.gene == x,
                "to.gene"]
        })

        .tarReg = function(nt = networkEdgeTable) {
            nt = nt[, 2:1]
            tar = unique(nt[, "to.gene"])
            res = lapply(stats::setNames(tar, tar), function(x) {
                nt[nt$to.gene == x, "from.gene"]
            })
            return(res)
        }
        tarReg = .tarReg(networkEdgeTable)
    } else {
        net = list()
        tarReg = list()
    }

    new("TopNetwork", net = net, netTable = networkEdgeTable,
        validReg = reg, tarReg = tarReg, directed = directed,
        networkConstruction = networkConstruction, percent = percent)
}


#' Enrich class
#'
#' The `Enrich` object is to store enrichment analysis results by either
#' `FET` method or `GSEA` method.
#'
#' @slot topResult data frame. The enrichment results that pass thresholds
#' (default threshold is 0.05).
#' @slot allResult data frame. The enrichment results by FET or GSEA methods.
#' @slot gene character vector indicating the genes used for enrichment
#' analysis.
#' @slot namedScores numeric vector, a vector of ranked scores (decendent),
#' the names of the scores are the genes to perform enrichment analysis.
#' Here the scores are p-value of each gene.
#' @slot type character indicating enrichment method, either 'FET' or 'GSEA'.
#' @export
setClass("Enrich", slots = c(topResult = "data.frame", allResult = "data.frame",
    gene = "vector", namedScores = "vector", type = "character"  # 'FET'/'SEA'
))



#' RegenrichSet class
#'
#' The `RegenrichSet` is the fundamental class that RegEnrich package
#' is working with.\cr
#' It contains two slots storing raw expression data
#' (@@rawData) and sample pheno data (@@phenoData), which remain not
#' changed in the whole RegEnrich analysis. \cr
#' In addition, it has two slots for storing input (@@paramsIn) and
#' output (@@paramsOut) parameters. \cr
#' Lastly, there are 6 slots (@@assayData, @@resDEA
#' , @@network, @@topNetP, @@resEnrich, @@resScore) storing different
#' parts of results during RegEnrich analysis.
#'
#' @slot rawData the raw expression data. This slot is used for
#' differential expression analysis, meanwhile generating @@assayData
#' slot.
#' @slot phenoData data frame indicating sample information. Each row
#' represent a sample and each column represent a feature of samples.
#' @slot assayData a slot for saving the expression data after filtering
#' (and after Variance Stabilizing Transformation, i.e. VST, if the
#' differential analysis method is 'Wald_DESeq2' or 'LRT_DESeq2'). This
#' slot is used for network inference and plotting.
#' @slot resDEA a slot for saving results after differential analysis.
#' After differential expression analysis, a \code{\link{DeaSet}}
#' object is assigned to this slot. Default is NULL.
#' @slot topNetP a slot for saving top network edges.
#' After regulator-target network inference, a \code{\link{TopNetwork-class}}
#' object is assigned to this slot, containing only top ranked edges
#' in the full network. Default is NULL.
#' @slot resEnrich a slot for saving enrichment analysis by Fisher's
#' exact test (FET) or gene set enrichment analysis (GSEA). It is a
#' \code{\link{Enrich-class}} object.
#' @slot resScore a slot for saving regulator ranking results, which is
#' a `regEnrichScore` (also `data.frame`) object. It contains five columns,
#' which are 'reg' (regulator), 'negLogPDEA' (-log10(p values of differential
#' expression analysis)), 'negLogPEnrich' (-log10(p values of enrichment
#' analysis),
#' 'logFC' (log2 fold changes), and 'score' (RegEnrich ranking score).
#' @slot paramsIn list. The parameters used in the whole RegEnrich
#' analysis. This slot can be updated by respecifying arguments in each step of
#' RegEnrich analysis.
#' @slot paramsOut a list of four elements: method (differential expression
#' method), network (regulator-target network construction method),
#' percent (what percentage of edges from the full network is used),
#' and enrichTest (enrichment method). By default, each element is NULL.
#' @slot network a slot for saving a full network.
#' Before performing network inference, this slot is NULL (default),
#' and after network inference, a \code{\link{TopNetwork-class}} object is
#' assigned to this slot.
#'
#' @export
#'
setClass("RegenrichSet", slots = c(rawData = "ANY", phenoData = "data.frame",
    assayData = "matrix", resDEA = "ANY", topNetP = "ANY", resEnrich = "ANY",
    resScore = "ANY", paramsIn = "list", paramsOut = "ANY", network = "ANY"),
    prototype = list(rawData = NULL, phenoData = data.frame(),
        assayData = matrix(), resDEA = NULL, topNetP = NULL,
        resEnrich = NULL, resScore = NULL, paramsIn = list(),
        paramsOut = list(method = NULL, network = NULL, percent = NULL,
            enrichTest = NULL), network = NULL))


#' RegenrichSet object creator
#'
#' This is `RegenrichSet` object creator function.
#' There are four types of parameters in this function.\cr
#' First, parameters to provide raw data and sample information;\cr
#' `expr` and `pData`.\cr\cr
#' Second, parameters to perform differential expression analysis;\cr
#' `method`, `minMeanExpr`, `design`, `reduced`, `contrast`,
#' `coef`, `name`, `fitType`, `sfType`, `betaPrior`, `minReplicatesForReplace`,
#' `useT`, `minmu`, `parallel`, `BPPARAM` (also for network inference),
#' `altHypothesis`, `listValues`, `cooksCutoff`, `independentFiltering`,
#' `alpha`, `filter`, `theta`, `filterFun`, `addMLE`, `blind`, `ndups`,
#' `spacing`, `block`, `correlation`, `weights`, `proportion`,
#' `stdev.coef.lim`, `trend`, `robust`, and `winsor.tail.p`.\cr\cr
#' Thrid, parameters to perform regulator-target network inference;\cr
#' `reg`, `networkConstruction`, `topNetPercent`, `directed`, `rowSample`,
#' `softPower`, `networkType`, `TOMDenom`, `RsquaredCut`, `edgeThreshold`,
#' `K`, `nbTrees`, `importanceMeasure`, `trace`,
#' `BPPARAM` (also for  differential expression analysis), and `minR`.\cr\cr
#' Fourth, parameters to perform enrichment analysis:\cr
#' `enrichTest`, `namedScoresCutoffs`, `minSize`, `maxSize`, `pvalueCutoff`,
#' `qvalueCutoff`, `regAltName`, `universe`, and `nperm`.\cr\cr
#'
#' @param expr matrix or data.frame, expression profile of a set of
#' genes or a set of proteins. If the \code{method = 'Wald_DESeq2' or
#' 'LRT_DESeq2'}
#' only non-negative integer matrix (read counts by RNA sequencing) is
#' accepted.
#' @param pData data.frame, sample phenotype data. The rows of pData must
#' correspond
#' to the columns of expr.
#' @param method either 'Wald_DESeq2', 'LRT_DESeq2', 'limma', or 'LRT_LM'
#' for the differential expression analysis.
#' \itemize{
#' \item When method = 'Wald_DESeq2', the Wald test in DESeq2 package is used;
#' \item When method = 'LRT_DESeq2', the likelihood ratio test (LRT) in DESeq2
#' package is used;
#' \item When method = 'limma', the `ls` method and empirical Bayes method in
#' limma package are used to calculate moderated t-statistics and differential
#' p-values;
#' \item When method = 'LRT_LM', a likelihood ratio test is performed for each
#' row of `expr` to compare two linear model specified by `design` and
#' `reduced` arguments. In this case, the fold changes are not calculated
#' but set to 0.
#' }
#' @param minMeanExpr numeric, the cutoff of gene average expression for
#' pre-filtering. The rows of `expr` with everage expression < minMeanExpr is
#' removed. The higher `minMeanExpr` is, the more genes are not included for
#' testing.
#' @param design either model formula or model matrix. For method =
#' 'LRT_DESeq2' or
#' 'LRT_LM', the design is the full model formula/matrix. For method =
#' 'limma',
#' and if design is a formula, the model matrix is constructed using
#' model.matrix(design, pData), so the name of each term in the design
#' formula must
#' be included in the column names of `pData`.
#' @param reduced The argument is used only when method = 'LRT_DESeq2' or
#' 'LRT_LM', it is a reduced formula/matrix to compare against.
#' If the design is a model matrix, `reduced` must also be a model matrix.
#' @param contrast The argument is used only when method = 'LRT_DESeq2',
#' 'Wald_DESeq2', or 'limma'. \cr
#' When method = 'LRT_DESeq2', or 'Wald_DESeq2', it specifies what comparison
#' to extract from the `DESeqDataSet` object to build a results table
#' (when method = 'LRT_DESeq2', this does not affect the value of `stat`,
#' `pvalue`, or `padj`). \cr
#' It can be one of following three formats:
#' \itemize{
#' \item a character vector with exactly three elements: the name of a
#' factor in the
#' design formula, the name of the numerator level for the fold change,
#' and the
#' name of the denominator level for the fold change;
#' \item a list of 1 or 2 character vector(s): the first element specifies
#' the names
#' of the fold changes for the numerator, and the second element (optional)
#' specifies the
#' names of the fold changes for the denominator. These names should be
#' elements
#' of \code{getResultsNames(design, pData)};
#' \item a numeric contrast vector with one element for each element in
#' \code{getResultsNames(design, pData)}.\cr
#' }
#'
#' When method = 'limma', It can be one of following two formats:
#' \itemize{
#' \item a numeric matrix with rows corresponding to coefficients in
#' design matrix and
#' columns containing contrasts;
#' \item a numeric vector if there is only one contrast. Each element of
#' the vector
#' corresponds to coefficients in design matrix. This is similar to the
#' third
#' format of contrast when method = 'LRT_DESeq2', or 'Wald_DESeq2'.
#' }
#'
#' @param coef The argument is used only when method = 'limma'. (Vector of)
#' column
#' number or column name specifying which coefficient or contrast of the
#' linear model
#' is of interest. Default is NULL.
#' @param name The argument is used only when method = 'LRT_DESeq2' or
#' 'Wald_DESeq2'.
#' the name of the individual effect (coefficient) for building a results
#' table.
#' Use this argument rather than contrast for continuous variables,
#' individual
#' effects or for individual interaction terms. The value provided to
#' name must
#' be an element of \code{getResultsNames(design, pData)}.
#'
#' @param fitType either 'parametric', 'local', or 'mean' for the type of
#' fitting
#' of dispersions to the mean intensity. This argument is used only when
#' method =
#' 'Wald_DESeq2' or 'LRT_DESeq2'. See \code{\link{DESeq}} from DESeq2
#' package
#' for more details. Default is 'parametric'.
#' @param sfType either 'ratio', 'poscounts', or 'iterate' for the type
#' of size
#' factor estimation. This argument is used only when method = either
#' 'Wald_DESeq2' or 'LRT_DESeq2'. See \code{\link{DESeq}} from DESeq2
#' package
#' for more details. Default is 'ratio'.
#' @param betaPrior This argument is used only when method = either
#' 'Wald_DESeq2' or 'LRT_DESeq2'. See \code{\link{DESeq}} from DESeq2
#' package
#' for more details.
#' @param minReplicatesForReplace This argument is used only when method
#' = either
#' 'Wald_DESeq2' or 'LRT_DESeq2'. See \code{\link{DESeq}} from DESeq2
#' package
#' for more details. Default is 7.
#' @param useT This argument is used only when method = either
#' 'Wald_DESeq2' or 'LRT_DESeq2'. See \code{\link{DESeq}} from DESeq2
#' package
#' for more details. Default is FALSE,
#' @param minmu This argument is used only when method = either
#' 'Wald_DESeq2' or 'LRT_DESeq2'. See \code{\link{DESeq}} from DESeq2
#' package
#' for more details. Default is 0.5.
#' @param parallel whether computing (only for differential analysis
#' with method = "Wald_DESeq2" or "LRT_DESeq2") is parallel (default
#' is FALSE).
#' @param BPPARAM parameters for parallel computing (default is
#' \code{bpparam()}).
#' @param altHypothesis = c('greaterAbs', 'lessAbs', 'greater', 'less').
#' This argument is used only when method = either
#' 'Wald_DESeq2' or 'LRT_DESeq2'. See \code{\link{results}} from DESeq2
#' package
#' for more details. Default is 'greaterAbs'.
#' @param listValues This argument is used only when method = either
#' 'Wald_DESeq2' or 'LRT_DESeq2'. See \code{\link{results}} from DESeq2
#' package
#' for more details. Default is c(1, -1),
#' @param cooksCutoff theshold on Cook's distance, such that if one or
#' more
#' samples for a row have a distance higher, the p-value for the row is
#' set to NA.
#' This argument is used only when method = either
#' 'Wald_DESeq2' or 'LRT_DESeq2'. See \code{\link{results}} from DESeq2
#' package
#' for more details.
#' @param independentFiltering logical, whether independent filtering
#' should be
#' applied automatically. This argument is used only when method = either
#' 'Wald_DESeq2' or 'LRT_DESeq2'. See \code{\link{results}} from DESeq2
#' package
#' for more details. Default is TRUE.
#' @param alpha the significance cutoff used for optimizing the independent
#' filtering.
#' This argument is used only when method = either
#' 'Wald_DESeq2' or 'LRT_DESeq2'. See \code{\link{results}} from DESeq2
#' package
#' for more details. Default is 0.1,
#' @param filter the vector of filter statistics over which the independent
#' filtering is optimized. By default the mean of normalized counts is used.
#' This argument is used only when method = either
#' 'Wald_DESeq2' or 'LRT_DESeq2'. See \code{\link{results}} from DESeq2
#' package
#' for more details.
#' @param theta the quantiles at which to assess the number of rejections
#' from
#' independent filtering. This argument is used only when method = either
#' 'Wald_DESeq2' or 'LRT_DESeq2'. See \code{\link{results}} from DESeq2
#' package
#' for more details.
#' @param filterFun an optional custom function for performing independent
#' filtering
#' and p-value adjustment. This argument is used only when method = either
#' 'Wald_DESeq2' or 'LRT_DESeq2'. See \code{\link{results}} from DESeq2
#' package
#' for more details.
#' @param addMLE if betaPrior=TRUE was used, whether the 'unshrunken' maximum
#' likelihood estimates (MLE) of log2 fold change should be added as a column
#' to the results table. This argument is used only when method = either
#' 'Wald_DESeq2' or 'LRT_DESeq2'. See \code{\link{results}} from DESeq2
#' package
#' for more details. Default is FALSE.
#' @param blind logical, whether to blind the transformation to the
#' experimental
#' design. This argument is used only when method = either
#' 'Wald_DESeq2' or 'LRT_DESeq2'. See \code{\link{vst}} from DESeq2 package
#' for
#' more details. Default is FALSE, which is different from the default of
#' vst function.
#'
#' @param ndups positive integer giving the number of times each distinct
#' probe is
#' printed on each array. This argument is used only when method = 'limma'.
#' See \code{\link{lmFit}} from limma package for more details. Default is 1.
#' @param spacing positive integer giving the spacing between duplicate
#' occurrences of
#' the same probe, spacing=1 for consecutive rows. This argument is used only
#' when method = 'limma'. See \code{\link{lmFit}} from limma package for
#' more details. Default is 1.
#' @param block vector or factor specifying a blocking variable on the arrays.
#' Has length equal to the number of arrays. Must be NULL if ndups > 2.
#' This argument is used only when method = 'limma'. See \code{\link{lmFit}}
#' from limma package for more details. Default is NULL.
#' @param correlation the inter-duplicate or inter-technical replicate
#' correlation.
#' The correlation value should be estimated using the
#' \code{\link{duplicateCorrelation}}
#' function. This argument is used only when method = 'limma'.
#' See \code{\link{lmFit}}
#' from limma package for more details.
#' @param weights non-negative precision weights. Can be a numeric matrix of
#' individual weights of same size as the object expression matrix, or a
#' numeric
#' vector of array weights with length equal to ncol of the expression matrix,
#' or a numeric vector of gene weights with length equal to nrow of the
#' expression
#' matrix. This argument is used only when method = 'limma' or 'LRT_LM'.
#' See \code{\link{lmFit}} from limma package for more details. Default
#' is NULL.
#' @param proportion numeric value between 0 and 1, assumed proportion of
#' genes which
#' are differentially expressed. This argument is used only when method =
#' 'limma'.
#' See \code{\link{eBayes}} from limma package for more details. Default is
#' 0.01.
#' @param stdev.coef.lim numeric vector of length 2, assumed lower and
#' upper limits
#' for the standard deviation of log2-fold-changes for differentially
#' expressed
#' genes. This argument is used only when method = 'limma'.
#' See \code{\link{eBayes}}
#' from limma package for more details. Default is c(0.1, 4).
#' @param trend logical, should an intensity-trend be allowed for the prior
#' variance?
#' This argument is used only when method = 'limma'. See \code{\link{eBayes}}
#' from limma package for more details. Default is FALSE, meaning that the
#' prior
#' variance is constant.
#' @param robust logical, should the estimation of df.prior and var.prior be
#' robustified against outlier sample variances? This argument is used only
#' when method = 'limma'. See \code{\link{eBayes}}
#' from limma package for more details. Default is FALSE.
#' @param winsor.tail.p numeric vector of length 1 or 2, giving left and right
#' tail proportions of x to Winsorize. Used only when method = 'limma' and
#' robust=TRUE. See \code{\link{eBayes}}
#' from limma package for more details. Default is c(0.05,0.1)
#'
#' @param reg a vector of regulator names (ID). By default, these are
#' transcription
#' (co-)factors defined by three literatures/databases, namely RegNet,
#' TRRUST, and Marbach2016. The type (for example ENSEMBL gene ID, Entrez
#' gene ID,
#' or gene symble/name) of names or IDs of these regulators must be the
#' same as the type of names or IDs in the regulator-target network.
#' @param networkConstruction the method to construct this network.
#' Possible can be:\cr
#' 'COEN', coexpression network;\cr
#' 'GRN', gene regulatory network by random forest;\cr
#' 'new' (default), meaning a network provided by user, rather than
#' infered based
#' on the expression data.\cr
#' @param topNetPercent numeric, what percentage of the top edges in the
#' full
#' network is ratained. Default is 5, meaning top 5\% of edges. This value
#' must
#' be between 0 and 100.
#' @param directed logical, whether the network is directed. Default is
#' FALSE.
#' @param rowSample logic, if TRUE, each row represents a sample.
#' Otherwise, each column represents a sample. Default is FALSE.
#' @param softPower numeric, a soft power to achieve scale free topology.
#' If not provided, the parameter will be picked automatically by
#' \code{\link{plotSoftPower}} function from the WGCNA package.
#' @param networkType network type. Allowed values are (unique abbreviations
#' of)
#' 'unsigned' (default), 'signed', 'signed hybrid'.
#' See \code{\link{adjacency}}.
#' @param TOMDenom a character string specifying the TOM variant to be used.
#' Recognized values are 'min' giving the standard TOM described in Zhang
#' and Horvath (2005), and 'mean' in which the min function in the
#' denominator is replaced by mean. The 'mean' may produce better results
#' but at this time should be considered experimental.
#' @param RsquaredCut desired minimum scale free topology fitting index R^2.
#' Default is 0.85.
#' @param edgeThreshold numeric, the threshold to remove the low weighted
#' edges, Default is NULL, which means no edges will be removed.
#'
#' @param K integer or character. The number of features in each tree,
#' can be either a integer number, `sqrt`, or `all`.
#' `sqrt` denotes sqrt(the number of `reg`), `all`
#' means the number of `reg`. Default is `sqrt`.
#' @param nbTrees integer. The number of trees. Default is 1000.
#' @param importanceMeasure character. importanceMeasure can be `\%IncMSE`
#' or `IncNodePurity`, corresponding to type = 1 and 2 in
#' \code{\link{importance}}
#' function, respectively. Default is `IncNodePurity`(decrease in node
#' impurity),
#' which is faster than `\%IncMSE` (decrease in accuracy).
#' @param trace logical. To show the progress or not (default).
#' @param minR numeric. The minimum correlation coefficient of
#' prediction is to
#' control model accuracy. Default is 0.3.
#'
#' @param enrichTest character, specifying the enrichment analysis method,
#' which
#' is either `FET` (Fisher's exact test) or `GSEA` (gene set enrichment
#' analysis).
#' @param namedScoresCutoffs numeric, the significance cutoff for the
#' differential
#' analysis p value. Default is 0.05.
#' @param minSize The minimum number (default 5) of target genes.
#' @param maxSize The maximum number (default 5000) of target genes.
#' @param pvalueCutoff numeric, the significance cutoff for adjusted
#' enrichment p value.
#' This is used for obtaining the `topResult` slot in the final `Enrich`
#' object. Default is 0.05.
#' @param qvalueCutoff numeric, the significance cutoff of enrichment
#' q-value.
#' Default is 0.2.
#' @param regAltName alternative name for regulator. Default is NULL.
#' @param universe a vector of charactors. Background target genes.
#' @param nperm integer, number of permutations. The minimial possible
#' nominal p-value is about 1/nperm. Default is 10000.
#' @return an object of RegenrichSet class.
#' @import WGCNA
#' @import BiocParallel
#' @include globals.R
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
RegenrichSet = function(expr, pData, method = c("Wald_DESeq2",
    "LRT_DESeq2", "limma", "LRT_LM"), minMeanExpr = NULL, design,
    reduced, contrast, coef = NULL, name, fitType = c("parametric",
        "local", "mean"), sfType = c("ratio", "poscounts", "iterate"),
    betaPrior, minReplicatesForReplace = 7, useT = FALSE, minmu = 0.5,
    parallel = FALSE, BPPARAM = bpparam(), altHypothesis = c("greaterAbs",
        "lessAbs", "greater", "less"), listValues = c(1, -1),
    cooksCutoff, independentFiltering = TRUE, alpha = 0.1, filter,
    theta, filterFun, addMLE = FALSE, blind = FALSE, ndups = 1,
    spacing = 1, block = NULL, correlation, weights = NULL, proportion = 0.01,
    stdev.coef.lim = c(0.1, 4), trend = FALSE, robust = FALSE,
    winsor.tail.p = c(0.05, 0.1), reg = TFs$TF_name,
    networkConstruction = c("COEN",
        "GRN", "new"), topNetPercent = 5, directed = FALSE, rowSample = FALSE,
    softPower = NULL, networkType = "unsigned", TOMDenom = "min",
    RsquaredCut = 0.85, edgeThreshold = NULL, K = "sqrt", nbTrees = 1000,
    importanceMeasure = "IncNodePurity", trace = FALSE,
    minR = 0.3, enrichTest = c("FET", "GSEA"), namedScoresCutoffs = 0.05,
    minSize = 5, maxSize = 5000, pvalueCutoff = 0.05, qvalueCutoff = 0.2,
    regAltName = NULL, universe = NULL, nperm = 10000) {

    if (missing(expr) || missing(pData)) {
        stop("Must provide both 'expr' and 'pData'")
    }

    if (ncol(expr) != nrow(pData)) {
        stop("The number of columns in 'expr' must be ",
            "the same as the number of rows in 'pData'")
        if (all(colnames(expr) != rownames(pData))) {
            stop("The columns of 'expr' must match the rows of 'pData'")
        }
    }

    if (!is.null(minMeanExpr) && is.numeric(minMeanExpr)) {
        exprHi = as.matrix(expr[rowMeans(expr) > minMeanExpr,
            , drop = FALSE])
    } else {
        exprHi = as.matrix(expr)
    }

    pData = as.data.frame(pData)

    resDEA = NULL
    network = topNetP = resEnrich = resScore = NULL
    paramsOut = list(method = NULL, network = NULL, enrichTest = NULL,
        percent = NULL)

    object = new(Class = "RegenrichSet", rawData = expr, assayData = exprHi,
        phenoData = pData, paramsIn = list(method = method,
                minMeanExpr = minMeanExpr,
            design = if (missing(design)) substitute() else design,
            reduced = if (missing(reduced)) substitute() else reduced,
            contrast = if (missing(contrast)) substitute() else contrast,
            coef = if (missing(coef)) substitute() else coef,
            name = if (missing(name)) substitute() else name,
            fitType = fitType, sfType = sfType,
            betaPrior = if (missing(betaPrior)) substitute() else betaPrior,
            minReplicatesForReplace = minReplicatesForReplace,
            useT = useT, minmu = minmu, parallel = parallel,
            BPPARAM = BPPARAM, altHypothesis = altHypothesis,
            listValues = listValues,
            cooksCutoff = if (missing(cooksCutoff)) {
                substitute()} else cooksCutoff,
            independentFiltering = independentFiltering, alpha = alpha,
            filter = if (missing(filter)) substitute() else filter,
            theta = if (missing(theta)) substitute() else theta,
            filterFun = if (missing(filterFun)) substitute() else filterFun,
            addMLE = addMLE, blind = blind, ndups = ndups, spacing = spacing,
            block = block,
            correlation = if (missing(correlation)) {
                substitute()} else correlation,
            weights = weights, proportion = proportion,
            stdev.coef.lim = stdev.coef.lim,
            trend = trend, robust = robust, winsor.tail.p = winsor.tail.p,
            reg = reg, rowSample = rowSample, softPower = softPower,
            networkConstruction = networkConstruction,
            networkType = networkType,
            TOMDenom = TOMDenom, RsquaredCut = RsquaredCut,
            edgeThreshold = edgeThreshold,
            K = K, nbTrees = nbTrees, importanceMeasure = importanceMeasure,
            trace = trace,
            minR = minR, topNetPercent = topNetPercent, directed = directed,
            enrichTest = enrichTest, namedScoresCutoffs = namedScoresCutoffs,
            minSize = minSize, maxSize = maxSize, pvalueCutoff = pvalueCutoff,
            qvalueCutoff = qvalueCutoff, regAltName = regAltName,
            universe = universe, nperm = nperm), resDEA = resDEA,
        network = network, topNetP = topNetP, resEnrich = resEnrich,
        resScore = resScore, paramsOut = paramsOut)
    return(object)
}

