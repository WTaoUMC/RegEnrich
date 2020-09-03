# Differential expression analysis
# @description DEA is a wraper function using DESeq2 and limma package
# to perform differential expression analysis.
# @param expr matrix or data.frame, expression profile of a set of
# genes or a set of proteins. If the \code{method = 'Wald_DESeq2' or
# 'LRT_DESeq2'} only non-negative integer matrix (read counts by RNA
# sequencing) is accepted.
# @param pData data.frame, sample phenotype data. The rows of pData must
# correspond to the columns of expr.
# @param method either 'Wald_DESeq2', 'LRT_DESeq2', 'limma', or 'LRT_LM' for
# the differential expression analysis.
# \itemize{
# \item When method = 'Wald_DESeq2', the Wald test in DESeq2 package is used;
# \item When method = 'LRT_DESeq2', the likelihood ratio test (LRT) in DESeq2
# package is used;
# \item When method = 'limma', the `ls` method and empirical Bayes method in
# limma package is used to calculate moderated t-statistics and differential
# p-values;
# \item When method = 'LRT_LM', a likelihood ratio test is performed for each
# row of `expr` to compare two linear model specified by `design` and
# `reduced` arguments. In this case, the fold changes are not calculated
# but set to 0.
# }
# @param minMeanExpr numeric, the cutoff of gene average expression for
# pre-filtering. The rows of `expr` with everage expression < minMeanExpr is
# removed. The higher `minMeanExpr` is, the more genes are not included for
# testing.
# @param design either model formula or model matrix.
# For method = 'LRT_DESeq2' or 'LRT_LM', the design is the full model
# formula/matrix. For method = 'limma',
# and if design is a formula, the model matrix is constructed using
# model.matrix(design, pData), so the name of each term in the design
# formula must
# be included in the column names of `pData`.
# @param reduced The argument is used only when method = 'LRT_DESeq2' or
# 'LRT_LM', it is a reduced formula/matrix to compare against.
# If the design is a model matrix, `reduced` must also be a model matrix.
# @param contrast The argument is used only when method = 'LRT_DESeq2',
# 'Wald_DESeq2', or 'limma'. \cr
# When method = 'LRT_DESeq2', or 'Wald_DESeq2', it specifies what comparison
# to extract from the `DESeqDataSet` object to build a results table
# (when method = 'LRT_DESeq2', this does not affect the value of `stat`,
# `pvalue`, or `padj`). \cr
# It can be one of following three formats:
# \itemize{
# \item a character vector with exactly three elements: the name of a
# factor in the design formula, the name of the numerator level for
# the fold change, and the name of the denominator level for the fold
# change;
# \item a list of 1 or 2 character vector(s): the first element specifies
# the names of the fold changes for the numerator, and the second
# element (optional) specifies the
# names of the fold changes for the denominator. These names should be
# elements
# of \code{getResultsNames(design, pData)};
# \item a numeric contrast vector with one element for each element in
# \code{getResultsNames(design, pData)}.\cr
# }
# When method = 'limma', It can be one of following two formats:
# \itemize{
# \item a numeric matrix with rows corresponding to coefficients in design
# matrix and columns containing contrasts;
# \item a numeric vector if there is only one contrast. Each element of the
# vector corresponds to coefficients in design matrix. This is similar to
# the third format of contrast when method = 'LRT_DESeq2', or 'Wald_DESeq2'.
# }
#
# @param coef The argument is used only when method = 'limma'. (Vector of)
# column number or column name specifying which coefficient or contrast of
# the linear model is of interest. Default is NULL.
# @param name The argument is used only when method = 'LRT_DESeq2' or
# 'Wald_DESeq2'. the name of the individual effect (coefficient) for
# building a results table.
# Use this argument rather than contrast for continuous variables, individual
# effects or for individual interaction terms. The value provided to name must
# be an element of \code{getResultsNames(design, pData)}.
#
# @param fitType either 'parametric', 'local', or 'mean' for the type of
# fitting of dispersions to the mean intensity. This argument is used only
# when method = 'Wald_DESeq2' or 'LRT_DESeq2'. See \code{\link{DESeq}} from
# DESeq2 package for more details. Default is 'parametric'.
# @param sfType either 'ratio', 'poscounts', or 'iterate' for the type of size
# factor estimation. This argument is used only when method = either
# 'Wald_DESeq2' or 'LRT_DESeq2'. See \code{\link{DESeq}} from DESeq2 package
# for more details. Default is 'ratio'.
# @param betaPrior This argument is used only when method = either
# 'Wald_DESeq2' or 'LRT_DESeq2'. See \code{\link{DESeq}} from DESeq2 package
# for more details.
# @param minReplicatesForReplace This argument is used only when
# method = either 'Wald_DESeq2' or 'LRT_DESeq2'. See \code{\link{DESeq}}
# from DESeq2 package
# for more details. Default is 7.
# @param useT This argument is used only when method = either
# 'Wald_DESeq2' or 'LRT_DESeq2'. See \code{\link{DESeq}} from DESeq2 package
# for more details. Default is FALSE,
# @param minmu This argument is used only when method = either
# 'Wald_DESeq2' or 'LRT_DESeq2'. See \code{\link{DESeq}} from DESeq2 package
# for more details. Default is 0.5.
# @param parallel whether computing is parallel (default is FALSE).
# @param BPPARAM parameters for parallel computing (default is
# \code{bpparam()}).
# @param altHypothesis = c('greaterAbs', 'lessAbs', 'greater', 'less').
# This argument is used only when method = either
# 'Wald_DESeq2' or 'LRT_DESeq2'. See \code{\link{results}} from DESeq2 package
# for more details. Default is 'greaterAbs'.
# @param listValues This argument is used only when method = either
# 'Wald_DESeq2' or 'LRT_DESeq2'. See \code{\link{results}} from DESeq2 package
# for more details. Default is c(1, -1),
# @param cooksCutoff theshold on Cook's distance, such that if one or more
# samples for a row have a distance higher, the p-value for the row is set
# to NA. This argument is used only when method = either
# 'Wald_DESeq2' or 'LRT_DESeq2'. See \code{\link{results}} from DESeq2 package
# for more details.
# @param independentFiltering logical, whether independent filtering should be
# applied automatically. This argument is used only when method = either
# 'Wald_DESeq2' or 'LRT_DESeq2'. See \code{\link{results}} from DESeq2 package
# for more details. Default is TRUE.
# @param alpha the significance cutoff used for optimizing the independent
# filtering.  This argument is used only when method = either
# 'Wald_DESeq2' or 'LRT_DESeq2'. See \code{\link{results}} from DESeq2 package
# for more details. Default is 0.1,
# @param filter the vector of filter statistics over which the independent
# filtering is optimized. By default the mean of normalized counts is used.
# This argument is used only when method = either
# 'Wald_DESeq2' or 'LRT_DESeq2'. See \code{\link{results}} from DESeq2 package
# for more details.
# @param theta the quantiles at which to assess the number of rejections from
# independent filtering. This argument is used only when method = either
# 'Wald_DESeq2' or 'LRT_DESeq2'. See \code{\link{results}} from DESeq2 package
# for more details.
# @param filterFun an optional custom function for performing independent
# filtering and p-value adjustment. This argument is used only when
# method = either 'Wald_DESeq2' or 'LRT_DESeq2'. See \code{\link{results}}
# from DESeq2 package for more details.
# @param addMLE if betaPrior=TRUE was used, whether the 'unshrunken' maximum
# likelihood estimates (MLE) of log2 fold change should be added as a column
# to the results table. This argument is used only when method = either
# 'Wald_DESeq2' or 'LRT_DESeq2'. See \code{\link{results}} from DESeq2 package
# for more details. Default is FALSE.
# @param blind logical, whether to blind the transformation to the
# experimental design. This argument is used only when method = either
# 'Wald_DESeq2' or 'LRT_DESeq2'. See \code{\link{vst}} from DESeq2 package for
# more details. Default is FALSE, which is different from the default of vst
# function.
#
# @param ndups positive integer giving the number of times each distinct
# probe is printed on each array. This argument is used only when
# method = 'limma'. See \code{\link{lmFit}} from limma package for more
# details. Default is 1.
# @param spacing positive integer giving the spacing between duplicate
# occurrences of the same probe, spacing=1 for consecutive rows.
# This argument is used only
# when method = 'limma'. See \code{\link{lmFit}} from limma package for
# more details. Default is 1.
# @param block vector or factor specifying a blocking variable on the arrays.
# Has length equal to the number of arrays. Must be NULL if ndups > 2.
# This argument is used only when method = 'limma'. See \code{\link{lmFit}}
# from limma package for more details. Default is NULL.
# @param correlation the inter-duplicate or inter-technical replicate
#  correlation. The correlation value should be estimated using
# the \code{\link{duplicateCorrelation}}
# function. This argument is used only when method = 'limma'.
# See \code{\link{lmFit}}
# from limma package for more details.
# @param weights non-negative precision weights. Can be a numeric matrix of
# individual weights of same size as the object expression matrix,
#  or a numeric
# vector of array weights with length equal to ncol of the expression matrix,
# or a numeric vector of gene weights with length equal to nrow of the
# expression
# matrix. This argument is used only when method = 'limma' or 'LRT_LM'.
# See \code{\link{lmFit}} from limma package for more details.
# Default is NULL.
# @param proportion numeric value between 0 and 1, assumed proportion of
# genes which
# are differentially expressed. This argument is used only when
# method = 'limma'.
# See \code{\link{eBayes}} from limma package for more details.
# Default is 0.01.
# @param stdev.coef.lim numeric vector of length 2,
# assumed lower and upper limits
# for the standard deviation of log2-fold-changes for differentially expressed
# genes. This argument is used only when method = 'limma'.
# See \code{\link{eBayes}}
# from limma package for more details. Default is c(0.1, 4).
# @param trend logical, should an intensity-trend be allowed for the
# prior variance?
# This argument is used only when method = 'limma'. See \code{\link{eBayes}}
# from limma package for more details. Default is FALSE, meaning that
# the prior variance is constant.
# @param robust logical, should the estimation of df.prior and var.prior be
# robustified against outlier sample variances? This argument is used only
# when method = 'limma'. See \code{\link{eBayes}}
# from limma package for more details. Default is FALSE.
# @param winsor.tail.p numeric vector of length 1 or 2, giving left and right
# tail proportions of x to Winsorize. Used only when method = 'limma' and
# robust=TRUE. See \code{\link{eBayes}}
# from limma package for more details. Default is c(0.05,0.1)
#' @include globals.R
#' @rawNamespace import(limma, except = plotMA)
#' @rawNamespace import(DESeq2, except = plotMA)
#' @import BiocParallel
# @importFrom SummarizedExperiment assay
#' @rawNamespace import(S4Vectors, except = c(cor, second, first, union,
#' intersect, setdiff, setequal, rename))
# @seealso \code{\link{lmFit}}, \code{\link{contrasts.fit}},
# \code{\link{eBayes}},  \code{\link{DESeq}}, and \code{\link{results}}
# @return a list of 3 elements: pFC (a data frame of gene, p value, and
# log2FC), fullRes (full version of differential analysis),
# fit (fitted model), and normExprHigh (expression data).
# @export
DEA = function(expr, pData, method = c("Wald_DESeq2", "LRT_DESeq2",
    "limma", "LRT_LM"), minMeanExpr = NULL, design, reduced,
    contrast, coef = NULL, name, fitType = c("parametric", "local",
        "mean"), sfType = c("ratio", "poscounts", "iterate"),
    betaPrior, minReplicatesForReplace = 7, useT = FALSE, minmu = 0.5,
    parallel = FALSE, BPPARAM = bpparam(), altHypothesis = c("greaterAbs",
        "lessAbs", "greater", "less"), listValues = c(1, -1),
    cooksCutoff, independentFiltering = TRUE, alpha = 0.1, filter,
    theta, filterFun, addMLE = FALSE, blind = FALSE, ndups = 1,
    spacing = 1, block = NULL, correlation, weights = NULL, proportion = 0.01,
    stdev.coef.lim = c(0.1, 4), trend = FALSE, robust = FALSE,
    winsor.tail.p = c(0.05, 0.1)) {
    stopifnot(is.matrix(expr) || is.data.frame(expr))
    stopifnot(colnames(expr) == rownames(pData))
    if (missing(design)) {
        stop("A model matrix or formula must ",
            "be assigned to 'design' argument")
    }
    stopifnot(is(design, "formula") || is.matrix(design))

    method = match.arg(method)
    if (method %in% c("LRT_DESeq2", "LRT_LM")) {
        if (missing(reduced)) {
            stop("When method is LRT_DESeq2 or LRT_LM, a ",
                "formula or model matrix must be",
                " assigned to 'reduced' argument.")
        }
        stopifnot(is(reduced, "formula") || is.matrix(reduced))
    } else if (method %in% c("limma", "Wald_DESeq2")) {
        if (!missing(reduced)) {
            message("When method is limma or Wald_DESeq2, ",
                "the 'reduced' argument is not used.\n")
        }
        if (method == "limma") {
            if (is(design, "formula")) {
                message("The method is limma, a formula is ",
                    "provided to the design,\n",
                    "so the model matrix is generated using ",
                    "model.matrix(design, data = pData).")
                design = stats::model.matrix(design, data = pData)
            }
        }
    }

    # pre-filtering low expression genes
    if (!is.null(minMeanExpr)) {
        expr = expr[rowMeans(expr) >= minMeanExpr, , drop = FALSE]
        if (nrow(expr) < 1)
            stop("No gene pass the minMeanExpr threshold.")
    }

    if (method == "Wald_DESeq2") {
        res = DEA_LRT_DESeq2(expr = as.matrix(expr), pData = pData,
            test = "Wald", full = design, contrast = contrast,
            name = name, fitType = fitType, sfType = sfType,
            betaPrior = betaPrior,
            minReplicatesForReplace = minReplicatesForReplace,
            useT = useT, minmu = minmu, parallel = parallel,
            BPPARAM = BPPARAM, altHypothesis = altHypothesis,
            listValues = listValues, cooksCutoff = cooksCutoff,
            independentFiltering = independentFiltering, alpha = alpha,
            filter = filter, theta = theta, filterFun = filterFun,
            addMLE = addMLE, blind = blind)
    } else if (method == "LRT_DESeq2") {
        res = DEA_LRT_DESeq2(expr = as.matrix(expr), pData = pData,
            test = "LRT", full = design, reduced = reduced,
            contrast = contrast,
            name = name, fitType = fitType, sfType = sfType,
            betaPrior = betaPrior,
            minReplicatesForReplace = minReplicatesForReplace,
            useT = useT, minmu = minmu, parallel = parallel,
            BPPARAM = BPPARAM, altHypothesis = altHypothesis,
            listValues = listValues, cooksCutoff = cooksCutoff,
            independentFiltering = independentFiltering, alpha = alpha,
            filter = filter, theta = theta, filterFun = filterFun,
            addMLE = addMLE, blind = blind)
    } else if (method == "limma") {
        res = DEA_limma(expr = as.matrix(expr), design = design,
            ndups = ndups, spacing = spacing, block = block,
            correlation = correlation, weights = weights, contrast = contrast,
            proportion = proportion, stdev.coef.lim = stdev.coef.lim,
            trend = trend, robust = robust, winsor.tail.p = winsor.tail.p,
            coef = coef, number = nrow(expr), genelist = rownames(expr),
            adjust.method = "BH", sort.by = "none", resort.by = NULL,
            p.value = 1, lfc = 0, confint = FALSE)
    } else if (method == "LRT_LM") {
        res = DEA_LRT_LM(expr = as.matrix(expr), pData = pData,
            full = design, reduced = reduced, weights = weights)
    } else {
        stop("Unknown 'method': ", method, ".\n'method' must be one of ",
            "\"LRT_DESeq2\", \"limma\", \"Wald_DESeq2\", and \"LRT_LM\".")
    }

    return(res)
}




####### 'LRT_DESeq2' and 'Wald_DESeq2'
DEA_LRT_DESeq2 = function(expr, pData, test = c("Wald", "LRT"),
    fitType = c("parametric", "local", "mean"), sfType = c("ratio",
        "poscounts", "iterate"), betaPrior, full, reduced,
    minReplicatesForReplace = 7,
    useT = FALSE, minmu = 0.5, parallel = FALSE, BPPARAM = bpparam(),
    contrast, name, altHypothesis = c("greaterAbs", "lessAbs",
        "greater", "less"), listValues = c(1, -1), cooksCutoff,
    independentFiltering = TRUE, alpha = 0.1, filter, theta,
    filterFun, addMLE = FALSE, blind = FALSE) {
    
    fitType = match.arg(fitType)
    expr = as.matrix(expr)
    object = DESeqDataSetFromMatrix(countData = expr, colData = pData,
        design = full)

    featureData = data.frame(gene = rownames(expr))
    mcols(object) = S4Vectors::DataFrame(mcols(object), featureData)

    object = DESeq(object = object, test = test, fitType = fitType,
        sfType = sfType, betaPrior = betaPrior, full = full,
        reduced = reduced, quiet = TRUE,
        minReplicatesForReplace = minReplicatesForReplace,
        modelMatrixType = "standard", useT = useT, minmu = minmu,
        parallel = parallel, BPPARAM = BPPARAM)
    res = results(object = object, contrast = contrast, name = name,
        lfcThreshold = 0, altHypothesis = altHypothesis,
        listValues = listValues,
        cooksCutoff = cooksCutoff, independentFiltering = independentFiltering,
        alpha = alpha, filter = filter, theta = theta, pAdjustMethod = "BH",
        filterFun = filterFun, format = "DataFrame", test = test,
        addMLE = addMLE, tidy = FALSE, parallel = parallel, BPPARAM = BPPARAM,
        minmu = minmu)

    res$pvalue[is.na(res$pvalue)] = 1

    pFC = data.frame(gene = rownames(res), p = res$pvalue,
        logFC = res$log2FoldChange, row.names = rownames(res))
    normExprHigh = assay(vst(object, blind = blind, nsub = min(1000,
        nrow(object)), fitType = fitType))
    return(list(pFC = pFC, fullRes = res, fit = NULL,
        normExprHigh = normExprHigh))
}


####### 'limma'
DEA_limma = function(expr, design = NULL, ndups = 1, spacing = 1,
    block = NULL, correlation, weights = NULL, contrast, proportion = 0.01,
    stdev.coef.lim = c(0.1, 4), trend = FALSE, robust = FALSE,
    winsor.tail.p = c(0.05, 0.1), coef = NULL, number = 10, genelist,
    adjust.method = "BH", sort.by = "none", resort.by = NULL,
    p.value = 1, lfc = 0, confint = FALSE) {
    # fit
    fit = lmFit(expr, design = design, ndups = ndups, spacing = spacing,
        block = block, correlation = correlation, weights = weights,
        method = "ls")
    # contrast
    if (!missing(contrast)) {
        fit = contrasts.fit(fit = fit, contrasts = contrast)
    }

    fit = eBayes(fit = fit, proportion = proportion,
        stdev.coef.lim = stdev.coef.lim,
        trend = trend, robust = robust, winsor.tail.p = winsor.tail.p)
    # topTable
    res = topTable(fit, coef = coef, number = number, genelist = genelist,
        adjust.method = adjust.method, sort.by = sort.by,
        resort.by = resort.by,
        p.value = p.value, lfc = lfc, confint = confint)

    pFC = data.frame(gene = rownames(res), p = res[, "P.Value",
        drop = TRUE], logFC = res[, "logFC", drop = TRUE],
        row.names = rownames(res))
    return(list(pFC = pFC, fullRes = res, fit = fit, normExprHigh = expr))
}


######### 'LRT_LM' ####### lm for matrix set.seed(123) weights =
######### matrix(runif(nrow(expr)*ncol(expr), min = 0, max = 1), m =
######### nrow(expr), n = ncol(expr)) DEA_LRT_LM(expr, pData, full,
######### reduced, weights) DEA_LRT_LM(expr, pData, full, reduced,
######### weights = NULL)
DEA_LRT_LM = function(expr, pData, full, reduced, weights = NULL) {
    if (!is.null(weights)) {
        stopifnot(is.matrix(weights))
        stopifnot(identical(dim(weights), dim(expr)))
    }

    if (is(full, "formula")) {
        stopifnot(all(all.vars(full) %in% colnames(pData)))
        x1 = stats::model.matrix(full, pData)
    } else {
        vars = all.vars(full)
        x1 = full
    }
    if (is(reduced, "formula")) {
        stopifnot(all(all.vars(reduced) %in% colnames(pData)))
        x2 = stats::model.matrix(reduced, pData)
    } else {
        x2 = reduced
    }

    z1 = if (is.null(weights)) {
        stats::lm.fit(x1, t(expr), offset = NULL, singular.ok = TRUE)
    } else {
        w = as.numeric(weights)
        stats::lm.wfit(x1, t(expr), w, offset = NULL, singular.ok = TRUE)
    }

    z2 = if (is.null(weights)) {
        stats::lm.fit(x2, t(expr), offset = NULL, singular.ok = TRUE)
    } else {
        w = as.numeric(weights)
        stats::lm.wfit(x2, t(expr), w, offset = NULL, singular.ok = TRUE)
    }

    resdf = c(z1$df.residual, z2$df.residual)
    resdev = cbind(colSums(stats::weighted.residuals(z1)^2, na.rm = TRUE),
        colSums(stats::weighted.residuals(z2)^2, na.rm = TRUE))
    df = -diff(resdf)
    sumOfSq = resdev[, 1] - resdev[, 2]

    scale = resdev[, 1]/resdf[1]
    val = sumOfSq/scale * sign(df)

    if (df == 0)
        val = NA

    val[!is.na(val) & val < 0] = NA
    val[is.na(val)] = NA
    pval = data.frame(p = stats::pchisq(val, abs(df), lower.tail = FALSE),
        lr = val, df = df)

    res = data.frame(baseMean = rowMeans(expr), logFC = 0, lr = pval[,
        "lr"], df = pval[, "df"], pvalue = pval[, "p"],
        padj = stats::p.adjust(pval[,
        "p"], "BH"), row.names = rownames(expr))
    pFC = data.frame(gene = rownames(expr), p = pval[, "p"],
        logFC = 0, row.names = rownames(expr))
    pFC$p[is.na(pFC$p)] = 1

    return(list(pFC = pFC, fullRes = res, fit = NULL, normExprHigh = expr))
}

#' Inference the name of results of DESeq analysis by a formula
#' (or model matrix) and sample information
#' @param design either a formula or a model matrix.
#' @param pData a data frame, showing the information of each sample.
#' If design is
#' a formula, the pData must be include the columns that identical
#' to the terms of
#' the design formula. If design is a model matrix, then pData is not used.
#' Default is NULL.
#' @return the names of contrast parameter (list of character format) that 
#' \code{\link{regenrich_diffExpr}} and \code{\link{results}}
#' function can use, and it is the same as the value that
#' \code{\link{resultsNames}} function returns.
#' @rawNamespace import(DESeq2, except = plotMA)
#' @export
#' @examples
#' # formula with intercept
#' design = ~condition
#' pData = data.frame(condition = factor(c('A', 'A', 'A', 'B', 'B', 'B'),
#'                                       c('A', 'B')))
#' getResultsNames(design, pData)
#'
#' # formula without intercept
#' design = ~0+condition
#' getResultsNames(design, pData)
#'
#' # formula with two terms
#' design = ~condition+treatment
#' pData = data.frame(condition = factor(rep(c('A', 'B'), each= 4),
#'                                       c('A', 'B')),
#'                    treatment = factor(rep_len(c('Ctrl', 'Treat'), 8),
#'                                       c('Ctrl', 'Treat')))
#' getResultsNames(design, pData)
#'
#' # formula with two terms and an interaction term
#' design = ~condition+treatment+condition:treatment
#' getResultsNames(design, pData)
#'
#' # design is a model matrix
#' pData = data.frame(condition = factor(rep(c('A', 'B'), each= 4),
#'                                       c('A', 'B')),
#'                    treatment = factor(rep_len(c('Ctrl', 'Treat'), 8),
#'                                       c('Ctrl', 'Treat')))
#' design = model.matrix(~condition+treatment, pData)
#' getResultsNames(design)
getResultsNames = function(design, pData = NULL) {
    stopifnot(is(design, "formula") || is(design, "matrix"))
    if (is(design, "formula")) {
        if (is.null(pData) || !is.data.frame(pData)) {
            stop("If design is a formula, pData must be a data frame")
        }
        if (!all(all.vars(design) %in% colnames(pData))) {
            stop("The terms in the design formula are not ",
                "all contained in the columns of pData.")
        }
        hasIntercept = attr(stats::terms(design), "intercept") ==
            1
        renameCols = hasIntercept
        modelMatrix = stats::model.matrix(design, data = pData)
        stopifnot(nrow(modelMatrix) == nrow(pData))
    } else {
        modelMatrix = design
        renameCols = FALSE
    }
    stopifnot(all(colSums(abs(modelMatrix)) > 0))
    modelMatrixNames = colnames(modelMatrix)

    modelMatrixNames[modelMatrixNames == "(Intercept)"] = "Intercept"
    modelMatrixNames = make.names(modelMatrixNames)

    if (renameCols) {
        convertNames = renameModelMatrixColumns(pData, design)
        convertNames = convertNames[convertNames$from %in% modelMatrixNames,
            , drop = FALSE]
        modelMatrixNames[match(convertNames$from,
            modelMatrixNames)] = convertNames$to
    }
    return(modelMatrixNames)
}
