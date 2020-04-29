
#' Compare the orders of two vectors
#' @description compare the orders of two vectors
#' @param name1 a vector with first order.
#' @param name2 a vector with anothoer second order.
#' @rawNamespace import(ggplot2, except = margin)
#' @return A plot of comparing two orders of vectors.
#' @examples
#' a = c('a1', 'a2', 'a5', 'a4')
#' b = c( 'a2', 'a5', 'a7', 'a4', 'a6')
#' plotOrders(a, b)
#' @export
plotOrders = function(name1, name2) {
    lab1 = deparse(substitute(name1))
    lab2 = deparse(substitute(name2))
    names = unique(c(name1, name2))
    m = length(name1)
    n = length(name2)
    x = match(name1, names)
    y = match(name2, names)
    xy = data.frame(order = c(seq_along(x), seq_along(y)), rank = c(x,
        y), group = c(rep(lab1, m), rep(lab2, n)))
    ggplot(data = xy, aes_string(x = "group", y = "order")) +
        geom_point() + geom_line(aes(group = rank)) + scale_y_reverse()+
        theme_bw()
}


#' Plot regulator and its targets expression
#' @param object a RegenrichSet object, to which at
#' least \code{\link{regenrich_diffExpr}} and \code{\link{regenrich_network}}
#' functions have been applied.
#' @param reg a regulator to plot.
#' @param n the maximun number of targets to plot.
#' @param scale logical, whether gene expression is z-score normalized.
#' @param tarCol the color of the lines for the targets of the regulator.
#' @param tarColAlpha numeric, ranging from 0 to 1, indicating transparancy
#' of target lines.
#' @param regCol the color of the line for the 'reg'.
#' @param xlab x label of plot.
#' @param ylab y label of plot.
#' @param ... other parameters in \code{\link{ggplot}} function.
#' @return a ggplot object.
#' @importFrom reshape2 melt
#' @export
#' @examples
#' # constructing a RegenrichSet object
#' colData = data.frame(patientID = paste0('Sample_', seq(50)),
#'                      week = rep(c('0', '1'), each = 25),
#'                      row.names = paste0('Sample_', seq(50)), 
#'                      stringsAsFactors = TRUE)
#' design = ~week
#' reduced = ~1
#' set.seed(123)
#' cnts = matrix(as.integer(rnbinom(n=1000*50, mu=100, size=1/0.1)), ncol=50,
#'               dimnames = list(paste0('gene', seq(1000)), rownames(colData)))
#' 
#' cnts[5,26:50] = cnts[5,26:50] + 50L # add reads to gene5 in some samples.
#' id = sample(31:1000, 20) # randomly select 20 rows, and assign reads.
#' cnts[id,] = vapply(cnts[5,], function(x){
#'   as.integer(rnbinom(n = 20, size = 1/0.02, mu = x))},
#'   FUN.VALUE = rep(1L, 20))
#'
#' object = RegenrichSet(expr = cnts,
#'                       colData = colData,
#'                       method = 'LRT_DESeq2', minMeanExpr = 0,
#'                       design = design, reduced = reduced, fitType = 'local',
#'                       networkConstruction = 'COEN',
#'                       enrichTest = 'FET',
#'                       reg = paste0('gene', seq(30)))
#'
# \donttest{
#' ## RegEnrich analysis
#' object = regenrich_diffExpr(object)
#' 
#' # Set a random softPower, otherwise it is difficult to achive a
#' # scale-free network because of a randomly generated count data.
#' object = regenrich_network(object, softPower = 3)
#' object = regenrich_enrich(object)
#' object = regenrich_rankScore(object)
#'
#' ## plot expression of a regulator and its targets.
#' plotRegTarExpr(object, reg = 'gene5')
#' plotRegTarExpr(object, reg = 'gene27')
# }
plotRegTarExpr = function(object, reg, n = 1000, scale = TRUE,
    tarCol = "black", tarColAlpha = 0.1, regCol = "#ffaa00",
    xlab = "Samples", ylab = "Z-scores", ...) {
    if (length(reg) != 1) {
        stop("The length of 'reg' must be one.")
    }

    if (ncol(object) < 2 || nrow(object) <
        2) {
        stop("Too little expression data.")
    } else {
        expr = assay(object)
    }

    if (isEmptyPFC(mcols(object))) {
        stop("Differential expression ",
            "analysis needs to be performed.")
    } else {
        pFC = as.data.frame(mcols(object))
    }

    if (is.null(object@topNetwork)) {
        stop("object@topNetwork is empty, network inference ", 
             "needs to be performed.")
    } else {
        topNet = object@topNetwork
    }

    .plotRegTarExpr(reg = reg, expr = expr, pFC = pFC, topNet = topNet,
        n = n, scale = scale, tarCol = alpha(tarCol, tarColAlpha),
        regCol = regCol, xlab = xlab, ylab = ylab, ...)
}

.plotRegTarExpr = function(reg, expr, pFC, topNet, n = 1000,
    scale = TRUE, tarCol = alpha("black", 0.1), regCol = "#ffaa00",
    xlab = "Samples", ylab = "Z-scores", ...) {
    if (length(reg) != 1) {
        stop("The length of reg must be one")
    }
    # stopifnot(class(topNet) %in% 'TopNetwork')
    stopifnot(is(topNet, "TopNetwork"))
    net = .net(topNet)
    stopifnot(reg %in% names(net))
    stopifnot("p" %in% colnames(pFC))
    stopifnot(n > 1)
    pReg = subset(pFC[net[[reg]], ], p < 0.05)
    
    gNames = rownames(sortDataframe(pReg, "p"))[seq(min(n, nrow(pReg)))]

    regG = c(reg, gNames)
    regG_id = regG %in% rownames(expr)
    if (!all(regG_id)) {
        regG = regG[regG_id]
    }
    expr_regTar = expr[regG, ]

    if (scale) {
        expr_regTar = t(scale(t(expr_regTar)))
    }

    # reshape::melt
    expr_regTar2 = melt(data.frame(expr_regTar,
        gene = factor(rownames(expr_regTar),
        rownames(expr_regTar))), id.vars = "gene")
    expr_regTar2 = sortDataframe(expr_regTar2, c("gene", "variable"))
    p = ggplot(data = expr_regTar2, aes_string(x = "variable",
        y = "value", group = "gene"), ...) + geom_line(aes(color = I(tarCol)),
        show.legend = FALSE) + geom_line(aes(color = I(regCol)),
        data = subset(expr_regTar2, expr_regTar2$gene == reg)) +
        xlab(xlab) + ylab(ylab) + ggtitle(reg) + theme_light() +
        theme(plot.title = element_text(hjust = 0.5),
            axis.text.x = element_text(angle = 90,
            hjust = 1, vjust = 0.5))
    return(p)
}



#' Plot soft power for WGCNA analysis
#' @description  Plot soft power and corresponding scale free topology
#' fitting index to find a proper soft power for WGCNA analysis.
#' @param expr Gene expression data, either a matrix or a data frame.
#' By default, each row represents a gene, each column represents a sample.
#' @param rowSample logic. If \code{TRUE}, each row represents a sample.
#' The default is \code{FALSE}.
#' @param weights, optional observation weights for \code{expr} to be used 
#' in correlation calculation. 
#' @param powerVector a vector of soft thresholding powers for which the scale
#' free topology fit indices are to be calculated.
#' @param RsquaredCut desired minimum scale free topology fitting index R^2.
#' The default is 0.85.
#' @param networkType character, network type. Allowed values are
#' (unique abbreviations of) "unsigned" (default), "signed", "signed hybrid".
#' See \code{\link{adjacency}}.
#' @param removeFirst, should the first bin be removed from the connectivity 
#' histogram? The default is FALSE.
#' @param nBreaks, number of bins in connectivity histograms. The default is
#' 10.
#' @param corFnc, correlation function to be used in adjacency calculation.
#' The default is the \code{cor} function in WGCNA.
#' @param corOptions, a named list of options to the correlation function 
#' specified in corFnc. The default is list(use = "p").
#' @param BPPARAM, a BiocParallelParam instance to determine the parameters of 
#' parallel computing, see \code{\link{bplapply}} for more details.
#' @return a list of three elements: \code{powerEstimate}, \code{fitIndices},
#' and \code{plot}.
#' \code{powerEstimate} is an estimate of an appropriate soft-thresholding 
#' power. \code{fitIndices} is a data frame containing the fit indices for 
#' scale free topology. The \code{plot} is a ggplot object.
#' @import WGCNA
#' @import BiocParallel
#' @examples
#' data(Lyme_GSE63085)
#' log2FPKM = log2(Lyme_GSE63085$FPKM + 1)
#' log2FPKMhi = log2FPKM[rowMeans(log2FPKM) >= 10^-3, , drop = FALSE]
#' log2FPKMhi = head(log2FPKMhi, 3000) # First 3000 genes for example
#' softP = plotSoftPower(log2FPKMhi, RsquaredCut = 0.85)
#' print(softP)
#' @export
#' 
plotSoftPower = function(expr, rowSample = FALSE,
                         weights = NULL, 
                         powerVector = c(seq(10), seq(12, 20, by=2)),
                         RsquaredCut = 0.85, networkType = "unsigned",
                         removeFirst = FALSE, nBreaks = 10, 
                         corFnc = WGCNA::cor, 
                         corOptions = list(use = "p"), 
                         BPPARAM = bpparam()){
    stopifnot(RsquaredCut < 1 || RsquaredCut > 0)
    
    if(!rowSample) {
        expr = t(expr)
        rowSample = !rowSample
    }
    # Call the network topology analysis function
    sft <- pickSoftThreshold2(
        expr, weights = weights, 
        RsquaredCut = RsquaredCut,
        powerVector = powerVector,
        networkType = networkType,
        removeFirst = removeFirst, 
        nBreaks = nBreaks, 
        corFnc = corFnc, 
        corOptions = corOptions, 
        BPPARAM = BPPARAM)
    # Scale-free topology fit index as a function
    # of the soft-thresholding power
    ylab = expression(paste("Scale Free Topology Model Fit (signed, ",
        R^2, ")"))

    dat = sft$fitIndices
    dat$textY = dat$mean.k./max(dat$mean.k.)
    sft$plot = ggplot(dat, aes_string(x = "Power")) +
        geom_text(aes_string(y = "SFT.R.sq", label = "Power"), colour = "red")+
        geom_text(aes_string(y = "textY", label = "Power"), colour = "blue")+
        geom_hline(yintercept = RsquaredCut, color = "red")+
        scale_y_continuous(breaks = seq(0, 1, length = 6),
            sec.axis = sec_axis(~.*max(dat$mean.k.),
            name = "Mean connectivity"))+
        labs(x = "Soft Threshold (power)", y = ylab)+
        coord_cartesian(ylim = c(0, 1))+
        theme_bw()+
        theme(axis.title.y.left = element_text(colour = "red"),
            axis.title.y.right = element_text(colour = "blue"))
    return(sft)
}

# A light version of pickSoftThreshold function from WGCNA
#' @importFrom stats median
#' @importFrom utils getFromNamespace
pickSoftThreshold2 = function (data, 
                               weights = NULL, RsquaredCut = 0.85, 
                               powerVector = c(seq(1, 10, by = 1), 
                                               seq(12, 20, by = 2)), 
                               removeFirst = FALSE, nBreaks = 10, 
                               corFnc = WGCNA::cor, 
                               corOptions = list(use = "p"), 
                               networkType = "unsigned", 
                               BPPARAM = bpparam()) {
    powerVector = sort(powerVector)
    .networkTypes = c("unsigned", "signed", "signed hybrid")
    intType = charmatch(networkType, .networkTypes)
    if (is.na(intType)) 
        stop(paste0("networkType must be one of: ", 
                    paste(.networkTypes, collapse = ", "), "."))
    nGenes = ncol(data)
    if (nGenes < 3) {
        stop("The input data data contain fewer than 3 rows (nodes).", 
             "\nThis would result in a trivial correlation network.")
    }
    colname1 = c("Power", "SFT.R.sq", "slope", "truncated R.sq", 
                 "mean(k)", "median(k)", "max(k)")
    corFnc = match.fun(corFnc)
    corOptions$x = data
    if (!is.null(weights)) {
        if (!isTRUE(all.equal(dim(data), dim(weights)))) 
            stop("When 'weights' are given, dimensions of 'data' and 'weights'",
                 " must be the same.")
        corOptions$weights.x = weights
    }
    corOptions$y = data
    if (!is.null(weights))
        corOptions$weights.y = weights
    corx = do.call(corFnc, corOptions)
    if (intType == 1) {
        corx = abs(corx)
    } else if (intType == 2) {
        corx = (1 + corx)/2
    } else if (intType == 3) {
        corx[corx < 0] = 0
    }
    if (any(is.na(corx))) warning("Some correlations are NA.")
    
    scaleFreeFitIndex = getFromNamespace("scaleFreeFitIndex", "WGCNA")
    corxPrev = matrix(1, nrow = nrow(corx), ncol = ncol(corx))
    powerSteps <- powerVector - c(0, head(powerVector, -1))
    uniquePowerSteps <- unique(powerSteps)
    corxPowers <- lapply(uniquePowerSteps, function(p) corx^p)
    names(corxPowers) <- uniquePowerSteps
    
    dat = list()
    
    for(i in seq_along(powerVector)){
        corxCur = corxPrev * corxPowers[[as.character(powerSteps[i])]]
        corxPrev = corxCur
        dat[[i]] = colSums(corxCur, na.rm = TRUE) - 1
    }
    
    datout = BiocParallel::bplapply(seq_along(powerVector), function(i){
        # corxCur <- corxPrev * corxPowers[[as.character(powerSteps[i])]]
        # assign("corxPrev", corxCur, envir = env)
        # khelp = colSums(corxCur, na.rm = TRUE) - 1
        khelp = dat[[i]]
        SFT1 = scaleFreeFitIndex(k = khelp, nBreaks = nBreaks, 
                                 removeFirst = removeFirst)
        c(powerVector[i], unlist(SFT1), 
          unlist(lapply(c(mean, median, max), function(f) {
              f(khelp, na.rm = TRUE)})))
    }, BPPARAM = BPPARAM)
    
    datout = do.call(rbind, datout)
    colnames(datout) = colname1
    
    ind1 = datout[, 2] > RsquaredCut
    indcut = ifelse(sum(ind1) > 0, min(c(1:length(ind1))[ind1]), NA)
    powerEstimate = powerVector[indcut][[1]]
    gc()
    list(powerEstimate = powerEstimate, fitIndices = data.frame(datout))
}


# function inside plot_Enrich
.plot_Enrich = function(object, reg = NULL, showCategory = 20,
                        regDescription = NULL, font.size = 12){
    if(is.null(regDescription)){
        regAll = results_enrich(object)@allResult[,1]
        DescriptionAll = results_enrich(object)@allResult[,1]
        regDescription = data.frame(reg = regAll, Description = DescriptionAll)
    }
    stopifnot(is.data.frame(regDescription) && ncol(regDescription) == 2)
    stopifnot(is.null(reg) || length(reg) == 1)
    object_enrich = results_enrich(object)
    if(is.null(object_enrich)){
        stop("FET/GSEA enrichment analysis has to be performed.")
    }
    if (object_enrich@type == "FET"){
        if (!is.null(reg)){
            message("'reg' argument is not used when ",
                    "the enrichment method is 'FET'.")
        }
        topResult = object_enrich@topResult
        stopifnot(nrow(topResult) > 0)
        df = topResult[seq(min(nrow(topResult), showCategory)), ]
        ratio = vapply(df$GeneRatio,
                       function(x) eval(parse(text = x)), FUN.VALUE = 1)
        idx = order(ratio, decreasing = FALSE)
        
        Description = regDescription[match(df$Description,
                                           regDescription[,1]), 2]
        df$Description <- factor(Description, levels = Description[idx])
        
        # replace Inf with the maximum finite value
        cB = -log10(df[["p.adjust"]])
        cB[is.infinite(cB)] = max(cB[is.finite(cB)])
        
        # color
        colorStr = paste0("-log10(p.adjust)")
        df$color = cB
        df$GeneRatio = ratio
        
        # plot
        pt = ggplot(df, aes_string(x = "ratio", y = "Description",
                                   size = "Count", color = "color")) +
            geom_point() +
            scale_color_gradient(name = colorStr, low = "blue", high = "red") +
            xlab("GeneRatio") + ylab("Regulator") +
            theme_dose(font.size)
    } else if (object_enrich@type ==  "GSEA"){
        namedScores = object_enrich@namedScores
        network = .net(results_topNet(object))
        
        if(is.null(reg)){
            message("reg is not provided, the most ",
                    "significantly enriched regulator is plotted.")
            reg = object_enrich@allResult$regulator[1]
        }
        
        Description = regDescription[match(reg, regDescription[,1]), 2][1]
        
        pt = fgsea::plotEnrichment(pathway = network[[reg]],
                                   stats = namedScores) +
            ggtitle(Description) +
            theme_dose(font.size)+
            theme(plot.title = element_text(hjust = 0.5))
    } else {
        stop("Unknown enrichment result.")
    }
    return(pt)
}


#' @rdname plot_Enrich
#' @export
setGeneric("plot_Enrich", function(object, ...) standardGeneric("plot_Enrich"))


#' Plot results of FET/GSEA enrichment analysis
#'
#' Plot FET/GSEA enrichment results. If the FET method is applied,
#' the top `showCategory` regulator will be plotted.
#' If the GSEA method is applied, the GSEA graph of regulator `reg`
#' will be plotted.
#'
#' @param object a \code{RegenrichSet} object.
#' @param reg The regulator to plot. This only works when the
#' GSEA enrichment method has used.
#' @param showCategory the number of regulator to plot.
#' @param regDescription NULL or a two-column data frame, in which
#' first column is the regulator IDs (for
#' example ENSEMBL IDs), and the second column is the description
#' of regulators (for example gene
#' name). Default is NULL, meaning both columns are the same regulator 
#' names/IDs in the network.
#' @param font.size font size of axis labels and axis tick mark
#' labels, default is 12.
#' @param ... other parameters.
#' @return a ggplot object of plotting FET or GSEA enrichment result.
#' @import fgsea
#' @import DOSE
#' @include regenrichClasses.R
#' @rdname plot_Enrich
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
#' object = RegenrichSet(expr = data1,
#'                       colData = colData,
#'                       method = "limma", minMeanExpr = 0,
#'                       design = design,
#'                       contrast = c(rep(0, ncol(design) - 1), 1),
#'                       networkConstruction = "COEN",
#'                       enrichTest = "FET")
#' # Differential expression analysis
#' object = regenrich_diffExpr(object)
#' # Network inference using "COEN" method
#' object = regenrich_network(object)
#' # Enrichment analysis by Fisher's exact test (FET)
#' object = regenrich_enrich(object)
#' # plot
#' plot_Enrich(object)
#'
#' # Enrichment analysis by Fisher's exact test (FET)
#' object = regenrich_enrich(object, enrichTest = "GSEA")
#' # plot
#' plot_Enrich(object)
setMethod("plot_Enrich", signature = "RegenrichSet", .plot_Enrich)





