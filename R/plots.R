
#' Compare the orders of two vectors
#' @description compare the orders of two vectors
#' @param name1 a vector with first order.
#' @param name2 a vector with anothoer second order.
#' @rawNamespace import(ggplot2, except = margin)
#'
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
#' \if{html}{\figure{plotRegTarExpr1.png}{options: width='70\%'
#' alt='Figure: plotRegTarExpr1.png'}}
#' \if{latex}{\figure{plotRegTarExpr1.pdf}{options: width=6cm}}
#' \if{html}{\figure{plotRegTarExpr2.png}{options: width='70\%'
#' alt='Figure: plotRegTarExpr2.png'}}
#' \if{latex}{\figure{plotRegTarExpr2.pdf}{options: width=6cm}}
#' @importFrom reshape2 melt
# @rawNamespace import(ggplot2, except = margin)
#' @export
#' @examples
#' # constructing a RegenrichSet object
#' pdata = data.frame(patientID = paste0('Sample_', seq_len(50)),
#'                    week = rep(c('0', '1'), each = 25),
#'                    row.names = paste0('Sample_', seq_len(50)))
#' design = ~week
#' reduced = ~1
#' set.seed(123)
#' cnts = matrix(rnbinom(n=1000*50, mu=100, size=1/0.1), ncol=50,
#'                dimnames = list(paste0('gene', seq_len(1000)),
#'                                rownames(pdata)))
#' cnts[5,26:50] = cnts[5,26:50] + 50 # add reads to gene5 in some samples.
#' id = sample(31:1000, 20) # randomly select 20 rows, and assign reads.
#' cnts[id,] = vapply(cnts[5,], function(x){
#'   rnbinom(n = 20, size = 1/0.02, mu = x)},
#'   FUN.VALUE = rep(1, 20))
#'
#' object = RegenrichSet(expr = cnts,
#'                       pData = pdata,
#'                       method = 'LRT_DESeq2', minMeanExpr = 0,
#'                       design = design, reduced = reduced, fitType = 'local',
#'                       networkConstruction = 'COEN',
#'                       enrichTest = 'FET',
#'                       reg = paste0('gene', seq_len(30)))
#'
#' ## RegEnrich analysis
#' object = regenrich_diffExpr(object)
#' # Set a random softPower, otherwise it is difficult to achive a
#' # scale-free network because of a randomly generated count data.
#' object = regenrich_network(object, softPower = 3)
#' object = regenrich_enrich(object)
#' object = regenrich_rankScore(object)
#' head(slot(object, 'resScore'))
#' tail(slot(object, 'resScore'))
#'
#' ## plot expression of a regulator and its targets.
#' plotRegTarExpr(object, reg = 'gene5')
#' plotRegTarExpr(object, reg = 'gene27')
plotRegTarExpr = function(object, reg, n = 1000, scale = TRUE,
    tarCol = "black", tarColAlpha = 0.1, regCol = "#ffaa00",
    xlab = "Samples", ylab = "Z-scores", ...) {
    if (length(reg) != 1) {
        stop("The length of 'reg' must be one.")
    }

    if (ncol(object@assayData) < 2 || nrow(object@assayData) <
        2) {
        stop("Too little expression data.")
    } else {
        expr = object@assayData
    }

    if (is.null(object@resDEA@pFC)) {
        stop("object@resDEA@pFC is NULL, differential expression ",
            "analysis needs to be performed.")
    } else {
        pFC = object@resDEA@pFC
    }

    if (is.null(object@topNetP)) {
        stop("object@topNetP is NULL, network inference needs to be performed.")
    } else {
        topNet = object@topNetP
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
    stopifnot(reg %in% names(topNet@net))
    stopifnot("p" %in% colnames(pFC))
    stopifnot(n > 1)
    pReg = subset(pFC[topNet@net[[reg]], ], p < 0.05)
    gNames = rownames(sortDataframe(pReg, "p"))[seq_len(min(n,
        nrow(pReg)))]

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
#' @param powerVector a vector of soft thresholding powers for which the scale
#' free topology fit indices are to be calculated.
#' @param RsquaredCut desired minimum scale free topology fitting index R^2.
#' The default is 0.85.
#' @param networkType character, network type. Allowed values are
#' (unique abbreviations of) "unsigned" (default), "signed", "signed hybrid".
#' See \code{\link{adjacency}}.
#' @param verbose integer level of verbosity. 0 (default) means silent, higher
#' values make the output progressively more and more verbose.
#' @param ... Parameters in \code{\link{pickSoftThreshold}} function.
#' @return a list of three elements: \code{powerEstimate}, \code{fitIndices},
#' and \code{plot}.
#' The details about \code{powerEstimate} and \code{fitIndices} are shown
#' in \code{\link{pickSoftThreshold}}. The \code{plot} is a ggplot object.
#' @seealso \code{\link{pickSoftThreshold}}
#' @import WGCNA
#' @examples
#' data(Lyme_GSE63085)
#' log2FPKM = log2(Lyme_GSE63085$FPKM + 1)
#' log2FPKMhi = log2FPKM[rowMeans(log2FPKM) >= 10^-3, , drop = FALSE]
#' softP = plotSoftPower(log2FPKMhi, RsquaredCut = 0.85)
#' @export
plotSoftPower = function(expr, rowSample = FALSE,
    powerVector = c(seq_len(10), seq(12, 20, by=2)),
    RsquaredCut = 0.85, networkType = "unsigned",
    verbose = 0, ...){
    stopifnot(RsquaredCut < 1 || RsquaredCut > 0)

    if(!rowSample) {
        expr = t(expr)
        rowSample = !rowSample
     }
    # Call the network topology analysis function
    # enableWGCNAThreads()
    # on.exit(disableWGCNAThreads())
    sft = pickSoftThreshold(expr, RsquaredCut = RsquaredCut,
         powerVector = powerVector,
         networkType = networkType,
         verbose = verbose, ...)
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


# function inside plot_Enrich
.plot_Enrich = function(object, showCategory = 20, reg = NULL,
     regDescription = data.frame(reg = results_enrich(object)@allResult[,1],
         Description = results_enrich(object)@allResult[,1]), font.size = 12){
    stopifnot(is.data.frame(regDescription) && ncol(regDescription) == 2)
    stopifnot(is.null(reg) || length(reg) != 1)
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
        df = topResult[seq_len(min(nrow(topResult), showCategory)), ]
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
        network = results_topNet(object)@net

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
#' @param showCategory the number of regulator to plot.
#' @param reg The regulator to plot. This only works when the
#' GSEA enrichment method has used.
#' @param regDescription a two-column data frame, in which
#' first column is the regulator ID (for
#' example ENSEMBL ID), and the second column is the description
#' of regulators (for example gene
#' name). Default is a data frame, in which the regulator in the
#' network were repeated in the
#' two columns.
#' @param font.size font size of axis labels and axis tick mark
#' labels, default is 12.
#' @param ... additional arguments.
#' @return a ggplot object of plotting FET or GSEA enrichment result.
#' @import fgsea
#' @import DOSE
#' @include regenrichClasses.R
#' @rdname plot_Enrich
#' @export
#' @examples
#' # library(RegEnrich)
#' # Initializing a "RegenrichSet" object
#' data = log2(Lyme_GSE63085$FPKM + 1)
#' x = apply(data, 1, sd)
#' pData = Lyme_GSE63085$sampleInfo
#' data1 = data[seq_len(2000), ]
#'
#' pData$week = as.factor(pData$week)
#' pData$patientID = as.factor(sub("(\\d+)-(\\d+)", "\\1_\\2",
#'                             pData$patientID))
#'
#' design = model.matrix(~0 + patientID + week, data = pData)
#' object = RegenrichSet(expr = data1,
#'                       pData = pData,
#'                       method = "limma", minMeanExpr = 0,
#'                       design = design,
#'                       contrast = c(rep(0, ncol(design) - 1), 1),
#'                       networkConstruction = "COEN",
#'                       enrichTest = "FET")
#'
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
#'
setMethod("plot_Enrich", signature = "RegenrichSet", .plot_Enrich)





