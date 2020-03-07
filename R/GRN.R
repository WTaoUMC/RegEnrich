#' Constructing Gene Regulatory Network (GRN)
#' @description \code{GRN} is a function to construct the GRN network
#' using random forest algorithm. This method was initially introduced in
#' GENIE3, but here, GRN support parallel computing. It can control the
#' model accuracy, and define the regulators in the network.
#' @param expr Gene expression data, either a matrix or a data frame.
#' By default (\code{rowSample = FALSE}), each row represents a gene,
#' each column represents a sample.
#' @param reg vector of charactors, representing gene regulators.
#' By default, these are transcription factors and co-factors,
#' defined by three literatures/databases, namely RegNet, TRRUST, and
#' Marbach2016.
#' @param rowSample logical If \code{TRUE}, each row represents a sample.
#' The default is \code{FALSE}.
#' @param K integer or character. The number of features in each tree,
#' can be either a integer number, \code{sqrt}, or \code{all}.
#' \code{sqrt} denotes sqrt(the number of \code{reg}), \code{all}
#' means the number of \code{reg}. The default is \code{sqrt}.
#' @param nbTrees integer. The number of trees. The default is 1000.
#' @param importanceMeasure character. The importance type in
#' \code{importance}. importanceMeasure can be \code{\%IncMSE}  or
#' \code{IncNodePurity}, corresponding to type = 1 and 2 in
#' \code{importance}. The default is \code{IncNodePurity}
#' (decrease in node impurity), which is faster than \code{\%IncMSE}
#' (decrease in accuracy).
#' @param trace logical. To show the progress or not (default).
#' @param BPPARAM parameters for parallel computing (default is
#' \code{bpparam()}).
#' @param maxMSE numeric. The maximum out-of-bag MSE is to control model
#' accuracy. The default is NULL, which means no filtering by this parameter.
#' @param minR numeric. The minimum correlation coefficient of prediction is to
#' control model accuracy. The default is 0.3.
#' @param ... the rest parameters in \code{\link{randomForest}} function.
#'
#' @return A list of 'weightHi' and 'performanceHi'. 'weightHi' is the edge
#' weights by high accurate model, and 'performanceHi' is the performances
#' of high accurate models.
#' @import randomForest
#' @import BiocParallel
#' @rawNamespace import(data.table, except = c(melt, dcast))
#' @include globals.R
# @export
GRN = function(expr, reg = TFs$TF_name, rowSample = FALSE, 
    K = "sqrt", nbTrees = 1000, importanceMeasure = "IncNodePurity", 
    trace = FALSE, BPPARAM = bpparam(), maxMSE = NULL, 
    minR = 0.3, ...) {
    stopifnot(is.null(maxMSE) || maxMSE > 0)
    stopifnot(is.null(minR) || (minR > 0 & minR < 1))
    
    inout = inOutput(expr = expr, reg = reg, rowSample = rowSample)
    inputMatrix = inout$inputMatrix
    outputMatrix = inout$outputMatrix
    
    net = grNet(inputMatrix, outputMatrix, K = K, nbTrees = nbTrees, 
        importanceMeasure = importanceMeasure, trace = trace, 
        BPPARAM = bpparam(), ...)
    weight = net$weight
    performance = net$performance
    
    # Only to use the good models
    performanceHi = performance
    if (!is.null(maxMSE)) {
        performanceHi = performanceHi[mse <= maxMSE]
    }
    if (!is.null(minR)) {
        performanceHi = performanceHi[r >= minR]
    }
    
    modelHi = performanceHi$gene
    if (length(modelHi) == 0) {
        stop("All models are not good enough, please check if 'maxMSE' ", 
            "and 'minR' were properly set, or try other ", 
            "network inference methods.")
    }
    weightHi = weight[to.gene %in% modelHi]
    return(list(weightHi = weightHi, performanceHi = performanceHi))
}

.local_grNet = function(targetGeneIdx, inputGeneNames, 
    outputGeneNames, inputMatrix, outputMatrix, K, 
    nbTrees, importanceMeasure, trace, nodesizeInArgs, 
    ...) {
    num.samples = nrow(outputMatrix)
    num.genes = ncol(outputMatrix)
    if (trace) {
        cat(paste("Computing gene ", targetGeneIdx, 
            "/", num.genes, "\n", sep = ""))
        utils::flush.console()
    }
    targetGeneName = outputGeneNames[targetGeneIdx]
    theseInputGeneNames = setdiff(inputGeneNames, targetGeneName)
    x = inputMatrix[, theseInputGeneNames]
    numInputGenes = length(theseInputGeneNames)
    y = outputMatrix[, targetGeneName]
    if (is(K, "numeric")) {
        mtry = K
    } else if (K == "sqrt") {
        mtry = round(sqrt(numInputGenes))
    } else if (K == "all") {
        mtry = numInputGenes
    } else {
        stop("Parameter K must be \"sqrt\", or \"all\", or an integer")
    }
    if (trace) {
        cat(paste("K = ", mtry, ", ", nbTrees, " trees\n\n", 
            sep = ""))
        utils::flush.console()
    }
    if (importanceMeasure == "%IncMSE") {
        y = y
        importance0 = TRUE
        type0 = 1
    } else {
        y = y/sd(y)
        importance0 = FALSE
        type0 = 2
    }
    if (nodesizeInArgs) {
        rf = randomForest(x, y, mtry = mtry, ntree = nbTrees, 
            importance = importance0, ...)
    } else {
        rf = randomForest(x, y, mtry = mtry, ntree = nbTrees, 
            importance = importance0, nodesize = 1, 
            ...)
    }
    im = importance(rf, type = type0)
    weight = data.table(from.gene = rownames(im), to.gene = rep(targetGeneName, 
        nrow(im)), imp = im/num.samples)
    y_hat = stats::predict(rf)
    return(list(weight = weight, mse = sum((y_hat - 
        y)^2)/length(y), r = as.numeric(cor(y_hat, 
        y)), p = stats::cor.test(y_hat, y, alternative = "greater")$p.value, 
        pVarExplaned = rf$rsq[length(rf$rsq)] * 100))
}


grNet = function(inputMatrix, outputMatrix, K = "sqrt", 
    nbTrees = 1000, importanceMeasure = "IncNodePurity", 
    trace = FALSE, BPPARAM = bpparam(), ...) {
    ## All of the matrixs are sample * gene inputMatrix
    ## is Gene Expression matrix of TFs outputMatrix is
    ## Gene Expression matrix of All genes (or TFs)
    
    # Report when parameter importanceMeasure is not
    # correctly spelled
    if (importanceMeasure != "IncNodePurity" && importanceMeasure != 
        "%IncMSE") {
        stop("Parameter importanceMeasure must be ", 
            "\"IncNodePurity\" or \"%IncMSE\"")
    }
    # Check if nodesize parameter is in the input
    # arguments
    args = list(...)
    nInArgs = "nodesize" %in% names(args)
    outputGeneNames = colnames(outputMatrix)
    inputGeneNames = colnames(inputMatrix)
    
    idx = stats::setNames(seq_len(ncol(outputMatrix)), 
        nm = outputGeneNames)
    
    tic = system.time({
        res = bplapply(idx, function(targetGeneIdx0) {
            .local_grNet(targetGeneIdx0, inputGeneNames, 
                outputGeneNames, inputMatrix, outputMatrix, 
                K, nbTrees, importanceMeasure, trace, 
                nodesizeInArgs = nInArgs, ...)
        }, BPPARAM = BPPARAM)
    })
    
    if (trace) {
        cat("GRN network inference costs", tic[3], 
            "seconds.\n")
    }
    weight = rbindlist(lapply(res, "[[", "weight"))
    colnames(weight) = c("from.gene", "to.gene", "weight")
    performance = data.table(gene = outputGeneNames, 
        mse = do.call("c", lapply(res, "[[", "mse")), 
        r = do.call("c", lapply(res, "[[", "r")), p = do.call("c", 
            lapply(res, "[[", "p")), pVarExplaned = do.call("c", 
            lapply(res, "[[", "pVarExplaned")))
    return(list(weight = weight, performance = performance))
}
