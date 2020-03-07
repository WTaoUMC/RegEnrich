# Obtain paramsIn slot from RegenrichSet object.
getParamsIn = function(object, arg = NULL) {
    stopifnot(is(object, "RegenrichSet"))
    if (!is.null(arg) && length(arg) != 1 && !is.character(arg)) {
        stop("arg can only be either NULL or character.")
    }
    if (is.null(arg)) {
        return(object@paramsIn)
    } else {
        return(object@paramsIn[[arg]])
    }
}

# Check if the names(argsInList) are all listed in paramsIn
# slot from RegenrichSet object.
checkParams = function(object, argsInList, mustInArgs = NULL) {
    argsName = names(argsInList)
    if (length(argsInList) > 0) {
        stopifnot(!is.null(argsName))
        if (!is.null(mustInArgs)) {
            indx = argsName %in% mustInArgs
            if (!all(argsName %in% names(object@paramsIn))) {
                stop("Unknown argument(s):\n", argsName[!argsName %in%
                    names(object@paramsIn)])
            }
            # arguments not in mustInArgs
            if (sum(!indx) > 0) {
                warning("Following argument(s) should not be respecified ",
                    "in the current function:\n", argsName[!indx])
            }
            # arguments in mustInArgs
            if (sum(indx) > 0) {
                argsInList = argsInList[indx]
            } else {
                argsInList = list()
            }
        }
    }
    if (length(argsInList) > 0) {
        object@paramsIn[names(argsInList)] = argsInList
    }
    return(object)
}

# sort data frame rows by its column data.
sortDataframe = function(x, by = x, decreasing = FALSE, returnID = FALSE) {
    stopifnot(is.data.frame(x))
    if (is.character(by)) {
        nm = by
    } else if (is.data.frame(by)) {
        nm = colnames(by)
    } else if (is.integer(by)) {
        nm = colnames(x)[by]
    } else {
        stop("Unknown class of 'by'")
    }
    stopifnot(all(nm %in% colnames(x)))
    cmd = paste0("with(x, order(", paste0(nm, collapse = ","),
        ", decreasing = decreasing))")
    id = eval(parse(text = cmd))
    y = x[id, ]
    if (returnID) {
        y = list(res = y, id = id)
    }
    return(y)
}


# Generate the input matrix and output matrix for network
# inference by random forest @description Standardize the
# inputMatrix and outputMatrix for \code{\link{grNet}}.
# @param expr Gene expression data, either a matrix or a data
# frame.  By default (\code{rowSample = FALSE}), each row
# represents a gene, each column represents a sample.  @param
# reg vector of charactors, representing gene regulators.  By
# default, these are transcription factors and co-factors,
# defined by three literatures/databases, namely RegNet,
# TRRUST, and Marbach2016.  @param rowSample logic.  If
# \code{TRUE}, each row represents a sample.  The default is
# \code{FALSE}.  @return A list of \code{inputMatrix}
# (expression of \code{reg}), \code{outputMatrix}
# (expression of all genes) and \code{validRegs} (the
# regulators exsist in \code{expr}).  @examples \dontrun{
# expr = matrix(rnorm(100*1000), nrow = 1000, ncol = 100,
# dimnames = list(paste0('G', seq_len(1000)), paste0('Samp',
# seq_len(100)))) set.seed(1234) TFs = paste0('G',
# sample(seq_len(1000),
# size = 50, replace = FALSE)) # rowSample = FALSE
# inOutput(expr, reg = TFs, rowSample = FALSE) # rowSample =
# TRUE inOutput(t(expr), reg = TFs, rowSample = TRUE) }
# @export
#' @include globals.R
inOutput = function(expr, reg = TFs$TF_name, rowSample = FALSE,
    trace = FALSE) {
    if (!rowSample) {
        outputMatrix = t(expr)
    } else {
        outputMatrix = expr
    }
    exprGenes = colnames(outputMatrix)
    exprSamp = rownames(outputMatrix)

    # only to use the regulators existing in both expr and reg.
    validRegs = reg[reg %in% exprGenes]
    if (length(validRegs) == 0) {
        stop("No valide regulators can be found. Please ",
            "change 'reg' or check gene ID.")
    }
    if (trace) {
        cat(length(validRegs), " regulators will be used. \n")
    }
    inputMatrix = outputMatrix[, validRegs, drop = FALSE]

    # inputMatrix is the gene expression matrix of regulators
    # (only reg) outputMatrix is the gene expression matrix of
    # all genes (including reg)
    return(list(inputMatrix = inputMatrix, outputMatrix = outputMatrix,
        validRegs = validRegs))
}


# derived from DESeq2:::renameModelMatrixColumns function
renameModelMatrixColumns = function (data, design){
    data = as.data.frame(data)
    designVars = all.vars(design)
    designVarsClass = vapply(designVars,
        function(v) is.factor(data[[v]]), FUN.VALUE = TRUE)
    factorVars = designVars[designVarsClass]
    colNamesFrom = make.names(do.call(c, lapply(factorVars,
        function(v) paste0(v, levels(data[[v]])[-1]))))
    colNamesTo = make.names(do.call(c, lapply(factorVars,
        function(v) paste0(v, "_", levels(data[[v]])[-1],
            "_vs_", levels(data[[v]])[1]))))
    data.frame(from = colNamesFrom, to = colNamesTo,
        stringsAsFactors = FALSE)
}


# Adjacency matrix to a data.frame of edges.
# @param mat adjacency matrix.
# @param mode Character, to specify the class of graph and which part of
# the matrix will be used. Possible values are: "directed" (default),
# "undirected", "upper", "lower".
# @param diag logic, whether to include the diagonal of the matrix.
# @return a data.frame of edge information. The first column is from node,
# the second column is to node, and the third is weight.
# @examples {
# \dontrun{
# mat = matrix(rnorm(4*4), nrow = 4,
#              dimnames = list(letters[seq_len(4)], LETTERS[seq_len(4)]))
# mat2Edge(mat, mode = "undirected", diag = TRUE)
# mat2Edge(mat, mode = "undirected", diag = FALSE)
# mat2Edge(mat, mode = "directed", diag = TRUE)
# mat2Edge(mat, mode = "upper", diag = TRUE)
# mat2Edge(mat, mode = "upper", diag = FALSE)
# }
# }
mat2Edge = function(mat, mode = c("directed", "undirected", "upper", "lower"),
                    diag = FALSE, removeEdgesBelowThisWeight = NULL){
    mode = match.arg(mode)
    rowN = nrow(mat)
    colN = ncol(mat)
    nameRow = rownames(mat)
    if(is.null(nameRow)) nameRow = seq_len(rowN)
    nameCol = colnames(mat)
    if(is.null(nameCol)) nameCol = seq_len(colN)

    if (mode == "directed"){
        id = !diag(!diag, rowN, colN)
    } else if (mode %in% c("undirected", "upper")){
        id = upper.tri(mat, diag = diag)
    } else if (mode == "lower"){
        id = lower.tri(mat, diag = diag)
    }

  if (!is.null(removeEdgesBelowThisWeight) &&
      is.numeric(removeEdgesBelowThisWeight)){
      id = id & (mat >= removeEdgesBelowThisWeight)
  }

  id = which(id, arr.ind = TRUE, useNames = TRUE)
  return(data.frame(from = nameRow[id[,1]],
      to = nameCol[id[,2]],
      weight = mat[id],
      stringsAsFactors = FALSE))
}
