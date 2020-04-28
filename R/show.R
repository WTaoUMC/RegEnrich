#' methods of generic function "show"
#' @param object one object of either \code{DeaSet}, \code{TopNetwork}, 
#' \code{Enrich}, \code{Score}, or \code{RegenrichSet} class.
#' @rdname methodsOfShow
#' @importMethodsFrom methods show
#' @return show returns an invisible original \code{object}.
#' @export
#' @examples
#' x = newScore(letters[1:5], 1:5, 1:5, -2:2, seq(2, 1, len = 5))
#' show(x)
#' 
setMethod(f = "show", signature = "DeaSet", function(object) {
  pFC = mcols(object)[, c("gene", "p", "logFC")]
  dim_rawData = dim(pFC)
  if (all(mcols(object)$p < 1e-13)) {
    message("Differential expression analysis needs to be performed.")
  }
  n_sig = sum(pFC$p < 0.05)
  
  cat(class(object), " object (total ", dim_rawData[1], " rows, of which ",
      n_sig, " with p < 0.05).\n", sep = "")
  print(pFC)
  return(invisible(object))
})


#' @rdname methodsOfShow
#' @export
# method in show generic for 'TopNetwork' object
setMethod(f = "show", signature = "TopNetwork", function(object) {
  dim_rawData = dim(object@elementset)
  
  cat(class(object), " object (", dim_rawData[1],
      " edges, networkConstruction: '",
      object@networkConstruction, "', percentage: ", object@percent,
      "%)\n", sep = "")
  show(object@elementset)
  return(invisible(object))
})


#' @rdname methodsOfShow
#' @export
#  method in show generic for 'Enrich' object
setMethod(f = "show", signature = "Enrich", definition = function(object) {
  enrichType = object@type
  
  if (enrichType == "FET") {
    topRes = object@topResult
    cat("Enrich object (FET method,", nrow(object@allResult),
        "regulators are used for enrichment, \n", nrow(topRes),
        "regulators pass the threshold)\n")
    if (nrow(topRes) > 0) {
      targ = paste0(apply(limma::strsplit2(object@topResult$geneID,
                                           split = "/")[, seq(2)], 1, paste0, collapse = "/"),
                    "...")
      topRes$geneID = targ
      
      print(topRes)
    }
  } else if (enrichType == "GSEA") {
    topRes = object@topResult
    cat("Enrich object (GSEA method,", nrow(object@allResult),
        "regulators are used for enrichment, \n", nrow(topRes),
        "regulators pass the threshold.\n")
    if (nrow(topRes) > 0) {
      targ = paste0(unlist(lapply(topRes$leadingEdge, function(x) {
        paste0(x[seq(2)], collapse = ",")
      })), ",...")
      topRes$leadingEdge = targ
      print(topRes)
    }
  }
  invisible(object)
})



#' @rdname methodsOfShow
#' @export
# Show Score object
setMethod("show", signature = "Score", function(object){
  print.tbl = utils::getFromNamespace("print.tbl", "tibble")
  print.tbl(S3Part(object))
  invisible(object)
})


# show RegenrichSet object
.showRegenrichSet = function(object) {
  dim_rawData = dim(object@assayRaw)
  dim_assayData = dim(assay(object))
  cat(class(object), " object \n assayData: ", dim_assayData[1],
      " rows, ", dim_assayData[2], " columns (filtered ", dim_rawData[1] -
        dim_assayData[1], " rows", sep = "")
  if (is.null(object@paramsIn$minMeanExpr)) {
    cat(")\n")
  } else {
    cat(" by average expression <= ", object@paramsIn$minMeanExpr,
        ")\n", sep = "")
  }
  
  cat("\n")
  pFC = mcols(object)
  dim_resDEA = if (isEmptyPFC(pFC)) {
    c(0, 0)
  } else dim(pFC)
  
  if (dim_resDEA[1] > 0) {
    p = pFC[, "p"]
    cat(" (1) ", sum(p < 0.05), " rows with differential ",
        "p-value < 0.05\n", sep = "")
  } else {
    message(" Differential expression analysis needs ",
            "to be performed.\n")
    return(invisible(object))
  }
  
  dim_network = if (isEmpty(object@network)) {
    c(0, 0)
  } else dim(object@network)
  if (dim_network[1] == 0) {
    message("\n Network inference needs to be performed, or ",
            "a 'TopNetwork' object needs to be provided.\n",
            sep = "")
    return(invisible(object))
  }
  
  # show top p network
  cat("\n")
  dim_topNetP = if (is.null(object@topNetwork)) {
    c(0, 0)
  } else dim(object@topNetwork)
  if (dim_topNetP[1] > 0) {
    cat(" (2) Top p% network info:\n")
    show(object@topNetwork)
  } else {
    message(" Top p% network needs to be retained.\n")
    return(invisible(object))
  }
  
  # show enrichment
  cat("\n")
  if (dim(object@resEnrich@allResult)[1] > 0) {
    cat(" (3) Enrichment info:\n")
    show(object@resEnrich)
  } else {
    message(" FET/GSEA enrichment needs to be performed.\n")
    return(invisible(object))
  }
  
  # show score
  cat("\n")
  dim_resScore = dim(object@resScore)
  if (dim_resScore[1] > 0) {
    cat(" (4) RegEnrich score:\n")
    show(object@resScore)
  } else {
    message(" 'regenrich_rankScore' needs to be performed.\n")
    return(invisible(object))
  }
}

#' @rdname methodsOfShow
#' @export
# method in show generic for 'RegenrichSet' object
setMethod(f = "show", signature = "RegenrichSet", .showRegenrichSet)



#' Print Score object
#' @param x a Score object.
#' @param ... optional arguments to print.
#' @export
#' @return print.Score returns the a \code{Score} object
#' @examples
#' x = newScore(letters[1:5], 1:5, 1:5, -2:2, seq(2, 1, len = 5))
#' print(x)
print.Score = function(x, ...){
  print.tbl = utils::getFromNamespace("print.tbl", "tibble")
  print.tbl(S3Part(x), ...)
  invisible(x)
}
