# method in show generic for 'DeaSet' object
setMethod(f = "show", signature = "DeaSet", function(object) {
    dim_rawData = dim(object@pFC)
    if (dim_rawData[1] < 1) {
        message("Differential expression analysis needs to be performed.")
    }
    n_sig = sum(object@pFC$p < 0.05)

    cat(class(object), " object (total ", dim_rawData[1], " rows, of which ",
        n_sig, " with p < 0.05).\n", sep = "")
    if (dim_rawData[1] > 6) {
        print(head(object@pFC, n = 3))
        cat("...\n")
        print.data.frame(tail(object@pFC, n = 3))
    } else {
        print(object@pFC)
    }
})

# method in show generic for 'TopNetwork' object
setMethod(f = "show", signature = "TopNetwork", function(object) {
    dim_rawData = dim(object@netTable)

    cat(class(object), " object (", dim_rawData[1],
        " edges, networkConstruction: '",
        object@networkConstruction, "', percentage: ", object@percent,
        "%)\n", sep = "")
    if (dim_rawData[1] > 6) {
        print(head(object@netTable, n = 3))
        cat("...\n")
        print.data.frame(tail(object@netTable, n = 3))
    } else {
        print(object@netTable)
    }
})


setMethod(f = "show", signature = "Enrich", definition = function(object) {
    enrichType = object@type

    if (enrichType == "FET") {
        topRes = object@topResult
        cat("Enrich object (FET method,", nrow(object@allResult),
            "regulators are used for enrichment, \n", nrow(topRes),
            "regulators pass the threshold)\n")
        if (nrow(topRes) > 0) {
            targ = paste0(apply(limma::strsplit2(object@topResult$geneID,
                split = "/")[, seq_len(2)], 1, paste0, collapse = "/"),
                "...")
            topRes$geneID = targ

            if (nrow(topRes) > 6) {
                print(head(topRes, 3))
                cat("...\n")
                print(tail(topRes, 3))
            } else {
                print(topRes)
            }
        }
    } else if (enrichType == "GSEA") {
        topRes = object@topResult
        cat("Enrich object (GSEA method,", nrow(object@allResult),
            "regulators are used for enrichment, \n", nrow(topRes),
            "regulators pass the threshold.\n")
        if (nrow(topRes) > 0) {
            targ = paste0(unlist(lapply(topRes$leadingEdge, function(x) {
                paste0(x[seq_len(2)], collapse = ",")
            })), ",...")
            topRes$leadingEdge = targ
            if (nrow(topRes) > 6) {
                print(head(topRes, 3))
                cat("...\n")
                print(tail(topRes, 3))
            } else {
                print(topRes)
            }
        }
    }
})

.showRegenrichSet = function(object) {
    dim_rawData = dim(object@rawData)
    dim_assayData = dim(object@assayData)
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
    dim_resDEA = if (is.null(object@resDEA)) {
        c(0, 0)
    } else dim(object@resDEA@pFC)
    if (dim_resDEA[1] > 0) {
        p = object@resDEA@pFC[, "p"]
        cat(" (1) ", sum(p < 0.05), " rows with differential ",
            "p-value < 0.05\n", sep = "")
    } else {
        message(" Differential expression analysis needs ",
            "to be performed.\n")
        return(invisible(NULL))
    }

    dim_network = if (is.null(object@network)) {
        c(0, 0)
    } else dim(object@network@netTable)
    if (dim_network[1] > 0) {
        # cat(' Whole network info:\n') show(object@network)
    } else {
        message("\n Network inference needs to be performed, or ",
            "a 'TopNetwork' object needs to be provided.\n",
            sep = "")
        return(invisible(NULL))
    }

    # show top p network
    cat("\n")
    dim_topNetP = if (is.null(object@topNetP)) {
        c(0, 0)
    } else dim(object@topNetP@netTable)
    if (dim_topNetP[1] > 0) {
        cat(" (2) Top p% network info:\n")
        show(object@topNetP)
    } else {
        message(" Top p% network needs to be retained.\n")
        return(invisible(NULL))
    }

    # show enrichment
    cat("\n")
    dim_resEnrich = if (is.null(object@resEnrich)) {
        c(0, 0)
    } else dim(object@resEnrich@allResult)
    if (dim_resEnrich[1] > 0) {
        cat(" (3) Enrichment info:\n")
        show(object@resEnrich)
    } else {
        message(" FET/GSEA enrichment needs to be performed.\n")
        return(invisible(NULL))
    }

    # show score
    cat("\n")
    dim_resScore = if (is.null(object@resScore)) {
        c(0, 0)
    } else dim(object@resScore)
    if (dim_resScore[1] > 0) {
        cat(" (4) RegEnrich score:\n")
        if (dim_resScore[1] > 10) {
            print(head(object@resScore, 10))
            cat("...", dim_resScore[1] - 10, "more rows ...\n")
        } else {
            print(object@resScore)
        }
    } else {
        message(" 'regenrich_rankScore' needs to be performed.\n")
        return(invisible(NULL))
    }
}
# method in show generic for 'RegenrichSet' object
setMethod(f = "show", signature = "RegenrichSet", .showRegenrichSet)

