#' Filter the top edges in the network.
#' @description Choose the top edges in the reulator-target network.
#'
#' @param network A data.frame of from.node (regulators),
#' to.node (target genes) and weight.
#' @param percent Numeric. The percentage (<= 10) of top edges, default is 5,
#' which means top 5\% of Reg-target network edges will be retained.
#' @param directed The type of the network: directed (\code{directed = TRUE})
#' and undirected
#' (\code{directed = FALSE}). Default is \code{TRUE}.
#' @param reg The gene regulators. The default is \code{NULL}.
#' @return A list of cutWeight (the cutoff edge weight),
#' net (the genes connected by
#' regulators), and validReg (valid regulators).
#' @rawNamespace import(data.table, except = c(melt, dcast))
# @export
topNet = function(network, percent = 5, directed = TRUE, reg = NULL) {
    stopifnot(identical(colnames(network), c("from.gene", "to.gene", 
        "weight")))
    stopifnot(percent <= 10)
    
    if (is.data.frame(network)) {
        network = as.data.table(network)
    }
    
    # from.gene are only from reg
    if (!is.null(reg)) {
        network = network[network$from.gene %in% reg]
    } else {
        reg = unique(as.character(network$from.gene))
    }
    
    if (!directed) {
        # reverse from.gene and to.gene
        networkRev = network
        networkRev$from.gene = network$to.gene
        networkRev$to.gene = network$from.gene
        networkRev = network[network$from.gene %in% reg]
        
        # combine edges of two directions
        network = rbind(network, networkRev)
        # remove duplicate
        id = !duplicated(network)
    }
    
    networkReg = subset(network, !duplicated(network))
    n = nrow(networkReg)
    if (n == 0) {
        stop("No edges were selected! Please check from.gene and reg")
    }
    
    # top percent % edges of the network
    cutWeight = stats::quantile(networkReg$weight, 1 - percent/100)
    networkRegTop = networkReg[networkReg$weight >= cutWeight]
    
    reg = as.character(unique(networkRegTop$from.gene))
    networkRegTop$from.gene = as.character(networkRegTop$from.gene)
    networkRegTop$to.gene = as.character(networkRegTop$to.gene)
    
    netTable = subset(networkRegTop, networkRegTop$from.gene %in% 
        reg)
    
    res = list(netTable = netTable)
    return(res)
}


