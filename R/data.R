#' Example RNAseq dataset [Human]
#'
#' Data from an RNA sequencing experiment on peripheral mononuclear blood
#' cells (PBMC) of Lyme disease patients against healthy controls.
#' It contains a gene expression (FPKM) table (data frame) and a
#' sample information table (data frame).
#'
#' @docType data
#'
#' @usage data(Lyme_GSE63085)
#'
#' @format A list of 2 elements: FPKM and sampleInfo. FPKM is the 'Fragments
#' Per Kilobase of transcript per Million mapped reads' data, which is a
#' 5000 (genes) * 52 (samples) data frame. sampleInfo is the information of
#' samples, which is 52 (samples) * 9 (features) data frame. The full version
#' of FPKM table contains 23615 rows, which can be downloaded from GEO
#' database.
#'
#' @keywords datasets
#'
#' @references Bouquet et al. (2016) mBio 7(1): e00100-16
#' (\href{https://www.ncbi.nlm.nih.gov/pubmed/26873097}{PubMed})
#'
#' @source
#' \href{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE63085}{URL}
#'
"Lyme_GSE63085"



#' Human gene regulators
#'
#' The transcrpiton factors and co-factors in humans are considered the
#' regulators in RegEnrich. And these regulators are obtained from
#' (Han et al. 2015; Marbach et al. 2016; and Liu et al. 2015).
#'
#' @docType data
#'
#' @usage data(TFs)
#'
#' @format An object of 2-column \code{data.frame}; The first column is
#' ENSEMBL ID of gene regulators. The second column is gene name of gene
#' regulators. The row name of this data frame is identical to the
#' ENSEMBL ID column.
#'
#' @keywords datasets
#'
#' @references Han et al. (2015) Scientific Reports, 5:11432
#' (\href{https://www.ncbi.nlm.nih.gov/pubmed/26066708}{PubMed}),
#' Liu et al. (2015) Database, bav095
#' (\href{https://www.ncbi.nlm.nih.gov/pubmed/26424082}{PubMed}),
#' Marbach et al. (2016) Nature Methods, 13(4):366-70
#' (\href{https://www.ncbi.nlm.nih.gov/pubmed/26424082}{PubMed}).
"TFs"

