setOldClass("hclust")

#' simprof_class
#'
#' @slot numgroups The number of groups which are found to bestatistically significant.. 
#' @slot significantclusters A list of length numgroups with each element containing the sample IDs (row/column numbers in the corresponding original data) that are in each significant cluster. 
#' @slot pval The merge component from the hclust results with an extra column of p-values. These p-values are for testing whether the two groups in that row are statistically different.
#' @slot hclust The S3 object returned by hclust.
#'
#' @return
#' S4 object.
#' @export
#' @examples
#' NULL
setClass("simprof",
         slots = c(numgroups = "integer",
                   significantclusters = "list",
                   pval = "matrix",
                   hclust = "hclust")
         )