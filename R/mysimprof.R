# Simprof from clustsig calculates the dist matrix internally.
# Ofcourse the Hellinger distance is not implemented.
# This is quite inconvinient to follow S Primpke approach.
# To avoid re-writing the code, we will implement a new method. 
# The 'hellinger' method.
mysimprof <- function (data, num.expected = 1000, num.simulated = 999, method.cluster = "average", 
                       method.distance = "euclidean", method.transform = "identity", 
                       alpha = 0.05, sample.orientation = "row", const = 0, 
                       silent = TRUE, increment = 100, undef.zero = TRUE, 
                       warn.braycurtis = TRUE, HellingerHack = NULL)
{
  if (!is.matrix(data)) 
    data <- as.matrix(data)
  rawdata <- data
  if (sample.orientation == "column") 
    rawdata <- t(rawdata)
  if (is.function(method.distance)) 
    rawdata.dist <- method.distance(rawdata)
  else if (method.distance == "braycurtis") {
    if (warn.braycurtis) {
      warning("This version of the Bray-Curtis index does not use standardization.", 
              call. = FALSE)
      warning("To use the standardized version, use \"actual-braycurtis\".", 
              call. = FALSE)
      warning("See the help documentation for more information.", 
              call. = FALSE)
    }
    rawdata.dist <- braycurtis(rawdata, const, undef.zero)
    if (!is.null(rownames(rawdata))) {
      attr(rawdata.dist, "Labels") <- rownames(rawdata)
    }
  }
  else if (method.distance == "czekanowski") {
    rawdata.dist <- czekanowski(rawdata, const, undef.zero)
    if (!is.null(rownames(rawdata))) {
      attr(rawdata.dist, "Labels") <- rownames(rawdata)
    }
  }
  else if (method.distance == "actual-braycurtis") {
    rawdata.dist <- braycurtis(preBCstandardization(rawdata), 
                               const, undef.zero)
    if (!is.null(rownames(rawdata))) {
      attr(rawdata.dist, "Labels") <- rownames(rawdata)
    }
  }
  else if (method.distance == "hellinger"){
    hclust.results <- HellingerHack
  }
  else {
    rawdata.dist <- dist(rawdata, method = method.distance)
  }
  if (!method.transform == "identity") 
    rawdata <- trans(rawdata, method.transform)
  
  if(method.distance != "hellinger"){
    hclust.results <- hclust(rawdata.dist, method = method.cluster)
  }
  pMatrix <- cbind(hclust.results$merge, matrix(data = NA, 
                                                nrow = nrow(hclust.results$merge), ncol = 1))
  simprof.results <- simprof.body(rawdata = rawdata, num.expected = num.expected, 
                                  num.simulated = num.simulated, method.cluster = method.cluster, 
                                  method.distance = method.distance, originaldata = rawdata, 
                                  alpha = alpha, clust.order = hclust.results$merge, startrow = nrow(hclust.results$merge), 
                                  pMatrix = pMatrix, const = const, silent = silent, increment = increment, 
                                  undef.zero)
  results <- list()
  results[["numgroups"]] <- length(simprof.results$samples)
  results[["pval"]] <- simprof.results$pval
  results[["hclust"]] <- hclust.results
  if (!is.null(rownames(rawdata))) {
    for (i in 1:length(simprof.results$samples)) {
      for (j in 1:length(simprof.results$samples[[i]])) {
        simprof.results$samples[[i]][j] <- rownames(rawdata)[as.numeric(simprof.results$samples[[i]][j])]
      }
    }
  }
  results[["significantclusters"]] <- simprof.results$samples
  return(results)
}
