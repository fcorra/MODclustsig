#' @include classes.R
NULL

#' simprof
#' 
#' A tool for determining the number of significant clusters produced using hclust() with the assumption of no a priori groups.
#'
#' @param data Input data in a matrix.
#' @param num.expected The number of similarity profiles to generate for creating the expected distribution of the data. This value should be large.
#' @param num.simulated The number of similarity profiles to generate for use in comparing the observed test statistic with its null distribution. This value should be large.
#' @param method.cluster The method of clustering to use with \code{\link{hclust}}. Standard values from \code{hclust} are \code{"ward"}, \code{"single"}, \code{"complete"}, \code{"average"}, \code{"mcquitty"}, \code{"median"} or \code{"centroid"}.
#' @param method.distance This value should be either an option to pass to the function \code{\link{dist}} (standard values are \code{"euclidean", "maximum", "manhattan", "canberra", "binary" or "minkowski"}), \code{"braycurtis"} or \code{"czekanowski"} for Czekanowski Dissimilarity (referred to as Bray-Curtis Disimilarity in some fields, particularly marine biology), or \code{"actual-braycurtis"} for the true Bray-Curtis Dissimilarity where the data are standardized before the dissimilarity is calculated. This value can also be any function which returns a \code{"dist"} object. In this version of clustsig the "hellinger" distance is also implemented. 
#' @param method.transform An option to specify a transform, if any, to be applied to the data. Possible values are \code{"identity"} (no transformation), \code{"squareroot"}, \code{"log"}, \code{"PA"}(Presence/Absence), or any numeric value (of type \code{"double"}). This transform is applied before the adjustment constant is applied, so choose a constant accordingly.
#' @param alpha The alpha level at which to reject the null hypothesis. If the null is rejected, the test continues and tests each sub-tree recursively until either all subtrees are exhausted by reaching the individual level or there are no significant distance. Due to the nature of multiple testing inherent in this process, care should be taken when choosing this alpha level.
#' @param sample.orientation The orientation of the data, either \code{"row"} or \code{"column"}. The practical effect of this is that the transpose will be examined if \code{"column"} is chosen.
#' @param const The value of the constant to be used in adjusting the Bray-Curtis Dissimilarity coefficient, if any is to be used. Any positive value of \code{"const"} will be appended as a new variable to each sample, acting as a sort of \dQuote{dummy species} (where that interpretation is appropriate).
#' @param silent A logical value indicating whether anything should be printed during the code execution. If \code{FALSE}, a message will be printed every \code{increment} (see below) number of times in the main looping procedure. This was implemented because the code can take a while to run due to many permutations and its recursive nature; however, for the same reason, many messages could be printed.
#' @param increment An integer value indicating, if \code{silent=FALSE}, one which iterations a message should be printed. (If the iteration number modulus \code{increment} equals 0, that number will be printed.)
#' @param undef.zero A logical value indicating whether undefined values arising from a denominator equal to 0 in the Bray-Curtis/Czekanowski Dissimilarity Indices should result in \code{NA} or 0. This defaults to \code{TRUE} so that NA values are replaced by 0. This default is to retain backward compatibility with the previous version of the package but may be changed in a future release.
#' @param warn.braycurtis A logical value indicating whether a warning should be printed when using the \code{"braycurtis"} option because of the naming confusion in some fields with the Czekanowski Dissimilarity Index. This defaults to \code{TRUE} but may change in future releases. For more information, see Yoshioka (2008) listed in the references.
#'
#' @description Simprof from clustsig calculates the dist matrix internally. The Hellinger distance is not implemented. This is quite inconvinient to follow S Primpke approach. To avoid writing a function each time, we will implement a new method. The 'hellinger' method.
#' 
#' @return
#' S4 object of class simprof. It has the following components:
#' \itemize{
#'   \item numgroups The number of groups which are found to bestatistically significant.
#'   \item significantclusters A list of length numgroups with each element containing the sample IDs (row/column numbers in the corresponding original data) that are in each significant cluster.
#'   \item pval The merge component from the hclust results with an extra column of p-values. These p-values are for testing whether the two groups in that row are statistically different.
#'   \item hclust An object of class hclust which is just the results of running hclust on the original data.
#' }
#' @references 
#' \itemize{
#'    \item Clarke, K.R., Somerfield, P.J., and Gorley, R.N., 2008. Testing of null hypotheses in exploratory community analyses similarity profiles and biota-environment linkage. J. Exp. Mar. Biol. Ecol. 366, 56--69.
#'    \item Yoshioka, P.M., 2008. Misidentification of the Bray-Curtis similarity index. Mar. Ecol. Prog. Ser. 368, 309--310.
#' }
#' @author Douglas Whitaker and Mary Christman.
#' 
#' @seealso \code{\link{hclust}}
#' 
#' @export
#'
#' @examples
#' \dontrun{
#' # Load the USArrests dataset included with R
#' # And use abbreviations of state names
#' # We leave out the third column because
#' # it is on a different scale
#' usarrests<-USArrests[,c(1,2,4)]
#' rownames(usarrests)<-state.abb
#' # Run simprof on the data
#' res <- simprof(data=usarrests, method.distance="braycurtis")
# Graph the result
#' pl.color <- simprof.plot(res)
#' }
simprof <- function(data, num.expected=1000, num.simulated=999, 
                    method.cluster="average", method.distance="euclidean", 
                    method.transform="identity", alpha=0.05, 
                    sample.orientation="row", const=0, 
                    silent=TRUE, increment=100, 
                    undef.zero=TRUE, warn.braycurtis=TRUE){
  
  # the basic way the data is passed is this:
  # simprof -> simprof.body -> diveDeep (splits into 'left' and 'right' subtrees -> simprof.body (same as before)
  # if you change anything (to pass it down the line), make sure you change the arguments of each of the above
  # TODO: pass arguments with ... instead of explicitly.
  
  if (!is.matrix(data)) # TODO: this could probably be re-worked to not be matrix-dependent
    data <- as.matrix(data) ### make it consistent handling of the data
  
  rawdata<-data # the user doesn't need to see "rawdata", but changing a lot of code isn't desirable right now
  ### the code is written for rows to be the samples/sites, so if it is the other way take the transpose
  if(sample.orientation == "column")
    rawdata<-t(rawdata)
  ### now we should either have a function here or something to pass to dist()	
  if (is.function(method.distance)) # allow arbitrary distance function choice
    rawdata.dist <- method.distance(rawdata)
  else if (method.distance == "braycurtis"){
    if (warn.braycurtis){
      warning("This version of the Bray-Curtis index does not use standardization.", call.=FALSE)
      warning("To use the standardized version, use \"actual-braycurtis\".", call.=FALSE)
      warning("See the help documentation for more information.", call.=FALSE)
    }
    rawdata.dist <- braycurtis(rawdata, const, undef.zero)
    # the next bit takes care of putting row names back onto the rawdata.dist
    if (!is.null(rownames(rawdata))){
      attr(rawdata.dist, "Labels") <- rownames(rawdata)
    }
  }
  else if (method.distance == "czekanowski"){
    rawdata.dist <- czekanowski(rawdata, const, undef.zero) # this is IDENTICAL to the braycurtis() function
    if (!is.null(rownames(rawdata))){
      attr(rawdata.dist, "Labels") <- rownames(rawdata)
    }
  }
  else if (method.distance == "actual-braycurtis"){
    rawdata.dist <- braycurtis(preBCstandardization(rawdata),const, undef.zero)
    if (!is.null(rownames(rawdata))){
      attr(rawdata.dist, "Labels") <- rownames(rawdata)
    }
  }
  else {
    rawdata.dist <- dist(rawdata, method=method.distance)
  }
  ### transforming the data
  if (!method.transform == "identity")
    rawdata <- trans(rawdata, method.transform)
  
  hclust.results <- hclust(rawdata.dist, method=method.cluster) ### Generate the cluster information
  pMatrix <- cbind(hclust.results$merge, matrix(data = NA, nrow=nrow(hclust.results$merge), ncol=1))
  
  ### Namely one that includes a vector indicating which significant cluster each sample belongs to
  simprof.results <- simprof.body(rawdata=rawdata, num.expected=num.expected, num.simulated=num.simulated, 
                                  method.cluster=method.cluster, method.distance=method.distance, 
                                  originaldata=rawdata, alpha=alpha, clust.order=hclust.results$merge, 
                                  startrow=nrow(hclust.results$merge), pMatrix=pMatrix, 
                                  const=const, silent=silent, increment=increment, undef.zero)
  
  results <- list()
  results[["numgroups"]] <- length(simprof.results$samples) # number of significant groups
  results[["pval"]] <- simprof.results$pval
  results[["hclust"]] <- hclust.results
  if(!is.null(rownames(rawdata))){
    for (i in 1:length(simprof.results$samples)){
      for(j in 1:length(simprof.results$samples[[i]])){
        simprof.results$samples[[i]][j] <- rownames(rawdata)[as.numeric(simprof.results$samples[[i]][j])]
      }
    }
  }
  results[["significantclusters"]] <- simprof.results$samples
  
  out <- new("simprof",
             numgroups = results[["numgroups"]],
             significantclusters = results[["significantclusters"]],
             pval = results[["pval"]],
             hclust = results[["hclust"]])
  
  return(out)
}