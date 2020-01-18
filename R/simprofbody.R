simprof.body <- function(rawdata, num.expected=1000, num.simulated=999, 
                         method.cluster="average", method.distance="euclidean", originaldata=NA,
                         alpha=0.05, clust.order=NA, startrow, pMatrix, currentsamples=NA, 
                         const=const, silent, increment, undef.zero){
  
  ### The below description is incorrect. Ignore and change.
  ### Find all of the basic data we'll need to form a test statistic
  ### To do this, first find all of the pairwise distances between sites/samples. Rank them.
  ### Second, permute each column separately. Find all of the new pairwise distances. Rank them.
  ### Third, repeat this num.expected (1000) times.
  ### Fourth, compute the sum of the absolute values of the distances for each rank between the original data and the expected profile.
  ### Fifth, repeat steps 2 num.simulated (999) times. This is the simulated data under the null for comparison.
  ### Sixth, compute the test statistic between each of these profiles and the expected.
  ### Seventh, compute the p-value given the original data and the num.simulated (999) simulated values.
  ### Finally, if the p-value is less than alpha (0.05), run simprof on the left and right sub-clusters.
  
  ### Real data simprof
  rawdata.simprof <- list(genSimilarityProfile(rawdata, method.distance, const, undef.zero))
  if (!dim(rawdata.simprof[[1]])[2]==1) # order screws up the list if it is only a column (and if it is only a single column, there is no need to order)
    rawdata.simprof[[1]] <- rawdata.simprof[[1]][,order(rawdata.simprof[[1]][2,])]
  ### maybe this is the slowdown point? 
  
  ### Expected Profile
  expectedprofile.simprof <- genProfile(rawdata, originaldata, num.expected, 
                                        method.distance, const, silent=silent, 
                                        increment=increment, type="Expected",
                                        undef.zero)
  expectedprofile.average <- computeAverage(expectedprofile.simprof, num.expected)
  
  ### Test Statistic
  teststatistic <- computeTestStatistic(rawdata.simprof[[1]], expectedprofile.average)
  
  ### Simulated Profile
  simulatedprofile.simprof <- genProfile(rawdata, originaldata, num.simulated, 
                                         method.distance, const, silent=silent, 
                                         increment=increment, type="Simulated",
                                         undef.zero)
  
  ### Comparison
  pval <- tsComparison(simulatedprofile.simprof, expectedprofile.average, num.simulated, teststatistic)
  
  findings <- c();
  findings[["pval"]] <- pTracker(pMatrix, startrow, pval)
  
  if (pval > alpha){
    if (!is.na(currentsamples[1])){ ### there are some warnings here... find out why
      ### warnings should stop if we reference a specific index
      ### none should be NA if currentsamples != NA
      # pop out and start working our way back
      findings[["samples"]] <- list(currentsamples)
      return(findings)
    }
    else{ # we failed to reject the first null: all samples are one group
      findings[["samples"]] <- list(c(1:nrow(rawdata)))
      return(findings)
      ### rework the above line (eventually) to be placed in an order determined by hclust()
      ### ... maybe 
    }
  }
  
  ### Need to cut down rawdata to only be the next set of samples
  
  ### See if we need to investigate subtrees
  if (pval <= alpha){ ## should this be strict inequality?
    
    simprof.left  <- diveDeep(rawdata=rawdata, num.expected=num.expected, num.simulated=num.simulated, 
                              method.cluster=method.cluster, method.distance=method.distance,
                              originaldata=originaldata, alpha=alpha, clust.order=clust.order, startrow=startrow, pMatrix=findings$pval,
                              side="LEFT", const=const, silent=silent, increment=increment)
    simprof.right <- diveDeep(rawdata=rawdata, num.expected=num.expected, num.simulated=num.simulated, 
                              method.cluster=method.cluster, method.distance=method.distance,
                              originaldata=originaldata, alpha=alpha, clust.order=clust.order, startrow=startrow, pMatrix=findings$pval,
                              side="RIGHT", const=const, silent=silent, increment=increment)
    
    findings[["samples"]] <- append(simprof.left$samples, simprof.right$samples)
    findings[["pval"]] <- pTrackerMerge(simprof.left$pval, simprof.right$pval)
    return(findings)
  }
}
