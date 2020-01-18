hellinger <- function(x){
  # first set all NA to 0
  x[is.na(x)] <- 0
  x <- sqrt(as.matrix(x))
  
  if(is.null(names(x))){
    rownames(x) <- seq(1, nrow(x), by =1)
  }
  
  # create matrix of 2 subsets based on rownumber
  # 1 first the diagonal with
  D <- cbind(matrix(rep(1:nrow(x),each=2),nrow=2),combn(1:nrow(x), 2))
  
  # create a dataframe with hellinger distances
  B <- data.frame(first=rownames(x)[D[1,]],
                  second=rownames(x)[D[2,]],
                  distance=apply(D, 2, function(y) HellDist(x[ y,]))
  )
  
  # reshape dataframe into a matrix with users on x and y axis
  B <- reshape(B, direction="wide", idvar="second", timevar="first")
  
  # convert wide table to distance table object
  d <- as.dist(B[,-1], diag = FALSE)
  attr(d, "Labels") <- B[, 1]
  return(d)
}
