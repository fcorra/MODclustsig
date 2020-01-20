hellinger <- function(x){
  
  x[is.na(x)] <- 0
  x <- sqrt(as.matrix(x))
  
  if(is.null(rownames(x))){
    rownames(x) <- seq(1, nrow(x), by =1)
  }
  
  # Calculating the Hellinger distance matrix
  D <- matrix(rep(NA), ncol = nrow(x), nrow = nrow(x))
  combination <- combn(1:nrow(x), 2)
  
  apply(combination, 2,
  function(tmp){
    D[tmp[1], tmp[2]] <<- HellDist(x[tmp, ])
  })

  # Adding the zero diagonal
  diag(D) <- 0
  
  # Formating the output to distance table
  d <- as.dist(t(D), diag = FALSE)
  attr(d, "Labels") <- rownames(x)
  d
}
