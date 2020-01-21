sam <- function(x){
  x <- tmp
  x[is.na(x)] <- 0
  
  k <- csam(x, x)
  k[upper.tri(k, diag = FALSE)] <- NA
  diag(k) <- 0
  
  d <- as.dist(k, diag = FALSE)
  attr(d, "Labels") <- rownames(x)
  d
}