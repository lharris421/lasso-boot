kl <- function(r1, r2) {
  # Estimate densities
  d1 <- density(r1)
  d2 <- density(r2)
  
  # KL divergence calculation needs probability, so convert density to probability
  p1 <- normal_density$y / sum(normal_density$y)
  p2 <- exponential_density$y / sum(exponential_density$y)
  
  # Compute KL divergence
  kl <- KL.divergence(p1, p2, k = 10)  # Using k=10 for KNN
  
  return(kl)
}