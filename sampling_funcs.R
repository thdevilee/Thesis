### function to take a sample that follows a checkerboard pattern. Here p represents the number of dimensions, n the number of observations per position on the checkerboard, n_clusts the number of clusters (use an uneven number, an even number gives a stripes pattern)
CheckerSamp <- function(n, n_clusts, p = 2, c = 1){
  data <- matrix(runif((n_clusts^2)*n * p), nrow = (n_clusts^2)*n, ncol = p) ### draw from uniform
  shift <- expand.grid(seq(0, n_clusts - 1), seq(0, n_clusts - 1)) ### compute shift
  shift <- do.call(rbind, replicate(n, shift, simplify = FALSE))*c ### replicate shift such that each cluster has n obs
  labs <- rep(c(FALSE, TRUE), ceiling(n_clusts^2/2))[1:(n_clusts^2)] ### generate labs
  labs <- rep(labs, n) ### replicate labs such that every obs has a label according to its cluster
  data <- data + shift ### compute shifted regressors
  data <- cbind(data, labs) ### combine regressors and labels into matrix
  names(data) <- c("x1", "x2", "y") ### rename columns
  indx <- order(data[["y"]])
  data <- data[indx, ]
  return(data)
}

LDASamp <- function(n_samples, p, effect){
  inv_vcov <- solve(diag(p)) ### inverse of the covariance matrix (mahalanobis distance) (must be identity in this case)
  e_vec <- rep(1,p) ### e vector (mahalanobis distance)
  c <- sqrt( effect / (n_samples * (e_vec %*% inv_vcov %*% e_vec)) ) ### compute constant as a function of the effect size (mahalanobis distance)
  shift <- matrix(effect*c, nrow=n_samples/2, ncol=p, byrow = TRUE) ### compute the shift for one class
  labs <- sample(rep(c(TRUE, FALSE), n_samples/2)) ### randomly sample labels (balanced)
  noise <- matrix(rnorm(n_samples*p), ncol = p, nrow = n_samples) ### generate regressor (independent standard normal)
  noise[labs,] <- noise[labs,] + shift ### compute shifted regressors 
  indx <- order(labs) ### sort labels (required for power_knn func)
  labs <- labs[indx] ### order labels
  noise <- noise[indx, ] ### order regressor in the same sequence
  output <- data.frame(noise, labs) ### concatenate regressors and labels into dataframe
  names(output) <- c(paste0("x", 1:p), "y") ### change column names
  return(output)
}