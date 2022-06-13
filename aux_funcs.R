EmpericalPower <- function(data){
  out <- apply(data, 2, function(x) mean(x < 0.05))
  return(out)
}

rename_folds <- function(folds){
  digits <- sort(unique(folds)) ### sort all unique fold ids
  for (i in 1:length(digits)){ ### for each element in unique sorted folds
    folds[folds == digits[i]] <- i ### rename to 1, .., k (k = length(unique sorted folds))
  }
  return(folds) 
}

is.int <- function(x) x%%1 == 0

create_folds <- function(n_folds, labs, balanced = TRUE, force = FALSE, ...){
  if(balanced){
    if (any(order(as.numeric(labs[["y"]])) != seq(1, nrow(labs)))) stop("Data (labels) is not sorted") ### checks if the above requirments are met
    n1 <- nrow(labs) - sum(as.numeric(labs[["y"]]))
    n2 <- sum(as.numeric(labs[["y"]]))
    reps1 <- n1/n_folds
    reps2 <- n2/n_folds
    if ((!is.int(reps1) | !is.int(reps2)) & !force) stop("Cannot balance folds")
    folds1 <- rep(seq(1, n_folds), ceiling(reps1))[1:n1] ### create n fold indexes for group 1
    folds2 <- rep(seq(1, n_folds), ceiling(reps2))[1:n2] ### create n fold indexes for group 2
    folds <- c(sample(folds1), sample(folds2)) ### randomly assign indexes to each observation in a balanced fashion
  } else{
    folds <- sample(rep(seq(1, n_folds), nrow(labs)/n_folds)) ### randomly assign folds if balanced is false
  }
  return(folds)
}