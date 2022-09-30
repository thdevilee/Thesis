comp_dist <- function(X){
  n <- nrow(X)
  nom_indx <- unlist(lapply(X, function(x) as.logical(is.factor(x) + is.character(x))))
  dist_nom <- dist_num <- matrix(0, nrow = n, ncol = n)
  
  if(ncol(X) - sum(as.numeric(nom_indx)) > 0){
    X_num <- scale(X[, !nom_indx])
    dist_num <- as.matrix(dist(X_num, diag = TRUE, upper = TRUE))^2
  }
  
  if(sum(as.numeric(nom_indx)) > 0){
    X_nom <- X[, nom_indx]
    for(i in 1:(n - 1)){
      for(j in (i + 1):n){
        dist_nom[i, j] <- sum(X_nom[i,] != X_nom[j,])
      }
    }
    dist_nom[lower.tri(dist_nom)] <- dist_nom[upper.tri(dist_nom)]
  }
  
  return(sqrt(dist_num + dist_nom))
} 

### function to compute knn ordering matrix for a given distance matrix
knn_mat <- function(dist_mat){
  n <- nrow(dist_mat) ### number of observations
  nn <- (n + 1) - t(apply(dist_mat, MARGIN = 1, FUN = rank)) ### ordering distance per observation
  return(nn)
}


InvPropWeight <- function(P, y){
  data <- data.frame(P, y = y)
  Z <- model.matrix(y ~ 1, data = data)
  X <- model.matrix(y ~ 0 + ., data = data)
  m <- ncol(Z)
  n <- nrow(Z)
  y <- residuals(lm(y ~ 1, data = data))
  sumyy <- sum(y*y)
  X <- X - Z %*% solve(crossprod(Z), crossprod(Z, X))
  csm <- colSums(X*X)
  X[,csm < max(csm)*1e-14] <- 0
  xy <- crossprod(X, y)
  S <- sum(xy * xy) / sumyy
  if (sumyy == 0) S <- 0
  lams <- eigen(crossprod(X), symmetric = TRUE, only.values=TRUE)$values
  if (length(lams) < n) lams <- c(lams, numeric(n-length(lams)))
  lams[1:(n-m)] <- lams[1:(n-m)] - S
  tr.term <- crossprod(X)
  var.num <- 2*sum(tr.term*tr.term)
  mu.num <- sum(lams) + (n-m)*S
  mu.den <- n-m
  var.den <- 2*(n-m)
  cov.term <- 2*mu.num
  VarS <- (mu.num^2/mu.den^2) * (var.num/(mu.num^2)  + var.den/(mu.den^2) - 2*cov.term/(mu.num*mu.den))
  weight <- 1/sqrt(VarS)
  if (is.nan(weight) | weight > 10^6) weight <- 0
  return(weight)
}

XXT <- function(ord_mat, weights = NA){
  weighting <- TRUE ### weighting true by default
  k <- dim(ord_mat)[1] ### number of observations/partition matrices
  if(all(is.na(weights))){ ### if no weights specified
    weighting <- FALSE ### set weighting to FALSE
    weights <- rep(1, k) ### set al weights to 1
  }
  if(weighting & k != length(weights)) stop("Incorrect dimensions weights") ### check if the dimensions of the weights match the number of partition matrices
  
  weights <- cumsum(rev(weights)^2) ### compute cumulative squared weights
  out <- matrix(data = NA, nrow = k, ncol = k) ### initialize output matrix
  diag(out) <- sum(weights) ### compute diagonal
  for(i in 1:(k - 1)){ ### compute all elements in the upper triangle (without diagonal) of the output matrix
    for(j in (i + 1):k){
      out[i, j] <- sum(weights[pmin(ord_mat[i,], ord_mat[j,])])
    }
  }
  out[lower.tri(out)] <- t(out)[lower.tri(out)] ### reflect upper triangle onto lower triangle
  return(out)
}


### All functions below require specification of  the regressors (regs) as data frame and label (labs) as data frame. Some functions also require the fold id (folds) as vector

### function to compute linear global test
linear.gt <- function(regs, labs){
  res <- gt(y~1, y~., data = data.frame(regs, y = labs)) ### linear global test
  return(res@result[1])  ### return p-value
}

### function to compute perform nested ridge
nestedcv.ridge <- function(folds, regs, labs){
  regs <- apply(regs, 2, as.numeric) ### makes sure all columns are numeric
  n_folds <- max(folds) ### number of folds
  res_outer <- numeric(n_folds) ### initialize output 
  for(i in 1:n_folds){ 
    log_indx <- folds != i ### create logical index for subset of the data
    regs_sub <- as.matrix(regs[log_indx, ]) ### select regressor and put in right format
    labs_sub <- labs[log_indx, , drop = TRUE] ### select labels
    folds_sub <- rename_folds(folds[log_indx]) ### rename folds to 1, ... k - 1 (required for cv.glment)
    res_cv <- glmnet::cv.glmnet(x = regs_sub, y = labs_sub, type.measure = "mse", foldid = folds_sub, alpha = 0) ### compute optimal lambda
    preds <- predict(res_cv, as.matrix(regs[folds == i,]), s = "lambda.min") >= 0.5 ### make predictions for the validation set from a model with the optimal lambda
    res_outer[i] <- mean(preds == labs[folds == i, , drop = TRUE]) ### compute average accuracy on outer validation set
  }
  return(weighted.mean(res_outer, table(folds)))
}

### function to performed nested knn based on a specified number of neighbours
nestedcv.knn <- function(folds, regs, labs){
  data <- data.frame(regs, labs)
  data[, "y"] <- factor(data[, "y"])
  n_folds <- max(folds) ### total number of folds
  res_outer <- numeric(n_folds) ### initialize outerloop output
  fold_id <- seq(1, n_folds) ### create fold numbers for innerloop
  n_trin <- (n_folds - 2)*min(table(folds))[[1]]
  neighbours <- round(seq(1, n_trin, length.out = round(sqrt(n_trin))))
  for(i in 1:n_folds){  ### for all folds
    res_inner <- matrix(NA, nrow = n_folds - 1, ncol = length(neighbours)) ### initialize innerloop output matrix
    for(j in 1:(n_folds - 1)){ ### for all folds - 1 (excluding one outer validations set)
      for(k in 1:length(neighbours)){ ### for every neighbour
        indx <- fold_id[fold_id != i]  ### exclude outer validation set
        model_inner <- caret::knn3(formula = y ~ ., data = data[(folds %in% indx[-j]), ], k = neighbours[k]) ### train on n_folds - 2
        preds_inner <- predict(model_inner, data[(folds == indx[j]),], type = "class")
        res_inner[j, k] <- mean(preds_inner != labs[folds == indx[j], , drop = TRUE]) ### validate on innerloop validation set
      }
    }
    res_inner <- apply(res_inner, MARGIN = 2, mean) ### compute average misclassification rate over the innerloop folds
    optimal_k <- neighbours[which.min(res_inner)]
    model_outer <- caret::knn3(formula = y ~ ., data = data[folds != i, ], k = optimal_k) ### train on n_folds - 1
    preds_outer <- predict(model_outer, data[folds == i, ], type = "class") ### train on n_folds - 1 for the k with the smallest misclassification rate
    res_outer[i] <- mean(preds_outer == labs[folds == i, , drop = TRUE])  ### compute average accuracy on outer validation set
  }
  return(weighted.mean(res_outer, table(folds)))
}

### function to perform global test on a partition matrix which corresponds to the (smallest) fraction of labels times the number of observations (rounded)
prop.knn_gt <- function(labs, ord_mat){
  prop <- mean(labs[["y"]]) ### compute fraction of labels
  prop <- min(prop, 1 - prop) ### select the smallest fraction ### compute ordering mat
  P <- ifelse(ord_mat > nrow(ord_mat) - prop*nrow(ord_mat), 1, 0) ### compute partition matrix corresponding to the smallest fraction of labels
  res <- gt(y ~ 1, y ~ ., data = cbind(P, labs))@result[1] ### compute gt statistic based on partition matrix
  return(res)
}

### function to perform the global test on the partition matrix of the validation set for the value of k which corresponds to the smallest p-value in the training set (splits data in two)
pval.knn_gt <- function(regs, labs, balanced = TRUE, force = FALSE){
  i <- 1 ### arbitrarily choose a validation fold
  folds <- create_folds(n_folds = 2, labs = labs, balanced = balanced, force = force) ### create folds such that the groups are balanced
  n_train <- sum(folds != i) ### compute the number of observations in each folds
  neighbours <- round(seq(2, n_train - 1, length.out = round(sqrt(n_train))))
  dist_mat <- comp_dist(regs[folds != i, ]) ### compute distance matrix
  ord_mat <- knn_mat(dist_mat) ### compute ordering mat
  p_vals <- numeric(length(neighbours)) ### initialize vector for p-values
  for (j in 1:length(neighbours)){ ### for all 2, ..., k neighbours
    P <- ifelse(ord_mat > (n_train - neighbours[j]), 1, 0) ### compute partition matrix and save in list
    p_vals[j] <- log(globaltest::gt(y ~ 1, y ~ ., data = cbind(P, labs[folds != i, ]))@result[1]) ### compute gt statistic based on partition matrix
  }
  optimal_k <- neighbours[which.min(p_vals)] ### select k which corresponds to smallest p-value
  n_test <- sum(folds == i) ### compute number of observations in validation set
  dist_mat <- comp_dist(regs[folds == i, ]) ### compute distance matrix
  ord_mat <- knn_mat(dist_mat) ### compute ordering mat
  P <- ifelse(ord_mat > (n_test - optimal_k), 1, 0) ### compute partition matrix for optimal k
  res <- gt(y ~ 1, y ~ ., data = cbind(P, labs[folds == i, ]))@result[1] ### perform global test
  return(res)
}

overallordinal.gt <- function(labs, ord_mat){
  X  <- apply(ord_mat, 2, rank)
  res <- globaltest::gt(y ~ 1, y ~ ., data = data.frame(X, y = labs))@result[1]
  return(res)
}




S_XXT <- function(XXt, y){
  out <- apply(y, 2, FUN = function(y, XXt) t(y) %*% XXt %*% y, XXt = XXt)
  return(drop(out))
}

GT_p <- function(ord_mat, y, n_perms, weights = NA){
  XXt <- XXT(ord_mat = ord_mat, weights = weights)
  S_t <- S_XXT(XXt = XXt, y = y)
  y_perm <- sapply(1:(n_perms - 1), function(x, y) y[sample(nrow(y)), ,drop = FALSE], y = y)
  S_p <- c(Inf, S_XXT(XXt, y = y_perm))
  res <- mean(S_t*(1-1e-14) <= S_p)
  return(res)
}

overall.gt <- function(labs, ord_mat, n_perms, invweighted = FALSE){
  n <- nrow(ord_mat) ### compute number of observations
  y <- matrix(residuals(lm(y ~ 1, data = labs)), nrow = nrow(labs))
  weights <- rep(NA, n)
  if(invweighted){
    for(i in 1:n){ ### for all 1, ..., k neighbours
      P <- ifelse(ord_mat > (n - i), 1, 0) ### compute partition matrix
      weights[i] <- InvPropWeight(P, labs)
    }
  }
  res <- GT_p(ord_mat = ord_mat, y = y, n_perms = n_perms)
  
  if(invweighted){ ### if inverse weights
    res <- cbind(res, GT_p(ord_mat = ord_mat, y = y, n_perms = n_perms, weights = weights)) #### perform global test on inversely weighted partition matrices
  }
  return(res)
}
