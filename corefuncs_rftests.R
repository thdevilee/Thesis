### function to generate a split given some regressor, partitions and hyperparameters
randomsplit <- function(regs, indx, skip, nodesize, leafsize = 1){
  unq_indx <- unique(indx) ### retrieve al unique indices
  unq_indx <- unq_indx[!unq_indx %in% skip] ### remove indices that should be skipped => not splittable
  splt_indx <- rep(NA, nrow(regs)) ### initialize empty split vector
  
  for (i in 1:length(unq_indx)){ ### for all indices
    sub_indx <- indx == unq_indx[i] ### logical from some index
    sub_regs <- regs[sub_indx, ] ### subset regressor
    n_sub <- sum(as.numeric(sub_indx)) ### number of observations with some index
    if(n_sub < 2*leafsize | n_sub < nodesize) next ### check if it possible to split some index according to leaf size and node size
    splt <- TRUE ### set split to true
    fls <- 0 ### zero fails
    while(splt){
      splt_indx[sub_indx] <- 1 ### set all indices to 1 for some index
      rcol <- sample(seq(1,ncol(sub_regs)), size = 1) ### sample a regressor
      rrow <- sample(1:nrow(sub_regs), size = 1) ### sample a observation
      ineq_vec <- parse(text = paste(sub_regs[rrow, rcol], sample(c("<", "<=", ">=", ">"), size = 1), sub_regs[, rcol])) 
      ### generate inequality (randomly generated sign) based on the randomly sampled regressor and observation for that same regressor 
      cond_indx <- sapply(ineq_vec, eval) ### evaluate the above expression
      splt_indx[sub_indx][cond_indx] <- 2 ### observations that satisfy this equality get assigned a 2 (instead of a one)
      if (var(splt_indx[sub_indx]) != 0) { ### check if a split occured
        splt <- FALSE ### if so quit the while loop
      } else {fls <- fls + 1} ### otherwise add 1 to the number of fails
      
      if (fls == 25){ ### if 25 fails occured
        spltb_reg <- regsplits(sub_regs, matrix(unq_indx[i], ncol = nrow(sub_regs)), leafsize = leafsize)
        ### check for the possibility in each of the regressor to split
        splt <- any(unlist(spltb_reg)) ### if any of the regressors can be split return back to the while loop
        if (!splt) splt_indx[sub_indx] <- NA ### if none of the regressor cant be split set al elements of the split vector to NA
        
      }
    }
  }
  return(list("splt_indx" = as.character(splt_indx), "skip" = unique(indx[is.na(splt_indx)])))
}

### function to split predictor space given some index and update the indices
updatesplit <- function(regs, indx, skip, nodesize, leafsize){
  out <- indx ### copy indices onto a new variable
  tmp <- randomsplit(regs = regs, indx = indx, skip = skip, nodesize = nodesize, leafsize = leafsize) ### generate random splits
  updt_indx <- tmp[["splt_indx"]] ### extract partitions
  log_indx <- !is.na(updt_indx) ### logical for which partitions are not "NA"
  out[log_indx] <- paste0(indx[log_indx], updt_indx[log_indx]) ### update indices for partitions that are not "NA"
  return(list("updt_indx" = out, "skip" = tmp[["skip"]]))
}

### function to check whether it is possible to split regressors
regsplits <- function(regs, indx, leafsize){
  indx_names <- unique(indx) ### extract unique index names
  out <- vector(mode = "list", length = length(indx_names)) ### initialize output list per unique index
  names(out) <- indx_names ### names each element according to a unique index
  
  for (i in 1:length(indx_names)){ ### for each index
    regs_sub <- regs[(indx == indx_names[i]), ] ### subset the regressor
    tmp <- vector(mode = "list", length = ncol(regs)) ### initialize output list per regressor
    
    for (j in 1:ncol(regs)){ ### for each regressor
      cumtable <- cumsum(table(regs_sub[, j])) ### cumpute number of observations per regressor value
      indx_cond <- cumtable > leafsize ### check which element from the above is larger than the leafsize
      tmp[[j]] <- (max(cumtable) - cumtable[indx_cond][1]) >= leafsize ### checks whether it is possible to generate a split to satisfies the leafsize
      
    }
    
    out[[indx_names[i]]] <- drop(do.call(rbind, tmp)) ### coerce possible splits per regressor in to a vector
    
  }
  return(out)
}

### function to generate a tree given the regressor and hyperparameters
gensplits <- function(regs, maxdepth = NA, maxnodes = NA, nodesize = NA, leafsize = 1){
  if (is.na(maxnodes)) maxnodes <- nrow(regs) ### set default maximum nodes unbounded when not specified
  if (is.na(maxdepth)) maxdepth <- nrow(regs) ### set default maximum depth unbounded when not specified
  if (is.na(nodesize)) nodesize <- 2 ### set default nodesize to 2 when not specified
  indx <- rep("1", nrow(regs)) ### initial partition ==> all the same partition
  out <- t(matrix(indx)) ### coerce into a row matrix
  out <- rbind(out, out) ### duplicate this row ==> matrix with two identical rows
  
  i <- 1 ### index
  psbl_split <- TRUE ### set possible to split to TRUE
  skip <- NA ### empty skip list
  fls <- 0 ### zero fails
  while (max(nchar(out[i, ])) - 1 < maxdepth & ### check whether the tree still satisfies the maximum depth
         length(unique(out[i, ])) < maxnodes & ### check whether the tree still satisfies the maximum number of nodes
         !all(table(out[i, ]) < 2*leafsize) & ### check whether the tree still satisfies the leaf size
         any(table(out[i, ] >= nodesize)) & ### check whether the tree still satisfies the node size
         psbl_split & ### possible to split
         !all(out[i, ] %in% skip) ### check if all indices are in the skip list
  ){
    
    res <- updatesplit(regs = regs, indx = out[i, ], skip = skip, nodesize = nodesize, leafsize = leafsize) ### generate possible split
    tmp <- res[["updt_indx"]] ### extract the candidate update for the index
    skip <- res[["skip"]] ### extract skip list
    
    ### check wheter the candidate update meets the requirements
    if(max(nchar(tmp)) - 1 <= maxdepth & ### check whether update satisfies the maximum depth
       length(unique(tmp)) <= maxnodes & ### check whether update still satisfies the maximum number of nodes
       all(table(tmp) >= leafsize) & ### check whether update satisfies the leaf size
       !isTRUE(all.equal(out[i, ], tmp))){ ### check whether an update occured in the first place
      
      fls <- 0 ### reset fails zo zero
      out[i + 1, ] <- tmp ### update index
      i <- i + 1 ### add 1 to the index
      out <- rbind(out[1:i, ], out[i, ]) ### duplicate the last row of the matrix
    } else fls <- fls + 1 ### add one to the number of fails
    
    if (fls == 25){ ### when there are 25 fails
      spltb_indx <- out[i, ] %in% skip ### extract indices that are not in the skip list
      spltb_reg <- regsplits(regs[!spltb_indx, ], out[i, !spltb_indx], leafsize = leafsize) ### per index check for the possibility to split
      psbl_split <- any(unlist(spltb_reg)) ### check wheter it possible to split any index
    }
  }
  return(out[i, ])
}

gentrees <- function(regs, n_tree, maxdepth = NA, maxnodes = NA, nodesize = NA, leafsize = 1){
  pckgs <- as.vector(lsf.str(.GlobalEnv)) ### required for parallel computation
  out <- foreach(j = 1:n_tree, .combine = "cbind", .export = pckgs) %do% { ### generate n_tree number of trees
    splts <- gensplits(regs = regs, maxdepth = maxdepth, maxnodes = maxnodes, nodesize = nodesize, leafsize = leafsize) ### grow a tree
    sapply(splts, FUN = function(x) as.numeric(substr(x, nchar(x), nchar(x)) == 2)) ### extract partition
  }
  out <- as.matrix(out) ### coerce result to matrix
  colnames(out) <- NULL ### remove colum names from out
  rownames(out) <- NULL### remove row names from out
  return(out)
}

### generates trees and compute test statistic (permutationa of y have to be supplied)
randomRF <- function(regs, y, y_perm, n_tree, maxdepth = NA, maxnodes = NA, nodesize = NA, leafsize = 1){
  X <- gentrees(regs = regs, n_tree = n_tree, leafsize = leafsize, maxdepth = maxdepth, maxnodes = maxnodes, nodesize = nodesize)
  ### generate forest based on some paramaters
  S_t <- S(X, y) ### compute true test statistic
  S_p <- c(Inf, S(X, y_perm)) ### compute permutation test statistic
  out <- mean(S_t*(1-1e-14) <= S_p) ### compute empirical p-value
  return(out)
}
### compute numerator GT statistic based on X and y
S <- function(X, y){
  Z <- matrix(1, nrow = nrow(X)) ### null hypothesis
  X <- X - Z %*% solve(crossprod(Z), crossprod(Z, X)) ### orthogonalize X with respect to H0
  tmp <- crossprod(X, y)^2 ### compute numerator GT statistic
  out <- apply(tmp, 2, sum)
  return(drop(out))
}


### function to compute linear global test
linear.gt <- function(regs, labs){
  res <- gt(y~1, y~., data = data.frame(regs, y = labs)) ### linear global test
  return(res@result[1])  ### return p-value
}



### function toperform nested ridge
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


### function to perform cross-validation with random forest
fixedcv.rf <- function(folds, regs, labs, mtry = "sqrt(p)", nodesize = 1, n_tree = 250){
  p <- ncol(regs) ### number of regressor
  mtry <- sapply(mtry, function(x, p) eval(parse(text = x)), p = p) ### evaluate mtry in terms of p
  
  
  n_folds <- max(folds) ### total number of folds
  res <- numeric(n_folds) ### initialize outerloop output
  for(i in 1:n_folds){  ### for all folds
    preds <- randomForest::randomForest(x = regs[folds != i,],  y = factor(labs[folds != i, , drop = TRUE]),
                                        xtest = regs[folds == i,], ytest = factor(labs[folds == i, , drop = TRUE]),
                                        mtry = mtry, nodesize = nodesize,
                                        ntree = n_tree)[["test"]][["predicted"]] ### train on n_folds - 1 and make predictions on the remaining one
    res[i] <- mean(preds == labs[folds == i, ]) ### compute accuracy
  }
  return(weighted.mean(res, table(folds))) ### return weighted accuracy
}



### function to perform the fixed rf GT
fixed.rf_gt <- function(regs, labs, leafsize = 1, n_perms = 1000, n_tree = 2500, maxdepth = 3, maxnodes = NA, nodesize = NA){
  y <- matrix(as.numeric(labs[, "y"])) ### coerce labs in to a matrix
  y_perm <- sapply(1:(n_perms - 1), function(x, y) y[sample(nrow(y)), ,drop = FALSE], y = y) ### generate permutations of y
  
  res <- randomRF(regs = regs, y = y, y_perm = y_perm, n_tree = n_tree, leafsize = leafsize, 
                  maxdepth = maxdepth, maxnodes = maxnodes, nodesize = nodesize) ### compute p-value
  return(res)
}


### function to perform the p-value rf GT
pval.rf_gt <- function(regs, labs, maxdepth = NA, n_perms = 1000, n_tree = 2500, balanced = TRUE, force = FALSE, leafsize = 1, maxnodes = NA, nodesize = NA){
  i <- 1 ### arbitrarily choose a validation fold
  folds <- create_folds(n_folds = 2, labs = labs, balanced = balanced, force = force) ### create folds such that the groups are balanced
  indx_train <- folds != i ### create training observation logical based on i
  n_train <- sum(folds != i) ### compute the number of observations in each folds
  if (all(is.na(maxdepth))) maxdepth <- c(seq(1, 7), 10) ### grid of hyperparameters for tree depth
  y_train <- matrix(as.numeric(labs[indx_train, ])) ### subset training labs and coerce labs in to a matrix
  p_vals <- numeric(length(maxdepth)) ### initialize vector for p-values
  
  y_perm_train <- sapply(1:(n_perms - 1), function(x, y) y[sample(nrow(y)), ,drop = FALSE], y = y_train) ### permutations of training y
  for (i in 1:length(maxdepth)){ ### for the grid of hyperparameters
    p_vals[i] <- randomRF(regs = regs[indx_train, ], y = y_train, y_perm = y_perm_train, n_tree = 125, leafsize = leafsize, 
                          maxdepth = maxdepth[i], maxnodes = maxnodes, nodesize = nodesize) ### compute p-value on training set
    
  }
  
  optimal_depth <- maxdepth[which.min(p_vals)] ### select tree depth which corresponds to smallest p-value
  y_test <- matrix(as.numeric(labs[!indx_train, ])) ### subset test labs and coerce labs in to a matrix
  y_perm_test <- sapply(1:(n_perms - 1), function(x, y) y[sample(nrow(y)), ,drop = FALSE], y = y_test) ### permutations of test y
  res <- randomRF(regs = regs[!indx_train, ], y = y_test, y_perm = y_perm_test, n_tree = n_tree, leafsize = leafsize, 
                  maxdepth = optimal_depth, maxnodes = maxnodes, nodesize = nodesize) ### compute p-value on test set
  return(res)
}


### function to perform the sampling overall RF GT
overall.Sgt <- function(regs, labs, n_tree = 2500, maxdepth = NA, n_perms = 1000, leafsize = 1, maxnodes = NA, nodesize = NA){
  pckgs <- as.vector(lsf.str(.GlobalEnv))
  n <- nrow(regs) ### number of observations
  y <- matrix(as.numeric(labs[["y"]])) ### coerce labs in to a matrix
  y_perm <- sapply(1:(n_perms - 1), function(x, y) y[sample(nrow(y)), ,drop = FALSE], y = y) ### permutations of y
  X <- gentrees(regs = regs, n_tree = 1, leafsize = leafsize, maxdepth = ceiling(runif(1, min = 0, max = 5)), maxnodes = maxnodes, nodesize = nodesize)
  for(i in 1:(n_tree-1)){ ### generate a tree based on a sampled max depth
    X <- cbind(X, gentrees(regs = regs, n_tree = 1, leafsize = leafsize, maxdepth = ceiling(runif(1, min = 0, max = 5)), maxnodes = maxnodes, nodesize = nodesize))
  }
  S_t <- S(X, y) ### compute true GT test statistic
  S_p <- c(Inf, S(X, y_perm)) ### compute permutation GT test statistic
  out <- mean(S_t*(1-1e-14) <= S_p) ### compute empirical p-value
  return(out)
}
