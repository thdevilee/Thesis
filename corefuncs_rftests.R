randomsplit <- function(regs, indx, skip, nodesize, leafsize = 1){
  unq_indx <- unique(indx)
  unq_indx <- unq_indx[!unq_indx %in% skip]
  splt_indx <- rep(NA, nrow(regs))
  
  for (i in 1:length(unq_indx)){
    sub_indx <- indx == unq_indx[i]
    sub_regs <- regs[sub_indx, ]
    n_sub <- sum(as.numeric(sub_indx))
    if(n_sub < 2*leafsize | n_sub < nodesize) next
    splt <- TRUE
    fls <- 0
    while(splt){
      splt_indx[sub_indx] <- 1
      rcol <- sample(seq(1,ncol(sub_regs)), size = 1)
      rrow <- sample(1:nrow(sub_regs), size = 1)
      ineq_vec <- parse(text = paste(sub_regs[rrow, rcol], sample(c("<", "<=", ">=", ">"), size = 1), sub_regs[, rcol]))
      cond_indx <- sapply(ineq_vec, eval)
      splt_indx[sub_indx][cond_indx] <- 2
      if (var(splt_indx[sub_indx]) != 0) {
        splt <- FALSE
      } else {fls <- fls + 1}
      
      if (fls == 25){
        spltb_reg <- regsplits(sub_regs, matrix(unq_indx[i], ncol = nrow(sub_regs)), leafsize = leafsize)
        splt <- any(unlist(spltb_reg))
        if (!splt) splt_indx[sub_indx] <- NA
        
      }
    }
  }
  return(list("splt_indx" = as.character(splt_indx), "skip" = unique(indx[is.na(splt_indx)])))
}


updatesplit <- function(regs, indx, skip, nodesize, leafsize){
  out <- indx
  tmp <- randomsplit(regs = regs, indx = indx, skip = skip, nodesize = nodesize, leafsize = leafsize)
  updt_indx <- tmp[["splt_indx"]]
  log_indx <- !is.na(updt_indx)
  out[log_indx] <- paste0(indx[log_indx], updt_indx[log_indx])
  return(list("updt_indx" = out, "skip" = tmp[["skip"]]))
}

regsplits <- function(regs, indx, leafsize){
  indx_names <- unique(indx)
  out <- vector(mode = "list", length = length(indx_names))
  names(out) <- indx_names
  
  for (i in 1:length(indx_names)){
    regs_sub <- regs[(indx == indx_names[i]), ]
    tmp <- vector(mode = "list", length = ncol(regs))
    
    for (j in 1:ncol(regs)){
      cumtable <- cumsum(table(regs_sub[, j]))
      indx_cond <- cumtable > leafsize
      tmp[[j]] <- (max(cumtable) - cumtable[indx_cond][1]) >= leafsize
      
    }
    
    out[[indx_names[i]]] <- drop(do.call(rbind, tmp))
    
  }
  return(out)
}

gensplits <- function(regs, maxdepth = NA, maxnodes = NA, nodesize = NA, leafsize = 1){
  if (is.na(maxnodes)) maxnodes <- nrow(regs)
  if (is.na(maxdepth)) maxdepth <- nrow(regs)
  if (is.na(nodesize)) nodesize <- 2
  indx <- rep("1", nrow(regs))
  out <- t(matrix(indx))
  out <- rbind(out, out)
  
  i <- 1
  psbl_split <- TRUE
  skip <- NA
  fls <- 0
  while (max(nchar(out[i, ])) - 1 < maxdepth &
         length(unique(out[i, ])) < maxnodes &
         !all(table(out[i, ]) < 2*leafsize) &
         any(table(out[i, ] >= nodesize)) &
         psbl_split &
         !all(out[i, ] %in% skip)
  ){
    
    res <- updatesplit(regs = regs, indx = out[i, ], skip = skip, nodesize = nodesize, leafsize = leafsize)
    tmp <- res[["updt_indx"]]
    skip <- res[["skip"]]
    
    if(max(nchar(tmp)) - 1 <= maxdepth &
       length(unique(tmp)) <= maxnodes &
       all(table(tmp) >= leafsize) &
       !isTRUE(all.equal(out[i, ], tmp))){
      
      fls <- 0
      out[i + 1, ] <- tmp
      i <- i + 1
      out <- rbind(out[1:i, ], out[i, ])
    } else fls <- fls + 1
    
    if (fls == 25){
      spltb_indx <- out[i, ] %in% skip
      spltb_reg <- regsplits(regs[!spltb_indx, ], out[i, !spltb_indx], leafsize = leafsize)
      psbl_split <- any(unlist(spltb_reg))
    }
  }
  return(out[i, ])
}

gentrees <- function(regs, n_tree, maxdepth = NA, maxnodes = NA, nodesize = NA, leafsize = 1){
  pckgs <- as.vector(lsf.str(.GlobalEnv))
  out <- foreach(j = 1:n_tree, .combine = "cbind", .export = pckgs) %do% { ################## change back
    splts <- gensplits(regs = regs, maxdepth = maxdepth, maxnodes = maxnodes, nodesize = nodesize, leafsize = leafsize)
    sapply(splts, FUN = function(x) as.numeric(substr(x, nchar(x), nchar(x)) == 2))
  }
  out <- as.matrix(out)
  colnames(out) <- NULL
  rownames(out) <- NULL
  return(out)
}

randomRF <- function(regs, y, y_perm, n_tree, maxdepth = NA, maxnodes = NA, nodesize = NA, leafsize = 1){
  X <- gentrees(regs = regs, n_tree = n_tree, leafsize = leafsize, maxdepth = maxdepth, maxnodes = maxnodes, nodesize = nodesize)
  S_t <- S(X, y)
  S_p <- c(Inf, S(X, y_perm))
  out <- mean(S_t*(1-1e-14) <= S_p)
  return(out)
}

S <- function(X, y){
  Z <- matrix(1, nrow = nrow(X))
  X <- X - Z %*% solve(crossprod(Z), crossprod(Z, X))
  tmp <- crossprod(X, y)^2
  out <- apply(tmp, 2, sum)
  return(drop(out))
}


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
  return(mean(res_outer))
}



fixedcv.rf <- function(folds, regs, labs, mtry = "sqrt(p)", nodesize = 1, n_tree = 250){
  p <- ncol(regs)
  mtry <- sapply(mtry, function(x, p) eval(parse(text = x)), p = p)
  
  
  n_folds <- max(folds) ### total number of folds
  res <- numeric(n_folds) ### initialize outerloop output
  for(i in 1:n_folds){  ### for all folds
    preds <- randomForest::randomForest(x = regs[folds != i,],  y = factor(labs[folds != i, , drop = TRUE]),
                                        xtest = regs[folds == i,], ytest = factor(labs[folds == i, , drop = TRUE]),
                                        mtry = mtry, nodesize = nodesize,
                                        ntree = n_tree)[["test"]][["predicted"]] ### train on n_folds - 2
    res[i] <- mean(preds == labs[folds == i, ])
  }
  return(mean(res))
}




fixed.rf_gt <- function(regs, labs, leafsize = 1, n_perms = 1000, n_tree = 2500, maxdepth = 3, maxnodes = NA, nodesize = NA){
  y <- matrix(as.numeric(labs[, "y"]))
  y_perm <- sapply(1:(n_perms - 1), function(x, y) y[sample(nrow(y)), ,drop = FALSE], y = y)
  
  res <- randomRF(regs = regs, y = y, y_perm = y_perm, n_tree = n_tree, leafsize = leafsize, 
                  maxdepth = maxdepth, maxnodes = maxnodes, nodesize = nodesize)
  return(res)
}



pval.rf_gt <- function(regs, labs, maxdepth = NA, n_perms = 1000, n_tree = 2500, balanced = TRUE, force = FALSE, leafsize = 1, maxnodes = NA, nodesize = NA){
  i <- 1 ### arbitrarily choose a validation fold
  folds <- create_folds(n_folds = 2, labs = labs, balanced = balanced, force = force) ### create folds such that the groups are balanced
  indx_train <- folds != i
  n_train <- sum(folds != i) ### compute the number of observations in each folds
  if (all(is.na(maxdepth))) maxdepth <- c(seq(1, 7), 10)
  y_train <- matrix(as.numeric(labs[indx_train, ]))
  p_vals <- numeric(length(maxdepth)) ### initialize vector for p-values
  
  y_perm_train <- sapply(1:(n_perms - 1), function(x, y) y[sample(nrow(y)), ,drop = FALSE], y = y_train)
  for (i in 1:length(maxdepth)){
    p_vals[i] <- randomRF(regs = regs[indx_train, ], y = y_train, y_perm = y_perm_train, n_tree = 125, leafsize = leafsize, 
                          maxdepth = maxdepth[i], maxnodes = maxnodes, nodesize = nodesize)
    
  }
  
  optimal_depth <- maxdepth[which.min(p_vals)] ### select k which corresponds to smallest p-value
  y_test <- matrix(as.numeric(labs[!indx_train, ]))
  y_perm_test <- sapply(1:(n_perms - 1), function(x, y) y[sample(nrow(y)), ,drop = FALSE], y = y_test)
  res <- randomRF(regs = regs[!indx_train, ], y = y_test, y_perm = y_perm_test, n_tree = n_tree, leafsize = leafsize, 
                  maxdepth = optimal_depth, maxnodes = maxnodes, nodesize = nodesize)
  return(res)
}



overall.Sgt <- function(regs, labs, n_tree = 2500, maxdepth = NA, n_perms = 1000, leafsize = 1, maxnodes = NA, nodesize = NA){
  pckgs <- as.vector(lsf.str(.GlobalEnv))
  n <- nrow(regs)
  y <- matrix(as.numeric(labs[["y"]]))
  y_perm <- sapply(1:(n_perms - 1), function(x, y) y[sample(nrow(y)), ,drop = FALSE], y = y)
  X <- gentrees(regs = regs, n_tree = 1, leafsize = leafsize, maxdepth = ceiling(runif(1, min = 0, max = 5)), maxnodes = maxnodes, nodesize = nodesize)
  for(i in 1:(n_tree-1)){
    X <- cbind(X, gentrees(regs = regs, n_tree = 1, leafsize = leafsize, maxdepth = ceiling(runif(1, min = 0, max = 5)), maxnodes = maxnodes, nodesize = nodesize))
  }
  S_t <- S(X, y)
  S_p <- c(Inf, S(X, y_perm))
  out <- mean(S_t*(1-1e-14) <= S_p)
  return(out)
}
