### computes power (alpha < 0.05) over columns of a matrix
EmpericalPower <- function(data){
  out <- apply(data, 2, function(x) mean(x < 0.05))
  return(out)
}

### required for the cv.ridge
rename_folds <- function(folds){
  digits <- sort(unique(folds)) ### sort all unique fold ids
  for (i in 1:length(digits)){ ### for each element in unique sorted folds
    folds[folds == digits[i]] <- i ### rename to 1, .., k (k = length(unique sorted folds))
  }
  return(folds) 
}

### function to check if some number is an integer
is.int <- function(x) x%%1 == 0

### function to create folds implicitly requires the group size per class has to be equal
### if balanced is true (works only for two classes), the labels have to be sorted (low to high) 
create_folds <- function(n_folds, labs, balanced = TRUE, force = FALSE, ...){
  if(balanced){
    if (any(order(as.numeric(labs[["y"]])) != seq(1, nrow(labs)))) stop("Data (labels) is not sorted") ### checks if the above requirments are met
    n1 <- nrow(labs) - sum(as.numeric(labs[["y"]])) ### number of observations in the group 2
    n2 <- sum(as.numeric(labs[["y"]])) ### number of observations group 2
    reps1 <- n1/n_folds ### number of observations per fold for group 1
    reps2 <- n2/n_folds ### number of observations per fold for group 2
    if ((!is.int(reps1) | !is.int(reps2)) & !force) stop("Cannot balance folds") ### checks whether it is possible to balance folds
    folds1 <- rep(seq(1, n_folds), ceiling(reps1))[1:n1] ### create n fold indexes for group 1
    folds2 <- rep(seq(1, n_folds), ceiling(reps2))[1:n2] ### create n fold indexes for group 2
    folds <- c(sample(folds1), sample(folds2)) ### randomly assign indexes to each observation in a balanced fashion
  } else{
    folds <- sample(rep(seq(1, n_folds), nrow(labs)/n_folds)) ### randomly assign folds if balanced is false
  }
  return(folds)
}

### function to generate power plot as function of one parameter
exprplot <- function(data, method = NA, x_names = NA, xlab = NA){
  mthd_spcfd <- !all(is.na(method))  ### check whether names of the methods are specified
  if(all(is.na(x_names))) x_names <- 1:nrow(data) ### check whether (x-axis) paramaters are specified otherwise use 1:nrow(data)
  
  data_wide <- data.frame(data, x_names) ### coerce data into dataframe
  colnames(data_wide)[ncol(data_wide)] <- xlab ### name of the x-axis
  data_long <- pivot_longer(data_wide, !last_col(), names_to = "method", values_to = "RF") ### transform data to long format
  out <- ggplot(data = data_long, aes(x = .data[[xlab]], y = RF, colour = method)) + geom_line(linetype = "dashed") +
    geom_point() + ylim(c(0, 1)) + theme_minimal() + ylab("Power") + theme(plot.title = element_text(face = "bold")) + 
    labs(color = "Test") ### construct power plot
  if (mthd_spcfd) out <- out + scale_colour_discrete(labels = method)
  return(out) ### return plot
}
