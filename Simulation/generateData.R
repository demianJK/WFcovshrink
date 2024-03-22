### this function simulates data in LF separately for the between and within level by specifying the population covariance matrices in lavaan syntax ###

generateData <- function(N, n, g, p, var_B, var_W, cor_B, cor_W) { 
  
  ## variances
  
  var_w_All <- c()
  var_b_All <- c()
  for (i in 1:p) {
    var_w_All[i] <- paste0("x", i, "~~", var_W, "*", "x", i)
    var_b_All[i] <- paste0("x", i, "~~", var_B, "*", "x", i)
  }
  var_w_All <- paste(var_w_All, collapse = "; ")
  var_b_All <- paste(var_b_All, collapse = "; ")
  
  ## covariances 
  
  # cov depends on corr (fix) and var (varying)
  # cov_xy  = corr_xy * var_x/y (because sd/var same for every variable, e.g., x and y)
  
  cov_w <- cor_W * var_W 
  cov_b <- cor_B * var_B
  covs_w <- c()
  covs_b <- c()
  count = 0
  for (i in 1:p) {
    for (j in 1:p) {
      if (i != j & j > i) {
        count <- count + 1
        covs_w[count] <- paste("x", i, "~~", cov_w, "*", "x", j, sep = "")
        covs_b[count] <- paste("x", i, "~~", cov_b, "*", "x", j, sep = "")
      }
    }
  }
  covs_w <- paste(covs_w, collapse = "; ")
  covs_b <- paste(covs_b, collapse = "; ")
  
  # put variances and covariances together (means are 0 per default)
  popModel_W <- paste(var_w_All, covs_w, sep = ";")
  popModel_B <- paste(var_b_All, covs_b, sep = ";")
  
  # generate data
  sample_B <- simulateData(popModel_B, sample.nobs = g, model.type = "lavaan")
  sample_W <- simulateData(popModel_W, sample.nobs = N, model.type = "lavaan")
  groups <- rep(1:g, each = n)
  data <- sample_W
  data[,] <- 0
  for (j in unique(groups)) {
    for (i in min(which(groups == j)):max(which(groups == j)))
      data[i,] <- sample_W[i,] + sample_B[j,]
  }
  data$groups <- as.factor(groups)
  data$persons <- rep(1:n, g)
  data <- cbind(data[,(p+1):(p+2)], data[,1:p])

  return(data)
}