### summarize results of Monte Carlo study

# (install and) load packages
packages <- c("dplyr",
              "xlsx",
              "DescTools"
              )
newPackages <- packages[!(packages %in% installed.packages()[, "Package"])]
if (length(newPackages)){install.packages(newPackages)}
lapply(packages, require, character.only = TRUE)

# set paths and load data
pathO <- "objects/" # path for aggregated infos
pathD <- "data/" # path for aggregated info from single repetitions
source("scripts/functions/projectInfo.R") # infos about simulation study
pI <- projectInfo()

SimConds <- readRDS(paste0(pathO, "SimConds.rds"))
ncond <- nrow(SimConds)
nrep <- pI$nrep
approaches <- pI$approaches
napproach <- length(approaches)

## Data -> ListStats ######################################################################################################################################################
# put all relevant info from single files into one list

critsList <- c("conv", # Did the model converge?
               "time", # How long did it take till the model converged?
               "W", # within-group model parameters
               "B", # between-group model parameters
               "errorMessages"
) 

# create list to save results of all replications in all conditions
ListStats <- setNames(as.list(approaches), approaches)
ListStats <- lapply(ListStats, function(x){
  setNames(vector("list", length(critsList)), critsList)
})
for (i in 1:napproach){ # initializing necessary
  for (j in 1:length(critsList)){
    ListStats[[i]][[j]] <- setNames(vector("list", ncond), 1:ncond)
    for (k in 1:ncond){
      ListStats[[i]][[j]][[k]] <- setNames(vector("list", nrep), 1:nrep)
    }
  }
}

for (rep in 1:nrep){    
  for (cond in 1:ncond){
    
    SingleList <- readRDS(paste0(pathD, "cond", cond, "_rep", rep, ".rds"))
    
    for (i in 1:napproach){
      
      # model properties
      ListStats[[approaches[i]]][["conv"]][[cond]][[rep]] <- SingleList[[approaches[i]]][["conv"]]
      ListStats[[approaches[i]]][["time"]][[cond]][[rep]] <- SingleList[[approaches[i]]][["time"]]
      ListStats[[approaches[i]]][["W"]][[cond]][[rep]] <- SingleList[[approaches[i]]][["W"]]
      ListStats[[approaches[i]]][["B"]][[cond]][[rep]] <- SingleList[[approaches[i]]][["B"]]
      ListStats[[approaches[i]]][["coverage_W"]][[cond]][[rep]] <- SingleList[[approaches[i]]][["coverage_W"]]
      ListStats[[approaches[i]]][["coverage_B"]][[cond]][[rep]] <- SingleList[[approaches[i]]][["coverage_B"]]
      ListStats[[approaches[i]]][["errorMessages"]][[cond]][[rep]] <- SingleList[[approaches[i]]][["errorMessages"]]
      if ("LF" %in% approaches[i]){
        ListStats[[approaches[i]]][["ICC_S"]][[cond]][[rep]] <- SingleList[[approaches[i]]][["ICC_S"]]
        ListStats[[approaches[i]]][["ICC_Sigma"]][[cond]][[rep]] <- SingleList[[approaches[i]]][["ICC_Sigma"]]
      } 
    }
    gc()
  }
  print(rep)
}

saveRDS(ListStats, file = paste0(pathO, "ListStats.rds"))



### ListStats --> TableStats  ######################################################################################################################################################
# summarize info from list into table

estim <- c("RMSE", "Bias")
critsTable <- c("conv", "time", "timesd",
                
                ## within-group
                # variance
                paste(estim, "var", "W", sep="_"),
                paste(paste0("rel", estim), "var", "W", sep="_"),
                # covariance
                paste(estim, "cov", "W", sep="_"),
                paste(paste0("rel", estim), "cov", "W", sep="_"),
                # overall
                paste(estim, "W", sep="_"),
                paste(paste0("rel", estim), "W", sep="_"),
                # negative variances
                "negVar_W",
                
                ## between-group
                # variance
                paste(estim, "var", "B", sep="_"),
                paste(paste0("rel", estim), "var", "B", sep="_"),
                # covariance
                paste(estim, "cov", "B", sep="_"),
                paste(paste0("rel", estim), "cov", "B", sep="_"),
                # overall
                paste(estim, "B", sep="_"),
                paste(paste0("rel", estim), "B", sep="_"),
                
                # negative variances
                "negVar_B",
                
                # Intraclass Correlation
                "ICC_Sigma_theta",
                "ICC_Sigma_theta_neg",
                "RMSE_ICC",
                "relRMSE_ICC",
                "Bias_ICC",
                "relBias_ICC"
)

colTable <- c()
for (i in 1:napproach){
  colTable <- append(colTable, paste(critsTable, approaches[i], sep="_"))
}
tmp <- as.data.frame(matrix(NA, ncol=length(colTable), nrow=ncond))
colnames(tmp) <- colTable
TableStats <- SimConds
TableStats <- cbind(TableStats, tmp)

for (cond in 1:ncond){
  
  # required for stat props
  p <- SimConds$p[cond]
  n <- SimConds$n[cond]
  popVar_B <- SimConds$ICC[cond]
  popVar_W <- (1 - SimConds$ICC[cond])
  popCov_B <- SimConds$cor_B[cond] * popVar_B
  popCov_W <- SimConds$cor_W[cond] * popVar_W
  
  # for indexing parameters
  pc <- p * (p + 1) / 2
  c <- pc - p
  # sorting of params in LF and WF: varW (p), covW (c), varB (p), covB (c)
  
  for (i in 1:napproach){
    
    # convergence
    tmpList <- unname(unlist(ListStats[[approaches[i]]][["conv"]][[cond]]))
    tmpList <- factor(tmpList, levels=c("TRUE", "FALSE"))
    tmp <- table(tmpList)
    TableStats[cond, paste("conv", approaches[i], sep="_")] <- tmp[1] / ( tmp[1] + tmp[2]) * 100
    
    if (TableStats[cond, paste("conv", approaches[i], sep="_")] > 0){
      
      # computation time
      tmp <- unname(unlist(ListStats[[approaches[i]]][["time"]][[cond]])) 
      TableStats[cond, paste("time", approaches[i], sep="_")] <- round(mean(tmp, na.rm=TRUE), 2) 
      TableStats[cond, paste("timesd", approaches[i], sep="_")] <- round(sd(tmp, na.rm=TRUE), 2) 
      
      ## Parameter Estimates
      
      ## within-group
      W <- na.omit(as.data.frame(do.call(rbind, ListStats[[approaches[i]]][["W"]][[cond]])))
      
      ## variance
      var_W <- W[,1:p]
      
      # negative variances 
      tv <- c()
      for (rep in 1:nrep){
        tv[rep] <- any(var_W[rep,] < 0) 
      }
      tv <- factor(tv, levels=c("TRUE", "FALSE"))
      tmp <- table(tv)
      TableStats[cond, paste("negVar_W", approaches[i], sep="_")] <- tmp[1] / ( tmp[1] + tmp[2]) * 100
      
      RMSE_var_W <- mean(sqrt(colMeans((var_W - popVar_W)^2)))
      TableStats[cond, paste("relRMSE_var_W", approaches[i], sep="_")] <- round( ((RMSE_var_W / popVar_W)* 100), 2) 
      
      Bias_var_W <- mean(colMeans((var_W - popVar_W)))
      TableStats[cond, paste("relBias_var_W", approaches[i], sep="_")] <- round( ((Bias_var_W / popVar_W)* 100), 2)
      
      ## covariance
      cov_W <- as.data.frame(W[,(p+1):(p+c)])
      
      RMSE_cov_W <- mean(sqrt(colMeans((cov_W - popCov_W)^2)))
      TableStats[cond, paste("relRMSE_cov_W", approaches[i], sep="_")] <- round( ((RMSE_cov_W / popCov_W)* 100), 2)
      
      Bias_cov_W <- mean(colMeans((cov_W - popCov_W)))
      TableStats[cond, paste("relBias_cov_W", approaches[i], sep="_")] <- round( ((Bias_cov_W / popCov_W)* 100), 2)
      
      ## overall (variance and covariance)

      TableStats[cond, paste("relRMSE_W", approaches[i], sep="_")] <- round( mean( c(((RMSE_var_W / popVar_W)* 100), ((RMSE_cov_W / popCov_W)* 100)) ), 2)
      TableStats[cond, paste("relBias_W", approaches[i], sep="_")] <- round( mean( c(((Bias_var_W / popVar_W)* 100), ((Bias_cov_W / popCov_W)* 100)) ), 2)
      
      
      ### between-group
      B <- na.omit(as.data.frame(do.call(rbind, ListStats[[approaches[i]]][["B"]][[cond]])))
      
      # variances
      var_B <- B[,1:p]
      
      # negative variances 
      tv <- c()
      for (rep in 1:nrep){
        tv[rep] <- any(var_B[rep,] < 0) 
      }
      tv <- factor(tv, levels=c("TRUE", "FALSE"))
      tmp <- table(tv)
      TableStats[cond, paste("negVar_B", approaches[i], sep="_")] <- tmp[1] / ( tmp[1] + tmp[2]) * 100
      
      RMSE_var_B <- mean(sqrt(colMeans((var_B - popVar_B)^2)))
      TableStats[cond, paste("relRMSE_var_B", approaches[i], sep="_")] <- round( ((RMSE_var_B / popVar_B)* 100), 2) 
      
      Bias_var_B <- mean(colMeans((var_B - popVar_B)))
      TableStats[cond, paste("relBias_var_B", approaches[i], sep="_")] <- round( ((Bias_var_B / popVar_B)* 100), 2)

      ## covariance
      cov_B <- as.data.frame(B[,(p+1):(p+c)])
      
      RMSE_cov_B <- mean(sqrt(colMeans((cov_B - popCov_B)^2)))
      TableStats[cond, paste("relRMSE_cov_B", approaches[i], sep="_")] <- round( ((RMSE_cov_B / popCov_B)* 100), 2)
      
      Bias_cov_B <- mean(colMeans((cov_B - popCov_B)))
      TableStats[cond, paste("relBias_cov_B", approaches[i], sep="_")] <- round( ((Bias_cov_B / popCov_B)* 100), 2)

      ## overall (variance and covariance)
      
      TableStats[cond, paste("relRMSE_B", approaches[i], sep="_")] <- round( mean( c(((RMSE_var_B / popVar_B)* 100), ((RMSE_cov_B / popCov_B)* 100)) ), 2)
      TableStats[cond, paste("relBias_B", approaches[i], sep="_")] <- round( mean( c(((Bias_var_B / popVar_B)* 100), ((Bias_cov_B / popCov_B)* 100)) ), 2)
      
      # ICC_Sigma_theta (based on model parameters of intercept-only model)
      
      ICC_Sigma_theta <- unlist((var_B/(var_W+var_B))) # vector
      tmp <- sum(ICC_Sigma_theta < 0)
      if (tmp > 0){
        TableStats[cond, paste("ICC_Sigma_theta_neg", approaches[i], sep="_")] <- round(tmp/length(ICC_Sigma_theta) * 100, 2)
        ICC_Sigma_theta[ICC_Sigma_theta < 0] <- 0
      }
      TableStats[cond, paste("ICC_Sigma_theta", approaches[i], sep="_")] <- round(mean(ICC_Sigma_theta, na.rm=TRUE), 2)

      ICC <- var_B/(var_W+var_B) # data frame
      
      RMSE_ICC <- mean(sqrt(colMeans((ICC - popVar_B)^2)))
      TableStats[cond, paste("relRMSE_ICC", approaches[i], sep="_")] <- round( ((RMSE_ICC / popVar_B)* 100), 2)
      
      Bias_ICC <- mean(colMeans(ICC - popVar_B))
      TableStats[cond, paste("relBias_ICC", approaches[i], sep="_")] <- round( ((Bias_ICC / popVar_B)* 100), 2)

    }  
  }
}

# aggregate over within- and between-group

for (i in 1:napproach){
  TableStats[paste("relRMSE", approaches[i], sep="_")] <- round(rowMeans(select(TableStats, paste(c("relRMSE_W", "relRMSE_B"), approaches[i], sep="_"))), 2)
  TableStats[paste("relBias", approaches[i], sep="_")] <- round(rowMeans(select(TableStats, paste(c("relBias_W", "relBias_B"), approaches[i], sep="_"))), 2)
}


write.xlsx( TableStats, file = paste0(pathO, "TableStats.xlsx"))



### TableStats --> dat  ######################################################################################################################################################
# rearrange summarized table into data frame that can be used for figures

ratio <- c()
for (i in 1:length(approaches)){
  if (grepl("LF", approaches[i])){
    ratio <- append(ratio, TableStats$ratio_LF)
  } else if (grepl("WF", approaches[i])){
    ratio <- append(ratio, TableStats$ratio_WF)
  }
}

dat <- data.frame(cond = rep(1:ncond, napproach),
                  n = rep(TableStats$n, napproach),
                  g = rep(TableStats$g, napproach),
                  p = rep(TableStats$p, napproach),
                  ratio  = ratio,
                  ICC = rep(TableStats$ICC, napproach),
                  cor_W = rep(TableStats$cor_W, napproach),
                  cor_B = rep(TableStats$cor_B, napproach),
                  approach = factor(rep(approaches, each=ncond), 
                                    level=c("LF", "WF",
                                            "WF_data_unbiased", "WF_S_biased", "WF_S_unbiased", 
                                            "WFcovreg_E", "WFcovreg_I", "WFcovreg_U"), ordered=TRUE),
                  conv = unlist(TableStats[, grepl("conv", colnames(TableStats))], use.names=FALSE),
                  time = c(unlist(TableStats[, grepl("time_", colnames(TableStats))], use.names=FALSE)),
                  ## type
                  # variance
                  relRMSE_var_B = unlist(TableStats[, paste0("relRMSE_var_B_", approaches)], use.names=FALSE),
                  relRMSE_var_W = unlist(TableStats[, paste0("relRMSE_var_W_", approaches)], use.names=FALSE),
                  relBias_var_B = unlist(TableStats[, paste0("relBias_var_B_", approaches)], use.names=FALSE),
                  relBias_var_W = unlist(TableStats[, paste0("relBias_var_W_", approaches)], use.names=FALSE),
                  negVar_B = unlist(TableStats[, paste0("negVar_B_", approaches)], use.names=FALSE),
                  negVar_W = unlist(TableStats[, paste0("negVar_W_", approaches)], use.names=FALSE),
                  # covariance
                  relRMSE_cov_B = unlist(TableStats[, paste0("relRMSE_cov_B_", approaches)], use.names=FALSE),
                  relRMSE_cov_W = unlist(TableStats[, paste0("relRMSE_cov_W_", approaches)], use.names=FALSE),
                  relBias_cov_B = unlist(TableStats[, paste0("relBias_cov_B_", approaches)], use.names=FALSE),
                  relBias_cov_W = unlist(TableStats[, paste0("relBias_cov_W_", approaches)], use.names=FALSE),
                  # level
                  relRMSE_B = unlist(TableStats[, paste0("relRMSE_B_", approaches)], use.names=FALSE),
                  relRMSE_W = unlist(TableStats[, paste0("relRMSE_W_", approaches)], use.names=FALSE),
                  relBias_B = unlist(TableStats[, paste0("relBias_B_", approaches)], use.names=FALSE),
                  relBias_W = unlist(TableStats[, paste0("relBias_W_", approaches)], use.names=FALSE),
                  # aggregated
                  relRMSE = unlist(TableStats[, paste0("relRMSE_", approaches)], use.names=FALSE),
                  relBias = unlist(TableStats[, paste0("relBias_", approaches)], use.names=FALSE),
                  # ICC
                  ICC_Sigma_theta = unlist(TableStats[, paste0("ICC_Sigma_theta_", approaches)], use.names=FALSE),
                  ICC_Sigma_theta_neg = unlist(TableStats[, paste0("ICC_Sigma_theta_neg_", approaches)], use.names=FALSE),
                  relRMSE_ICC = unlist(TableStats[, paste0("relRMSE_ICC_", approaches)], use.names=FALSE),
                  relBias_ICC = unlist(TableStats[, paste0("relBias_ICC_", approaches)], use.names=FALSE)
)


saveRDS(dat, paste0(pathO, "dat.rds"))
