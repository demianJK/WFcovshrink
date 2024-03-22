projectInfo <- function(){
  out <- list(

  approaches = c("LF", "WF", # WF is WF_data_biased
                 "WF_data_unbiased", "WF_S_biased", "WF_S_unbiased", # controls
                 "WFcovreg_I", "WFcovreg_E", "WFcovreg_U"), # WF with shrinkage estimation
  nrep = 1000,
  p = c(2, 5, 10),
  ICC = c(0.05, 0.25),
  cor_B = c(0.1, 0.3),
  cor_W = c(0.1, 0.3),
  g = c(4, 10, 30, 50, 100),
  n = c(2, 5, 10)
  )

  if ( !file.exists("objects/SimConds.rds") ){ # create Sim Table if it does not exist
    final <- expand.grid(out[3:length(out)])
    final$N <- final$g * final$n # total sample size
    final$pn <- final$p * final$n # nb of "observed" variables (i.e., columns) in WF approach
    final <- final[order(final$p),] # sort with increasing p
    cond <- 1:nrow(final) # give simulation conditions nb
    table <- cbind(cond, final)
    xlsx::write.xlsx( table, file="objects/SimConds.xlsx" )
    saveRDS(table, file = "objects/SimConds.rds" ) 
  }
  
  return(out)
  
}


