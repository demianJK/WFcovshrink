# multiple scripts in individual sessions are run for using multiple cores
# copy "run_1.R" times the number of cores you want to use and change the number in the name and in the respective script (it's used for the random seed below)
script <- 1

# set paths
pathS <- "scripts/functions/" # path for local functions 
pathO <- "objects/" # path for aggregated infos
pathD <- "data/" # path for aggregated info from single runs

# source local functions
source(paste0(pathS, "generateData.R"))
source(paste0(pathS, "analyseData.R"))
source(paste0(pathS, "projectInfo.R"))
pI <- projectInfo()

# (install and) load packages
packages <- c("lavaan", # for SEM
              "ShrinkCovMat", # for shrinkage estimation
              "tidyr", # for reformating into WF
              "DescTools" # for %()%
)
newPackages <-
  packages[!(packages %in% installed.packages()[, "Package"])]
if (length(newPackages)){install.packages(newPackages)}
lapply(packages, require, character.only = TRUE)

# get simulation design
pI <- projectInfo()
SimConds <- readRDS(paste0(pathO, "SimConds.rds"))
ncond <- nrow(SimConds)

# use multiple cores through initializing multiple R sessions
ncores <- 20
nrep <- pI$nrep
scriptReps <- ((script-1)*nrep/ncores+1):(script*nrep/ncores)
set.seed(as.numeric(paste0(script, nrep, ncores, collapse="")))

# Do you want to center the data?
center <- FALSE


approaches <- pI$approaches
napproach <- length(approaches)


for (rep in scriptReps){
  for (cond in 1:ncond){
    
    ListStats <- list()
    
    
    ### 0) preparation ###
    
    # sample characteristics
    n <- SimConds$n[cond]
    p <- SimConds$p[cond]
    g <- SimConds$g[cond]
    
    # population characteristics
    popVar_B <- SimConds$ICC[cond]
    popVar_W <- (1 - SimConds$ICC[cond])
    popCov_B <- SimConds$cor_B[cond] * popVar_B
    popCov_W <- SimConds$cor_W[cond] * popVar_W
    
    # for indexing parameters
    pc <- p * (p + 1) / 2 # number of unique covariance matrix elements (in LF)
    c <- pc-p # numb of unique covariances
    popParams_W <- c(rep(popVar_W, p), rep(popCov_W, c))
    popParams_B <- c(rep(popVar_B, p), rep(popCov_B, c))
    
    
    ### 1) generate data ####

    data <- generateData(N=SimConds$N[cond], n=n, g=g, p=p, var_B=SimConds$ICC[cond], var_W=(1-SimConds$ICC[cond]), cor_B=SimConds$cor_B[cond], cor_W=SimConds$cor_W[cond])

    # prepare data
    data_WF <- pivot_wider(data, names_from = "persons", values_from = 3:(p+2))
    if (center == TRUE){  # grand-mean centering
      data_LF_centered <- as.data.frame(scale(data[,3:(p+2)], scale = F))
      data_LF <- cbind(data[,1:2], data_LF_centered)
      data_WF_centered <- as.data.frame(scale(data_WF[,2:(p*n+1)], scale = F))
      data_WF <- cbind(data_WF[,1], data_WF_centered)
    } else {
      data_LF <- data
    }
    
    
    ### 2) analyse data ####
    
    fitList <- analyseData(data_LF=data_LF,
                           data_WF=data_WF,
                           p=p, n=n, g=g, center=center)
    
    
    
    ### 3) extract relevant information ###
    
    for (i in 1:napproach){ 
      fit <- fitList[[(paste0("fit_", approaches[i]))]]

      # a) if not converged
      ListStats[[approaches[i]]][["conv"]] <- FALSE
      ListStats[[approaches[i]]][["time"]]<- NA
      ListStats[[approaches[i]]][["W"]]<- NA
      ListStats[[approaches[i]]][["B"]] <- NA
      ListStats[[approaches[i]]][["coverage_W"]] <- NA
      ListStats[[approaches[i]]][["coverage_B"]] <- NA
      ListStats[[approaches[i]]][["errorMessages"]] <- NA
      
      # b) if converged
      if(!is.atomic(fit) && fit@Fit@converged){ # if exists and optimizer says converged
        ListStats[[approaches[i]]][["conv"]] <- fit@Fit@converged
        ListStats[[approaches[i]]][["time"]] <- unname(fit@timing[["total"]]) # in s
        if ("WFcovreg_I" %in% approaches[i]){ 
          ListStats[[approaches[i]]][["time"]] <- round(unname(fit@timing[["total"]]) + fitList$WFcovreg_I$covreg_time, 2) # in s
        } else if ("WFcovreg_E" %in% approaches[i]){ 
          ListStats[[approaches[i]]][["time"]] <- round(unname(fit@timing[["total"]]) + fitList$WFcovreg_E$covreg_time, 2) # in s
        } else if ("WFcovreg_U" %in% approaches[i]){ 
          ListStats[[approaches[i]]][["time"]]<- round(unname(fit@timing[["total"]]) + fitList$WFcovreg_U$covreg_time, 2) # in s
        }
        if ("LF" %in% approaches[i]){
          
          ## parameter estimates
          ListStats[[approaches[i]]][["W"]] <- fit@Fit@x[1:pc]
          ListStats[[approaches[i]]][["B"]] <- fit@Fit@x[(pc + 1):(2 * pc)] # the following elements would be the between level intercepts
          
        } else if (grepl("WF", approaches[i])){
          
          ## parameter estimates
          # within-group: n*p times
          idxW_var <- c() # all equals (n) consecutively
          for (j in 1:p){
            idxW_var <- append(idxW_var, (j-1)*n+1)
          }
          endVarW <- p*n
          idxW_cov <- (endVarW+1):(endVarW+c) # all uniques, then equal (n) again
          idxW_all <- c(idxW_var, idxW_cov) 
          # between-group: p times
          startB_var <- tail(idxW_all,1)+1+c*(n-1)
          idxB_var <- startB_var:(startB_var+(p-1))
          startB_cov <- tail(idxB_var,1)+1
          idxB_cov <- startB_cov:(startB_cov+(c-1))
          idxB_all <- c(idxB_var, idxB_cov)
          # the following elements would be the between level (factor) intercepts (means)
          ListStats[[approaches[i]]][["W"]] <- fit@Fit@x[idxW_all]  
          ListStats[[approaches[i]]][["B"]] <- fit@Fit@x[idxB_all] 

        }
      } else if (!is.null(attr(fit, "condition"))) { # save error note 
        ListStats[[approaches[i]]][["errorMessages"]] <- attr(fit, "condition")[[1]] # only error, not function call
      }
    }
    # save information for every repetition of every condition in single RDS
    saveRDS(ListStats, file = paste0(pathD, "cond", cond, "_rep", rep, ".rds"))
    
  }
  rm(list = c("n", "g", "p", "popVar_B", "popCov_B", "popVar_W", "popCov_W", "pc", "c", "popParams_B", "popParams_W", "data", "data_LF", "data_WF", "fit"))
  gc()
  print(rep)
}
