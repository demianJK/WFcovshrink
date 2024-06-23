### this function analyses data in the LF and WF approaches with an intercept-only model ###

analyseData <- function(data_LF, data_WF, p, n, g, center){
  
  ### LF approach
  
  ## 1) set model specification (same for within and between)
  
  # variances
  Var <- c()
  for (i in 1:p) {
    Var[i] <- paste0("x", i, "~~", "x", i)
  }
  Var <- paste(Var, collapse="; ")
  
  # covariances
  Cov <- c()
  count=0
  for(i in 1:p){   
    for(j in 1:p){
      if(i != j & j > i){
        count = count+1
        Cov[count] <- paste("x", i, "~~", "x", j, sep="") 
      }
    }
  }
  Cov <- paste(Cov, collapse="; ")
  
  model_LF_W <- paste0("Level: 1\n", paste(Var, Cov, sep=";"))
  model_LF_B <- paste0("Level: 2\n", paste(Var, Cov, sep=";"))
  model_LF <- paste(model_LF_W, model_LF_B, sep="\n")
  
  # 2) fit model
  
  fit_LF <- try(sem(model = model_LF, # for multilevel SEM, use sem(), multilevel meanstructure: within ints=0, between ints estimated, see https://www.lavaan.ugent.be/tutorial/multilevel.html
                    data = data_LF,
                    cluster = "groups"), 
                silent = TRUE)
  
  
  ### WF approach
  
  ## 1) set model specification (different for within and between)
  
  ## within (p * n)
  
  # variances (p * n)
  tmp2 <- c()
  resid_w <- c() # residual variances (p * n - equality among n blocks)
  tmp3 <- c()
  for (j in 1:p){
    for (i in 1:n){
      tmp2[i] <- paste0("x", j,  ".", i)
      tmp3[i] <- paste0(tmp2[i], "~~Vx", j, "_w*", tmp2[i])
    }
    resid_w[j] <- paste(tmp3, collapse="; ")
  }
  resid_w <- paste(resid_w, collapse="; ")
  
  # covariances (p-pc=c*n with c-wise equality constraints)
  resid_cov <- c() # manifest correlations (p * n - equality among n blocks)
  count <- 0
  for (i in 1:n){ # n-unit-wise ordering
    for(j in 1:p){   
      for(m in 1:p){
        if(j != m & m > j){
          count <- count + 1
          resid_cov[count] <- paste0("x", j,  ".", i, "~~", "Cx", j, m, "_w*", "x", m,  ".", i)
        }
      } 
    }
  }
  resid_cov <- paste(resid_cov, collapse="; ")
  
  # "If variables are both at the within- and the between-level, the intercepts at the within-level should be fixed to zero." (p.701f, Barendse & Rosseel, 2020)
  fac_int_w <- c()
  tmp <- c()
  count <- 0
  for (j in 1:p){
    for (i in 1:n){
      count <- count + 1
      tmp[count] <- paste0("x", j,  ".", i, "~0*1")
    }
  }
  fac_int_w <- paste(tmp, collapse = "; ")
  
  model_WF_W <- paste(resid_w, resid_cov, fac_int_w, sep = "; ")
  
  ## between (p)
  
  fac_b <- c() # latent factor = random intercepts (p - loadings fixed to 1)
  tmp <- c()
  for (j in 1:p){
    for (i in 1:n){
      tmp[i] <- paste0("1*x", j,  ".", i)
    }
    fac_b[j] <- paste0("fx", j, "=~", paste(tmp, collapse="+"))
  }
  fac_b <- paste(fac_b, collapse="; ")
  
  # variances
  fac_var_b <- c() # factor variance (p)
  fac_int_b <- c() # to estimate means
  for (j in 1:p){
    fac_var_b[j] <- paste0("fx", j, "~~fx", j)
    fac_int_b[j] <- paste0("fx", j, "~1")
  }
  fac_var_b <- paste(fac_var_b, collapse="; ")
  fac_int_b <- paste(fac_int_b, collapse="; ")
  
  # covariances
  fac_cov_b <- c() # correlations between factors
  count <- 0
  for(j in 1:p){   
    for(m in 1:p){
      if(j != m & m > j){
        count <- count + 1
        fac_cov_b[count] <- paste0("fx", j, "~~", "fx", m)
      }
    } 
  }
  fac_cov_b <- paste(fac_cov_b, collapse = "; ")
  
  model_WF_B <- paste(fac_b, fac_var_b, fac_cov_b, fac_int_b, sep="; ")
  
  model_WF <- paste(model_WF_W, model_WF_B, sep="; ")
  
  ## 2) fit model
  
  # Does it make a difference whether we  ...
  # (1) supply data (a, b) or S (c, d)?
  # (2) use biased (a, c) or unbiased S (b, d)?
  
  # (a) with data and biased S (/g; default)
  fit_WF <- try(sem(model = model_WF, # Barendse & Rossel (2020) use sem() as well (see Appendix A, p.718)
                    #likelihood = "normal", # S, SE, and test statistics based on /(g)
                    data = data_WF),
                    silent = TRUE)
  
  # (b) with data and unbiased S (/(g-1))
  fit_WF_data_unbiased <- try(sem(model = model_WF, # Barendse & Rossel (2020) use sem() as well (see Appendix A, p.718)
                              likelihood = "wishart", # S, SE, and test statistics based on /(g-1)
                              data = data_WF),
                          silent = TRUE)

  S_WF_biased <- cov(data_WF[,2:(p*n+1)]) * (g-1)/g # biased, normal theory based S (/g); in contrast to unbiased, wishart distributed S (/(g-1))

  # (c) supplying biased S (/g) instead of data
  fit_WF_S_biased <- try(sem(model = model_WF,
                           sample.cov = S_WF_biased,
                           sample.cov.rescale = FALSE, # rescale sample.cov with (g-1/g)?
                           #likelihood = "normal", # S, SE, and test statistics based on /(g)
                           sample.mean = colMeans(data_WF[,2:(p*n+1)]),
                           sample.nobs = g),
                           silent = TRUE)

  # (d) supplying unbiased S (/(g-1))
  fit_WF_S_unbiased <- try(sem(model = model_WF,
                             sample.cov = S_WF_biased,
                             sample.cov.rescale = TRUE, # rescale sample.cov with (g-1/g)?
                             likelihood = "wishart", # S, SE, and test statistics based on /(g-1)
                             sample.mean = colMeans(data_WF[,2:(p*n+1)]),
                             sample.nobs = g),
                         silent = TRUE)
  
  
  ### WF approach with shrinkage estimaton
  
  # preparation
  data_WF_trans <- t(data_WF[,-1]) # transpose because shrinkcovmat(data, ..) expects rows of the data matrix data correspond to variables and the columns to subjects
  means <- colMeans(data_WF[,-1])
  varNames <- colnames(data_WF[,-1])
  
  ## 1) with identity matrix
  
  WFcovreg_I_start <- Sys.time()
  WFcovreg_I <- ShrinkCovMat::shrinkcovmat.identity(data=data_WF_trans, centered=center)
  WFcovreg_I_end <- Sys.time()
  WFcovreg_I$covreg_time <- as.numeric(difftime(WFcovreg_I_end, WFcovreg_I_start, units = "secs"))
  
  # name regularized matrix
  rownames(WFcovreg_I$Sigmahat) <- varNames
  colnames(WFcovreg_I$Sigmahat) <- varNames
  
  fit_WFcovreg_I <- try(lavaan::sem(model_WF,
                                    sample.cov = WFcovreg_I$Sigmahat, # note that covreg is based on unbiased and not biased S (as lavaan)
                                    sample.cov.rescale = FALSE, # rescale sample.cov with (g-1/g)?
                                    sample.nobs = g,
                                    sample.mean = means),
                        silent=TRUE)
  
  
  ## 2) with equal target matrix
  
  WFcovreg_E_start <- Sys.time()
  WFcovreg_E <- ShrinkCovMat::shrinkcovmat.equal(data=data_WF_trans, centered=center)
  WFcovreg_E_end <- Sys.time()
  WFcovreg_E$covreg_time <- as.numeric(difftime(WFcovreg_E_end, WFcovreg_E_start, units = "secs"))
  
  # name regularized matrix
  rownames(WFcovreg_E$Sigmahat) <- varNames
  colnames(WFcovreg_E$Sigmahat) <- varNames
  
  fit_WFcovreg_E <- try(lavaan::sem(model_WF,
                                     sample.cov = WFcovreg_E$Sigmahat, # note that covreg is based on unbiased and not biased S (as lavaan)
                                     sample.cov.rescale = FALSE, # rescale sample.cov with (g-1/g)?
                                     sample.nobs = g,
                                     sample.mean = means),
                         silent=TRUE)
  
  
  ## 3) with unequal target matrix
  
  WFcovreg_U_start <- Sys.time()
  WFcovreg_U <- ShrinkCovMat::shrinkcovmat.unequal(data=data_WF_trans, centered=center)
  WFcovreg_U_end <- Sys.time()
  WFcovreg_U$covreg_time <- as.numeric(difftime(WFcovreg_U_end, WFcovreg_U_start, units = "secs"))
  
  # name regularized matrix
  rownames(WFcovreg_U$Sigmahat) <- varNames
  colnames(WFcovreg_U$Sigmahat) <- varNames
  
  fit_WFcovreg_U <- try(lavaan::sem(model_WF,
                                    sample.cov = WFcovreg_U$Sigmahat, # note that covreg is based on unbiased and not biased S (as lavaan)
                                    sample.cov.rescale = FALSE, # rescale sample.cov with (g-1/g)?
                                    sample.nobs = g,
                                    sample.mean = means),
                        silent=TRUE)
  
  return(list(fit_LF=fit_LF, 
              fit_WF=fit_WF,
              fit_WF_data_unbiased=fit_WF_data_unbiased,
              fit_WF_S_biased=fit_WF_S_biased,
              fit_WF_S_unbiased=fit_WF_S_unbiased,
              fit_WFcovreg_I=fit_WFcovreg_I,
              fit_WFcovreg_E=fit_WFcovreg_E,
              fit_WFcovreg_U=fit_WFcovreg_U,
              WFcovreg_I=WFcovreg_I,
              WFcovreg_E=WFcovreg_E,
              WFcovreg_U=WFcovreg_U
  ) ) 
}

# Barendse, M. T., & Rosseel, Y. (2020). Multilevel Modeling in the ‘Wide Format’ Approach with Discrete Data: A Solution for Small Cluster Sizes. Structural Equation Modeling: A Multidisciplinary Journal, 27(5), 696–721. https://doi.org/10.1080/10705511.2019.1689366
