### Example Code 
# - for data, sample covariance matrix, and two-level random-intercept model in Fig. 1
# - for shrinkage estimate of the sample covariance matrix, and two-level random-intercept model in Fig.2 

# Note that the code only works for p=2 and n=2.
# If you want to examine other settings, check out the
# the code for the simulation study on Github:
# - Gitlink-


## (0) preparation ########################################################

# load required packages
library(lavaan) # for model estimation
library(tidyr) # for reformating LF to WF with pivot_wider()
library(ShrinkCovMat) # for shrinkage estimation 
# (code runs in ShrinkCovMat 1.4.0 which is the latest on CRAN, Feb 21st 2024)

# set random number seed to obtain example data
set.seed(4395) 


## (1) Population Characteristics #########################################

# We use the lavaan syntax to set the population models.
popModel_B <- "x1~~0.05*x1; x2~~0.05*x2; x1~~0.015*x2" # between level
popModel_W <- "x1~~0.95*x1; x2~~0.95*x2; x1~~0.285*x2" # within level
# means are zero by default

# We have two variables x1 and x2.
p <- 2
# The variances for both variables are the same at each level.
# The variance at the between level is 0.05.
# The variance at the within level is 0.95.
# Thus, the ICC=0.25.
# The correlation of the two variables is the same at both levels (.3).
# The covariances differ.
# Transform the correlation formula to get the covariances.
# corr_x1x2 = cov_x1x2 / (sd_x1 * var_x2)    
# | * (sd_x1 * sd_x2) and sd_x1 = sd_x2 thus | * var_x1
# corr_x1x2 * var_x1 = cov_x1x2

## (2) Sample Characteristics #############################################

g <- 4 # number of groups (you can change this for exploration)
n <- 2 # group size (balanced data)
N <- g * n # total sample size

# the data sampling is done in long format (LF)
sample_B <- simulateData(popModel_B, sample.nobs = g, 
                         model.type = "lavaan") # between level 
sample_W <- simulateData(popModel_W, sample.nobs = N, # within level
                         model.type = "lavaan")
groups <- rep(1:g, each = n) # group numbers ("j" in Fig. 1)
LF_T <- sample_W # create data frame with the same dimensions
LF_T[,] <- 0 # .. and clear all entries
for (j in unique(groups)) { # merge the sampled data from both levels
  for (i in min(which(groups == j)):max(which(groups == j)))
    LF_T[i,] <- sample_W[i,] + sample_B[j,]
}
LF_T$persons <- rep(1:n, g) # unit numbers ("i" in Fig. 1)
LF_T$groups <- as.factor(groups)
LF_T <- cbind(LF_T[,(p+1):(p+2)], LF_T[,1:p]) # rearrange columns
# LF-T is the total data matrix in long format (LF) ..
# .. which we'll need to reformat to wide format (WF) ..
# .. but first, note that we rounded the data for Fig. 1 and 2
LF_T[,3:(p+2)] <- round(LF_T[,3:(3+p-1)], 0)
WF_T <- pivot_wider(LF_T, names_from = "persons", values_from = 3:4)

# shrinkage estimate S^*_E with equal target Matrix vI_p
# note that unbiased S is employed

WF_T_trans <- t(WF_T[,-1]) # transpose because ShrinkCovMat(data, ..) expects 
# that rows correspond to variables and columns to observations

# estimate S*_E (note that the approach uses unbiased S_WF-T)
WFcovshrink_E <- ShrinkCovMat::shrinkcovmat.equal(data=WF_T_trans, 
                                               centered=FALSE)
round(WFcovshrink_E$Sigmasample, 2) # unbiased S, round for Fig. 2
round(WFcovshrink_E$Target, 2) # vI_p, round for Fig. 2
round(WFcovshrink_E$Sigmahat, 2) # S^*_E, round for Fig. 2
round(WFcovshrink_E$lambdahat, 2) # lambda_E, round for Fig. 2

colnames(WFcovshrink_E$Sigmahat) <- colnames(S_WF_T)
rownames(WFcovshrink_E$Sigmahat) <- rownames(S_WF_T)




## (3) estimate models ####################################################

model_WF <- # Level: 1 (unique factors)
            paste0(
            "x1_1~~Vx1_w*x1_1; x1_2~~Vx1_w*x1_2; x2_1~~Vx2_w*x2_1;   
             x2_2~~Vx2_w*x2_2; x1_1~~Cx12_w*x2_1; x1_2~~Cx12_w*x2_2;",
             # these are the desired within variances and covariances
             # Vx1_w, Vx2_w, and Vx12_w are equality constraints 
            # Level: 2 (common factors)
             "x1_1~0*1; x1_2~0*1; ; x2_1~0*1; x2_2~0*1;",
             # if level-2 variables are aggregates of level-1 variables, 
             # intercepts at level-1 have to be fixed to 0
             "fx1=~1*x1_1+1*x1_2; fx2=~1*x2_1+1*x2_2;", 
             # measurement model with factor loadings set to 1
             "fx1~~fx1; fx2~~fx2; fx1~~fx2;", 
             # these are the desired between variances and covariances
             "fx1~1; fx2~1")  # between means
fit_WF <- sem(model = model_WF, 
              data = WF_T)
summary(fit_WF)

# check sample covariance matrix S_WF-T (biased, normal theory based S)
round(lavInspect(fit_WF, "sampstat")$cov, 2)

fit_WFcovshrink_E <- sem(model_WF,
                                  sample.cov = WFcovshrink_E$Sigmahat, 
                                  sample.cov.rescale = FALSE, 
                                  # rescale sample.cov with (g-1/g)?
                                  sample.nobs = g,
                                  sample.mean = colMeans(WF_T[,-1]))
summary(fit_WFcovshrink_E)

# check shrinkage estimate S*_E (based on unbiased S)
round(S_E,lavInspect(fit_WFcovshrink_E, "sampstat")$cov, 2)

## (4) estimate ICCs ######################################################
# ICC in population: 0.05 (see "Population Characteristics")
# the ICCs are estimated by the parameters of the model-implied matrices

## WF
# both ICC=0, because the between variances of x1 and x2 are negative
fit_WF@Fit@x[7]/(fit_WF@Fit@x[7]+fit_WF@Fit@x[1])
fit_WF@Fit@x[8]/(fit_WF@Fit@x[8]+fit_WF@Fit@x[3])

## WFcovshrink_E
# both ICC=0, because the between variances of x1 and x2 are negative
fit_WFcovshrink_E@Fit@x[7]/(fit_WFcovshrink_E@Fit@x[7]+fit_WFcovshrink_E@Fit@x[1])
fit_WFcovshrink_E@Fit@x[8]/(fit_WFcovshrink_E@Fit@x[8]+fit_WFcovshrink_E@Fit@x[3])

