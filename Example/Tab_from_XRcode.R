################################################################################
## Table 1: Model Parameter Estimates of Two-Level Intercept-Only Model for Example Data Set}

# load packages
library(dplyr)
library(huxtable)

# indexes for WF
idxW_var <- c() # all equals (n) consecutively
for (j in 1:p){
  idxW_var <- append(idxW_var, (j-1)*n+1)
}
endVarW <- p*n
idxW_cov <- (endVarW+1):(endVarW+c) # all uniques, then equal (n) again
idxW_all <- c(idxW_var, idxW_cov) 
startB_var <- tail(idxW_all,1)+1+c*(n-1)
idxB_var <- startB_var:(startB_var+(p-1))
startB_cov <- tail(idxB_var,1)+1
idxB_cov <- startB_cov:(startB_cov+(c-1))
idxB_all <- c(idxB_var, idxB_cov)

# order: var_x1, var_x2, cov_x1x2

# WF
theta_WF_W <- fit_WF@Fit@x[idxW_all]
theta_WF_B <- fit_WF@Fit@x[idxB_all]

# WFcovshrink_E
theta_WFcovshrink_E_W <- fit_WFcovshrink_E@Fit@x[idxW_all]  
theta_WFcovshrink_E_B <- fit_WFcovshrink_E@Fit@x[idxB_all] 

table <- rbind(
  c(theta_WF_B, theta_WF_W),
  c(theta_WFcovshrink_E_B, theta_WFcovshrink_E_W))
table <- as.data.frame(table)
table <- cbind(c(#"LF", 
  "WF", "WFcovshrink_E"), table)
colnames(table) <- c("",
  "var_x1_B", "var_x2_B", "cov_x1x2_B",
  "var_x1_W", "var_x2_W", "cov_x1x2_W")

# set nb of decimals in table
ndec <- 2 # number of decimals

table <- as_hux(table)
number_format(table) <- ndec
table <- insert_row(table, "Approaches", "Between", "", "", "Within", "", "", after = 0) 
table <- merge_cells(table, 1, 2:4)
table <- merge_cells(table, 1, 5:7)
table <- set_align(table, 1, 1:ncol(table), "center") # all headers centered
table <- set_header_rows(table, 1:2, TRUE)
table <- map_align(table, by_cols(".")) # align at decimal
table <- set_align(table, 1:nrow(table), 1, "left") # first col left-aligned
table <- theme_article(table)
table <- style_headers(table, bold=FALSE)
print_screen(table)

# quick_latex(table, file="MS_tab1.tex")
