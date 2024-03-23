## Code for Fig. 3 and following

# load packages
library(ggplot2)
library(dplyr) 
library(patchwork) # arrange plots
library(cowplot) # get_legend
library(tidyr)

# load data

pathF <- "/figures"
dat <- readRDS("objects/dat.rds")

# data preparation
dat$n_named <- paste0("n == ", dat$n)
dat$n_named <- factor(dat$n_named, ordered=TRUE, levels=unique(dat$n_named))
dat$p_named <- paste0("p == ", dat$p)
dat$p_named <- factor(dat$p_named, ordered=TRUE, levels=unique(dat$p_named))
dat$ICC_named <- paste0("ICC == ", dat$ICC) 
dat$ICC_named <- factor(dat$ICC_named, ordered=TRUE, levels=unique(dat$ICC_named))
dat$cor_W_named <- factor(dat$cor_W, levels = c("0.1", "0.3"), 
                          labels = c("rho[W]==.10", "rho[W]==.30"))
dat$cor_B_named <- factor(dat$cor_B, levels = c("0.1", "0.3"), 
                          labels = c("rho[B]==.10", "rho[B]==.30"))

# data frame with all approaches
allNames <- c("LF", "WF", "WFcovreg_E", "WFcovreg_I", "WFcovreg_U")
allNames_new <- c("LF", "WF", "WFcovshrink(E)", "WFcovshrink(I)", "WFcovshrink(U)")
datAll <- filter(dat, approach %in% allNames)
allCols <- c("#31a354", "#2c7fb8", "#810f7c", "#8856a7", "#c51b8a")
# select only conditions where all models converged
idxConv <- filter(filter(datAll, approach == "WF"), conv > 0)$cond
datAllconv <- filter(datAll, cond %in% idxConv)


# data frame with all standard WF approaches
WFNames <- c("WF", "WF_data_unbiased", "WF_S_biased", "WF_S_unbiased")
datWF <- filter(dat, approach %in% WFNames)
datWF$Input <- c(rep("data", 2*length(unique(dat$cond))), rep("S", 2*length(unique(dat$cond))))
datWF$Estimator <- c(rep(c("biased", "unbiased"), each=length(unique(dat$cond)), times=2))
WFCols <- c("#2c7fb8", "#b2e2e2") # biased, unbiased


# cols:rows
# LF-B: p / g
# WF-T: (p*n) / g
limCR <- c(0, max(dat$ratio))
breaksCR <- sort(unique(dat$ratio))
labelsCR <- c(rep("", 15), 1, rep("", 5), 5, "", 10, "", 25)

# figure specification
jitterWidth <- 3
pointSize <- 3



### Fig. 3: Convergence Rates by Sample Characteristics

a <- ggplot(datAll, aes(y=conv, x=g, col=approach)) + 
  stat_summary(geom = "line", linetype="dashed") +
  stat_summary(fun="mean") +
  scale_y_continuous(name="Convergence Rate (%)", expand=c(0.05,0.05), limits=c(0,100)) +
  scale_x_continuous(breaks=unique(dat$g)) +
  theme_minimal() + 
  scale_colour_manual(name="Approach", labels=allNames_new, values=allCols) +
  theme(text = element_text(family="serif"), panel.grid.minor = element_blank(), legend.pos = "none",
        panel.border = element_rect(color = "grey",  fill = NA, linewidth = 0.5)) +
  labs(title="A")

b <- ggplot(datAll, aes(y=conv, x=ratio, col=approach)) + 
  stat_summary(geom = "line", linetype="dashed") +
  stat_summary(fun="mean") +
  facet_grid(cols=vars(approach)) +
  scale_y_continuous(name="Convergence Rate (%)", expand=c(0.05,0.05)) +
  scale_x_continuous(breaks=breaksCR, labels = labelsCR, name=expression(cols:rows)) +
  geom_vline(xintercept=1, color="grey", size=0.2) +
  theme_minimal() + 
  scale_colour_manual(name="Approach", labels=allNames_new, values=allCols) +
  theme(text = element_text(family="serif"), panel.grid.minor = element_blank(), 
        strip.text.x = element_blank(), 
        panel.border = element_rect(color = "grey",  fill = NA, linewidth = 0.5)) +
  labs(title="B")

a + b + plot_layout(width = c(0.5, 1), guides='collect') & theme(legend.position = "bottom")

# ggsave(
#   filename = "Fig_3.jpeg",
#   path = pathF,
#   device = "jpeg",
#   plot = last_plot(),
#   width = 900*3, # values used for png multiplied with 3 (bc for png, dpi=97)
#   height = 300*3,
#   units = "px",
#   dpi = 300
# )


### Fig. 4: Relative Bias of Parameter Estimates

a <-
  ggplot(datAllconv, aes(x=g, y=relBias_var_B, colour=approach)) + 
  geom_rect(aes(xmin = Inf, xmax = -Inf, ymin = 0, ymax = Inf), alpha = 0.02, col="#fbc8c7", fill="#fbc8c7") +
  geom_rect(aes(xmin = Inf, xmax = -Inf, ymin = -Inf, ymax = 0), alpha = 0.02, col="#99e5e8", fill="#99e5e8") +
  stat_summary(geom = "line", linetype="dashed", alpha=0.5) +
  stat_summary(fun="mean", alpha=0.5) +
  scale_colour_manual(name="Approach", labels=allNames_new, values=allCols) +
  scale_y_continuous(name="Relative Bias (%)", expand=c(0.05,0.05), breaks=seq(-300, 0, 50)
  ) +
  scale_x_continuous(breaks=unique(datAllconv$g), labels=unique(datAllconv$g), limits=range(datAllconv$g), expand=c(0.05,0.05)) +
  coord_cartesian(ylim = c(-310, 0)) +
  theme_minimal() + 
  theme(text = element_text(family="serif"), panel.grid.minor = element_blank(), 
        axis.text.x = element_text(margin=margin(10,0,0,0)), legend.pos = "none",
        panel.border = element_rect(color = "grey",  fill = NA, size = 0.5)) +
  labs(title="Between-Group Variance")

b <-
  ggplot(datAllconv, aes(x=g, y=relBias_cov_B, colour=approach)) + 
  geom_rect(aes(xmin = Inf, xmax = -Inf, ymin = 0, ymax = Inf), alpha = 0.02, col="#fbc8c7", fill="#fbc8c7") +
  geom_rect(aes(xmin = Inf, xmax = -Inf, ymin = -Inf, ymax = 0), alpha = 0.02, col="#99e5e8", fill="#99e5e8") +
  stat_summary(geom = "line", linetype="dashed", alpha=0.5) +
  stat_summary(fun="mean", alpha=0.5) +
  scale_colour_manual(name="Approach", labels=allNames_new, values=allCols) +
  scale_y_continuous(name="Relative Bias (%)", expand=c(0.05,0.05), breaks=seq(-300, 0, 50)
  ) +
  scale_x_continuous(breaks=unique(datAllconv$g), labels=unique(datAllconv$g), limits=range(datAllconv$g), expand=c(0.05,0.05)) +
  coord_cartesian(ylim = c(-310, 0)) +
  theme_minimal() + 
  theme(text = element_text(family="serif"), panel.grid.minor = element_blank(), 
        axis.text.x = element_text(margin=margin(10,0,0,0)), legend.pos = "none",
        panel.border = element_rect(color = "grey",  fill = NA, size = 0.5)) +
  labs(title="Between-Group Covariance")

c <- ggplot(datAllconv, aes(x=g, y=relBias_var_W, colour=approach)) + 
  geom_rect(aes(xmin = Inf, xmax = -Inf, ymin = 0, ymax = Inf), alpha = 0.02, col="#fbc8c7", fill="#fbc8c7") +
  geom_rect(aes(xmin = Inf, xmax = -Inf, ymin = -Inf, ymax = 0), alpha = 0.02, col="#99e5e8", fill="#99e5e8") +
  stat_summary(geom = "line", linetype="dashed", alpha=0.5) +
  stat_summary(fun="mean", alpha=0.5) +
  scale_colour_manual(name="Approach", labels=allNames_new, values=allCols) +
  scale_y_continuous(name="Relative Bias (%)", expand=c(0.05,0.05), breaks=seq(0, 75, 25)
  ) +
  scale_x_continuous(breaks=unique(datAllconv$g), labels=unique(datAllconv$g), limits=range(datAllconv$g), expand=c(0.05,0.05)) +
  coord_cartesian(ylim = c(0, 75)) +
  theme_minimal() + 
  theme(text = element_text(family="serif"), panel.grid.minor = element_blank(), 
        axis.text.x = element_text(margin=margin(10,0,0,0)), legend.pos = "none",
        panel.border = element_rect(color = "grey",  fill = NA, size = 0.5)) +
  labs(title="Within-Group Variance")

d <- ggplot(datAllconv, aes(x=g, y=relBias_cov_W, colour=approach)) + 
  geom_rect(aes(xmin = Inf, xmax = -Inf, ymin = 0, ymax = Inf), alpha = 0.02, col="#fbc8c7", fill="#fbc8c7") +
  geom_rect(aes(xmin = Inf, xmax = -Inf, ymin = -Inf, ymax = 0), alpha = 0.02, col="#99e5e8", fill="#99e5e8") +
  stat_summary(geom = "line", linetype="dashed", alpha=0.5) +
  stat_summary(fun="mean", alpha=0.5) +
  scale_colour_manual(name="Approach", labels=allNames_new, values=allCols) +
  scale_y_continuous(name="Relative Bias (%)", expand=c(0.05,0.05), breaks=seq(-75, 0, 25)
  ) +
  scale_x_continuous(breaks=unique(datAllconv$g), labels=unique(datAllconv$g), limits=range(datAllconv$g), expand=c(0.05,0.05)) +
  coord_cartesian(ylim = c(-80, 1)) +
  theme_minimal() + 
  theme(text = element_text(family="serif"), panel.grid.minor = element_blank(), 
        axis.text.x = element_text(margin=margin(10,0,0,0)), legend.pos = "none",
        panel.border = element_rect(color = "grey",  fill = NA, size = 0.5)) +
  labs(title="Within-Group Covariance")

e <-
  ggplot(datAllconv, aes(x=g, y=relBias_ICC, colour=approach)) + 
  geom_rect(aes(xmin = Inf, xmax = -Inf, ymin = 0, ymax = Inf), alpha = 0.02, col="#fbc8c7", fill="#fbc8c7") +
  geom_rect(aes(xmin = Inf, xmax = -Inf, ymin = -Inf, ymax = 0), alpha = 0.02, col="#99e5e8", fill="#99e5e8") +
  stat_summary(geom = "line", linetype="dashed", alpha=0.5) +
  stat_summary(fun="mean", alpha=0.5) +
  scale_colour_manual(name="Approach", labels=allNames_new, values=allCols) +
  scale_y_continuous(name="Relative Bias (%)", expand=c(0.05,0.05), breaks=seq(-300, 0, 50)
  ) +
  scale_x_continuous(breaks=unique(datAllconv$g), labels=unique(datAllconv$g), limits=range(datAllconv$g), expand=c(0.05,0.05)) +
  coord_cartesian(ylim = c(-310, 0)) +
  theme_minimal() + 
  theme(text = element_text(family="serif"), panel.grid.minor = element_blank(), 
        axis.text.x = element_text(margin=margin(10,0,0,0)), legend.pos = "bottom", legend.direction = "vertical", legend.box = "horizontal",
        panel.border = element_rect(color = "grey",  fill = NA, size = 0.5)) +
  labs(title="Estimates of ICC")

legAll <- get_legend(e)

a + b + c + d + e + legAll + plot_layout(ncol=2, guides='collect') & theme(legend.position = "none")

# ggsave(
#   filename = "Fig_4.jpeg",
#   path = pathF,
#   device = "jpeg",
#   plot = last_plot(),
#   width = 500*3,
#   height = 700*3,
#   units = "px",
#   dpi = 300
# )


### Fig. 5: Overall Estimation Accuracy by Sample Characteristics

a <-
  ggplot(datAllconv, aes(x=g, y=relRMSE_B, colour=approach)) + 
  facet_grid(p_named ~ n_named, labeller=label_parsed) +
  stat_summary(geom = "line", linetype="dashed", alpha=0.5) +
  stat_summary(fun="mean", alpha=0.5) +
  scale_colour_manual(name="Approach", labels=allNames_new, values=allCols) +
  scale_y_continuous(name="Relative RMSE (%)", expand=c(0.01,0.01)) +
  scale_x_continuous(breaks=unique(datAllconv$g), labels=unique(datAllconv$g), limits=range(datAllconv$g), expand=c(0.05,0.05)) +
  coord_cartesian(ylim = c(0, 2500)) +
  theme_minimal() + 
  theme(text = element_text(family="serif"), panel.grid.minor = element_blank(), 
        axis.text.x = element_text(margin=margin(10,0,0,0)), legend.pos = "none",
        panel.border = element_rect(color = "grey",  fill = NA, size = 0.5)) +
  labs(title="Between-Group Parameters")

b <- ggplot(datAllconv, aes(x=g, y=relRMSE_W, colour=approach)) + 
  facet_grid(p_named ~ n_named, labeller=label_parsed) +
  stat_summary(geom = "line", linetype="dashed", alpha=0.5) +
  stat_summary(fun="mean", alpha=0.5) +
  scale_colour_manual(name="Approach", labels=allNames_new, values=allCols) +
  scale_y_continuous(name="Relative RMSE (%)", expand=c(0.01,0.01)) +
  scale_x_continuous(breaks=unique(datAllconv$g), labels=unique(datAllconv$g), limits=range(datAllconv$g), expand=c(0.05,0.05)) +
  coord_cartesian(ylim = c(0, 250)) +
  theme_minimal() + 
  theme(text = element_text(family="serif"), panel.grid.minor = element_blank(), 
        axis.text.x = element_text(margin=margin(10,0,0,0)), legend.pos = "none",
        panel.border = element_rect(color = "grey",  fill = NA, size = 0.5)) +
  labs(title="Within-Group Parameters")

c <- ggplot(datAllconv, aes(x=g, y=relRMSE_ICC, colour=approach)) + 
  facet_grid(p_named ~ n_named, labeller=label_parsed) +
  stat_summary(geom = "line", linetype="dashed", alpha=0.5) +
  stat_summary(fun="mean", alpha=0.5) +
  scale_colour_manual(name="Approach", labels=allNames_new, values=allCols) +
  scale_y_continuous(name="Relative RMSE (%)", expand=c(0.01,0.01), breaks=seq(0,750,250)) +
  scale_x_continuous(breaks=unique(datAllconv$g), labels=unique(datAllconv$g), limits=range(datAllconv$g), expand=c(0.05,0.05)) +
  coord_cartesian(ylim = c(0, 755)) +
  theme_minimal() + 
  theme(text = element_text(family="serif"), panel.grid.minor = element_blank(), 
        axis.text.x = element_text(margin=margin(10,0,0,0)), legend.pos = "none",
        panel.border = element_rect(color = "grey",  fill = NA, size = 0.5)) +
  labs(title="Estimates of ICC")

a + b + c + plot_layout(ncol=1, guides='collect') & theme(legend.position = "bottom")

# ggsave(
#   filename = "Fig_5.jpeg",
#   path = pathF,
#   device = "jpeg",
#   plot = last_plot(),
#   width = 600*3,
#   height = 1200*3,
#   units = "px",
#   dpi = 300
# )

# for information in article:
datMeans <- filter(datAllconv, approach=="LF", n==5)
aggregate(datMeans$relRMSE_ICC, list(datMeans$g), FUN=mean) 


### Fig. 6: Overall Estimation Accuracy by Population Characteristics

a <-
  ggplot(datAllconv, aes(x=ICC, y=relRMSE_B, colour=approach)) + 
  facet_grid(cor_B_named ~ cor_W_named, labeller=label_parsed) +
  stat_summary(geom = "line", linetype="dashed", alpha=0.5) +
  stat_summary(fun="mean", alpha=0.5) +
  scale_colour_manual(name="Approach", labels=allNames_new, values=allCols) +
  scale_y_continuous(name="Relative RMSE (%)", expand=c(0.01,0.01)) +#,
  scale_x_continuous(breaks=unique(datAllconv$ICC), labels=unique(datAllconv$ICC), limits=range(datAllconv$ICC), expand=c(0.05,0.05)) +
  theme_minimal() + 
  coord_cartesian(ylim = c(0, 1200)) + # then scales="free_y" deprecated
  theme(text = element_text(family="serif"), panel.grid.minor = element_blank(), 
        axis.text.x = element_text(margin=margin(10,0,0,0)), legend.pos = "none",
        panel.border = element_rect(color = "grey",  fill = NA, size = 0.5)) +
  labs(title="Between-Group Parameters")

b <- 
  ggplot(datAllconv, aes(x=ICC, y=relRMSE_W, colour=approach)) + 
  facet_grid(cor_B_named ~ cor_W_named, labeller=label_parsed) +
  stat_summary(geom = "line", linetype="dashed", alpha=0.5) +
  stat_summary(fun="mean", alpha=0.5) +
  scale_colour_manual(name="Approach", labels=allNames_new, values=allCols) +
  scale_y_continuous(name="Relative RMSE (%)", expand=c(0.01,0.01), breaks=seq(0,90, 30)) +
  scale_x_continuous(breaks=unique(datAllconv$ICC), labels=unique(datAllconv$ICC), limits=range(datAllconv$ICC), expand=c(0.05,0.05)) +
  theme_minimal() + 
  coord_cartesian(ylim = c(0, 90)) + # then scales="free_y" deprecated
  theme(text = element_text(family="serif"), panel.grid.minor = element_blank(), 
        axis.text.x = element_text(margin=margin(10,0,0,0)), legend.pos = "bottom",
        panel.border = element_rect(color = "grey",  fill = NA, size = 0.5)) +
  labs(title="Within-Group Parameters")

c <- 
  ggplot(datAllconv, aes(x=ICC, y=relRMSE_ICC, colour=approach)) + 
  facet_grid(cor_B_named ~ cor_W_named, labeller=label_parsed) +
  stat_summary(geom = "line", linetype="dashed", alpha=0.5) +
  stat_summary(fun="mean", alpha=0.5) +
  scale_colour_manual(name="Approach", labels=allNames_new, values=allCols) +
  scale_y_continuous(name="Relative RMSE (%)", expand=c(0.01,0.01)) +
  scale_x_continuous(breaks=unique(datAllconv$ICC), labels=unique(datAllconv$ICC), limits=range(datAllconv$ICC), expand=c(0.05,0.05)) +
  theme_minimal() + 
  coord_cartesian(ylim = c(0, 300)) + # then scales="free_y" deprecated
  theme(text = element_text(family="serif"), panel.grid.minor = element_blank(), 
        axis.text.x = element_text(margin=margin(10,0,0,0)), legend.pos = "bottom",
        panel.border = element_rect(color = "grey",  fill = NA, size = 0.5)) +
  labs(title="Estimates of ICC")

a + b + c + plot_layout(ncol=1, guides='collect') & theme(legend.position = "bottom")

# ggsave(
#   filename = "Fig_6.jpeg",
#   path = pathF,
#   device = "jpeg",
#   plot = last_plot(),
#   width = 600*3,
#   height = 900*3,
#   units = "px",
#   dpi = 300
# )




##### Appendix

### Fig. B1: The Standard WF Approach and Its Different Input and Estimator Possibilities

a <- ggplot(datWF, aes(y=conv, x=approach, col=Estimator, pch=Input)) + 
  stat_summary(fun="mean") +
  scale_y_continuous(name="Convergence Rate (%)", expand=c(0,0)) +
  scale_x_discrete(labels=NULL, name=NULL) +
  theme_minimal() + 
  scale_colour_manual(name="Estimator for S", values=WFCols) +
  theme(text = element_text(family="serif"), panel.grid.minor = element_blank(), legend.pos = "none",
        panel.border = element_rect(color = "grey",  fill = NA, linewidth = 0.5)) +
  coord_cartesian(ylim = c(0, 100)) +
  labs(title="A")

b <- ggplot(datWF, aes(x=approach, y=relRMSE_B, col=Estimator, pch=Input)) + 
  stat_summary(fun="mean") +
  scale_colour_manual(name="Estimator for S", values=WFCols) +
  scale_y_continuous(name="Relative RMSE (%)", expand=c(0,0), breaks=seq(0,500,100)) +
  scale_x_discrete(labels=NULL, name=NULL) +
  geom_text(label="hat(theta)[between]", aes(x=3.5, y=100), col="black", parse=TRUE) +
  theme_minimal() + 
  theme(text = element_text(family="serif"), panel.grid.minor = element_blank(), 
        axis.text.x = element_text(margin=margin(10,0,0,0)), legend.pos = "none",
        panel.border = element_rect(color = "grey",  fill = NA, size = 0.5)) +
  coord_cartesian(ylim = c(0, 500)) +
  labs(title="B")

c <- ggplot(datWF, aes(x=approach, y=relRMSE_W, col=Estimator, pch=Input)) + 
  stat_summary(fun="mean") +
  scale_colour_manual(name="Estimator for S", values=WFCols) +
  scale_y_continuous(name="Relative RMSE (%)", expand=c(0,0)) +
  scale_x_discrete(labels=NULL, name=NULL) +
  geom_text(label="hat(theta)[within]", aes(x=3.5, y=20), col="black", parse=TRUE) +
  theme_minimal() + 
  theme(text = element_text(family="serif"), panel.grid.minor = element_blank(), 
        axis.text.x = element_text(margin=margin(10,0,0,0)), legend.pos = "bottom",
        panel.border = element_rect(color = "grey",  fill = NA, size = 0.5)) +
  coord_cartesian(ylim = c(0, 100)) +
  labs(title="C")

a + b + c + plot_layout(ncol=3, guides='collect') & theme(legend.position = "bottom")

# ggsave(
#   filename = "Fig_B1.jpeg",
#   path = pathF,
#   device = "jpeg",
#   plot = last_plot(),
#   width = 550*3,
#   height = 200*3,
#   units = "px",
#   dpi = 300
# )



### Fig. B2: Computation Time by Sample Characteristics

a <- ggplot(datAllconv, aes(x=approach, y=time, colour=approach), alpha=0.9) + 
  stat_summary(geom = "line", linetype="dashed") +
  stat_summary(fun="mean") +
  scale_colour_manual(name="Approach", labels=allNames_new, values=allCols) +
  scale_y_continuous(name="Computation Time (s)", expand=c(0,0), 
                     breaks=seq(0,15,5)) + 
  scale_x_discrete(name=NULL, labels=NULL, expand=c(0.1,0.1)) +
  coord_cartesian(ylim = c(0, 15)) + # error when supplied in scale_y_cont see https://stackoverflow.com/questions/40415423/ggplot-stat-summary-when-axis-is-limited
  theme_minimal() + 
  theme(text = element_text(family="serif"), panel.grid.minor = element_blank(), 
        axis.text.x = element_text(margin=margin(10,0,0,0)), #legend.pos = "none",
        panel.border = element_rect(color = "grey",  fill = NA, size = 0.5)) +
  labs(title="A")

# for information in article:
aggregate(datAllconv$time, list(datAllconv$approach), FUN="mean", na.rm = TRUE)

b <- ggplot(datAllconv, aes(x=p, y=time, colour=approach)) + 
  stat_summary(geom = "line", linetype="dashed") +
  stat_summary(fun="mean") +
  scale_colour_manual(name="Approach", labels=allNames_new, values=allCols) +
  scale_y_continuous(name="Computation Time (s)", expand=c(0,0), breaks=seq(0,30,5)) +
  scale_x_continuous(expand=c(0.05,0.05), breaks=unique(dat$p)) +
  coord_cartesian(ylim = c(0, 15)) +
  theme_minimal() + 
  theme(text = element_text(family="serif"), panel.grid.minor = element_blank(), 
        axis.text.x = element_text(margin=margin(10,0,0,0)), #legend.pos = "none",
        panel.border = element_rect(color = "grey",  fill = NA, size = 0.5)) +
  labs(title="B")

c <- ggplot(datAllconv, aes(x=n, y=time, colour=approach)) + 
  stat_summary(geom = "line", linetype="dashed") +
  stat_summary(fun="mean") +
  scale_colour_manual(name="Approach", labels=allNames_new, values=allCols) +
  scale_y_continuous(name="Computation Time (s)", expand=c(0,0), breaks=seq(0,30,5)) +
  scale_x_continuous(expand=c(0.05,0.05), breaks=unique(dat$n)) +
  coord_cartesian(ylim = c(0, 15)) +
  theme_minimal() + 
  theme(text = element_text(family="serif"), panel.grid.minor = element_blank(), 
        axis.text.x = element_text(margin=margin(10,0,0,0)), 
        panel.border = element_rect(color = "grey",  fill = NA, size = 0.5)) +
  labs(title="C")

a + b + c + plot_layout(ncol=3, guides='collect', widths = c(0.5, 1, 1)) & theme(legend.position = "bottom")

# ggsave(
#   filename = "Fig_B2.jpeg",
#   path = pathF,
#   device = "jpeg",
#   plot = last_plot(),
#   width = 800*3,
#   height = 300*3,
#   units = "px",
#   dpi = 300
# )


### Fig. B3: Negatively Estimated Variances at the Between-Group Level and ICC

a <- ggplot(datAllconv, aes(x=g, y=negVar_B, colour=approach)) + 
  facet_grid(rows=vars(ICC_named), cols=vars(n_named), labeller=label_parsed) +
  stat_summary(geom = "line", linetype="dashed") +
  stat_summary(fun="mean") +
  scale_colour_manual(name="Approach", labels=allNames_new, values=allCols) +
  scale_y_continuous(name="Any Negative Variance in Model (%)", expand=c(0.01,0.01), limits=c(0,100)) +
  scale_x_continuous(breaks=unique(datAllconv$g), labels=unique(datAllconv$g), limits=range(datAllconv$g), expand=c(0.05,0.05)) +
  theme_minimal() + 
  theme(text = element_text(family="serif"), panel.grid.minor = element_blank(), 
        axis.text.x = element_text(margin=margin(10,0,0,0)), legend.pos = "none",
        panel.border = element_rect(color = "grey",  fill = NA, size = 0.5)) +
  labs(title="A")

b <- ggplot(datAllconv, aes(x=g, y=ICC_Sigma_theta_neg, colour=approach)) + 
  facet_grid(rows=vars(ICC_named), cols=vars(n_named), labeller=label_parsed) +
  stat_summary(geom = "line", linetype="dashed") +
  stat_summary(fun="mean") +
  scale_colour_manual(name="Approach", labels=allNames_new, values=allCols) +
  scale_y_continuous(name="Negative Estimates of ICC (%)", expand=c(0.01,0.01), limits=c(0,100)) +
  scale_x_continuous(breaks=unique(datAllconv$g), labels=unique(datAllconv$g), limits=range(datAllconv$g), expand=c(0.05,0.05)) +
  theme_minimal() + 
  theme(text = element_text(family="serif"), panel.grid.minor = element_blank(), 
        axis.text.x = element_text(margin=margin(10,0,0,0)), legend.pos = "bottom",
        panel.border = element_rect(color = "grey",  fill = NA, size = 0.5)) +
  labs(title="B")

a + b + plot_layout(ncol=2, guides='collect') & theme(legend.position = "bottom")

# ggsave(
#   filename = "Fig_B3.jpeg",
#   path = pathF,
#   device = "jpeg",
#   plot = last_plot(),
#   width = 800*3,
#   height = 400*3,
#   units = "px",
#   dpi = 300
# )
