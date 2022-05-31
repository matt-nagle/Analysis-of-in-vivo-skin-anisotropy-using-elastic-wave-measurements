library(bpnreg)
library(circular)
library(rstanarm)
library(ggplot2)
library(ggpubr)
library(egg)

#####====================== Import Datasets ======================#####

setwd("C:/Users/mattn/OneDrive - University College Dublin/Documents/00. Research/99. Publications/Analysis of in-vivo skin anisotropy using elastic wave measurements")


all_fit_ellipse_params <- read.csv("Ellipse Params - Natural and Stretched Config.csv", header = TRUE)

# Rescale the area for numerical stability
all_fit_ellipse_params$Rescaled_Area <- all_fit_ellipse_params$Area/1000





#####====================== Bayesian multivariate outcome regression ======================#####

##### Model for Effect of Age Gender & Config #####

bayesianmultivar <- stan_mvmer(formula = list(log(Eccentricity) ~ 1 + Age + Gender + Configuration + (1 | Subject_ID),
                                Rescaled_Area ~ 1  + Age + Gender + Configuration + (1 | Subject_ID),
                                Semi.minor_Axis ~ 1  + Age + Gender + Configuration + (1 | Subject_ID)),
                               data = all_fit_ellipse_params,
                               iter = 20000,
                               seed = 42)

# Examine model summary
summary(bayesianmultivar)

# Extract the chains of posterior estimates
sims3 <- as.matrix(bayesianmultivar)
sims3_df <- data.frame(sims3)


## Create ggplot variables for the posterior distributions of the coefficients of interest
ecc_age <- ggplot(data = sims3_df, mapping = aes(x = y1.Age)) +
  geom_vline(xintercept = bayesianmultivar$stan_summary["y1|Age", "2.5%"], size = 1, color = "black") +
  geom_vline(xintercept = bayesianmultivar$stan_summary["y1|Age", "97.5%"], size = 1, color = "black") +
  annotate("rect", 
           xmin = bayesianmultivar$stan_summary["y1|Age", "2.5%"], 
           xmax = bayesianmultivar$stan_summary["y1|Age", "97.5%"], 
           ymin = 0, ymax = Inf, fill = "grey", alpha = 0.5) +
  geom_histogram(aes(y = ..density..), bins = 30, alpha = 0.5, position = "identity", colour = "black") +
  geom_density(lwd = 1) +
  ggtitle("Age") +
  ylab("Posterior Density") +
  xlab("Age Coefficient") +
  #ylim(0,5) +
  theme_light()

ecc_male <- ggplot(data = sims3_df, mapping = aes(x = y1.GenderMale)) +
  geom_vline(xintercept = bayesianmultivar$stan_summary["y1|GenderMale", "2.5%"], size = 1, color = "black") +
  geom_vline(xintercept = bayesianmultivar$stan_summary["y1|GenderMale", "97.5%"], size = 1, color = "black") +
  annotate("rect", 
           xmin = bayesianmultivar$stan_summary["y1|GenderMale", "2.5%"], 
           xmax = bayesianmultivar$stan_summary["y1|GenderMale", "97.5%"], 
           ymin = 0, ymax = Inf, fill = "grey", alpha = 0.5) +
  geom_histogram(aes(y = ..density..), bins = 30, alpha = 0.5, position = "identity", colour = "black") +
  geom_density(lwd = 1) +
  ggtitle("Gender") +
  ylab("Posterior Density") +
  xlab("Male Coefficient") +
  geom_vline(xintercept = 0, linetype = "dashed", size = 1.25, color = "red") +
  #ylim(0,5) +
  theme_light()

ecc_config <- ggplot(data = sims3_df, mapping = aes(x = y1.ConfigurationStretched)) +
  geom_vline(xintercept = bayesianmultivar$stan_summary["y1|ConfigurationStretched", "2.5%"], size = 1, color = "black") +
  geom_vline(xintercept = bayesianmultivar$stan_summary["y1|ConfigurationStretched", "97.5%"], size = 1, color = "black") +
  annotate("rect", 
           xmin = bayesianmultivar$stan_summary["y1|ConfigurationStretched", "2.5%"], 
           xmax = bayesianmultivar$stan_summary["y1|ConfigurationStretched", "97.5%"], 
           ymin = 0, ymax = Inf, fill = "grey", alpha = 0.5) +
  geom_histogram(aes(y = ..density..), bins = 30, alpha = 0.5, position = "identity", colour = "black") +
  geom_density(lwd = 1) +
  ggtitle("Configuration") +
  ylab("Posterior Density") +
  xlab("Stretched Coefficient") +
  geom_vline(xintercept = 0, linetype = "dashed", size = 1.25, color = "red") +
  #ylim(0,5) +
  theme_light()





area_age <- ggplot(data = sims3_df, mapping = aes(x = y2.Age)) +
  geom_vline(xintercept = bayesianmultivar$stan_summary["y2|Age", "2.5%"], size = 1, color = "black") +
  geom_vline(xintercept = bayesianmultivar$stan_summary["y2|Age", "97.5%"], size = 1, color = "black") +
  annotate("rect", 
           xmin = bayesianmultivar$stan_summary["y2|Age", "2.5%"], 
           xmax = bayesianmultivar$stan_summary["y2|Age", "97.5%"], 
           ymin = 0, ymax = Inf, fill = "grey", alpha = 0.5) +
  geom_histogram(aes(y = ..density..), bins = 30, alpha = 0.5, position = "identity", colour = "black") +
  geom_density(lwd = 1) +
  ggtitle("Age") +
  ylab("Posterior Density") +
  xlab("Age Coefficient") +
  geom_vline(xintercept = 0, linetype = "dashed", size = 1.25, color = "red") +
  #ylim(0,1) +
  theme_light()

area_male <- ggplot(data = sims3_df, mapping = aes(x = y2.GenderMale)) +
  geom_vline(xintercept = bayesianmultivar$stan_summary["y2|GenderMale", "2.5%"], size = 1, color = "black") +
  geom_vline(xintercept = bayesianmultivar$stan_summary["y2|GenderMale", "97.5%"], size = 1, color = "black") +
  annotate("rect", 
           xmin = bayesianmultivar$stan_summary["y2|GenderMale", "2.5%"], 
           xmax = bayesianmultivar$stan_summary["y2|GenderMale", "97.5%"], 
           ymin = 0, ymax = Inf, fill = "grey", alpha = 0.5) +
  geom_histogram(aes(y = ..density..), bins = 30, alpha = 0.5, position = "identity", colour = "black") +
  geom_density(lwd = 1) +
  ggtitle("Gender") +
  ylab("Posterior Density") +
  xlab("Male Coefficient") +
  geom_vline(xintercept = 0, linetype = "dashed", size = 1.25, color = "red") +
  #ylim(0,1) +
  theme_light()

area_config <- ggplot(data = sims3_df, mapping = aes(x = y2.ConfigurationStretched)) +
  geom_vline(xintercept = bayesianmultivar$stan_summary["y2|ConfigurationStretched", "2.5%"], size = 1, color = "black") +
  geom_vline(xintercept = bayesianmultivar$stan_summary["y2|ConfigurationStretched", "97.5%"], size = 1, color = "black") +
  annotate("rect", 
           xmin = bayesianmultivar$stan_summary["y2|ConfigurationStretched", "2.5%"], 
           xmax = bayesianmultivar$stan_summary["y2|ConfigurationStretched", "97.5%"], 
           ymin = 0, ymax = Inf, fill = "grey", alpha = 0.5) +
  geom_histogram(aes(y = ..density..), bins = 30, alpha = 0.5, position = "identity", colour = "black") +
  geom_density(lwd = 1) +
  ggtitle("Configuration") +
  ylab("Posterior Density") +
  xlab("Stretched Coefficient") +
  #ylim(0,1) +
  theme_light()





semi_min_age <- ggplot(data = sims3_df, mapping = aes(x = y3.Age)) +
  geom_vline(xintercept = bayesianmultivar$stan_summary["y3|Age", "2.5%"], size = 1, color = "black") +
  geom_vline(xintercept = bayesianmultivar$stan_summary["y3|Age", "97.5%"], size = 1, color = "black") +
  annotate("rect",
           xmin = bayesianmultivar$stan_summary["y3|Age", "2.5%"],
           xmax = bayesianmultivar$stan_summary["y3|Age", "97.5%"],
           ymin = 0, ymax = Inf, fill = "grey", alpha = 0.5) +
  geom_histogram(aes(y = ..density..), bins = 30, alpha = 0.5, position = "identity", colour = "black") +
  geom_density(lwd = 1) +
  ggtitle("Age") +
  ylab("Posterior Density") +
  xlab("Age Coefficient") +
  geom_vline(xintercept = 0, linetype = "dashed", size = 1.25, color = "red") +
  #ylim(0,1) +
  theme_light()

semi_min_male <- ggplot(data = sims3_df, mapping = aes(x = y3.GenderMale)) +
  geom_vline(xintercept = bayesianmultivar$stan_summary["y3|GenderMale", "2.5%"], size = 1, color = "black") +
  geom_vline(xintercept = bayesianmultivar$stan_summary["y3|GenderMale", "97.5%"], size = 1, color = "black") +
  annotate("rect",
           xmin = bayesianmultivar$stan_summary["y3|GenderMale", "2.5%"],
           xmax = bayesianmultivar$stan_summary["y3|GenderMale", "97.5%"],
           ymin = 0, ymax = Inf, fill = "grey", alpha = 0.5) +
  geom_histogram(aes(y = ..density..), bins = 30, alpha = 0.5, position = "identity", colour = "black") +
  geom_density(lwd = 1) +
  ggtitle("Gender") +
  ylab("Posterior Density") +
  xlab("Male Coefficient") +
  geom_vline(xintercept = 0, linetype = "dashed", size = 1.25, color = "red") +
  #ylim(0,1) +
  theme_light()

semi_min_config <- ggplot(data = sims3_df, mapping = aes(x = y3.ConfigurationStretched)) +
  geom_vline(xintercept = bayesianmultivar$stan_summary["y3|ConfigurationStretched", "2.5%"], size = 1, color = "black") +
  geom_vline(xintercept = bayesianmultivar$stan_summary["y3|ConfigurationStretched", "97.5%"], size = 1, color = "black") +
  annotate("rect", 
           xmin = bayesianmultivar$stan_summary["y3|ConfigurationStretched", "2.5%"], 
           xmax = bayesianmultivar$stan_summary["y3|ConfigurationStretched", "97.5%"], 
           ymin = 0, ymax = Inf, fill = "grey", alpha = 0.5) +
  geom_histogram(aes(y = ..density..), bins = 30, alpha = 0.5, position = "identity", colour = "black") +
  geom_density(lwd = 1) +
  ggtitle("Configuration") +
  ylab("Posterior Density") +
  xlab("Stretched Coefficient") +
  #ylim(0,1) +
  theme_light()



# create a separate row of plots for each of the output variables and annotate as such
plot_ecc_coeff_row <- ggarrange(ecc_age, 
                                ecc_male + theme(axis.title.y = element_blank()), 
                                ecc_config + theme(axis.title.y = element_blank()),
                                nrow = 1)

annotate_figure(plot_ecc_coeff_row, 
                left = text_grob("Log Eccentricity Output\n", face = "bold", 
                                 size = 18, rot = 90),
                fig.lab.pos = "top.left")

ggsave(file = "Figures/figure6a.eps", units = "px", width = 2426, height = 912, device=cairo_ps)



plot_area_coeff_row <- ggarrange(area_age, 
                                area_male + theme(axis.title.y = element_blank()), 
                                area_config + theme(axis.title.y = element_blank()),
                                nrow = 1)

annotate_figure(plot_area_coeff_row, 
                left = text_grob("Area Output\n", face = "bold", 
                                 size = 18, rot = 90),
                fig.lab.pos = "top.left")

ggsave(file = "Figures/figure6b.eps", units = "px", width = 2426, height = 912, device=cairo_ps)




plot_semi_min_coeff_row <- ggarrange(semi_min_age, 
                                     semi_min_male + theme(axis.title.y = element_blank()), 
                                     semi_min_config + theme(axis.title.y = element_blank()),
                                     nrow = 1)

annotate_figure(plot_semi_min_coeff_row, 
                left = text_grob("Semi-minor Axis Output\n", face = "bold", 
                                 size = 18, rot = 90),
                fig.lab.pos = "top.left")

ggsave(file = "Figures/figure6c.eps", units = "px", width = 2426, height = 912, device=cairo_ps)














#####====================== Circular ANNOVA style model for Angle vs Config ======================#####

fit.config <- bpnr(pred.I = Angle_to_semi.maj_Rads ~ 1 + Configuration, 
                   data = all_fit_ellipse_params,
                   its = 10000, burn = 100, n.lag = 3, seed = 101)

fit.config

# Follwing the discussion in Cremers:
# Calculate the highest posterior density interval (HPDI) for both the Natural and Stretched Config

Natural <- atan2(fit.config$beta2[,1], fit.config$beta1[,1])
Natural <- circular(Natural, type = "angles", units = "radians")
hpd_est_circ(Natural) #0.95 HPDI for Natural Config

Stretched <- atan2(fit.config$beta2[,1] + fit.config$beta2[,2], fit.config$beta1[,1] + fit.config$beta1[,2])
Stretched <- circular(Stretched, type = "angles", units = "radians")
hpd_est_circ(Stretched) #0.95 HPDI for Stretched Config



Natural_df <- data.frame(cbind(Natural))
Natural_df$Configuration <- "Natural"
colnames(Natural_df) <- c("HPD", "Configuration")

Stretched_df <- data.frame(cbind(Stretched))
Stretched_df$Configuration <- "Stretched"
colnames(Stretched_df) <- c("HPD", "Configuration")

all_HPD_df <- rbind(Natural_df, Stretched_df)

ggplot(data = all_HPD_df, mapping = aes(x = HPD, group = Configuration, fill = Configuration)) +
  geom_vline(xintercept = fit.config$circ.coef.means[1,4], size = 1, color = "#F8766D") +
  geom_vline(xintercept = fit.config$circ.coef.means[1,5], size = 1, color = "#F8766D") +
  annotate("rect", 
           xmin = fit.config$circ.coef.means[1,4], 
           xmax = fit.config$circ.coef.means[1,5], 
           ymin = 0, ymax = Inf, fill = "#F8766D", alpha = 0.2) +
  geom_histogram(aes(y = ..density..), bins = 50, alpha = 0.5, position = "identity", colour = "black") +
  geom_vline(xintercept = fit.config$circ.coef.means[2,4], size = 1, color = "#00BFC4") +
  geom_vline(xintercept = fit.config$circ.coef.means[2,5], size = 1, color = "#00BFC4") +
  annotate("rect", 
           xmin = fit.config$circ.coef.means[2,4], 
           xmax = fit.config$circ.coef.means[2,5], 
           ymin = 0, ymax = Inf, fill = "#00BFC4", alpha = 0.2) +
  geom_histogram(aes(y = ..density..), bins = 50, alpha = 0.5, position = "identity", colour = "black") +
  ggtitle("Natural vs Stretched Configuration HPD Interval") +
  xlab("Angle (Rad)") +
  ylab("Posterior Density") +
  theme_light() +
  theme(legend.position = "bottom")

ggsave(file = "Figures/figure9.eps", units = "px", width = 1724, height = 1294, device=cairo_ps)
