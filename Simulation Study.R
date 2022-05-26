library(conicfit)
library(ggplot2)
library(gridExtra)
library(ggforce)


setwd("C:/Users/mattn/OneDrive - University College Dublin/Documents/00. Research/99. Publications/Analysis of in-vivo skin properties using elastic wave measurements")


##### ======================= Simulation Study with Various Eccentricities ======================= #####

n_sim <- 10000

eccentricities <- c(0.3, 0.5, 0.7, 0.9)
sigmas <- c(0.1, 1, 10, 20, 25, 30)

set.seed(1)
for(e in eccentricities){
  for(s in sigmas){
    print(paste("e = ", e, ",  s = ", s))
    
    ## Create dataframe to store results
    
    sim_results <- data.frame(matrix(ncol = 51, nrow = n_sim))
    colnames(sim_results) <- c("Simulation_#", "Max_RRT", "Min_RRT", "Max_Min_Ratio",
                               "Centre_x", "Centre_y", "Axis_1", "Axis_2", "Semi-major_Axis", "Semi-minor_Axis", 
                               "Angle_Rads", "Angle_Deg", "Angle_to_semi-maj_Rads", "Angle_to_semi-maj_Deg", 
                               "Eccentricity", 1:36)
    for(i in 1:n_sim){
      
      semi_min_axis <- 160 # Let's use the semi min axis length of 160 (avg from reviscometer data)
      angles <- seq(0,350,10)
      angles_rad <- deg2rad(angles)
      
      # Create core shape in polar coordinates
      base <- sqrt((semi_min_axis^2)/(1 - (e)^2*cos(angles_rad)^2))
      
      # Add random noise to core shape
      sim_data <- base + rnorm(36, mean = 0, sd = s)
      
      sim_results[i, "Simulation_#"] <-  i
      sim_results[i, "Max_RRT"] <- max(sim_data)
      sim_results[i, "Min_RRT"] <- min(sim_data)
      sim_results[i, "Max_Min_Ratio"] <- max(sim_data)/min(sim_data)
      sim_results[i, 16:51] <- sim_data
      
      Angle <- seq(0,350,10)
      
      # Convert Measurements from Polar to Cartesian
      x <- sim_data*cos(Angle*(pi/180)) 
      y <- sim_data*sin(Angle*(pi/180))
      
      subject_RRT_cart <- data.frame(x = x, y = y)
      
      
      # The following returns algebraic parameters for general conic equation 
      # f(x,y) = ax^2 + bxy + cy^2 + dx + ey + f = 0
      ellipse_fit <- EllipseDirectFit(subject_RRT_cart)
      # Converts algebraic parameters to geometric
      ellipse_fit_geo <- AtoG(ellipse_fit)$ParG 
      #calculates xy coords for fit ellipse
      ellipse_xy <- calculateEllipse(ellipse_fit_geo[1],ellipse_fit_geo[2],ellipse_fit_geo[3],
                                     ellipse_fit_geo[4],(180/pi)*ellipse_fit_geo[5], steps = 100)
      

      if (ellipse_fit_geo[3] < ellipse_fit_geo[4]) {
        if (ellipse_fit_geo[5] >= pi/2){
          angle_to_semi_maj <- ellipse_fit_geo[5] - pi/2
        } else {
          angle_to_semi_maj <- ellipse_fit_geo[5] + pi/2
        }
        ellipse_fit_geo[5] + pi/2
      } else {
        angle_to_semi_maj <- ellipse_fit_geo[5]
      }
      
      sim_results[i,5:15] <- c(ellipse_fit_geo[1],ellipse_fit_geo[2], # Centre x, Centre y
                               ellipse_fit_geo[3],ellipse_fit_geo[4], # Axis 1, Axis 2
                               max(c(ellipse_fit_geo[3],ellipse_fit_geo[4])), # Semi-maj Axis
                               min(c(ellipse_fit_geo[3],ellipse_fit_geo[4])), # Semi-min Axis
                               ellipse_fit_geo[5], # Angle (Rads)
                               ellipse_fit_geo[5]*(180/pi), # Angle (Deg)
                               angle_to_semi_maj, # Angle to the semi-major axis (Rads)
                               angle_to_semi_maj*(180/pi), # Angle to the semi-major axis (Degs)
                               sqrt(1 - (min(c(ellipse_fit_geo[3],ellipse_fit_geo[4]))/max(c(ellipse_fit_geo[3],ellipse_fit_geo[4])))^2)) # Eccentricity
    }
    
    # create unique dataset for each eccentricity and noise combination
    assign(paste("sim_results", s, "e", e, sep=""), sim_results)
  }
}



sim_results0.1e0.3$Sigma <- 0.1
sim_results0.1e0.5$Sigma <- 0.1
sim_results0.1e0.7$Sigma <- 0.1
sim_results0.1e0.9$Sigma <- 0.1

sim_results1e0.3$Sigma <- 1
sim_results1e0.5$Sigma <- 1
sim_results1e0.7$Sigma <- 1
sim_results1e0.9$Sigma <- 1

sim_results10e0.3$Sigma <- 10
sim_results10e0.5$Sigma <- 10
sim_results10e0.7$Sigma <- 10
sim_results10e0.9$Sigma <- 10

sim_results20e0.3$Sigma <- 20
sim_results20e0.5$Sigma <- 20
sim_results20e0.7$Sigma <- 20
sim_results20e0.9$Sigma <- 20

sim_results25e0.3$Sigma <- 25
sim_results25e0.5$Sigma <- 25
sim_results25e0.7$Sigma <- 25
sim_results25e0.9$Sigma <- 25

sim_results30e0.3$Sigma <- 30
sim_results30e0.5$Sigma <- 30
sim_results30e0.7$Sigma <- 30
sim_results30e0.9$Sigma <- 30

all_sim_resultse0.3 <- rbind(sim_results0.1e0.3, sim_results1e0.3, sim_results10e0.3,
                             sim_results20e0.3, sim_results25e0.3, sim_results30e0.3)
all_sim_resultse0.5 <- rbind(sim_results0.1e0.5, sim_results1e0.5, sim_results10e0.5,
                             sim_results20e0.5, sim_results25e0.5, sim_results30e0.5)
all_sim_resultse0.7 <- rbind(sim_results0.1e0.7, sim_results1e0.7, sim_results10e0.7,
                             sim_results20e0.7, sim_results25e0.7, sim_results30e0.7)
all_sim_resultse0.9 <- rbind(sim_results0.1e0.9, sim_results1e0.9, sim_results10e0.9,
                             sim_results20e0.9, sim_results25e0.9, sim_results30e0.9)



all_sim_resultse0.5$TrueEccentricity <- 0.5
all_sim_resultse0.7$TrueEccentricity <- 0.7
all_sim_resultse0.9$TrueEccentricity <- 0.9

all_sim_resultse0.5e0.7e0.9 <- rbind(all_sim_resultse0.5, all_sim_resultse0.7, all_sim_resultse0.9)
all_sim_resultse0.5e0.7e0.9$Sigma <- as.factor(all_sim_resultse0.5e0.7e0.9$Sigma)
all_sim_resultse0.5e0.7e0.9$TrueEccentricity <- as.factor(all_sim_resultse0.5e0.7e0.9$TrueEccentricity)




##### ======================= Creating Figures for paper =======================#####

## Example Dataset where Eccentricity outperforms Anisotropic Ratio
s <- 5
e <- 0.7
semi_min_axis <- 160 # Let's use the semi min axis length of 160 (avg from reviscometer data)
semi_maj_axis <- semi_min_axis/sqrt(1 - e^2)
angles <- seq(0,350,10)
angles_rad <- deg2rad(angles)
base <- sqrt((semi_min_axis^2)/(1 - (e)^2*cos(angles_rad)^2))

set.seed(40)
sim_data <- base + rnorm(36, mean = 0, sd = s)

sim_data[2] <- sim_data[2] + 70
sim_data[12] <- sim_data[12] - 50

Angle <- seq(0,350,10)

x <- sim_data*cos(Angle*(pi/180)) # Convert Measurements from Polar to Cartesian
y <- sim_data*sin(Angle*(pi/180)) # Convert Measurements from Polar to Cartesian

x_outliers <- sim_data[c(2,12)]*cos(Angle[c(2,12)]*(pi/180))
y_outliers <- sim_data[c(2,12)]*sin(Angle[c(2,12)]*(pi/180))

subject_RRT_cart <- data.frame(x = x, y = y)
outliers_cart <- data.frame(x = x_outliers, y = y_outliers)


# The following returns algebraic parameters for general conic equation 
# f(x,y) = ax^2 + bxy + cy^2 + dx + ey + f = 0
ellipse_fit <- EllipseDirectFit(subject_RRT_cart)
ellipse_fit_geo <- AtoG(ellipse_fit)$ParG # Converts algebraic parameters to geometric
ellipse_xy <- calculateEllipse(ellipse_fit_geo[1],ellipse_fit_geo[2],ellipse_fit_geo[3],
                               ellipse_fit_geo[4],(180/pi)*ellipse_fit_geo[5], steps = 100) #calculates xy coords for fit ellipse

ellipse_xy_df <- data.frame(x = ellipse_xy[,1], y = ellipse_xy[,2])

Anisotropic_Ratio <- max(sim_data)/min(sim_data)
Eccentricity <- sqrt(1 - (min(c(ellipse_fit_geo[3],ellipse_fit_geo[4]))/max(c(ellipse_fit_geo[3],ellipse_fit_geo[4])))^2) # Eccentricity


ggplot(data = subject_RRT_cart, mapping = aes(x = x, y = y)) + 
  geom_point(size = 3) + 
  geom_point(data = outliers_cart, fill = "grey", color = "black", size = 5, pch = 21) +
  geom_ellipse(aes(x0 = 0, y0 = 0, a = semi_maj_axis, b = semi_min_axis, angle = 0), color = "red", lwd = 1, lty = "dashed") +
  geom_ellipse(aes(x0 = ellipse_fit_geo[1], y0 = ellipse_fit_geo[2], a = max(c(ellipse_fit_geo[3],ellipse_fit_geo[4])), b = min(c(ellipse_fit_geo[3],ellipse_fit_geo[4])), angle = ellipse_fit_geo[5]), color = "blue", lwd = 1) +
  geom_hline(yintercept = 0, color = "black", size = 1) +
  geom_vline(xintercept = 0, color = "black", size = 1) +
  ggtitle("True Eccentricity = 0.7, True AR = 1.4003") + 
  ylim(-210, 210) +
  xlim(-300, 300) +
  theme_light()

ggsave(file = "Figures/figure4-noannotations.eps", units = "px", width = 1714, height = 1274)








## Example of Simulation Data with Noise
e <- 0.7
semi_min_axis <- 160 # Let's use the semi min axis length of 160 (avg from reviscometer data)
semi_maj_axis <- semi_min_axis/sqrt(1 - e^2)

Angle <- seq(0,350,10)

sim_data <- all_sim_resultse0.5e0.7e0.9[all_sim_resultse0.5e0.7e0.9$`Simulation_#` == 1 &
                                          all_sim_resultse0.5e0.7e0.9$TrueEccentricity == 0.7 & 
                                          all_sim_resultse0.5e0.7e0.9$Sigma == 1, 16:51]

x <- t(sim_data)*cos(Angle*(pi/180)) # Convert Measurements from Polar to Cartesian
y <- t(sim_data)*sin(Angle*(pi/180)) # Convert Measurements from Polar to Cartesian

subject_RRT_cart <- data.frame(x = x, y = y)


ggplot(data = subject_RRT_cart, mapping = aes(x = x, y = y)) + 
  geom_point(size = 3) + 
  geom_ellipse(aes(x0 = 0, y0 = 0, a = semi_maj_axis, b = semi_min_axis, angle = 0), color = "red", lwd = 1, lty = "dashed") +
  geom_hline(yintercept = 0, color = "black", size = 1) +
  geom_vline(xintercept = 0, color = "black", size = 1) +
  ggtitle("True Eccentricity = 0.7, Sigma = 1") + 
  ylim(-210, 210) +
  xlim(-300, 300) +
  theme_light()

ggsave(file = "Figures/figure3a.eps", units = "px", width = 1724, height = 1284)




sim_data <- all_sim_resultse0.5e0.7e0.9[all_sim_resultse0.5e0.7e0.9$`Simulation_#` == 1 &
                                          all_sim_resultse0.5e0.7e0.9$TrueEccentricity == 0.7 & 
                                          all_sim_resultse0.5e0.7e0.9$Sigma == 10, 16:51]

x <- t(sim_data)*cos(Angle*(pi/180)) # Convert Measurements from Polar to Cartesian
y <- t(sim_data)*sin(Angle*(pi/180)) # Convert Measurements from Polar to Cartesian

subject_RRT_cart <- data.frame(x = x, y = y)


ggplot(data = subject_RRT_cart, mapping = aes(x = x, y = y)) + 
  geom_point(size = 3) + 
  geom_ellipse(aes(x0 = 0, y0 = 0, a = semi_maj_axis, b = semi_min_axis, angle = 0), color = "red", lwd = 1, lty = "dashed") +
  geom_hline(yintercept = 0, color = "black", size = 1) +
  geom_vline(xintercept = 0, color = "black", size = 1) +
  ggtitle("True Eccentricity = 0.7, Sigma = 10") + 
  ylim(-210, 210) +
  xlim(-300, 300) +
  theme_light()

ggsave(file = "Figures/figure3b.eps", units = "px", width = 1724, height = 1284)





sim_data <- all_sim_resultse0.5e0.7e0.9[all_sim_resultse0.5e0.7e0.9$`Simulation_#` == 1 &
                                          all_sim_resultse0.5e0.7e0.9$TrueEccentricity == 0.7 & 
                                          all_sim_resultse0.5e0.7e0.9$Sigma == 20, 16:51]

x <- t(sim_data)*cos(Angle*(pi/180)) # Convert Measurements from Polar to Cartesian
y <- t(sim_data)*sin(Angle*(pi/180)) # Convert Measurements from Polar to Cartesian

subject_RRT_cart <- data.frame(x = x, y = y)


ggplot(data = subject_RRT_cart, mapping = aes(x = x, y = y)) + 
  geom_point(size = 3) + 
  geom_ellipse(aes(x0 = 0, y0 = 0, a = semi_maj_axis, b = semi_min_axis, angle = 0), color = "red", lwd = 1, lty = "dashed") +
  geom_hline(yintercept = 0, color = "black", size = 1) +
  geom_vline(xintercept = 0, color = "black", size = 1) +
  ggtitle("True Eccentricity = 0.7, Sigma = 20") + 
  ylim(-210, 210) +
  xlim(-300, 300) +
  theme_light()

ggsave(file = "Figures/figure3c.eps", units = "px", width = 1724, height = 1284)



sim_data <- all_sim_resultse0.5e0.7e0.9[all_sim_resultse0.5e0.7e0.9$`Simulation_#` == 1 &
                                          all_sim_resultse0.5e0.7e0.9$TrueEccentricity == 0.7 & 
                                          all_sim_resultse0.5e0.7e0.9$Sigma == 30, 16:51]

x <- t(sim_data)*cos(Angle*(pi/180)) # Convert Measurements from Polar to Cartesian
y <- t(sim_data)*sin(Angle*(pi/180)) # Convert Measurements from Polar to Cartesian

subject_RRT_cart <- data.frame(x = x, y = y)


ggplot(data = subject_RRT_cart, mapping = aes(x = x, y = y)) + 
  geom_point(size = 3) + 
  geom_ellipse(aes(x0 = 0, y0 = 0, a = semi_maj_axis, b = semi_min_axis, angle = 0), color = "red", lwd = 1, lty = "dashed") +
  geom_hline(yintercept = 0, color = "black", size = 1) +
  geom_vline(xintercept = 0, color = "black", size = 1) +
  ggtitle("True Eccentricity = 0.7, Sigma = 30") + 
  ylim(-210, 210) +
  xlim(-300, 300) +
  theme_light()

ggsave(file = "Figures/figure3d.eps", units = "px", width = 1724, height = 1284)



















## Boxplot of Simulation Study Results

# Create user-defined function, which extracts legends from ggplots
extract_legend <- function(my_ggp) {
  step1 <- ggplot_gtable(ggplot_build(my_ggp))
  step2 <- which(sapply(step1$grobs, function(x) x$name) == "guide-box")
  step3 <- step1$grobs[[step2]]
  return(step3)
}

p1_legend <- ggplot(data = all_sim_resultse0.5e0.7e0.9, mapping = aes(x = Sigma, y = Eccentricity, fill = TrueEccentricity)) + 
  geom_boxplot() + 
  geom_hline(yintercept = 0.5, linetype = "dashed", col = "#F8766D", size = 0.75) +
  geom_hline(yintercept = 0.7, linetype = "dashed", col = "#00BA38", size = 0.75) +
  geom_hline(yintercept = 0.9, linetype = "dashed", col = "#619CFF", size = 0.75) +
  scale_x_discrete(limits=c("1", "10", "20", "30")) +
  ggtitle("Simulation Results for Eccentricity") +
  theme_light() +
  theme(legend.position = "bottom") +
  scale_fill_manual(name="True Eccentricity of Sim. Data:", values=c("#F8766D", "#00BA38", "#619CFF")) +
  ylim(0,1)

p1 <- ggplot(data = all_sim_resultse0.5e0.7e0.9, mapping = aes(x = Sigma, y = Eccentricity, fill = TrueEccentricity)) + 
  geom_boxplot() + 
  geom_hline(yintercept = 0.5, linetype = "dashed", col = "#F8766D", size = 0.75) +
  geom_hline(yintercept = 0.7, linetype = "dashed", col = "#00BA38", size = 0.75) +
  geom_hline(yintercept = 0.9, linetype = "dashed", col = "#619CFF", size = 0.75) +
  scale_x_discrete(limits=c("1", "10", "20", "30")) +
  ggtitle("Simulation Results for Eccentricity") +
  theme_light() +
  theme(legend.position = "none") +
  scale_fill_manual(name="True Eccentricity of Sim. Data:", values=c("#F8766D", "#00BA38", "#619CFF")) +
  ylim(0,1)

p2 <- ggplot(data = all_sim_resultse0.5e0.7e0.9, mapping = aes(x = Sigma, y = Max_Min_Ratio, fill = TrueEccentricity)) + 
  geom_boxplot() + 
  geom_hline(yintercept = (2*sqrt(3))/3, linetype = "dashed", col = "#F8766D", size = 0.75) +
  geom_hline(yintercept = (10*sqrt(51))/51, linetype = "dashed", col = "#00BA38", size = 0.75) +
  geom_hline(yintercept = (10*sqrt(19))/19, linetype = "dashed", col = "#619CFF", size = 0.75) +
  scale_y_log10() +
  scale_x_discrete(limits=c("1", "10", "20", "30")) +
  ggtitle("Simulation Results for AR") +
  ylab("Anisotropic Ratio") +
  theme_light() +
  theme(legend.position = "none") + 
  annotation_logticks(sides = "l")  

shared_legend <- extract_legend(p1_legend)

grid.arrange(arrangeGrob(p1, p2, ncol = 2), shared_legend, nrow = 2, heights = c(10,1))

g <- arrangeGrob(p1, p2, ncol = 2)

ggsave(file = "Figures/figure5.eps", g, units = "px", width = 2586, height = 1926)
















## Worst example of each measure

# Worst Ratio Example
e <- 0.7
semi_min_axis <- 160 # Let's use the semi min axis length of 160 (avg from reviscometer data)
semi_maj_axis <- semi_min_axis/sqrt(1 - e^2)

sim_results20e0.7[sim_results20e0.7$Max_Min_Ratio == max(sim_results20e0.7$Max_Min_Ratio),]
worst_ratio_polar <- as.numeric(sim_results20e0.7[sim_results20e0.7$Max_Min_Ratio == max(sim_results20e0.7$Max_Min_Ratio),16:51])

centre_x_fit <- sim_results20e0.7[sim_results20e0.7$Max_Min_Ratio == max(sim_results20e0.7$Max_Min_Ratio),"Centre_x"]
centre_y_fit <- sim_results20e0.7[sim_results20e0.7$Max_Min_Ratio == max(sim_results20e0.7$Max_Min_Ratio),"Centre_y"]
semi_maj_axis_fit <- sim_results20e0.7[sim_results20e0.7$Max_Min_Ratio == max(sim_results20e0.7$Max_Min_Ratio),"Semi-major_Axis"]
semi_min_axis_fit <- sim_results20e0.7[sim_results20e0.7$Max_Min_Ratio == max(sim_results20e0.7$Max_Min_Ratio),"Semi-minor_Axis"]
angle_fit <- sim_results20e0.7[sim_results20e0.7$Max_Min_Ratio == max(sim_results20e0.7$Max_Min_Ratio),"Angle_Rads"]
worst_ecc_polar <- as.numeric(sim_results20e0.7[sim_results20e0.7$Eccentricity == min(sim_results20e0.7$Eccentricity),16:51])


Angle <- seq(0,350,10)

x <- worst_ratio_polar*cos(Angle*(pi/180)) # Convert Measurements from Polar to Cartesian
y <- worst_ratio_polar*sin(Angle*(pi/180)) # Convert Measurements from Polar to Cartesian

worst_ratio_cart <- data.frame(x = x, y = y)
worst_ratio_cart$dist_to_origin <- sqrt(worst_ratio_cart$x^2 + worst_ratio_cart$y^2)
worst_ratio_cart$min_max <- 0
worst_ratio_cart[worst_ratio_cart$dist_to_origin %in% c(max(worst_ratio_cart$dist_to_origin), 
                                                        min(worst_ratio_cart$dist_to_origin)),"min_max"] <- 1

ggplot(data = worst_ratio_cart, mapping = aes(x = x, y = y)) + 
  geom_point(size = 3) + 
  geom_point(data = worst_ratio_cart[worst_ratio_cart$min_max == 1,], mapping = aes(x = x, y = y), size = 6, col = "blue") +
  geom_ellipse(aes(x0 = 0, y0 = 0, a = semi_maj_axis, b = semi_min_axis, angle = 0), color = "red", lwd = 1, lty = "dashed") +
  geom_ellipse(aes(x0 = centre_x_fit, y0 = centre_y_fit, a = semi_maj_axis_fit, b = semi_min_axis_fit, angle = angle_fit), color = "blue", lwd = 1) +
  geom_hline(yintercept = 0, color = "black", size = 1) +
  geom_vline(xintercept = 0, color = "black", size = 1) +
  ggtitle("Poor AR Performance, e = 0.7, sigma = 20") + 
  theme_light()

ggsave(file = "Figures/figureA1b-noannotations.eps", units = "px", width = 1714, height = 1250)


# Worst Ellipse Example
e <- 0.7
semi_min_axis <- 160 # Let's use the semi min axis length of 160 (avg from reviscometer data)
semi_maj_axis <- semi_min_axis/sqrt(1 - e^2)

sim_results20e0.7[sim_results20e0.7$Eccentricity == min(sim_results20e0.7$Eccentricity),]
centre_x_fit <- sim_results20e0.7[sim_results20e0.7$Eccentricity == min(sim_results20e0.7$Eccentricity),"Centre_x"]
centre_y_fit <- sim_results20e0.7[sim_results20e0.7$Eccentricity == min(sim_results20e0.7$Eccentricity),"Centre_y"]
semi_maj_axis_fit <- sim_results20e0.7[sim_results20e0.7$Eccentricity == min(sim_results20e0.7$Eccentricity),"Semi-major_Axis"]
semi_min_axis_fit <- sim_results20e0.7[sim_results20e0.7$Eccentricity == min(sim_results20e0.7$Eccentricity),"Semi-minor_Axis"]
angle_fit <- sim_results20e0.7[sim_results20e0.7$Eccentricity == min(sim_results20e0.7$Eccentricity),"Angle_Rads"]
worst_ecc_polar <- as.numeric(sim_results20e0.7[sim_results20e0.7$Eccentricity == min(sim_results20e0.7$Eccentricity),16:51])

Angle <- seq(0,350,10)

x <- worst_ecc_polar*cos(Angle*(pi/180)) # Convert Measurements from Polar to Cartesian
y <- worst_ecc_polar*sin(Angle*(pi/180)) # Convert Measurements from Polar to Cartesian

worst_ecc_cart <- data.frame(x = x, y = y)
worst_ecc_cart$dist_to_origin <- sqrt(worst_ecc_cart$x^2 + worst_ecc_cart$y^2)
worst_ecc_cart$min_max <- 0
worst_ecc_cart[worst_ecc_cart$dist_to_origin %in% c(max(worst_ecc_cart$dist_to_origin),
                                                    min(worst_ecc_cart$dist_to_origin)),"min_max"] <- 1

ggplot(data = worst_ecc_cart, mapping = aes(x = x, y = y)) + 
  geom_point(size = 3) + 
  geom_point(data = worst_ecc_cart[worst_ecc_cart$min_max == 1,], mapping = aes(x = x, y = y), size = 6, col = "blue") +
  geom_ellipse(aes(x0 = 0, y0 = 0, a = semi_maj_axis, b = semi_min_axis, angle = 0), color = "red", lwd = 1, lty = "dashed") +
  geom_ellipse(aes(x0 = centre_x_fit, y0 = centre_y_fit, a = semi_maj_axis_fit, b = semi_min_axis_fit, angle = angle_fit), color = "blue", lwd = 1) +
  geom_hline(yintercept = 0, color = "black", size = 1) +
  geom_vline(xintercept = 0, color = "black", size = 1) +
  ggtitle("Poor Eccentricity Performance, e = 0.7, sigma = 20") + 
  theme_light()

ggsave(file = "Figures/figureA1a-noannotations.eps", units = "px", width = 1714, height = 1250)









