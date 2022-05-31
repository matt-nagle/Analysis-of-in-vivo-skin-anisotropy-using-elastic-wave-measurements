library(conicfit)
library(ggplot2)


## Import Averaged Datasets from the Reviscometer

setwd("C:/Users/mattn/OneDrive - University College Dublin/Documents/00. Research/99. Publications/Analysis of in-vivo skin anisotropy using elastic wave measurements")
revis_data_orig_units <- read.csv("Reviscometer Avg Data - RRT Units.csv", header = TRUE)
revis_data_stretch_orig_units <- read.csv("Reviscometer Avg Data - Stretched - RRT Units.csv", header = TRUE)

revis_data_orig_units$Gender <- as.factor(revis_data_orig_units$Gender)
revis_data_stretch_orig_units$Gender <- as.factor(revis_data_stretch_orig_units$Gender)



## Create Dataframes for the parameters of the ellipse fit to the averaged Reviscometer Dataset

fit_ellipse_params <- data.frame(matrix(ncol = 16, nrow = nrow(revis_data_orig_units)))
colnames(fit_ellipse_params) <- c("Subject_ID", "Gender", "Age", "Centre_x", "Centre_y", "Axis_1", "Axis_2", "Semi-major_Axis", "Semi-minor_Axis", "Angle_Rads", "Angle_Deg", "Angle_to_semi-maj_Rads", "Angle_to_semi-maj_Deg", "Eccentricity", "RSS", "RSS_Area_Norm") 

fit_ellipse_params$Subject_ID <- revis_data_orig_units$Subject_Id
fit_ellipse_params$Gender <- revis_data_orig_units$Gender
fit_ellipse_params$Age <- revis_data_orig_units$Age


fit_stretch_ellipse_params <- data.frame(matrix(ncol = 16, nrow = nrow(revis_data_stretch_orig_units)))
colnames(fit_stretch_ellipse_params) <- c("Subject_ID", "Gender", "Age", "Centre_x", "Centre_y", "Axis_1", "Axis_2", "Semi-major_Axis", "Semi-minor_Axis", "Angle_Rads", "Angle_Deg", "Angle_to_semi-maj_Rads", "Angle_to_semi-maj_Deg", "Eccentricity", "RSS", "RSS_Area_Norm") 

fit_stretch_ellipse_params$Subject_ID <- revis_data_stretch_orig_units$Subject_Id
fit_stretch_ellipse_params$Gender <- revis_data_stretch_orig_units$Gender
fit_stretch_ellipse_params$Age <- revis_data_stretch_orig_units$Age




#####========= Define a function to plot an ellipse for a given Subject =========#####

plot_ellipse <- function(subject, plot_stretched = TRUE, plot_maj_min_axes = TRUE, 
                         plot_revis_data = TRUE, display_ellipse = TRUE, display_legend = TRUE){
  Angle <- seq(0,350,10)
  
  ## Natural Configuration ##
  # Convert Measurements from Polar to Cartesian
  x <- as.numeric(revis_data_orig_units[revis_data_orig_units$Subject_Id == subject, 4:39]*cos(Angle*(pi/180)))
  y <- as.numeric(revis_data_orig_units[revis_data_orig_units$Subject_Id == subject, 4:39]*sin(Angle*(pi/180)))
  
  subject_RRT_cart <- data.frame(x = x, y = y)
  
  
  # Isolate ellipse params for the chosen subject
  ellipse_params <- fit_ellipse_params[fit_ellipse_params$Subject_ID == subject, 
                                       c("Centre_x", "Centre_y", "Axis_1", "Axis_2", "Angle_Rads", "Angle_to_semi-maj_Rads")]
  
  # Calculates xy coords for fit ellipse
  ellipse_xy <- calculateEllipse(ellipse_params[[1]],ellipse_params[[2]],ellipse_params[[3]],
                                 ellipse_params[[4]],(180/pi)*ellipse_params[[5]], steps = 100)
  
  # Find the x and y axis limits for the plot
  max_axis <- max(c(ellipse_params[[3]],ellipse_params[[4]]))
  min_axis <- min(c(ellipse_params[[3]],ellipse_params[[4]]))
  
  global_max_axis <- max_axis
  
  
  ## Stretched Configuration ##
  if(plot_stretched){
    # Convert Measurements from Polar to Cartesian
    x_stretch <- as.numeric(revis_data_stretch_orig_units[revis_data_stretch_orig_units$Subject_Id == subject, 4:39]*cos(Angle*(pi/180)))
    y_stretch <- as.numeric(revis_data_stretch_orig_units[revis_data_stretch_orig_units$Subject_Id == subject, 4:39]*sin(Angle*(pi/180)))
    
    subject_RRT_cart_stretch <- data.frame(x = x_stretch, y = y_stretch)
    
    
    ellipse_params_stretch <- fit_stretch_ellipse_params[fit_stretch_ellipse_params$Subject_ID == subject, 
                                                 c("Centre_x", "Centre_y", "Axis_1", "Axis_2", "Angle_Rads", "Angle_to_semi-maj_Rads")]
    
    ellipse_xy_stretch <- calculateEllipse(ellipse_params_stretch[[1]],ellipse_params_stretch[[2]],ellipse_params_stretch[[3]],
                                           ellipse_params_stretch[[4]],(180/pi)*ellipse_params_stretch[[5]], steps = 100) #calculates xy coords for fit ellipse
    
    max_axis_stretch <- max(c(ellipse_params_stretch[[3]],ellipse_params_stretch[[4]]))
    min_axis_stretch <- min(c(ellipse_params_stretch[[3]],ellipse_params_stretch[[4]]))
    
    global_max_axis <- max(c(max_axis, max_axis_stretch))
  }
  
  
  if(plot_revis_data){
    # Plot the Reviscometer Data
    plot(subject_RRT_cart,
         xlim=c(ellipse_params[[1]] - global_max_axis*1.05,ellipse_params[[1]] + global_max_axis*1.05),
         ylim=c(ellipse_params[[2]] - global_max_axis*1.05,ellipse_params[[2]] + global_max_axis*1.05),
         main=paste("Subject ",  subject,
                    ", Age = ", revis_data_orig_units[revis_data_orig_units$Subject_Id == subject, "Age"],
                    ", Gender = ", revis_data_orig_units[revis_data_orig_units$Subject_Id == subject, "Gender"],
                    sep=""),
         xlab = "x coordinate (RRT)",
         ylab = "y coordinate (RRT)",
         pch = 20,
         cex = 1.5,
         asp = 1,
         col = "black");par(new=TRUE)

    # Add stretched Reviscometer Data
    if(plot_stretched){
      points(subject_RRT_cart_stretch,
             xlab = "",
             ylab = "",
             pch = 20,
             cex = 1.5,
             asp = 1,
             col = "red");par(new=TRUE)
    }
  }

  if(display_ellipse){
    # Plot the Fit Ellipse
    plot(ellipse_xy[,1],
         ellipse_xy[,2],
         xlim=c(ellipse_params[[1]] - global_max_axis*1.05,ellipse_params[[1]] + global_max_axis*1.05),
         ylim=c(ellipse_params[[2]] - global_max_axis*1.05,ellipse_params[[2]] + global_max_axis*1.05),
         xlab = "",
         ylab = "",
         type='l',
         col='black',
         lwd=2,
         asp = 1);par(new=TRUE)
    
    # Add stretched Reviscometer Ellipse
    if(plot_stretched){
      plot(ellipse_xy_stretch[,1],
           ellipse_xy_stretch[,2],
           xlim=c(ellipse_params[[1]] - global_max_axis*1.05,ellipse_params[[1]] + global_max_axis*1.05),
           ylim=c(ellipse_params[[2]] - global_max_axis*1.05,ellipse_params[[2]] + global_max_axis*1.05),
           xlab = "",
           ylab = "",
           type='l',
           col='red',
           lwd=2,
           asp = 1);par(new=TRUE)
    }
  }


  if(plot_maj_min_axes){
    # Plot the Centre Point
    plot(ellipse_params[[1]], ellipse_params[[2]],
         xlim=c(ellipse_params[[1]] - global_max_axis*1.05,ellipse_params[[1]] + global_max_axis*1.05),
         ylim=c(ellipse_params[[2]] - global_max_axis*1.05,ellipse_params[[2]] + global_max_axis*1.05),
         xlab = "", ylab = "",
         pch=16,
         asp = 1);par(new=TRUE)
    
    # Add stretched Centre Point
    if(plot_stretched){
      plot(ellipse_params_stretch[[1]], ellipse_params_stretch[[2]],
           xlim=c(ellipse_params[[1]] - global_max_axis*1.05,ellipse_params[[1]] + global_max_axis*1.05),
           ylim=c(ellipse_params[[2]] - global_max_axis*1.05,ellipse_params[[2]] + global_max_axis*1.05),
           xlab = "", ylab = "",
           pch=16,
           asp = 1);par(new=TRUE)
    }
    
    # Plot the Semi-major and Semi-minor Axes
    segments(x0 = ellipse_params[[1]], y0 = ellipse_params[[2]],
             x1 = max_axis*cos(ellipse_params[[6]]) + ellipse_params[[1]],
             y1 = max_axis*sin(ellipse_params[[6]]) + ellipse_params[[2]],
             lwd = 2, col = "blue")
    segments(x0 = ellipse_params[[1]], y0 = ellipse_params[[2]],
             x1 = min_axis*cos(ellipse_params[[6]] + pi/2) + ellipse_params[[1]],
             y1 = min_axis*sin(ellipse_params[[6]] + pi/2) + ellipse_params[[2]],
             lwd = 2, col = "green")
    
    # Add stretched Semi-major and Semi-minor Axes
    if(plot_stretched){
      segments(x0 = ellipse_params_stretch[[1]], y0 = ellipse_params_stretch[[2]],
               x1 = max_axis_stretch*cos(ellipse_params_stretch[[6]]) + ellipse_params_stretch[[1]],
               y1 = max_axis_stretch*sin(ellipse_params_stretch[[6]]) + ellipse_params_stretch[[2]],
               lwd = 2, col = "blue")
      segments(x0 = ellipse_params_stretch[[1]], y0 = ellipse_params_stretch[[2]],
               x1 = min_axis_stretch*cos(ellipse_params_stretch[[6]] + pi/2) + ellipse_params_stretch[[1]],
               y1 = min_axis_stretch*sin(ellipse_params_stretch[[6]] + pi/2) + ellipse_params_stretch[[2]],
               lwd = 2, col = "green")
    }
  }

  
  # Plot the x and y axes
  abline(v = 0)
  abline(h = 0)
  
  if(display_legend){
    legend("topright",
           legend=c("Natural Config", "Stretched Config", "Semi-Maj", "Semi-Min"),
           col = c("black", "red", "blue", "green"), lwd=2)
  }
  
}




#####=============== Fit an Ellipse to Each Subject Measurement and Record the Ellipse Parameters ===============#####

for (subject in 1:78){
  Angle <- seq(0,350,10)
  
  ## Natural Configuration ##
  
  # Convert Measurements from Polar to Cartesian
  x <- as.numeric(revis_data_orig_units[revis_data_orig_units$Subject_Id == subject, 4:39]*cos(Angle*(pi/180)))
  y <- as.numeric(revis_data_orig_units[revis_data_orig_units$Subject_Id == subject, 4:39]*sin(Angle*(pi/180))) 
  
  subject_RRT_cart <- data.frame(x = x, y = y)
  
  
  # The following returns algebraic parameters for general conic equation 
  # f(x,y) = ax^2 + bxy + cy^2 + dx + ey + f = 0
  ellipse_fit <- EllipseDirectFit(subject_RRT_cart)
  # Converts algebraic parameters to geometric
  ellipse_fit_geo <- AtoG(ellipse_fit)$ParG 
  #calculates xy coords for fit ellipse
  ellipse_xy <- calculateEllipse(ellipse_fit_geo[1],ellipse_fit_geo[2],ellipse_fit_geo[3],
                                 ellipse_fit_geo[4],(180/pi)*ellipse_fit_geo[5], steps = 100) 
  
  
  # How well does the fit ellipse match the points?
  ellipse_res <- Residuals.ellipse(data.matrix(subject_RRT_cart), matrix(c(ellipse_fit_geo), ncol=1));par(new=FALSE)
  
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
  
  
  # Add natural config ellipse params to the dataframe
  fit_ellipse_params[fit_ellipse_params$Subject_ID == subject,4:16] <- c(ellipse_fit_geo[1],ellipse_fit_geo[2], # Centre x, Centre y
                                                                         ellipse_fit_geo[3],ellipse_fit_geo[4], # Axis 1, Axis 2
                                                                         max(c(ellipse_fit_geo[3],ellipse_fit_geo[4])), # Semi-maj Axis
                                                                         min(c(ellipse_fit_geo[3],ellipse_fit_geo[4])), # Semi-min Axis
                                                                         ellipse_fit_geo[5], # Angle (Rads)
                                                                         ellipse_fit_geo[5]*(180/pi), # Angle (Deg)
                                                                         angle_to_semi_maj, # Angle to the semi-major axis (Rads)
                                                                         angle_to_semi_maj*(180/pi), # Angle to the semi-major axis (Degs)
                                                                         sqrt(1 - (min(c(ellipse_fit_geo[3],ellipse_fit_geo[4]))/max(c(ellipse_fit_geo[3],ellipse_fit_geo[4])))^2), # Eccentricity
                                                                         ellipse_res$RSS, # Residual Sum of Squares
                                                                         ellipse_res$RSS/(pi*ellipse_fit_geo[3]*ellipse_fit_geo[4])) # Scaled RSS
  
  
  ## Stretched Configuration ##
  # Convert Measurements from Polar to Cartesian
  x <- as.numeric(revis_data_stretch_orig_units[revis_data_stretch_orig_units$Subject_Id == subject, 4:39]*cos(Angle*(pi/180))) 
  y <- as.numeric(revis_data_stretch_orig_units[revis_data_stretch_orig_units$Subject_Id == subject, 4:39]*sin(Angle*(pi/180)))
  
  subject_RRT_cart <- data.frame(x = x, y = y)
  
  
  
  # The following returns algebraic parameters for general conic equation 
  # f(x,y) = ax^2 + bxy + cy^2 + dx + ey + f = 0
  ellipse_fit <- EllipseDirectFit(subject_RRT_cart)
  # Converts algebraic parameters to geometric
  ellipse_fit_geo <- AtoG(ellipse_fit)$ParG 
  #calculates xy coords for fit ellipse
  ellipse_xy <- calculateEllipse(ellipse_fit_geo[1],ellipse_fit_geo[2],ellipse_fit_geo[3],
                                 ellipse_fit_geo[4],(180/pi)*ellipse_fit_geo[5], steps = 100) 
  
  
  # How well does the fit ellipse match the points?
  ellipse_res <- Residuals.ellipse(data.matrix(subject_RRT_cart), matrix(c(ellipse_fit_geo), ncol=1));par(new=FALSE)
  
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
  
  
  # Add stretched config ellipse params to the dataframe
  fit_stretch_ellipse_params[fit_stretch_ellipse_params$Subject_ID == subject,4:16] <- c(ellipse_fit_geo[1],ellipse_fit_geo[2], # Centre x, Centre y
                                                                         ellipse_fit_geo[3],ellipse_fit_geo[4], # Axis 1, Axis 2
                                                                         max(c(ellipse_fit_geo[3],ellipse_fit_geo[4])), # Semi-maj Axis
                                                                         min(c(ellipse_fit_geo[3],ellipse_fit_geo[4])), # Semi-min Axis
                                                                         ellipse_fit_geo[5], # Angle (Rads)
                                                                         ellipse_fit_geo[5]*(180/pi), # Angle (Deg)
                                                                         angle_to_semi_maj, # Angle to the semi-major axis (Rads)
                                                                         angle_to_semi_maj*(180/pi), # Angle to the semi-major axis (Degs)
                                                                         sqrt(1 - (min(c(ellipse_fit_geo[3],ellipse_fit_geo[4]))/max(c(ellipse_fit_geo[3],ellipse_fit_geo[4])))^2), # Eccentricity
                                                                         ellipse_res$RSS, # Residual Sum of Squares
                                                                         ellipse_res$RSS/(pi*ellipse_fit_geo[3]*ellipse_fit_geo[4])) # Scaled RSS
  
  # Plot each fit ellipse
  par(mfrow=c(1,1))
  plot_ellipse(subject, plot_stretched = TRUE, plot_maj_min_axes = FALSE, plot_revis_data = TRUE, display_ellipse = TRUE)
}



## Create a combined ellipse dataset and write to csv for use in model building

fit_ellipse_params$Area <- pi*fit_ellipse_params$`Semi-major_Axis`*fit_ellipse_params$`Semi-minor_Axis`
fit_stretch_ellipse_params$Area <- pi*fit_stretch_ellipse_params$`Semi-major_Axis`*fit_stretch_ellipse_params$`Semi-minor_Axis`

fit_ellipse_params$Configuration <- "Natural"
fit_stretch_ellipse_params$Configuration <- "Stretched"

all_fit_ellipse_params <- rbind(fit_ellipse_params, fit_stretch_ellipse_params)

write.csv(all_fit_ellipse_params, "Ellipse Params - Natural and Stretched Config.csv")






#####=============== Visualisations for Publication ===============#####

## Example of an ellipse fit

plot_ellipse(29, plot_stretched = FALSE, plot_maj_min_axes = FALSE, 
             plot_revis_data = TRUE, display_ellipse = FALSE, display_legend = FALSE)

setEPS()
postscript("Figures/figure2a.eps", width = 6.7134, height = 5)
plot_ellipse(29, plot_stretched = FALSE, plot_maj_min_axes = FALSE, 
             plot_revis_data = TRUE, display_ellipse = FALSE, display_legend = FALSE)
dev.off()


plot_ellipse(29, plot_stretched = FALSE, plot_maj_min_axes = TRUE, 
             plot_revis_data = TRUE, display_ellipse = TRUE, display_legend = FALSE)
legend("topright",
       legend=c("Natural Config", "Semi-Maj", "Semi-Min"),
       col = c("black", "blue", "green"), lwd=2)

setEPS()
postscript("Figures/figure2b-noannotations.eps", width = 6.7134, height = 5)
plot_ellipse(29, plot_stretched = FALSE, plot_maj_min_axes = TRUE, 
             plot_revis_data = TRUE, display_ellipse = TRUE, display_legend = FALSE)
legend("topright",
       legend=c("Natural Config", "Semi-Maj", "Semi-Min"),
       col = c("black", "blue", "green"), lwd=2)
dev.off()



## Example of Natural vs Stretched Ellipse fit

plot_ellipse(29, plot_stretched = TRUE, plot_maj_min_axes = TRUE, plot_revis_data = TRUE, display_ellipse = TRUE)

setEPS()
postscript("Figures/figure8.eps", width = 6.7134, height = 5)
plot_ellipse(29, plot_stretched = TRUE, plot_maj_min_axes = TRUE, plot_revis_data = TRUE, display_ellipse = TRUE)
dev.off()





## Eccentricity of Natural vs Stretched Ellipses

p <- ggplot(data = all_fit_ellipse_params, aes(x = Configuration, y = Eccentricity, group = Subject_ID))

p + geom_line() + 
  stat_smooth(aes(group = 1), color = "darkred", fill = "brown1", size = 2) + 
  stat_summary(aes(group = 1),geom = "point", fun = mean, shape = 19, size = 5, color = "darkred") +
  ylim(0,1) +
  scale_x_discrete(expand = c(0.1,0.1)) +
  ggtitle("Eccentricity vs Configuration") +
  theme_light()

ggsave(file = "Figures/figureB2c.eps", units = "px", width = 1724, height = 1258, device=cairo_ps)



## Semi-min Axis of Natural vs Stretched Ellipses

p <- ggplot(data = all_fit_ellipse_params, aes(x = Configuration, y = `Semi-minor_Axis`, group = Subject_ID))

p + geom_line() + 
  stat_smooth(aes(group = 1), color = "darkred", fill = "brown1", size = 2) + 
  stat_summary(aes(group = 1),geom = "point", fun = mean, shape = 19, size = 5, color = "darkred") +
  scale_x_discrete(expand = c(0.1,0.1)) +
  ggtitle("Semi-minor Axis vs Configuration") +
  ylab("Semi-minor Axis") +
  theme_light()

ggsave(file = "Figures/figureB2b.eps", units = "px", width = 1724, height = 1258, device=cairo_ps)




## Angle of Natural vs Stretched Ellipses

p <- ggplot(data = all_fit_ellipse_params, aes(x = Configuration, y = `Angle_to_semi-maj_Deg`, group = Subject_ID))

p + geom_line() + 
  stat_smooth(aes(group = 1), color = "darkred", fill = "brown1", size = 2) + 
  stat_summary(aes(group = 1),geom = "point", fun = mean, shape = 19, size = 5, color = "darkred") +
  scale_x_discrete(expand = c(0.1,0.1)) +
  ggtitle("Angle vs Configuration") +
  ylab("Angle") +
  theme_light()

ggsave(file = "Figures/figureB2a.eps", units = "px", width = 1724, height = 1258, device=cairo_ps)






## Log Model of Eccentricity vs Age

Ecc_Age_log <- lm(formula = Eccentricity ~ 1 + log10(Age),
                  data = fit_ellipse_params)

summary(Ecc_Age_log)


fit_log_model_points <- data.frame(cbind(revis_data_orig_units$Age, Ecc_Age_log$fitted.values))
colnames(fit_log_model_points) <- c("Age", "fitted.values")

ggplot(data = fit_ellipse_params, mapping = aes(x = Age, y = Eccentricity)) + 
  geom_point() + 
  geom_line(data = fit_log_model_points, mapping = aes(x = Age, y = fitted.values), col = "red", lwd = 1) +
  ggtitle("Eccentricity vs Age") + 
  theme_light()

ggsave(file = "Figures/figure7.eps", units = "px", width = 1724, height = 1284)