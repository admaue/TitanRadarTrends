### Analyzing backscatter trends in Titan's radar-bright fluvial features
### by A. Maue

# loading libraries
library(readxl)      
library(openxlsx)    
library(writexl)     
library(dplyr)     
library(corrplot)   
library(minpack.lm) 
library(nlstools) 
library(finalfit) 
library(interactions) 
library(tidyverse) 
library(GGally)    

# define variables
sheets <- getSheetNames("E:/Titan/R/data_table_downstreammeasurements.xlsx") 
num_sheets <- length(sheets)  

# preallocate vectors
w_m <- numeric(num_sheets)  
w_b <- numeric(num_sheets)
w_r2 <- numeric(num_sheets)
w_mean <- numeric(num_sheets)
w_sd <- numeric(num_sheets)
z_m <- numeric(num_sheets)
z_b <- numeric(num_sheets)
z_r2 <- numeric(num_sheets)
z_mean <- numeric(num_sheets)
z_sd <- numeric(num_sheets)
s0_m <- numeric(num_sheets)
s0_b <- numeric(num_sheets)
s0_r2 <- numeric(num_sheets)
s0_model_se <- numeric(num_sheets)
s0_mean <- numeric(num_sheets)
s0_sd <- numeric(num_sheets)
clat_mean <- numeric(num_sheets)
clon_mean <- numeric(num_sheets)

# read data, perform regression fits, and store results
for (i in seq_along(sheets)) {
  dat <- read_excel("E:/Titan/R/data_table_downstreammeasurements.xlsx", sheet = sheets[i])
  
  # linear model for width vs. distance
  w_model <- lm(width_km ~ x_km, data = dat)
  w_m[i] <- coef(w_model)[2]  # Slope
  w_b[i] <- coef(w_model)[1]  # Intercept
  w_r2[i] <- summary(w_model)$r.squared  # R² value
  w_mean[i] <- mean(dat$width_km, na.rm = TRUE)
  w_sd[i] <- sd(dat$width_km, na.rm = TRUE)
  
  # linear model for elevation vs. distance
  z_model <- lm(z_m ~ x_km, data = dat)
  z_m[i] <- coef(z_model)[2]
  z_b[i] <- coef(z_model)[1]
  z_r2[i] <- summary(z_model)$r.squared
  z_mean[i] <- mean(dat$z_m, na.rm = TRUE)
  z_sd[i] <- sd(dat$z_m, na.rm = TRUE)
  
  # linear model for s0_mean vs. distance (weighted)
  s0_se <- dat$s0_sd / sqrt(dat$N_s0)  # Compute standard error
  s0_model <- lm(s0_mean ~ x_km, weights = 1/s0_se, data = dat)
  s0_m[i] <- coef(s0_model)[2]
  s0_b[i] <- coef(s0_model)[1]
  s0_r2[i] <- summary(s0_model)$r.squared
  s0_model_se[i] <- sqrt(deviance(s0_model) / df.residual(s0_model))
  s0_mean[i] <- mean(dat$s0_mean, na.rm = TRUE)
  s0_sd[i] <- sd(dat$s0_mean, na.rm = TRUE)
  
  # mean latitude and longitude
  clat_mean[i] <- mean(dat$clat, na.rm = TRUE)
  clon_mean[i] <- mean(dat$clon, na.rm = TRUE)
}

# store results storage as a data frame
fits <- data.frame(
  sheet = sheets,
  clat_mean = numeric(length(sheets)),
  clon_mean = numeric(length(sheets)),
  w_m = numeric(length(sheets)), w_b = numeric(length(sheets)), w_r2 = numeric(length(sheets)), 
  w_mean = numeric(length(sheets)), w_sd = numeric(length(sheets)), 
  z_m = numeric(length(sheets)), z_b = numeric(length(sheets)), z_r2 = numeric(length(sheets)), 
  z_mean = numeric(length(sheets)), z_sd = numeric(length(sheets)), 
  s0_m = numeric(length(sheets)), s0_b = numeric(length(sheets)), s0_r2 = numeric(length(sheets)), 
  s0_model_se = numeric(length(sheets)), s0_mean = numeric(length(sheets)), s0_sd = numeric(length(sheets))
)

# loop through sheets and populate the data frame
for (i in seq_along(sheets)) {
  dat <- read_excel("E:/Titan/R/data_table_downstreammeasurements.xlsx", sheet = sheets[i])
  
  # linear models
  w_model <- lm(width_km ~ x_km, data = dat)
  z_model <- lm(z_m ~ x_km, data = dat)
  s0_se <- dat$s0_sd / sqrt(dat$N_s0)
  s0_model <- lm(s0_mean ~ x_km, weights = 1/s0_se, data = dat)
  
  # store results
  fits[i, ] <- list(
    sheets[i], mean(dat$clat, na.rm = TRUE), mean(dat$clon, na.rm = TRUE),
    coef(w_model)[2], coef(w_model)[1], summary(w_model)$r.squared,
    mean(dat$width_km, na.rm = TRUE), sd(dat$width_km, na.rm = TRUE),
    coef(z_model)[2], coef(z_model)[1], summary(z_model)$r.squared,
    mean(dat$z_m, na.rm = TRUE), sd(dat$z_m, na.rm = TRUE),
    coef(s0_model)[2], coef(s0_model)[1], summary(s0_model)$r.squared,
    sqrt(deviance(s0_model) / df.residual(s0_model)), mean(dat$s0_mean, na.rm = TRUE), sd(dat$s0_mean, na.rm = TRUE)
  )
}

# export results
write_xlsx(fits, "E:/Titan/R/fits_downstreammeasurements.xlsx")

# generate and save plots
plot_dir <- "E:/Titan/R/plots/feature_plots/"
for (i in seq_along(sheets)) {
  dat <- read_excel("E:/Titan/R/data_table_downstreammeasurements.xlsx", sheet = sheets[i])
  
  # compute standard error and fit linear model
  s0_se <- dat$s0_sd / sqrt(dat$N_s0)
  s0_model <- lm(s0_mean ~ x_km, weights = 1/s0_se, data = dat)
  s0_slope <- coef(s0_model)[2]
  
  # dynamic y-axis limits
  y_min <- min(dat$s0_mean - (0.5 * mean(dat$s0_sd, na.rm = TRUE)), na.rm = TRUE)
  y_max <- max(dat$s0_mean + (0.5 * mean(dat$s0_sd, na.rm = TRUE)), na.rm = TRUE)
  
  # generate plot
  plot_file <- paste0(plot_dir, sheets[i], ".png")
  png(file = plot_file, width = 450, height = 450)
  plot(
    dat$x_km, dat$s0_mean, main = sheets[i],
    ylab = "NRCS", xlab = "x (km)", cex.lab = 2.5, cex.axis = 2, cex.main = 3, cex = 2,
    mgp = c(2.4, 0.8, 0), ylim = c(y_min, y_max)
  )
  abline(s0_model, col = ifelse(s0_slope > 0, "blue", "red"), lwd = 2)
  text_x <- max(dat$x_km, na.rm = TRUE) * 0.5
  text_y <- y_min
  text_label <- paste0(formatC(s0_slope, digits = 3), " 1/km")
  text(text_x, text_y, text_label, col = ifelse(s0_slope > 0, "blue", "red"), cex = 2)
  dev.off()
}
