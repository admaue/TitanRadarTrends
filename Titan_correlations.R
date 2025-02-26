### Correlating variables measured in Titan's radar-bright fluvial features
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
library(stats)
library(purrr)

# define variables
dat <- read_excel("E:/Titan/R/data_table_featureproperties_r.xlsx")
latitude <- dat$clat
length <- dat$clength_km
sinuosity <- dat$sinuosity
s0_trend <- dat$s0_slope 
s0_r2 <- dat$s0_r2
s0_mean <- dat$s0_mean
s0_sd <- dat$s0_sd
s0_model_se <- dat$s0_model_se
elevation <- dat$z_mean_m
slope <- dat$z_slope 
width <- dat$width_mean_km
width_trend <- dat$w_slope  
downstream_CL <- dat$downstream_confidence
downstream_CLabc <- dat$CL_ABC

# derive additional variables
abs_lat <- abs(latitude)
abs_slope <- abs(slope)
abs_width_trend <- abs(width_trend)
abs_s0_trend <- abs(s0_trend)
log_length <- log10(length)
log_abs_s0_trend <- log(abs(s0_trend))
log_abs_slope <- log(abs(slope))
log_abs_width_trend <- log(abs(width_trend))

panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor) {
  usr <- par("usr")
  on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- abs(cor(x, y, use = "complete.obs"))
  txt <- format(c(r, 0.123456789), digits = digits)[1]
  txt <- paste(prefix, txt, sep = "")
  if (missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
  text(0.5, 0.5, txt, cex =  cex.cor * (1 + r) / 2)
}


## plotting variables
# scatter plot matrix with basic variables
pairs(cbind(s0_trend,s0_mean,abs_lat,length,sinuosity,elevation,slope,width,width_trend), 
      lower.panel=panel.lm,upper.panel = panel.cor)

# scatter plot matrix with derived variables
pairs(cbind(log_abs_s0_trend,s0_mean,abs_lat,length,sinuosity,elevation,abs_slope,width,log_abs_width_trend), 
      lower.panel=panel.lm,upper.panel = panel.cor)

# set up boxplot
dev.new(height = 4.3, width= 7.2, noRStudioGD = TRUE)
par(mfrow = c(1, 3), mar = c(4.5, 5, 2, 1))  # 3 rows, adjusted margins
plot_box_strip <- function(x, y, xlab, ylab, ylim) {
  boxplot(y ~ x, 
          xlab = xlab, 
          ylab = ylab, 
          cex.axis = 1.5, 
          cex.lab = 1, 
          ylim = ylim, 
          outline = FALSE)
  stripchart(y ~ x, 
             method = "jitter", 
             pch = 1, 
             cex = 0.75, 
             col = 'black', 
             vertical = TRUE, 
             add = TRUE)
  abline(h = 0, col = "red", lty = 2)  # Horizontal reference line
}

# generate the three plots
plot_box_strip(downstream_CL, s0_trend, 
               xlab = "Flow Direction Confidence Level", 
               ylab = "Radar Backscatter Trend (1/km)", 
               ylim = c(-0.013, 0.013))
plot_box_strip(downstream_CLabc, s0_trend, 
               xlab = "Confidence Level from A, B, and C", 
               ylab = "Radar Backscatter Trend (1/km)", 
               ylim = c(-0.013, 0.013))
plot_box_strip(downstream_CLabc, slope, 
               xlab = "Confidence Level from A, B, and C", 
               ylab = "Slope (m/km)", 
               ylim = c(-3.6, 3.6))


## multivariate linear regression models
# remove missing values
dat2 <- na.omit(data.frame(latitude, length, sinuosity, s0_trend, s0_mean, s0_sd, 
                           slope, elevation, width, width_trend, downstream_CL, 
                           downstream_CLabc, s0_r2, s0_model_se))

# transform variables
dat2$slope <- dat2$slope / 1000  # Convert slope to km scale
dat2$log_length <- log10(dat2$length)
dat2$abs_s0_trend <- abs(dat2$s0_trend)
dat2$abs_slope <- abs(dat2$slope)
dat2$abs_latitude <- abs(dat2$latitude)
dat2$abs_width_trend <- abs(dat2$width_trend)
dat2$log_abs_width_trend <- log10(dat2$abs_width_trend)

# stepwise regression model
fit_stepwise_model <- function(response, predictors, weights = NULL) {
  null_model <- lm(as.formula(paste(response, "~ 1")), data = dat2, weights = weights)
  full_model <- lm(as.formula(paste(response, "~", paste(predictors, collapse = " + "))), 
                   data = dat2, weights = weights)
  best_model <- step(null_model, scope = list(upper = full_model), direction = "forward", test = "F")
  return(best_model)
}

# define predictor variables for different models
predictors_full <- c("abs_latitude", "latitude", "log_length", "s0_mean", "elevation", 
                     "abs_slope", "slope", "width", "abs_width_trend", "log_abs_width_trend",
                     "abs_latitude:elevation", "abs_width_trend:elevation", "log_abs_width_trend:elevation",
                     "width:abs_width_trend", "s0_mean:elevation", "width:elevation")
predictors_some <- c("latitude", "log_length", "s0_mean", "elevation", 
                     "abs_slope", "width", "log_abs_width_trend")

# fit models with weights (1/s0_model_se)
models0.best <- fit_stepwise_model("s0_trend", predictors_full, weights = 1 / dat2$s0_model_se)
summary(models0.best)
models0.best_abs <- fit_stepwise_model("abs_s0_trend", predictors_some, weights = 1 / dat2$s0_model_se)
summary(models0.best_abs)

# fit models without weights
models0.best_no_wt <- fit_stepwise_model("s0_trend", predictors_full)
summary(models0.best_no_wt)


## relationship between s0 trend and feature slope for various subgroups
# plot slope vs. radar backscatter trend
plot_slope_vs_s0 <- function(sheet_name, equation_text) {
  slopedata <- read_excel("E:/Titan/R/data_slope_relationship.xlsx", sheet = sheet_name)
  
  # extract and transform variables
  slope <- slopedata$z_slope / 1000
  s0_trend <- slopedata$s0_slope
  s0_model_se <- slopedata$s0_model_se
  
  # plot the data
  dev.new(width = 5, height = 5, noRStudioGD = TRUE)
  plot(s0_trend ~ slope, xlab = "Slope (m/m)", ylab = "Radar Backscatter Trend (1/km)",
       cex.lab = 1.3, cex.axis = 1.2, cex = 1.5,
       ylim = c(-0.015, 0.015),
       xlim = range(slope, na.rm = TRUE) * 1.1) 
  arrows(x0 = slope, y0 = s0_trend - ((1 / s0_model_se) / 500),
         x1 = slope, y1 = s0_trend + ((1 / s0_model_se) / 500),
         code = 3, angle = 90, length = 0.05, col = "grey")
  
  # linear model with weighted regression
  models0.slope <- lm(s0_trend ~ slope, weights = 1 / s0_model_se)
  abline(models0.slope, col = "red", lwd = 2)
  usr <- par("usr")
  text(usr[1] + 0.02 * (usr[2] - usr[1]), -0.0135, equation_text, col = "red", cex = 1.,adj = c(0, 0.5))
  
  # model summary
  print(summary(models0.slope))
}

# define regression equations
sheets <- list(
  "all"       = expression("s0 (dB) = -1.02 slope - 0.001 (R² = 0.18)"),
  "R2_02"     = expression("s0 (dB) = -1.80 slope - 0.003 (R² = 0.31)"),
  "R2_02pos"  = expression("s0 (dB) = -1.36 slope - 0.002 (R² = 0.05)"),
  "R2_02flip" = expression("s0 (dB) = -2.47 slope - 0.003 (R² = 0.13)")
)

# generate plots
for (sheet in names(sheets)) {
  plot_slope_vs_s0(sheet, sheets[[sheet]])
}


## write spreadsheet of all dataframes
# define file path and read sheets
file_path <- "E:/Titan/R/data_table_downstreammeasurements.xlsx"
sheet_names <- excel_sheets(file_path)
alldata <- map_dfr(sheet_names, ~ read_excel(file_path, sheet = .x) %>% mutate(sheet_name = .x))

# save the combined dataframe to a new excel file
write_xlsx(alldata, "E:/Titan/R/downstreammeasurements_onetable.xlsx")

