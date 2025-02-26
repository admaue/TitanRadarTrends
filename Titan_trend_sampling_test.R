### Testing subsamples of data from all mapped radar-bright fluvial features
### by A. Maue

# loading libraries
library(xlsx)
library(readxl)
library(ggplot2)
library(dplyr)

# read the combined data
alldata <- read_excel("E:/Titan/R/downstreammeasurements_onetable.xlsx")
z_km <- alldata$z_m/1000

# log-transformed backscatter (decibels)
s0_db <- 10 * log10(alldata$s0_mean)

# basic scatter plot with linear regression
dev.new(width=7.5, height=6, noRStudioGD = TRUE)
plot(z_km, s0_db, 
     ylab="Mean Radar Backscatter (dB)", xlab="Elevation (km)",
     cex.lab=1.5, cex.axis=1.5, pch=1)

# weighted linear model
wt <- 1 / alldata$s0_sd
lm1_wt <- lm(s0_db ~ z_km, weights=wt)
summary(lm1_wt)
abline(lm1_wt, col="red", lwd=2)
eq_text <- bquote(sigma[0] ~ "(dB) =" ~ .(round(coef(lm1_wt)[2], 2)) * z + .(round(coef(lm1_wt)[1], 2)) ~ " (R"^2 ~ "=" ~ .(round(summary(lm1_wt)$r.squared, 2)) ~ ")")
usr <- par("usr")  # Get plot limits (xmin, xmax, ymin, ymax)
text(usr[1] + 0.02 * (usr[2] - usr[1]), usr[4] - 0.05 * (usr[4] - usr[3]), 
     labels = eq_text, col = "red", cex = 1.3, adj = 0)


## random sampling from full dataset
# read data
dat <- read_excel("E:/Titan/R/data_distance_estimate.xlsx") %>%
  na.omit() %>%
  arrange(d_est)

# plot backscatter as a function of estimated distance
plot(s0_mean ~ d_est, data = dat, 
     xlab = "Estimated Distance (km)", 
     ylab = "Mean Radar Backscatter (dB)", 
     pch = 1) 

# define sampling and iterations
sample_lengths <- c(50, 100, 200, 300, 400, 500)
N_samples = 40
N_iterations = 100
results_list <- list() # list to store results for each sample length

# data frame to store basic stats for each sample length
samples_stats <- data.frame(
  length = numeric(length(sample_lengths)),
  m_mean = numeric(length(sample_lengths)),
  m_sd = numeric(length(sample_lengths))
)

# loop over each sample length
for (sample_length in sample_lengths) {
  
  # preallocate results dataframe
  samples_trends <- data.frame(
    m = numeric(N_iterations),
    b = numeric(N_iterations),
    r2 = numeric(N_iterations)
  )
  
  # sampling loop
  for (i in 1:N_iterations) {
    # select random start point within valid range
    section_start <- runif(1, min = min(dat$d_est), max = max(dat$d_est) - sample_length)
    
    # find closest indices for start and end
    section_start_i <- which.min(abs(dat$d_est - section_start))
    section_end_i <- which.min(abs(dat$d_est - (section_start + sample_length)))
    
    # extract section
    section <- dat[section_start_i:section_end_i, ]
    
    # proceed only if we have enough samples
    if (nrow(section) > N_samples) {
      samples <- section[sample(nrow(section), N_samples), ]
      
      # plot sampled data
      plot(s0_mean ~ d_est, data = samples)
      
      # fit linear model
      samples_model <- lm(s0_mean ~ d_est, data = samples)
      
      # store results
      samples_trends$m[i] <- coef(samples_model)[2]
      samples_trends$b[i] <- coef(samples_model)[1]
      samples_trends$r2[i] <- summary(samples_model)$r.squared
    } else {
      next  # skip this iteration if not enough data points
    }
  }
  
  # store results
  results_list[[paste(sample_length, "km")]] <- samples_trends
  
  # store m mean and SD for each sampling length
  m_mean <- mean(samples_trends$m, na.rm = TRUE)
  m_sd <- sd(samples_trends$m, na.rm = TRUE)
  row_index <- which(samples_stats$length == 0)[1]
  samples_stats[row_index, ] <- c(sample_length, m_mean, m_sd)
}

# write all results to one xlsx
write.xlsx(results_list, file = "E:/Titan/R/trend_sampling_test.xlsx")


## boxplot of results for different sampling lengths
dev.new(width = 3.2, height = 5, noRStudioGD = TRUE)
m_vals <- unlist(lapply(results_list, function(x) x$m))
l_vals <- rep(sample_lengths, each = N_iterations)
par(mar = c(4.2, 4.2, .5, .5))
boxplot(m_vals ~ factor(l_vals), 
        xlab = "Window Size (km)", 
        ylab = "Radar Backscatter Trend (1/km)", 
        cex.axis = 1.5, 
        cex.lab = 1.5, 
        ylim = c(-0.013, 0.013), 
        xaxt = "n")
axis(1, at = 1:length(sample_lengths), labels = sample_lengths, cex.axis = 1.3, las = 3)
abline(0, 0)
abline(-0.001227, 0, col = 'red')

