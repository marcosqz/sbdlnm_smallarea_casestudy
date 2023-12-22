################################################################################
# IMPORTANT: THE MORTALITY DATA provided are simulated FROM MODEL 4 IN THE
# PAPER (SB-DLNM WITH TIME-SERIES DESIGN) due to confidentiality restrictions of 
# the original data set. Original data may be requested to the data owner (The 
# Barcelona Public Health Agency) who should share them in similar terms than 
# those applying to this study. The supplied mortality, with exactly the same 
# structure than the original data set, allows reproducing the code provided so 
# that interested readers may run it as an example of use. However, WinBUGS 
# results from our original analyses are additionally supplied (input/result_paper
# folder) so that readers should be able to reproduce exactly our analyses after 
# those WinBUGS calls (model summaries and plots).
################################################################################

################################################################################
# In this R project we implement Bayesian and Spatial Bayesian Distributed 
# Lag Non-Linear Models (B-DLNM and SB-DLNM) for the case Study of short-term 
# associations between temperature and mortality in the city of Barcelona.
################################################################################

################################################################################
# CODE 3: PLOT SB-DLNMs
# Handling with output from SB-DLNMs: Construct and plot 3D surface associations, 
# overall cumulative associations, and mapping spatial patterns in R.
################################################################################

# Load libraries

library(dlnm)
library(leaflet)
library(sf)
library(lubridate)

# Load necessary files

load("input/daily_data.RData")
shapefile_bcn <- read_sf("input/shapefile_bcn.shp")
# WARNING: WE LOAD THE OUTPUT FROM MODEL 3 EXECUTED IN THE PAPER WITH THE 
# REAL MORTALITY DATA. IT CAN BE REPLACED FOR ANY OF THE OTHER OUTPUTS OBTAINED
# WITH REAL DATA OR OBTAINED WITH THE PREDICTED DATA FROM CODE 02_run_sbdlm.R
load("input/result_paper/final_simsmatrix_model3_spatial_casecrossover.RData")
load("output/dlnm_configuration.RData")

# DATA PREPARATION

# Subset data to only summer months of 2007 to 2016

data <- subset(data, month(date) %in% 6:9)
data <- subset(data, year(date) %in% 2007:2016)

# Define the percentiles of temperature to be calculated 
percentiles <- c(seq(0, 1, by = 0.1), 
                 seq(2, 98, by = 1), 
                 seq(99, 100, by = 0.1)) /100

# Create the temperature values used in the DLNM models 
list_temp <- lapply(1:dlnm_var$n_reg, function(i_reg) {
  
  label_reg <- formatC(i_reg, width = 2, flag = "0")
  temp <- subset(data, region == label_reg, 
                 select = c("temp", paste0("lag", 1:dlnm_var$max_lag)))
  temp_knots <- quantile(temp[["temp"]], dlnm_var$var_prc, na.rm = TRUE)
  temp_boundary <- range(temp, na.rm = TRUE)
  x_temp <- quantile(temp[["temp"]], percentiles, na.rm = TRUE)
  
  return(list(temp_knots = temp_knots, # knots of the exposure-response function  
              temp_boundary = temp_boundary, # range of temperatures
              x_temp = x_temp)) # temperatures in which risk are calculated
  
})

temp_knots <- sapply(list_temp, function(x) x[["temp_knots"]])
temp_boundary <- sapply(list_temp, function(x) x[["temp_boundary"]])
x_temp <- sapply(list_temp, function(x) x[["x_temp"]])
rm(list_temp)

# Create a list with the basis for exposure and the basis for lags 
# in each neighbourhood

basis_all <- lapply(1:dlnm_var$n_reg, function(i_reg) {
  
  f.temp_knots <- temp_knots[,i_reg]
  f.temp_boundary <- temp_boundary[,i_reg]
  f.x_temp <- x_temp[,i_reg]
  
  # basis temperatures
  Q <- onebasis(f.x_temp, fun = dlnm_var$var_fun, knots = f.temp_knots, 
                Boundary.knots = f.temp_boundary)
  
  # basis lags
  C <- onebasis(0:dlnm_var$max_lag, fun = "ns", 
                knots = logknots(dlnm_var$max_lag, dlnm_var$lagnk), 
                intercept = TRUE)
  
  return(list(basis_exp = Q, basis_lag = C))
  
})

# Create a list with the cross-basis in each neighbourhood

cb <- lapply(1:dlnm_var$n_reg, function(i_reg) {
  
  cb <- crossbasis(matrix(rep(x_temp[,i_reg], dlnm_var$max_lag + 1), 
                          ncol = dlnm_var$max_lag + 1),
                   argvar = list(fun = dlnm_var$var_fun,
                                 knots = temp_knots[,i_reg],
                                 Boundary.knots = temp_boundary[,i_reg]),
                   arglag = list(fun = "ns",
                                 knots = logknots(dlnm_var$max_lag, dlnm_var$lagnk),
                                 intercept = TRUE))
  
  return(cb)
  
})

Center_RR <- function(f.rr, f.cen, f.temp){
  
  cen <- f.temp[which.min(abs(f.temp - f.cen))]
  rr <- apply(f.rr, 2, function(x) x - x[f.temp == cen])
  rr
  
}

#--------------------------------------------------------------
# A) PANEL 3-DIMENSIONAL TEMPERATURE-LAG-RISK ASSOCIATION
# FOR REGION 1 / ITERATION 1 (SAME AS FIRST PANEL FIGURE 1 IN THE PAPER)
#--------------------------------------------------------------

# Specify the neighborhood of Barcelona (any integer from 1 to 73)
i_reg <- 1 

basis_exp <- basis_all[[i_reg]]$basis_exp
basis_lag <- basis_all[[i_reg]]$basis_lag
f.x_temp <- x_temp[,i_reg]

# Extract the ensemble of sample coefficients of the crossbasis
beta_reg <- 
  winbugs_res[, grepl(
    paste0("^beta\\[", i_reg,","), colnames(winbugs_res))]

# Initialize the array to store the RR for each temperature and lag combination
cp <- array(NA, dim = c(nrow(beta_reg), length(f.x_temp), dlnm_var$max_lag + 1))

for(i_iter in 1:nrow(beta_reg)) {
  for(i_temp in 1:length(f.x_temp)) {
    for(i_lag in 1:(dlnm_var$max_lag + 1)) {
      cp[i_iter, i_temp, i_lag] <- sum(
        (t(basis_exp[i_temp,,drop = FALSE]) %*% basis_lag[i_lag,, drop = FALSE]) *
          matrix(beta_reg[i_iter,], ncol = ncol(basis_lag), byrow = TRUE))
    }
  }
}
rm(i_iter, i_temp, i_lag)

# Position of the temperature percentile where we center the associations
i_cen <- percentiles == 0.3

# New array to store the centered association
cp_cen <- array(NA, dim = c(nrow(beta_reg), length(f.x_temp), dlnm_var$max_lag + 1))

for(i in 1:length(f.x_temp)){
  cp_cen[,i,] <- cp[,i,]-cp[,i_cen,]
}; rm(i, cp)

# Choose an iteration for plotting (e.g., the first iteration)
cp_plot <- exp(cp_cen[1,,])

# Plot three-dimensional associations
pdf("plot/figureA_exposure_lag_response_surface.pdf", 
    width = 5, height = 5)
par(mex = 0.3)
association <- 
  persp(x = f.x_temp, y = 0:dlnm_var$max_lag, z = cp_plot, border = NA,
        ticktype = "detailed", theta = 240, phi = 50,
        xlab = "\nTemperature (ºC)", ylab = "\nLag (days)", zlab = "\nRelative risk", 
        col = "lightskyblue", zlim = c(min(cp_plot), max(cp_plot)), 
        ltheta = 170, shade = 0.75, r = sqrt(3), d = 5, 
        cex.lab = 0.8, cex.axis = 0.8, main = "A) Exposure-lag-response surface")

# Add lines for risk at lag 1, risk at percentile 99, and centring temperature
lines(trans3d(x = f.x_temp, y = 1,
              z = cp_plot[,2],
              pmat = association), lwd = 2, col = "red")
lines(trans3d(x = f.x_temp[percentiles == 0.99], y = 0:dlnm_var$max_lag,
              z = cp_plot[percentiles == 0.99,],
              pmat = association), lwd = 2, col = "orange")
lines(trans3d(x = f.x_temp[i_cen], y = 0:dlnm_var$max_lag,
              z = cp_plot[i_cen,],
              pmat = association), lwd = 2, col = "black")

dev.off()

# Clean up variables
rm(i_reg, basis_exp, basis_lag, f.x_temp, beta_reg, cp_cen, i_cen, association, 
   cp_plot)

#--------------------------------------------------------------
### B) OVERALL CUMULATIVE TEMPERATURE-MORTALITY ASSOCIATION FOR REGION 1
### (SAME AS FIGURE 2 IN THE PAPER)
#--------------------------------------------------------------

# Specify the neighborhood of Barcelona (any integer from 1 to 73)
i_reg <- 1 

# Extract the ensemble of sample coefficients of the crossbasis
beta_reg <- 
  winbugs_res[, grepl(
    paste0("^beta\\[", i_reg,","), colnames(winbugs_res))]

# Calculate directly the RR of the overall cumulative temperature-mortality 
# associations for region i_reg
rr <- apply(beta_reg, 1, function(x) {
  sapply(1:length(x_temp[,i_reg]), function(i) cb[[i_reg]][i,] %*% x)
})

# Keep the exposure variable for the region
x_plot <- x_temp[, 1]

# Center the curve
i_cen <- percentiles == 0.3
rr_plot <- t(apply(rr, 1, function(x) x - rr[i_cen,]))

# Plot the overall cumulative temperature-mortality association
pdf("plot/figureB_overall_cumulative-temperature_mortality_association.pdf", 
    width = 6, height = 5)
plot(x_plot, rep(1, length(x_plot)), type = "n", log = "y",
     xlim = range(x_plot), ylim = c(0.5, 5),
     xlab = "Temperature (ºC)", ylab = "Relative risk",
     main = "B) Overall cumulative temperature-mortality association",
     axes = FALSE)
axis(1, at = seq(15, 40, by = 5))
axis(2, at = c(0.5, seq(1, 5, by = 1)), labels = c(0.5, seq(1, 5, by = 1)))

# Draw the ensemble of overall accumulated associations
for(i in 1:ncol(rr_plot)){
  lines(x_plot, exp(rr_plot[,i]), lwd = 0.2, col = "lightgrey")
}; rm(i)
abline(h = 1)

# Calculate and draw the summary associations
rr_plot_summary <- apply(rr_plot, 1, quantile, c(0.025, 0.5, 0.975))
lines(x_plot, exp(rr_plot_summary[2,]), lwd = 3)
lines(x_plot, exp(rr_plot_summary[1,]), lwd = 2, lty = 2)
lines(x_plot, exp(rr_plot_summary[3,]), lwd = 2, lty = 2)

dev.off()

# Clean up variables
rm(i_reg, beta_reg, rr, x_plot, i_cen, rr_plot, rr_plot_summary)

#--------------------------------------------------------------
### C) MAP OF THE RISKS AT PERCENTILE 99
### (SAME AS FIGURE 4C IN THE PAPER)
#--------------------------------------------------------------

# Calculate directly the RR of the overall cumulative temperature-mortality 
# associations for all regions
rr <- lapply(1:dlnm_var$n_reg, function(i_reg) {
  
  # Extract all the iterations of the coefficients of the crossbasis
  beta_reg <- winbugs_res[, 
   grepl(paste0("^beta\\[", i_reg,","), colnames(winbugs_res))]
  
  # The RR in each temperature x is the sum of the product of x transformed
  # through the crossbasis function and the coefficients of the crossbasis
  rr <- apply(beta_reg, 1, function(x) {
    sapply(1:length(x_temp[,i_reg]), function(i) cb[[i_reg]][i,] %*% x)
  })
  
  rr
  
})

# Create a function for centring the risk in each region
Center_RR <- function(f.rr, f.cen, f.temp){
  
  cen <- f.temp[which.min(abs(f.temp - f.cen))]
  rr <- apply(f.rr, 2, function(x) x - x[f.temp == cen])
  rr
  
}

# Center the relative risk in each region and extract the point estimate at 
# percentile 99
rr_plot <- sapply(1:dlnm_var$n_reg, function(i_reg) {
  x_plot <- x_temp[, i_reg]
  cen_plot <- x_plot[percentiles == 0.3]
  rr_plot <- Center_RR(f.rr = rr[[i_reg]], 
                       f.cen = cen_plot,
                       f.temp = x_plot)
  
  # Point estimate as the median of the values at percentile 99
  median(exp(rr_plot[percentiles == 0.99,]))
  
})

# Set the minimum and maximum RR in the plots
rr_max <- 2
rr_min <- exp(-log(rr_max))

# Pallete of colours for the maps
pal <- colorNumeric(palette = rev(
  c("#A90C38", "#C52A40", "#E24848", "#F16B61", "#F89183", "#FEB6A8", "#FEDAD3",
    "#FFFFFF", "#D3E5F2", "#A8CCE5", "#88B4D5", "#6D9CC3", "#5585B1", "#416F9C", 
    "#2E5A87")), domain = c(log(rr_min), log(rr_max)), reverse = FALSE)

# Plot the map with the risks
pdf("plot/figureC_map_relative_risks.pdf", 
    width = 3, height = 5)

par(mar = c(1, 1, 1, 1), omi = c(0, 0, 0.4, 0))
plot(shapefile_bcn$geometry, col = pal(log(rr_plot)), 
     main = "C) Map relative risks")

dev.off()

# Clean up variables
rm(rr, Center_RR, rr_plot, rr_max, rr_min, pal)