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
# CODE 2: RUN B-DLNMs AND SB-DLNMs
# Model specification and execution: independent B-DLNMs and SB-DLNMs in R and
# WinBUGS with simulated mortality data.
################################################################################

# Load libraries

library(R2WinBUGS)
library(pbugs) # remotes::install_github("fisabio/pbugs")

# Load data

load("output/data_casecrossover.RData")
load("output/crossbasis_casecrossover.RData")

load("output/data_timeseries.RData")
load("output/crossbasis_timeseries.RData")
load("output/trend_timeseries.RData")
load("output/seasonality_timeseries.RData")

load("output/list_neighbours.RData")
load("output/dlnm_configuration.RData")

# Define execution models 
execute_model <- list(
  "model1" = TRUE, # Independent B-DLNM with case-crossover design
  "model2" = TRUE, # Independent B-DLNM with time-series design
  "model3" = TRUE, # SB-DLNM with case-crossover design
  "model4" = TRUE  # SB-DLNM with time-series design
)

n.iter <- list(
  "model1" = 10,
  "model2" = 10,
  "model3" = 10,
  "model4" = 10
  )
n.chains <- 3 

# The paper's results were obtained with a specific number of iterations and chains.
# By default, we have set fewer iterations and chains here for quick execution,
# but readers can adapt these values according to their preferences.

# For reference, the original settings used in the paper were:
# n.iter <- list(
#   "model1" = 150000,
#   "model2" = 100000,
#   "model3" = 100000,
#   "model4" = 50000
# )
# n.chains <- 12

#...............................................................................
### MODEL 1 (Independent B-DLNM – case-crossover design) ####
#...............................................................................

if(execute_model$model1) {
  
  winbugs_model <- function(){
    for(i in 1:n_data) {
      y[i] ~ dpois(lambda[i])
      log(lambda[i]) <- 
        p.alpha[strata[i]] +
        inprod2(beta[region[i],], w[i,])
    }
    
    for(i in 1:n_strata) {
      p.alpha[i] ~ dnorm(0, 0.001)
    }
    
    for(i in 1:n_reg){
      for(j in 1:n_coef) {
        beta[i, j] ~ dnorm(0, 0.001)
      }
    }
  }
  
  winbugs_data <-
    list(y = data_cco$mort,
         strata = as.numeric(factor(data_cco$strata)),
         region = as.numeric(data_cco$region),
         w = cb_cco,
         n_data = nrow(data_cco),
         n_strata = max(as.numeric(factor(data_cco$strata))),
         n_reg = dlnm_var$n_reg,
         n_coef = dlnm_var$n_coef)
  
  winbugs_inits <- function(){list(
    p.alpha = rnorm(winbugs_data$n_strata, 0, sd = 0.5),
    beta = matrix(rnorm(winbugs_data$n_reg * winbugs_data$n_coef, 0, sd = 1),
                  nrow = winbugs_data$n_reg, ncol = winbugs_data$n_coef)
  )}
  
  winbugs_param <- c("p.alpha", "beta")
  
  winbugs_res <- pbugs(data = winbugs_data, inits = winbugs_inits, 
                       param = winbugs_param, model = winbugs_model,
                       n.iter = n.iter$model1, n.burn = n.iter$model1/10, 
                       n.chains = n.chains, DIC = TRUE, debug = FALSE, 
                       bugs.directory = paste0(getwd(), "/WinBUGS14/"))
  
  winbugs_res <- winbugs_res$sims.matrix
  # save(winbugs_res,
  #      file = "output/predicted_simsmatrix_model1_independent_casecrossover.RData")
  rm(winbugs_model, winbugs_data, winbugs_inits, winbugs_param, winbugs_res)
  
}

#...............................................................................
### MODEL 2 (Independent B-DLNM – time-series design) ####
#...............................................................................

if(execute_model$model2) {
  
  winbugs_model <- function(){
    for(i in 1:n_data) {
      y[i] ~ dpois(lambda[i])
      log(lambda[i]) <- alpha[region[i], o[i]] +
        inprod2(beta[region[i],], w[i,]) +
        inprod2(gamma[region[i],], t[i,]) +
        inprod2(delta[region[i], year[i],], s[i,])
    }
    
    for(i in 1:n_reg){
      for(j in 1:n_coef) {
        beta[i, j] ~ dnorm(0, 0.001)
      }
      
      for(j in 1:7) { 
        alpha[i, j] ~ dnorm(0, 0.001)
      }
      
      for(j in 1:n_trend){
        gamma[i, j] ~ dnorm(0, 0.001)
      }
      
      for(j in 1:n_year) {
        for(k in 1:n_seas) {
          delta[i, j, k] ~ dnorm(0, 0.001)
        }
      }
    }
  }
  
  winbugs_data <-
    list(y = data_ts$mort,
         region = as.numeric(data_ts$region),
         o = data_ts$day_of_week,
         w = cb_ts,
         t = trend,
         s = seas,
         year = data_ts$year - min(data_ts$year) + 1, # We are assuming the same years in all regions and that all years are consecutive
         n_data = nrow(data_ts),
         n_reg = dlnm_var$n_reg,
         n_coef = dlnm_var$n_coef,
         n_trend = dlnm_var$df_trend,
         n_seas = dlnm_var$df_seas,
         n_year = length(unique(data_ts$year))
    )
  
  winbugs_inits <- function(){list(
    alpha = matrix(rnorm(winbugs_data$n_reg * 7, 0, sd = 0.5),
                   nrow = winbugs_data$n_reg, ncol = 7),
    beta = matrix(rnorm(winbugs_data$n_reg * winbugs_data$n_coef, 0, sd = 1),
                  nrow = winbugs_data$n_reg, ncol = winbugs_data$n_coef),
    gamma = matrix(rnorm(winbugs_data$n_reg * winbugs_data$n_trend, 0, sd = 0.5),
                   nrow = winbugs_data$n_reg, ncol = winbugs_data$n_trend),
    delta = array(rnorm(winbugs_data$n_reg * winbugs_data$n_year * winbugs_data$n_seas, 0, sd = 0.5),
                  dim = c(winbugs_data$n_reg, winbugs_data$n_year, winbugs_data$n_seas))
  )}
  
  winbugs_param <- c("alpha", "beta", "gamma", "delta")
  
  winbugs_res <- pbugs(data = winbugs_data, inits = winbugs_inits, 
                       param = winbugs_param, model = winbugs_model,
                       n.iter = n.iter$model2, n.burn = n.iter$model2/10, 
                       n.chains = n.chains, DIC = TRUE, debug = FALSE, 
                       bugs.directory = paste0(getwd(), "/WinBUGS14/"))
  
  winbugs_res <- winbugs_res$sims.matrix
  # save(winbugs_res, 
  #      file = "output/predicted_simsmatrix_model2_independent_timseseries.RData")
  rm(winbugs_model, winbugs_data, winbugs_inits, winbugs_param, winbugs_res)
  
}

#...............................................................................
### MODEL 3 (SB-DLNM – case-crossover design) ####
#...............................................................................

if(execute_model$model3) { 
  
  winbugs_model <- function(){
    for(i in 1:n_data) {
      y[i] ~ dpois(lambda[i])
      log(lambda[i]) <- 
        p.alpha[strata[i]] +
        inprod2(beta[region[i],], w[i,])
    }
    
    for(i in 1:n_reg){
      for(j in 1:n_coef) {
        beta[i, j] <- mu[j] + sigma[j] * beta_ast[i, j]
        beta_ast[i, j] ~ dnorm(mean_leroux[i, j], prec_leroux[i])
      }
      prec_leroux[i] <- (1 - rho + rho * n[i])
    }
    
    for(j in 1:n_coef) {
      for(i in 1:n_reg) {
        mean_leroux[i, j] <- (rho / prec_leroux[i]) * 
          sum(beta_ast_r[(index_r[i] + 1):index_r[i + 1], j]) # sum the beta_ast of neighbouring regions of i
      }
      for(k in 1:sum_n) {
        beta_ast_r[k, j] <- beta_ast[r_prime[k], j] # r_prime is the list of neighbors of r
      }
    }
    
    for(j in 1:n_coef) {
      mu[j] ~ dunif(-100, 100)
      sigma[j] ~ dunif(0, 10)
    }
    
    rho ~ dunif(0, 1)
    
    for(i_strata in 1:n_strata) {
      p.alpha[i_strata] ~ dnorm(0, 0.001)
    }
  }
  
  winbugs_data <-
    list(y = data_cco$mort,
         strata = as.numeric(factor(data_cco$strata)),
         region = as.numeric(data_cco$region),
         w = cb_cco,
         n_data = nrow(data_cco),
         n_strata = max(as.numeric(factor(data_cco$strata))),
         n_reg = dlnm_var$n_reg,
         n_coef = dlnm_var$n_coef,
         n = list_neig$num, # number of neighbours of each region
         sum_n = sum(list_neig$num), # total number of neighbouring relationships
         r_prime = list_neig$adj, # vector of neighbouring relationships
         index_r = c(0, cumsum(list_neig$num)) # indexes defining the neighbouring relationships in the vector
    )
  
  winbugs_inits <- function(){list(
    p.alpha = rnorm(winbugs_data$n_strata, 0, sd = 0.5),
    mu = rnorm(winbugs_data$n_coef, 0, 0.1),
    rho = runif(1),
    beta_ast = matrix(rnorm(winbugs_data$n_reg * winbugs_data$n_coef, 0, 1), 
                      ncol = winbugs_data$n_coef),
    sigma = runif(winbugs_data$n_coef, 0, 0.2)
  )}
  
  winbugs_param <- c("p.alpha", "beta", "mu", "rho", "beta_ast", "sigma")
  
  winbugs_res <- pbugs(data = winbugs_data, inits = winbugs_inits, 
                       param = winbugs_param, model = winbugs_model,
                       n.iter = n.iter$model3, n.burn = n.iter$model3/10, 
                       n.chains = n.chains, DIC = TRUE, debug = FALSE, 
                       bugs.directory = paste0(getwd(), "/WinBUGS14/"))
  
  winbugs_res <- winbugs_res$sims.matrix
  # save(winbugs_res, 
  #      file = "output/predicted_simsmatrix_model3_spatial_casecrossover.RData")
  rm(winbugs_model, winbugs_data, winbugs_inits, winbugs_param, winbugs_res)
  
}

#...............................................................................
### MODEL 4 (SB-DLNM – time-series design) ####
#...............................................................................

if(execute_model$model4) { 
  
  winbugs_model <- function(){
    for(i in 1:n_data) {
      y[i] ~ dpois(lambda[i])
      log(lambda[i]) <- alpha[region[i], o[i]] +
        inprod2(beta[region[i],], w[i,]) +
        inprod2(gamma[region[i],], t[i,]) +
        inprod2(delta[region[i], year[i],], s[i,])
    }
    
    for(i in 1:n_reg){
      for(j in 1:n_coef) {
        beta[i, j] <- mu[j] + sigma[j] * beta_ast[i, j]
        beta_ast[i, j] ~ dnorm(mean_leroux[i, j], prec_leroux[i])
      }
      prec_leroux[i] <- (1 - rho + rho * n[i])
    }
    
    for(j in 1:n_coef) {
      for(i in 1:n_reg) {
        mean_leroux[i, j] <- (rho / prec_leroux[i]) * 
          sum(beta_ast_r[(index_r[i] + 1):index_r[i + 1], j]) # sum the beta_ast of neighbouring regions of i
      }
      for(k in 1:sum_n) {
        beta_ast_r[k, j] <- beta_ast[r_prime[k], j] # r_prime is the list of neighbors of r
      }
    }
    
    for(j in 1:n_coef) {
      mu[j] ~ dunif(-100, 100)
      sigma[j] ~ dunif(0, 10)
    }
    
    rho ~ dunif(0, 1)
    
    for(i in 1:n_reg){
      
      for(j in 1:7) { 
        alpha[i, j] ~ dnorm(0, 0.001)
      }
      
      for(j in 1:n_trend){
        gamma[i, j] ~ dnorm(0, 0.001)
      }
      
      for(j in 1:n_year) {
        for(k in 1:n_seas) {
          delta[i, j, k] ~ dnorm(0, 0.001)
        }
      }
    }
  }
  
  winbugs_data <-
    list(y = data_ts$mort,
         region = as.numeric(data_ts$region),
         o = data_ts$day_of_week,
         w = cb_ts,
         t = trend,
         s = seas,
         year = data_ts$year - min(data_ts$year) + 1, # We are assuming the same years in all regions and that all years are consecutive
         n_data = nrow(data_ts),
         n_reg = dlnm_var$n_reg,
         n_coef = dlnm_var$n_coef,
         n_trend = dlnm_var$df_trend,
         n_seas = dlnm_var$df_seas,
         n_year = length(unique(data_ts$year)),
         n = list_neig$num, # number of neighbours of each region
         sum_n = sum(list_neig$num), # total number of neighbouring relationships
         r_prime = list_neig$adj, # vector of neighbouring relationships
         index_r = c(0, cumsum(list_neig$num)) # indexes defining the neighbouring relationships in the vector
    )
  
  winbugs_inits <- function(){list(
    alpha = matrix(rnorm(winbugs_data$n_reg * 7, 0, sd = 0.5),
                   nrow = winbugs_data$n_reg, ncol = 7),
    gamma = matrix(rnorm(winbugs_data$n_reg * winbugs_data$n_trend, 0, sd = 0.5),
                   nrow = winbugs_data$n_reg, ncol = winbugs_data$n_trend),
    delta = array(rnorm(winbugs_data$n_reg * winbugs_data$n_year * winbugs_data$n_seas, 0, sd = 0.5),
                  dim = c(winbugs_data$n_reg, winbugs_data$n_year, winbugs_data$n_seas)),
    mu = rnorm(winbugs_data$n_coef, 0, 0.1),
    rho = runif(1),
    beta_ast = matrix(rnorm(winbugs_data$n_reg * winbugs_data$n_coef, 0, 1), 
                      ncol = winbugs_data$n_coef),
    sigma = runif(winbugs_data$n_coef, 0, 0.2)
  )}
  
  winbugs_param <- c("alpha", "beta", "gamma", "delta", "mu", "rho", "beta_ast", "sigma")
  
  winbugs_res <- pbugs(data = winbugs_data, inits = winbugs_inits, 
                       param = winbugs_param, model = winbugs_model,
                       n.iter = n.iter$model4, n.burn = n.iter$model4/10, 
                       n.chains = n.chains, DIC = TRUE, debug = FALSE, 
                       bugs.directory = paste0(getwd(), "/WinBUGS14/"))
  
  winbugs_res <- winbugs_res$sims.matrix
  # save(winbugs_res, 
  #      file = "output/predicted_simsmatrix_model4_spatial_timeseries.RData")
  rm(winbugs_model, winbugs_data, winbugs_inits, winbugs_param, winbugs_res)
  
}