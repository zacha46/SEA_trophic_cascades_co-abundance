#### Re-vamping the simulation code with additional feedback and context from humans (not just ChatGPT!)
### And aiming for a simple model with building complexity 
## and still working with chat, but much more critically. 

#### Remember, the key focus of the simulation is to ensure we can recover the species interaction parameter
### And ensure that value remains robust under different common biases in camera trap data 
## or at least understand how it varies structurally. 

## Zachary Amir, Z.Amir@uq.edu.au
## Code created: August 11th, 2025
## Last updated: 

## start fresh 
rm(list = ls())

# load library
library(tidyverse)     ## For lots of functions 
# library(jagsUI)        ## for running and inspecting bayes models
library(rjags)         ## For running and inspecting bayes models
# library(coda)          ## for manipulating bayes models, but loaded w/ rjags

## Set seed for reproducibiliy
set.seed(2025)


#### Matt H's code ####

# basic idea: simulate data by writing your model in reverse

# metadata
j <- 100  # number of sites
k <- 3  # number of surveys per site
S <- 2  # number of species

# parameters
lambda_a <- rgamma(S, 3, 1)  # baseline abundance of species
lambda_b <- rnorm(1, mean = -1)  # SIV -> effect of species 1 on species 2
p <- rbeta(S, 1, 1)  # detection probabilities

# site-level abundance
N <- matrix(0, j, S)
N[, 1] <- rpois(j, lambda_a[1])  # true abundance of species 1
log_lambda <- rep(log(lambda_a[2]), j)  # linear predictor of species 2
for (i in 1:j) {
  if (N[i]) {
    log_lambda[i] <- log_lambda[i] + lambda_b * log(N[i, 1])  # increment abundance of species 1
  }
}
N[, 2] <- rpois(j, exp(log_lambda))  # true abundance of species 2

# observation process/data for both species 
y <- array(0, c(j, k, S))
for (i in 1:j) {
  for (s in 1:S) {
    y[i, , s] <- rbinom(k, N[i, s], p[s])
  }
}

## See if poission GLM can recover similar coefficent as lambda_b
summary(glm(N[,2] ~ N[,1], family = "poisson")) # still there and significant! 

## keep it clean for the next sim
rm(list = ls())


#### Simulate simple co-abundance model ##### 

# Terminology: 
# Minimal 2-species N-mixture (no ZIP)
# Labels: a0[1:2], b0[1:2], a5
# No random effects or extra covariates yet

# --- dimensions ---
nsites <- 60   # j
nreps  <- 3    # k

## inverse logit function 
inv_logit <- function(x) 1/(1+exp(-x))

# --- simulate one dataset (simple, no ZIP) ---
sim_one <- function(nsites, nreps, true_siv, log_dom = FALSE){
  ## simulate covariates 
  a0 <- rnorm(2, 0, 1)             # state intercepts
  a5 <- rnorm(1, true_siv, 0.5)    # interaction: sp2 depends on sp1
  b0 <- rnorm(2, 0, 1)             # detection intercepts
  
  ## create iZIP by incorporating occupancy prob w/ sites 
  # Here: first half sites have dom present, sub present at random
  Z_dom <- c(rep(1, nsites/2), rep(0, nsites/2))
  Z_sub <- rbinom(nsites, 1, 0.7)  # tweak rule as needed
  
  # expected abundance of dominant species 
  lambda_dom <- exp(a0[2])
  
  # Draw dominant counts first w/ iZIP 
  N_dom <- rpois(nsites, lambda_dom) * Z_dom
  
  # Subordinate mean depends on latent N.dom 
  if(log_dom){
    # can use log for stability
    lambda_sub <- exp(a0[1] + a5 * log(N_dom + 1))
  }else{
    lambda_sub <- exp(a0[1] + a5 * N_dom)
  }

   # Draw subordinate counts
  N_sub <- rpois(nsites, lambda_sub) * Z_sub
  
  ## Detection is same for both species, w/ different intercepts
  p_dom <- inv_logit(b0[2])
  p_sub <- inv_logit(b0[1])
  # create the count matrix based on true abundance and detection prob 
  y_dom <- matrix(0L, nsites, nreps)
  y_sub <- matrix(0L, nsites, nreps)
  for(j in 1:nsites){
    for(k in 1:nreps){
      y_dom[j,k] <- rbinom(1, N_dom[j], p_dom)
      y_sub[j,k] <- rbinom(1, N_sub[j], p_sub)
    }
  }
  
  list(
    data = list(nsites=nsites, nreps=nreps, y.dom=y_dom, y.sub=y_sub, Z.dom = Z_dom, Z.sub = Z_sub),
    truth = list(a0=a0, a5=a5, b0=b0, Z_dom = Z_dom, Z_sub = Z_sub,
                 lambda_dom=lambda_dom, lambda_sub=lambda_sub)
  )
} # end function

## Use the function to simulate a dataset 
sim <- sim_one(nsites, nreps, true_siv = -2, log_dom = F)

# --- JAGS model (matches the simulator) ---
jags_model <- "
model{
  # Priors
  for(i in 1:2){
    a0[i] ~ dnorm(0, 1)
    b0[i] ~ dnorm(0, 1)
  }
  a5 ~ dnorm(0, 1)

 # State process
  for(j in 1:nsites){
    # Dominant
    log(lambda.dom[j]) <- a0[2]
    N.dom[j] ~ dpois(lambda.dom[j] * Z.dom[j])

    # Subordinate depends on **latent N.dom** (can use log for stability, but not now)
    log(lambda.sub[j]) <- a0[1] + a5 * N.dom[j]
    N.sub[j] ~ dpois(lambda.sub[j] * Z.sub[j])
  }

  # Detection
  for(j in 1:nsites){
    for(k in 1:nreps){
      # dominant species detection formula
      logit(p.dom[j,k]) <- b0[2]
      y.dom[j,k] ~ dbin(p.dom[j,k], N.dom[j])
      # and fill in replicated matrix
      y.rep.dom[j,k] ~ dbin(p.dom[j,k], N.dom[j])
      
      # subordinate species detection formula 
      logit(p.sub[j,k]) <- b0[1]
      y.sub[j,k] ~ dbin(p.sub[j,k], N.sub[j])
      # and fill in replicated matrix
      y.rep.sub[j,k] ~ dbin(p.sub[j,k], N.sub[j])
      
      # Chi-square discrepancy
      exp_dom[j,k] <- N.dom[j] * p.dom[j,k]
      exp_sub[j,k] <- N.sub[j] * p.sub[j,k]
      
      E.dom[j,k] <- pow((y.dom[j,k] - exp_dom[j,k]), 2) / (exp_dom[j,k] + 0.5)
      E.rep.dom[j,k] <- pow((y.rep.dom[j,k] - exp_dom[j,k]), 2) / (exp_dom[j,k] + 0.5)
      
      E.sub[j,k] <- pow((y.sub[j,k] - exp_sub[j,k]), 2) / (exp_sub[j,k] + 0.5)
      E.rep.sub[j,k] <- pow((y.rep.sub[j,k] - exp_sub[j,k]), 2) / (exp_sub[j,k] + 0.5)
      
    } # end per nreps 
  } # end per nsites 
  
  ## summarize the ppc values 
  fit.dom <- sum(E.dom[,])
  fit.rep.dom <- sum(E.rep.dom[,])
  
  fit.sub <- sum(E.sub[,])
  fit.rep.sub <- sum(E.rep.sub[,])
}
" # end model 

# --- safe initial values so y <= N holds at start ---
make_inits <- function(data_list){
  # max values observed in matricies 
  max_y_dom <- apply(data_list$y.dom, 1, max)
  max_y_sub <- apply(data_list$y.sub, 1, max)
  
  # ensure N starts >= observed maxima, and >=1 so log(N+eps) is finite
  Ndom_init <- pmax(as.integer(max_y_dom + 1L), 1L)
  Nsub_init <- pmax(as.integer(max_y_sub + 1L), 1L)
  # Force zero where Z==0
  Ndom_init[data_list$Z.dom == 0] <- 0L
  Nsub_init[data_list$Z.sub == 0] <- 0L
  
  ## NA vector for iZIP
  # Z.dom = as.vector(rep(as.numeric(NA), length.out = length(data_list$Z.dom)))
  # Z.sub = as.vector(rep(as.numeric(NA), length.out = length(data_list$Z.sub)))
  
  # bundle in a list 
  list(
    a0 = rnorm(2, 0, 0.2),
    b0 = rnorm(2, 0, 0.2),
    a5 = rnorm(1, 0, 0.2),
    N.dom = Ndom_init, #pmax(as.integer(max_y_dom + 1L), 1L),
    N.sub = Nsub_init #pmax(as.integer(max_y_sub + 1L), 1L)
    # Z.dom = Z.dom,
    # Z.sub = Z.sub
  )
}

data_list <- sim$data
inits_list <- list(make_inits(data_list), make_inits(data_list))

# --- fit ---
jm <- jags.model(textConnection(jags_model),
                 data=data_list,
                 inits=inits_list,
                 n.chains=2, n.adapt=500)
update(jm, 1000)
vars <- c("a0","a5","b0",
          "fit.dom","fit.rep.dom",
          "fit.sub","fit.rep.sub")
samp <- coda.samples(jm, variable.names=vars, n.iter=2000, thin=2)

# ---Inspect results  ---
# extract MCMC matrix
mcmc_mat <- do.call(rbind, lapply(samp, as.matrix))

# a5 summaries
a5_draws <- mcmc_mat[, "a5"]
ci <- quantile(a5_draws, c(0.025, 0.975))
post_mean <- mean(a5_draws)
rhat <- tryCatch(gelman.diag(samp)$psrf["a5","Point est."], error=function(e) NA_real_)

# Bayes p-values for PPC
bayes_p_dom <- mean(mcmc_mat[, "fit.rep.dom"] > mcmc_mat[, "fit.dom"])
bayes_p_sub <- mean(mcmc_mat[, "fit.rep.sub"] > mcmc_mat[, "fit.sub"])

# c-hat (overdispersion) = mean(fit) / df
# df = total counts minus parameters; here approx nsites*nreps - #params
df <- (nsites * nreps) - length(c("a0","a5","b0"))
chat_dom <- mean(mcmc_mat[, "fit.dom"]) / df
chat_sub <- mean(mcmc_mat[, "fit.sub"]) / df

## display results
cat("\nTruth a5:", sim$truth$a5,
    "\nPost mean a5:", round(post_mean,3),
    "\n95% CI:", round(ci,3),
    "\nRhat(a5):", round(rhat,3),
    "\nBayes p (dom):", round(bayes_p_dom,3),
    "\nBayes p (sub):", round(bayes_p_sub,3),
    "\nc-hat (dom):", round(chat_dom,3),
    "\nc-hat (sub):", round(chat_sub,3), "\n")



# summ <- summary(samp)
# a5_draws <- do.call(rbind, lapply(samp, as.matrix))[, "a5"]
# ci <- quantile(a5_draws, c(0.025, 0.975))
# post_mean <- mean(a5_draws)
# rhat <- tryCatch(gelman.diag(samp)$psrf["a5","Point est."], error=function(e) NA_real_)
# 
# cat("\nTruth a5:", sim$truth$a5,
#     "\nPost mean a5:", round(post_mean,3),
#     "\n95% CI:", round(ci,3),
#     "\nRhat(a5):", round(rhat,3), "\n")



