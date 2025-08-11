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
library(jagsUI)

## Set seed for reproduciblity
set.seed(2025)

##### Function for simulating data  ##### 

# load dimensions for testing + terminiology 
nsites  <- 60    # j, rows in the matrix
nreps   <- 3     # k, columns in the matrix
narea   <- 4     # number of landscapes for RE
nyear   <- 3     # number of years for RE 
nsource <- 2     # number of data providers for RE

## inverse logit function, used inside function 
inv_logit <- function(x) 1/(1+exp(-x))

# Load the function
coabundance_simulator <- function(nsites, nreps, narea, nyear, nsource, 
                                  true_siv, log_dom = FALSE,
                                  area = NULL, Zdom_area = NULL, Zsub_area = NULL,
                                  # --- bias toggles & strengths ---
                                  bias_doublecount = FALSE,   dc_rate = 0.4,               # 1) double counting, with a customizable rate 
                                  bias_sp_autocorr = FALSE,   rho_sp = 0.3, p_cross = 0.4, # 2) spillover across sites, 
                                  # rho = how much neighbor contributes, and p = the chance critter is seen at neighbor site 
                                  bias_unmeas_state = FALSE,  u_state_sd = 0.6,            # 3) unmeasured state cov, w/ customizable rate
                                  bias_unmeas_det = FALSE,    u_det_sd = 0.6){             # 4) unmeasured det cov, w/ customizable rate 
  ### Assign areas if not supplied
  if (is.null(area)) {
    q <- rep(floor(nsites / narea), narea)
    # if there are less less levels than sites 
    if (sum(q) < nsites){ 
      # add 1
      q[seq_len(nsites - sum(q))] <- q[seq_len(nsites - sum(q))] + 1L
    }
    # create area groupings 
    area <- rep(seq_len(narea), times = q)
  } # end null condition
  
  ## Validate Z vectors supplied to the function
  if (is.null(Zdom_area) || is.null(Zsub_area)) {
    # default to present in all areas if not provided
    Zdom_area <- rep(1L, narea)
    Zsub_area <- rep(1L, narea)
  } # end null z condition
  ## Make sure Z is kosher and equals area groupings in length
  stopifnot(length(Zdom_area) == narea, length(Zsub_area) == narea)
  Zdom_area <- as.integer(Zdom_area); Zsub_area <- as.integer(Zsub_area)
  
  ## Assign year and source random effects
  year = sample(1:nyear, size = nsites, replace = TRUE)
  source = sample(1:nsource, size = nsites, replace = TRUE)
  
  ## simulate key covariates 
  a0 <- rnorm(2, 0, 1)             # state intercepts
  a5 <- rnorm(1, true_siv, 0.5)    # interaction: sp2 depends on sp1
  b0 <- rnorm(2, 0, 1)             # detection intercepts
  ## Additional covaraites 
  a1 <- rnorm(2, 0, 1)             # flii effects
  a2 <- rnorm(2, 0, 1)             # hfp effects
  a3 <- rnorm(2, 0, 1)             # elev effects
  a4 <- rnorm(2, 0, 1)             # community detections 
  ## Detection covariate 
  b1 <- rnorm(2, 0, 1)             # sampling effort observation covariate 
  sd.p <- runif(2, 0, 1)           # OD RE per species 
  ## bias covariates
  u_state <- rnorm(nsites, 0, u_state_sd)   # site latent driver

  # Site covariates (scale for stability)
  flii <- scale(rnorm(nsites, 0, 1))[,1]
  hfp  <- scale(rnorm(nsites, 0, 1))[,1]
  elev <- scale(rnorm(nsites, 0, 1))[,1]
  comm_det <- scale(log(runif(nsites, 1, 100) + 1))[,1]
  
  # Detection coviariate 
  # raw_effort <- matrix(runif(nsites * nreps, min = 0.5, max = 1.5), nsites, nreps)
  raw_effort = matrix(rlnorm(nsites * nreps, meanlog=0, sdlog=0.6), nsites, nreps) # add wider spread 
  cams <- scale(raw_effort)
  
  ## generate the random intercepts using moderate SDs to avoid issues 
  sigma_a6 <- 0.5 ; sigma_a7 <- 0.5   # dom/sub ~ area
  a6_area  <- rnorm(narea, 0, sigma_a6)
  a7_area  <- rnorm(narea, 0, sigma_a7)
  
  ## Do the same for the two other random effects: 
  # year
  sigma_a8 <- 0.5 ; sigma_a9 <- 0.5
  a8_year  <- rnorm(nyear, 0, sigma_a8) # subordinate 
  a9_year  <- rnorm(nyear, 0, sigma_a9) # dominant 
  # source 
  sigma_b2 <- 0.5 ; sigma_b3 <- 0.5
  b2_source <- rnorm(nsource, 0, sigma_b2) # subordinate
  b3_source <- rnorm(nsource, 0, sigma_b3) # dominant 
  
  #
  ##
  ### State formula 
  
  # expected abundance of dominant species 
  lambda_dom <- exp(a0[2] + a1[2]*flii + a2[2]*hfp + a3[2]*elev + a4[2]*comm_det + a7_area[area] + a9_year[year] + if(bias_unmeas_state){u_state}else{0}) 
  # Draw dominant counts first w/ iZIP 
  N_dom <- rpois(nsites, lambda_dom) * Zdom_area[area]
  
  # Subordinate mean depends on latent N.dom 
  if(log_dom){
    # can use log for stability
    lp_sub <- a0[1] + a1[1]*flii + a2[1]*hfp + a3[1]*elev + a4[1]*comm_det + a5*log(N_dom + 1) + a6_area[area] + a8_year[year] + if(bias_unmeas_state){u_state}else{0}
  }else{
    ## Center N_dom by calculating the Z-score: value - mean / sd 
    N_dom_c <- (N_dom - mean(N_dom))/sd(N_dom)
    lp_sub <- a0[1] + a1[1]*flii + a2[1]*hfp + a3[1]*elev + a4[1]*comm_det + a5*N_dom_c + a6_area[area] + a8_year[year] + if(bias_unmeas_state){u_state}else{0}
  } # end log condition
  ## clamp the linear predictor to a more narrow range to prevent overflow
  lp_sub <- pmax(pmin(lp_sub, 20), -20)   
  lambda_sub <- exp(lp_sub)
  
  # Draw subordinate counts
  N_sub <- rpois(nsites, lambda_sub) * Zsub_area[area]
  
  ### Prepare for the spatial autocorrelation bias 
  # create a simple neighbor map: within each area, neighbors are the previous/next site indices
  neighbors <- vector("list", nsites)
  idx_by_area <- split(seq_len(nsites), area)
  for (grp in idx_by_area) {
    for (pos in seq_along(grp)) {
      nb <- integer(0)
      if (pos > 1) nb <- c(nb, grp[pos - 1])
      if (pos < length(grp)) nb <- c(nb, grp[pos + 1])
      neighbors[[grp[pos]]] <- nb
    }
  }
  
  #
  ##
  ### Detection formula 
  
  # create the count matrix based on true abundance and detection prob 
  y_dom <- matrix(0L, nsites, nreps)
  y_sub <- matrix(0L, nsites, nreps)
  # specify the OD RE outside of the loop 
  eps_dom <- rnorm(1, 0, sd.p[2])
  eps_sub <- rnorm(1, 0, sd.p[1])
  
  ## Create a bias for unobserved covariate in detection formula
  if (bias_unmeas_det) {
    u_det <- matrix(rnorm(nsites*nreps, 0, u_det_sd), nsites, nreps)
  } else {
    # defaults to zero if not used. 
    u_det <- matrix(0, nsites, nreps)
  }
  
  # for each site 
  for(j in 1:nsites){
    # and for each rep 
    for(k in 1:nreps){
      
      ## Detection linear predictor is same for both species, w/ different intercepts/coefficents
      lp_dom = b0[2] + b1[2]*cams + eps_dom + b2_source[source[j]] + u_det[j,k]
      lp_sub = b0[1] + b1[1]*cams + eps_sub + b3_source[source[j]] + u_det[j,k]
      
      ## take the logit transformation of the linear predictors to calculate det prob 
      p_dom <- inv_logit(max(min(lp_dom, 250), -250))
      p_sub <- inv_logit(max(min(lp_sub, 250), -250))
      
      ## now combine det prob with latent abundance to inform the count history matrix using a binomial distribution
      y_base_dom <- rbinom(1, N_dom[j], p_dom)
      y_base_sub <- rbinom(1, N_sub[j], p_sub)
      
      ### Add a double counting bias
      if (bias_doublecount) {
        # extra "double-counted" opportunities scale with abundance
        dc_dom <- rpois(1, lambda = dc_rate * max(N_dom[j], 0))
        dc_sub <- rpois(1, lambda = dc_rate * max(N_sub[j], 0))
        # draw the extra counts 
        y_ext_dom <- rbinom(1, dc_dom, p_dom)
        y_ext_sub <- rbinom(1, dc_sub, p_sub)
        # combine with base counts
        y_dom[j,k] <- y_base_dom + y_ext_dom
        y_sub[j,k] <- y_base_sub + y_ext_sub
      } else {
        # but if false, dont add any double counting! 
        y_dom[j,k] <- y_base_dom
        y_sub[j,k] <- y_base_sub
      }
      
      ### Add a spatial autocorrelation bias 
      if (bias_sp_autocorr) {
        # neighbor abundance pool (dominant & subordinate)
        Nnb_dom <- if (length(neighbors[[j]])>0) round(rho_sp * mean(N_dom[neighbors[[j]]])) else 0L
        Nnb_sub <- if (length(neighbors[[j]])>0) round(rho_sp * mean(N_sub[neighbors[[j]]])) else 0L
        # spillover detection from neighbors (can exceed local N)
        y_dom[j,k] <- y_dom[j,k] + rbinom(1, max(Nnb_dom,0), p_cross * p_dom)
        y_sub[j,k] <- y_sub[j,k] + rbinom(1, max(Nnb_sub,0), p_cross * p_sub)
      }
      
    } # end per rep
  } # end per site
  
  ## save all the data in a list 
  list(
    data = list(nsites=nsites, nreps=nreps, 
                y.dom=y_dom, y.sub=y_sub, 
                Zdom_area=Zdom_area, Zsub_area=Zsub_area, 
                flii=as.numeric(flii), hfp=as.numeric(hfp), elev=as.numeric(elev), comm_det = as.numeric(comm_det),
                narea = narea, area = as.integer(area),
                cams = cams,
                nyear = nyear, year = as.numeric(year), 
                nsource = nsource, source = as.numeric(source),
                N.dom = N_dom, N.sub = N_sub),
    # and also keep the true values for comparison later 
    truth = list(a0=a0, a5=a5, b0=b0,
                 Zdom_area=Zdom_area, Zsub_area=Zsub_area, sd.p = sd.p,
                 lambda_dom=lambda_dom, lambda_sub=lambda_sub)
  )
} # end function

#### Test the function! 

## Define key values for the function
nsites  <- 60; nreps <- 3
narea   <- 4
nyear   <- 2
nsource <- 2
area   <- rep(1:narea, times=c(20,15,15,10))
Zdom_area <- c(1, 1, 0, 1)  # area 3 absent for dominant
Zsub_area <- c(0, 1, 1, 1)  # area 1 absent for subordinate

## Use the function to simulate a dataset 
sim <- coabundance_simulator(nsites, nreps, narea, nyear, nsource, 
                             true_siv = 1, log_dom = T, area = area, 
                             Zdom_area = Zdom_area, Zsub_area = Zsub_area,
                             bias_doublecount = FALSE, bias_sp_autocorr = FALSE, 
                             bias_unmeas_state = FALSE, bias_unmeas_det = FALSE)


##### Load a JAGS model and inital values #####
jags_model <- "
model{
  # Regular priors for both species
  for(i in 1:2){
    a0[i] ~ dnorm(0, 1)   # state intercept 
    b0[i] ~ dnorm(0, 1)   # det intercept 
    a1[i] ~ dnorm(0, 1)   # flii
    a2[i] ~ dnorm(0, 1)   # hfp
    a3[i] ~ dnorm(0, 1)   # elev
    a4[i] ~ dnorm(0, 1)   # comm_det
    b1[i] ~ dnorm(0, 1)   # sampling effort 
    # OD params for observation model
    sd.p[i] ~ dunif(0, 1) 
    tau.p[i] <- pow(sd.p[i], -2)
  }
  # SIV prior 
  a5 ~ dnorm(0, 1)
  
  # landscape RE hyper priors 
  sigma.a6 ~ dunif(0, 5)
  var.a6 <- 1 / pow(sigma.a6, 2)
  sigma.a7 ~ dunif(0, 5)
  var.a7 <- 1 / pow(sigma.a7, 2)
  
  for (k in 1:narea){
    a6[k] ~ dnorm(0, var.a6)   # sub ~ area
    a7[k] ~ dnorm(0, var.a7)   # dom ~ area
  }
  
  # Year RE hyper priors 
  sigma.a8 ~ dunif(0, 5)
  var.a8 <- 1/pow(sigma.a8, 2)
  sigma.a9 ~ dunif(0, 5)
  var.a9 <- 1/pow(sigma.a9, 2)
  
  for (k in 1:nyear){
    a8[k] ~ dnorm(0, var.a6)   # sub ~ year
    a9[k] ~ dnorm(0, var.a7)   # dom ~ year
  }
  
  # Source RE hyper prior 
  sigma.b2 ~ dunif(0, 5)
  var.b2 <- 1/pow(sigma.b2, 2)
  sigma.b3 ~ dunif(0, 5)
  var.b3 <- 1/pow(sigma.b3, 2)
  
  for (k in 1:nsource){
    b2[k] ~ dnorm(0, var.b2)   # sub ~ source
    b3[k] ~ dnorm(0, var.b3)   # dom ~ source
  }

 # State process
  for(j in 1:nsites){
    # Dominant
    log(lambda.dom[j]) <- a0[2] + a1[2]*flii[j] + a2[2]*hfp[j] + a3[2]*elev[j] + a4[2]*comm_det[j] + a7[area[j]] + a9[year[j]]
    N.dom[j] ~ dpois(lambda.dom[j] * Zdom_area[area[j]])

    # Subordinate depends on **latent N.dom** (can use log for stability, but not now)
    log(lambda.sub[j]) <- a0[1] + a5 * log(N.dom[j] + 1.0E-6) + a1[1]*flii[j] + a2[1]*hfp[j] + a3[1]*elev[j] + a4[1]*comm_det[j] + a6[area[j]] + a8[year[j]]
    N.sub[j] ~ dpois(lambda.sub[j] * Zsub_area[area[j]])
  }

  # Detection
  for(j in 1:nsites){
    for(k in 1:nreps){
    
      # dominant species detection formula
      lp.dom[j,k] <- b0[2] + b1[2]*cams[j,k] + eps.p.dom[j,k] + b2[source[j]]
      # implement a stable logit transformation
      p.dom[j,k] <- 1 / (1 + exp(-max(-250, min(250, lp.dom[j,k])))) 
      # calculate det prob
      y.dom[j,k] ~ dbin(p.dom[j,k], N.dom[j])
      # and the ODRE  
      eps.p.dom[j,k] ~ dnorm(0, tau.p[2])
      
      # and fill in replicated matrix
      y.rep.dom[j,k] ~ dbin(p.dom[j,k], N.dom[j])
      
      # subordinate species detection formula 
      lp.sub[j,k] <- b0[1] + b1[1]*cams[j,k] + eps.p.sub[j,k] + b3[source[j]]
      # implement a stable logit transformation
      p.sub[j,k] <- 1 / (1 + exp(-max(-250, min(250, lp.sub[j,k])))) 
      # calculate det prob
      y.sub[j,k] ~ dbin(p.sub[j,k], N.sub[j])
      # and the ODRE
      eps.p.sub[j,k] ~ dnorm(0, tau.p[1])
      
      # and fill in replicated matrix
      y.rep.sub[j,k] ~ dbin(p.sub[j,k], N.sub[j])
      
      # Chi-square discrepancy
      exp_dom[j,k] <- N.dom[j] * p.dom[j,k]
      exp_sub[j,k] <- N.sub[j] * p.sub[j,k]
      
      E.dom[j,k]     <- pow((y.dom[j,k]     - exp_dom[j,k]), 2) / (exp_dom[j,k] + 0.5)
      E.rep.dom[j,k] <- pow((y.rep.dom[j,k] - exp_dom[j,k]), 2) / (exp_dom[j,k] + 0.5)

      E.sub[j,k]     <- pow((y.sub[j,k]     - exp_sub[j,k]), 2) / (exp_sub[j,k] + 0.5)
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
  Ndom_init[data_list$Zdom_area[data_list$area] == 0] <- 0L
  Nsub_init[data_list$Zsub_area[data_list$area] == 0] <- 0L
  
  ## NA vector for iZIP
  # Z.dom = as.vector(rep(as.numeric(NA), length.out = length(data_list$Z.dom)))
  # Z.sub = as.vector(rep(as.numeric(NA), length.out = length(data_list$Z.sub)))
  
  # bundle in a list 
  list(
    # key params
    a0 = rnorm(2, 0, 0.2),
    b0 = rnorm(2, 0, 0.2),
    a5 = rnorm(1, 0, 0.2),
    # site covs 
    a1 = rnorm(2, 0, 0.2),
    a2 = rnorm(2, 0, 0.2),
    a3 = rnorm(2, 0, 0.2),
    a4 = rnorm(2, 0, 0.2),
    # det covs 
    b1 = rnorm(2, 0, 0.2),
    sd.p = runif(2, 0, 0.3),
    # abundance inital values 
    N.dom = Ndom_init, #pmax(as.integer(max_y_dom + 1L), 1L),
    N.sub = Nsub_init #pmax(as.integer(max_y_sub + 1L), 1L)
  )
}

inits_list <- list(make_inits(sim$data), make_inits(sim$data))

##### Fit the model and inspect diagnostics #####

# # load DIC to track deviance in the model 
# load.module("dic") 
# 
# # Create the model object and compile it 
# jm <- jags.model(textConnection(jags_model),
#                  data=sim$data,
#                  inits=inits_list,
#                  n.chains=2, n.adapt=500)
# # Burn in 1000 iterations (i.e. update w/out saving )
# update(jm, 1000)
# # select which variables to monitor 
# vars <- c("a0","a1","a2","a3","a4","a5","b0","b1","sd.p",
#           "fit.dom","fit.rep.dom","fit.sub","fit.rep.sub","deviance")
# # Draw the MCMC samples 
# samp <- coda.samples(jm, variable.names = vars, n.iter=2000, thin=2)

### Can also run jagsUI style?
# 1) Save the model string to a real file
model_path <- tempfile(fileext = ".jags")
writeLines(jags_model, con = model_path)

# 2) Run jagsUI with the file path
params <- c("a0","a1","a2","a3","a4","a5","b0","b1","sd.p", "N.dom", "N.sub",
            "fit.dom","fit.rep.dom","fit.sub","fit.rep.sub") # no 'deviance' unless DIC=TRUE
## remove simulated abundances from list and save in a separate DF 
abund = data.frame("N.dom_simulated" = sim$data$N.dom,
                   "N.sub_simulated" = sim$data$N.sub)
sim$data$N.dom = NULL
sim$data$N.sub = NULL

# 3) call the model 
mod <- jagsUI::jags(data   = sim$data,
                    inits  = inits_list,
                    parameters.to.save = params,
                    model.file = model_path,
                    n.chains = 2,
                    n.adapt  = 1000,
                    n.burnin = 1000,
                    n.iter   = 3000,  # total iterations INCLUDING burn-in
                    n.thin   = 50,
                    parallel = TRUE,   # TRUE on HPC
                    DIC      = FALSE)
# print(mod, 2)

## Grab info about a5
coeff = data.frame(mod$summary[,c(1:3,7:10)])
post_mean = coeff$mean[rownames(coeff) == "a5"]
lower = round(coeff$X2.5.[rownames(coeff) == "a5"], 3) 
upper = round(coeff$X97.5.[rownames(coeff) == "a5"], 3)
ci = paste(lower, upper, collapse = " - ")
rhat = coeff$Rhat[rownames(coeff) == "a5"]

# Bayes p-values for PPC
bayes_p_dom <- mean(mod$sims.list$fit.rep.dom > mod$sims.list$fit.dom)
bayes_p_sub <- mean(mod$sims.list$fit.rep.sub > mod$sims.list$fit.sub) 

# Pull posterior draws
a_sub = mod$sims.list$fit.sub
b_sub = mod$sims.list$fit.rep.sub
a_dom = mod$sims.list$fit.dom
b_sub = mod$sims.list$fit.rep.dom

# c-hat calculation 
chat_dom <- mean(a_dom / b_dom)
chat_sub <- mean(a_sub / b_sub)

## display results as text 
cat("\nTrue SIV:", sim$truth$a5,
    "\nPosterior mean SIV:", round(post_mean,3),
    "\n95% CI:", ci, #round(ci,3),
    "\nRhat(SIV):", round(rhat,3),
    "\nBayes p-value (dom):", round(bayes_p_dom,3),
    "\nBayes p-value (sub):", round(bayes_p_sub,3),
    "\nc-hat (dom):", round(chat_dom,3),
    "\nc-hat (sub):", round(chat_sub,3), "\n")

## Also compare true vs simulated latent abundance 

# Create empty df to fill in estimates abundance of both species
est.dat = data.frame(matrix(NA, nrow = 0, ncol = 7))
names(est.dat) = c("Sub_abundance", "Dom_abundance",
                   "lower_sub", "upper_sub",
                   "lower_dom", "upper_dom",
                   # "Sampling_Unit", "Landscape",
                   "sim_test")
# Then fill it in!
est.dat[1:length(colMeans(mod$sims.list$N.sub)), 1] = colMeans(mod$sims.list$N.sub)
est.dat[1:length(colMeans(mod$sims.list$N.dom)), 2] = colMeans(mod$sims.list$N.dom)
est.dat[1:length(colMeans(mod$sims.list$N.sub)), 3] = apply(mod$sims.list$N.sub, 2, quantile, probs = 0.025)
est.dat[1:length(colMeans(mod$sims.list$N.sub)), 4] = apply(mod$sims.list$N.sub, 2, quantile, probs = 0.975)
est.dat[1:length(colMeans(mod$sims.list$N.dom)), 5] = apply(mod$sims.list$N.dom, 2, quantile, probs = 0.025)
est.dat[1:length(colMeans(mod$sims.list$N.dom)), 6] = apply(mod$sims.list$N.dom, 2, quantile, probs = 0.975)
est.dat[1:length(colMeans(mod$sims.list$N.sub)), 7] = "test"

## Now add in the true simulated abundance values that we removed earlier 
est.dat = cbind(abund, est.dat)

## options to calculate the differences
# option 1, coverage w/in credible interval
coverage_sub <- mean(est.dat$N.sub_simulated >= est.dat$lower_sub & est.dat$N.sub_simulated <= est.dat$upper_sub)
coverage_dom <- mean(est.dat$N.dom_simulated >= est.dat$lower_dom & est.dat$N.dom_simulated <= est.dat$upper_dom)
# ideally coverage is greater than .90%
if(any(coverage_sub > .9 & coverage_dom > .9)){
  print("The estimated abundance values provide acceptable coverage within 95% CI of the estimated values")
}else{
  print("The estimated abundance values DO NOT provide acceptable coverage within 95% CI of the estimated values")
}

# option 2, root mean square deviation -> how far are the posteriors from the truth?
# RMSE summarizes both the bias and the variance in abundance estimates --> shouldnt be greater than 5
rmse_sub <- sqrt(mean((est.dat$N.sub_simulated - est.dat$Sub_abundance)^2))
rmse_dom <- sqrt(mean((est.dat$N.dom_simulated - est.dat$Dom_abundance)^2))
# lower rmse is better!
if(any(rmse_sub > 5 & rmse_sub > 5)){
  print("The estimated abundance values are NOT close to the simulated truth.")
}else{
  print("The estimated abundance values are close to the simulated truth.")
}


## clean everything up! (except functions)
rm(abund, coeff, dic.out, est.dat, inits_list, jm, M, mcmc_mat, mod, samp, sim, source, 
   a_dom, a_sub, a5_draws, area, b_dom, b_sub, bayes_p_dom, bayes_p_sub, chat_dom, chat_sub, ci, 
   coverage_dom, coverage_sub, Dbar, df_eff, jags_model, lower, model_path, narea, nobs, 
   nreps, nsites, nsource, nyear, params, pD, post_mean, post_means, rhat, rmse_dom, rmse_sub, upper, 
   vars, years, Zdom_area, Zsub_area)


#### Simulate different scenarios #### 

## the idea here is that we get a range of SIVs and different biases
a5_values = c(-2, -1, 0, 1, 2)
biases = c("none", "double_count", "spatial", 
           "unmeasured_state", "unmeasured_detection","all_biases")

## Define key values for the function
nsites  <- 500
nreps   <- 5
narea   <- 12
nyear   <- 3
nsource <- 4

## Even split across areas (simple & transparent)
sites_per_area <- rep(floor(nsites / narea), narea)
if (sum(sites_per_area) < nsites) {
  # add the remainder to the first few areas
  sites_per_area[seq_len(nsites - sum(sites_per_area))] <- sites_per_area[seq_len(nsites - sum(sites_per_area))] + 1L
}
area <- rep(seq_len(narea), times = sites_per_area)

## Area-level iZIP: simple fixed pattern (transparent)
# Example: dominant absent in 4 areas; subordinate absent in 3 areas
Zdom_area <- rep(1L, narea); Zdom_area[c(1,3,7,10)] <- 0L  
Zsub_area <- rep(1L, narea); Zsub_area[c(2,5,8)]       <- 0L


## run a loop! 
data_list = list() # store data bundles here
for(i in 1:length(a5_values)){
  # grab one 
  a5 = a5_values[i]
  # store temp values here
  temp = list()
  for(l in 1:length(biases)){
    # condition against each bias
    if(biases[l] == "none"){
      # switch all biases to false 
      sim = coabundance_simulator(nsites, nreps, narea, nyear, nsource, 
                                  true_siv = a5, log_dom = T, area = area, 
                                  Zdom_area = Zdom_area, Zsub_area = Zsub_area,
                                  bias_doublecount = FALSE, bias_sp_autocorr = FALSE, 
                                  bias_unmeas_state = FALSE, bias_unmeas_det = FALSE)
    } # end none 
    if(biases[l] == "double_count"){
      # switch relevant bias to true 
      sim = coabundance_simulator(nsites, nreps, narea, nyear, nsource, 
                                  true_siv = a5, log_dom = T, area = area, 
                                  Zdom_area = Zdom_area, Zsub_area = Zsub_area,
                                  bias_doublecount = TRUE, bias_sp_autocorr = FALSE, 
                                  bias_unmeas_state = FALSE, bias_unmeas_det = FALSE)
    } # end double count 
    if(biases[l] == "spatial"){
      # switch relevant bias to true 
      sim = coabundance_simulator(nsites, nreps, narea, nyear, nsource, 
                                  true_siv = a5, log_dom = T, area = area, 
                                  Zdom_area = Zdom_area, Zsub_area = Zsub_area,
                                  bias_doublecount = FALSE, bias_sp_autocorr = TRUE, 
                                  bias_unmeas_state = FALSE, bias_unmeas_det = FALSE)
    } # end spatial
    if(biases[l] == "unmeasured_state"){
      # switch relevant bias to true 
      sim = coabundance_simulator(nsites, nreps, narea, nyear, nsource, 
                                  true_siv = a5, log_dom = T, area = area, 
                                  Zdom_area = Zdom_area, Zsub_area = Zsub_area,
                                  bias_doublecount = FALSE, bias_sp_autocorr = FALSE, 
                                  bias_unmeas_state = TRUE, bias_unmeas_det = FALSE)
    } # end state
    if(biases[l] == "unmeasured_detection"){
      # switch relevant bias to true 
      sim = coabundance_simulator(nsites, nreps, narea, nyear, nsource, 
                                  true_siv = a5, log_dom = T, area = area, 
                                  Zdom_area = Zdom_area, Zsub_area = Zsub_area,
                                  bias_doublecount = TRUE, bias_sp_autocorr = FALSE, 
                                  bias_unmeas_state = FALSE, bias_unmeas_det = TRUE)
    } # end det
    if(biases[l] == "all_biases"){
      # switch relevant bias to true 
      sim = coabundance_simulator(nsites, nreps, narea, nyear, nsource, 
                                  true_siv = a5, log_dom = T, area = area, 
                                  Zdom_area = Zdom_area, Zsub_area = Zsub_area,
                                  bias_doublecount = TRUE, bias_sp_autocorr = FALSE, 
                                  bias_unmeas_state = FALSE, bias_unmeas_det = FALSE)
    } # end all 
    # save in the temp list 
    temp[[l]] = sim
    names(temp)[l] = paste("a5", a5, "bias", biases[l], sep = "_")
  } # end per bias 
  # save in the full list 
  data_list[[i]] = temp
}
rm(sim, temp, a5, area, i, l, narea, nreps, nsites, nsource, nyear, 
   sites_per_area, Zdom_area, Zsub_area)

## flatten all values into a single list
final_list = flatten(data_list)
# check that they are all here
names(final_list) # all good! 

## grab the date
day<-substr(Sys.Date(),9, 10)
month<-substr(Sys.Date(),6,7)
year<-substr(Sys.Date(),1,4)
date = paste(year, month, day, sep = "")
rm(day, month, year)

# make a path 
path = paste("~/Dropbox/Zach PhD/Ch3 Trophic release project/SEA_TC_GitHub_data_storage/data/step2_output_CoA_bundles/Bundled_data_for_bayes_co-abudance_mods_SIMULATION_", length(final_list), "_species_pairs_", 
             date, ".RDS", sep = "")

## save this as a RDS object
saveRDS(final_list, path)
  
