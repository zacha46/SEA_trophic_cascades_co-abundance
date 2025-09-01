### Co-abundance models - simulation study 

# Import libraries
library(ggplot2)
library(tidyverse)
library(jagsUI)

# Scenarios

# these are just ideas. We may need to structure main scenarios around common cam trap pitfalls such as double counting, overdispersion, autocorrelation, etc

# Scenario 1 - Top Down (these are actually implied by the direction of SIV)
# 1a - with strong covariate effects
# - varying intensity of SIV
# 1b - with weak covariate effects
# - varying intensity of SIV

# Scenario 2 - Bottom Up
# 2a - with strong covariate effects
# - varying intensity of SIV
# 2b - with weak covariate effects
# - varying intensity of SIV

# Scenario 3 - No Effects
# 3a - with strong covariate effects
# 3b - with weak covariate effects

# I think all 'expected' variation should just be across all scenarios e.g. random effects, some overdispersion is all fine to be included throughout

# Getting started with a framework

# For later inspection, set random seed
set.seed(2026)

# Bundled all previous work into a function. We construct count histories here
simulate_coabundance_matrix <- function(n_sites = 100, # Number of sampling units
                                        n_landscapes = 3, # Number of landscapes
                                        sites_per_landscape = c(40, 40, 20), # Allocation of SUs in landscapes
                                        range_land = matrix(c(1, 1, 0,
                                                             1, 0, 0), 
                                                            byrow = T,
                                                            nrow=2), # Species ranges per landscape
                                        n_species = 2, # Number of species - always 2 for co-abundance models
                                        n_rep = 10, # Number of sampling occasions or 'visits'
                                        psi = c(1, 1), # Base probability for a species to inhabit a site. If 1, defaults to the full range
                                        b0 = c(log(2), log(2)), # State formula intercept
                                        bFLLI = c(0.8, 0.8), # Forest Integrity Beta; good for both species
                                        bHFP = c(-0.7, -0.3), # HFP Beta; worse for predators, but bad for both
                                        bELEV = c(0.5, -0.2), # Elevation Beta; good for predators, a bit bad for prey usually
                                        bSIV = 0.5, # Species Interaction Value. Negative means Dom suppresses sub (top down) and vice versa is bottom up.
                                        p = c(0.6, 0.4),
                                        sd_landscapes = 0.3,
                                        unmeasured_SD = 0.1){
  
  # bring in some controlling logic
  stopifnot(n_sites == sum(sites_per_landscape))
  stopifnot(n_landscapes == length(sites_per_landscape))
  # Create the landscape attribute based on # sites
  landscape <- rep(1:n_landscapes, times = sites_per_landscape)
  
  # Introduce landscapes where species are extirpated (for IZIP modelling)
  is_in_range <- matrix(NA, nrow = n_sites, ncol = n_species)
  is_in_range[,1] <- range_land[1,landscape]
  is_in_range[,2] <- range_land[2,landscape]
  
  # Simulate presence and absence of species in each landscape
  z <- matrix(0, nrow = n_sites, ncol = n_species)
  for (s in 1:n_species) {
    z[, s] <- rbinom(n_sites, 1, psi[s]) * is_in_range[,s]
  }
  
  ##### Generate covariates
  
  # Landscapes random effects
  b_landscapes <- rnorm(n_landscapes, mean = 0, sd=sd_landscapes) 
  
  # A forest integrity - type covariate; good for both species
  # Forest Landscape Integrity Index (0–10 scaled) - the actual covariate values
  flli <- as.numeric(scale(pmin(pmax(runif(n_sites, 2, 9) + rnorm(n_sites, 0, 0.7), 0), 10))) # site-level spread
  
  # A hunting/poaching type covariate, maybe worse for the predator (e.g. tiger)
  # Human Footprint Index (0–50 scaled, but right-skewed)
  hfp <- as.numeric(scale(100 * rbeta(n_sites, shape1 = 3, shape2 = 6)))
  
  # A elevation type covariate - mixed results. good for predator bad for prey (just an example)
  # Elevation: mostly lowland, with random hilltops/uplands
  elev <- as.numeric(scale(ifelse(runif(n_sites) < 0.8,
                                  runif(n_sites, 0, 300),                    # 80% lowland
                                  ifelse(runif(n_sites) < 0.7,
                                         runif(n_sites, 300, 800),           # 14% mid-elev
                                         runif(n_sites, 800, 2000)))))         # 6% upland
  
  ##### Generate Lambdas and Ns
  lambda <- matrix(NA, n_sites, n_species)
  N <- matrix(0, n_sites, n_species)
  
  # Dominant - always sp 1
  lambda[,1] <- exp(b0[1] + bFLLI[1] * flli + bHFP[1] * hfp + bELEV[1] * elev + b_landscapes[landscape] + rnorm(n_sites, 0, unmeasured_SD)) 
  for(i in 1:n_sites){
    N[i,1] <- rpois(1, lambda[i,1] * z[i,1])
  }
  
  # Subordinate - always sp 2
  abu.dom <- log1p(N[,1])
  lambda[,2] <- exp(b0[2] + bFLLI[2] * flli + bHFP[2] * hfp + bELEV[2] * elev + bSIV * abu.dom + b_landscapes[landscape] + rnorm(n_sites, 0, unmeasured_SD))
  for(i in 1:n_sites){
    N[i,2] <- rpois(1, lambda[i,2] * z[i,2])
  }
  
  ##### Generate detections and repeated visits
  # Generate the fill detection matrix
  y <- array(NA, c(n_sites, n_rep, n_species))
  for (s in 1:n_species) {
    for (i in 1:n_sites) {
      for (j in 1:n_rep) {
        y[i,j,s] <- rbinom(1, N[i,s], p[s])
      }
    }
  }
  
  ##### Bundle up data for modelling
  sim <- list()
  
  sim$data <- list(
    n_sites = n_sites,
    nreps = n_rep,
    y.dom = y[,,1],
    y.sub = y[,,2],
    Z.dom = z[,1],
    Z.sub = z[,2],
    flii = flli,
    hfp = hfp,
    elev = elev,
    narea = n_landscapes,
    area = landscape
  )
  
  sim$true <- list(
    a0 = b0,
    a1 = bFLLI,
    a2 = bHFP,
    a3 = bELEV,
    a5 = bSIV,
    sigma.a6 = sd_landscapes,
    lambda.dom = lambda[,1],
    lambda.sub = lambda[,2],
    N.dom = N[,1],
    N.sub = N[,2]
  )
  
  return(sim)
} # Species specific detection probabilties

# Assign our count history matrix
sim <- simulate_coabundance_matrix()

##### Fit the model!

# Zach's co-abundance model altered for a simpler initial run.

# removed detection covariates and OD for now. 

Co_abundance_simple <- "simulations_vJMS/Co-abundance Simple.txt"
writeLines("
model{

  # Priors for both species
  for (i in 1:2) {
      
    # Intercept
    a0[i] ~ dnorm(0, 0.01)
      
    # FLII
    a1[i] ~ dnorm(0, 0.01)
      
    # HFP
    a2[i] ~ dnorm(0, 0.01)
    
    # Elev
    a3[i] ~ dnorm(0, 0.01)
  
  
    # Det intercept
    #b0[i] ~ dnorm(0, 0.01)
    
    # Cams
    #b2[i] ~ dnorm(0, 0.01)
      
    ## OD params for observation model
    #tau.p[i] <- pow(sd.p[i], -2) 
    #sd.p[i] ~ dunif(0, 1)  #not sure how to define variance here... leaving it the same as K&R
      
    }
    
  ## Species interaction
  a5 ~ dnorm(0, 0.01)
    
  # Landscape RE hyper prior --> define it's variance
  sigma.a6 ~ dunif(0,5)
  var.a6 <- 1/(sigma.a6*sigma.a6) 
    
  #sigma.a7 ~ dunif(0,5)
  #var.a7 <- 1/(sigma.a7*sigma.a7)
    
  for (k in 1:narea) {
      
    a6[k] ~ dnorm(0,var.a6)
    #a7[k] ~ dnorm(0,var.a7)
      
  }
  
  # year RE hyper prior --> define it's variance
  #sigma.a8 ~ dunif(0,5)
  #var.a8 <- 1/(sigma.a8*sigma.a8)

  #sigma.a9 ~ dunif(0,5)
  #var.a9 <- 1/(sigma.a9*sigma.a9)

  #for (k in 1:nyear) {

    #a8[k] ~ dnorm(0,var.a8)
    #a9[k] ~ dnorm(0,var.a9)

  #}
  
    # source RE hyper prior --> define it's variance
  #sigma.b3 ~ dunif(0,5)
  #var.b3 <- 1/(sigma.b3*sigma.b3)

  #sigma.b4 ~ dunif(0,5)
  #var.b4 <- 1/(sigma.b4*sigma.b4)

  #for (k in 1:nsource) {

    #b3[k] ~ dnorm(0,var.b3)
    #b4[k] ~ dnorm(0,var.b4)

 # }
 
 # Prior for P fixed
 alpha.p.sub ~ dnorm(0, 0.001)
 alpha.p.dom ~ dnorm(0, 0.001)
    
  # Likelihood
  # Ecological model for true abundance per site
  for (j in 1:n_sites) {
    
    # Abundance of Subordinate Species w/ iZIP
    N.sub[j] ~ dpois(lambda.sub[j] * Z.sub[j])

      log(lambda.sub[j]) <- a0[2] + a1[2]*flii[j] + a2[2]*hfp[j] + a3[2]*elev[j] + a5 * log(1+N.dom[j]) + a6[area[j]] 
    
    # Abundance of Dominant Species w/ iZIP 
    N.dom[j] ~ dpois(lambda.dom[j] * Z.dom[j])
    
      log(lambda.dom[j]) <- a0[1] + a1[1]*flii[j] + a2[1]*hfp[j] + a3[1]*elev[j] + a6[area[j]] 
                      
                      
    # Observation model for counts per replicated observation with OD params 
    for (k in 1:nreps) {
    
      ## Subordinate species
      y.sub[j,k] ~ dbin(p.sub[j,k], N.sub[j])
      
        p.sub[j,k] <- 1 / (1 + exp(-lp.lim.sub[j,k]))
        
          lp.lim.sub[j,k]<- min(250, max(-250, lp.sub[j,k])) #stabilize logit
          
             lp.sub[j,k] <- alpha.p.sub 
        
              #eps.p.sub[j,k] ~ dnorm(0, tau.p[1])
          
      
      ## Dominant species
      y.dom[j,k] ~ dbin(p.dom[j,k], N.dom[j])
      
        p.dom[j,k] <- 1 / (1 + exp(-lp.lim.dom[j,k]))
        
          lp.lim.dom[j,k]<- min(250, max(-250, lp.dom[j,k])) #stabilize logit
          
            lp.dom[j,k] <- alpha.p.dom
        
              #eps.p.dom[j,k] ~ dnorm(0, tau.p[2])
         
         
    ### PPC- Subordinate
      
    ## Expected count at site j, occasion k
    exp.sub[j,k]<- N.sub[j] * p.sub[j,k]
      
    ## Discrepency from Real Data
    # small denominator added to prevent dividing by zero
    E.sub[j,k] <- pow((y.sub[j,k] - exp.sub[j,k]), 2) / (exp.sub[j,k] + 0.5) 
      
    ## Simulate new sub counts from model 
    y.rep.sub[j,k] ~ dbin(p.sub[j,k], N.sub[j])
      
    ## Discrepency from Simulation
    E.rep.sub[j,k] <- pow((y.rep.sub[j,k] - exp.sub[j,k]), 2) / (exp.sub[j,k] + 0.5)
      
      
      
    ### PPC- Dominant
      
    ## Expected count at site j, occasion k
    exp.dom[j,k]<- N.dom[j] * p.dom[j,k]
      
    ## Discrepency from Real Data
    # small denominator added to prevent dividing by zero
    E.dom[j,k] <- pow((y.dom[j,k] - exp.dom[j,k]), 2) / (exp.dom[j,k] + 0.5) 
      
    ## Simulate new dom counts from model 
    y.rep.dom[j,k] ~ dbin(p.dom[j,k], N.dom[j])
      
    ## Discrepency from Simulation
    E.rep.dom[j,k] <- pow((y.rep.dom[j,k] - exp.dom[j,k]), 2) / (exp.dom[j,k] + 0.5)
      

    } #k
  } #j

  ## Derived Parameters
  
  #Chi-Square Test Statistic- Subordinate
  fit.sub = sum(E.sub[,])
  fit.rep.sub = sum(E.rep.sub[,])

  #Chi-Square Test Statistic- Dominant
  fit.dom = sum(E.dom[,])
  fit.rep.dom = sum(E.rep.dom[,])

} 
", con = Co_abundance_simple) # end model 

# Initial values
# --- safe initial values so y <= N holds at start ---
make_inits <- function(data_list){
  
  # max values observed in matricies 
  max_y_dom <- apply(data_list$y.dom, 1, max)
  max_y_sub <- apply(data_list$y.sub, 1, max)
  
  # ensure N starts >= observed maxima, and >=1 so log(N+eps) is finite
  Ndom_init <- pmax(as.integer(max_y_dom + 1L), 1L)
  Nsub_init <- pmax(as.integer(max_y_sub + 1L), 1L)
  
  # Force zero where Z==0
  Ndom_init <- data_list$Z.dom * Ndom_init
  Nsub_init <- data_list$Z.sub * Nsub_init
  
  # bundle in a list 
  return(list(
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
    alpha.p.sub = plogis(0),
    alpha.p.dom = plogis(0),
    # abundance inital values 
    N.dom = Ndom_init,
    N.sub = Nsub_init 
  ))
}

inits_list <- list(make_inits(sim$data), 
                   make_inits(sim$data),
                   make_inits(sim$data))

# 2) Run jagsUI with the file path
params <- c("a0","b0", "a5", "a1","a2","a3",#"a4","a5","b0","b1","sd.p", "N.dom", "N.sub",
            "sigma.a6", "lambda.sub", "lambda.dom", "N.dom", "N.sub") # no 'deviance' unless DIC=TRUE

# MCMC settings
ni = 6000; nb = 1000; nt = 10; nc = 3

## set parallell processing power
options(mc.cores = nc)

# 3) call the model 
mod <- jagsUI::jags(data   = sim$data,
                    inits  = inits_list,
                    parameters.to.save = params,
                    model.file = Co_abundance_simple,
                    n.chains = nc,
                    n.burnin = nb,
                    n.iter   = ni,  
                    n.thin   = nt,
                    parallel = TRUE,  
                    DIC      = FALSE)

##### Post-fitting analyses!

##### Plotting key parameters
extract_true_and_jags <- function(true_list, jags_out, param_names) {
  
  # --- Extract true values ---
  true_df <- do.call(rbind, lapply(param_names, function(nm) {
    if (!nm %in% names(true_list)) {
      stop("Parameter ", nm, " not found in the list")
    }
    
    vals <- true_list[[nm]]
    
    if (length(vals) == 2) {
      names_out <- paste0(nm, c("_dom", "_sub"))
    } else {
      names_out <- nm
    }
    
    data.frame(parameter = names_out, true = vals, stringsAsFactors = FALSE)
  }))
  
  rownames(true_df) <- NULL
  
  # --- Extract JAGS posterior summaries ---
  summary_table <- jags_out$summary
  jags_df <- do.call(rbind, lapply(true_df$parameter, function(pn) {
    if (grepl("_dom$", pn)) {
      jags_name <- sub("_dom$", "[1]", pn)
    } else if (grepl("_sub$", pn)) {
      jags_name <- sub("_sub$", "[2]", pn)
    } else {
      jags_name <- pn
    }
    
    if (!jags_name %in% rownames(summary_table)) {
      stop("Parameter ", jags_name, " not found in JAGS output")
    }
    
    vals <- summary_table[jags_name, c("2.5%", "50%", "97.5%")]
    data.frame(
      parameter = pn,
      lower = vals["2.5%"],
      median = vals["50%"],
      upper = vals["97.5%"],
      stringsAsFactors = FALSE
    )
  }))
  
  rownames(jags_df) <- NULL
  
  # --- Merge true and estimated ---
  merged_df <- merge(true_df, jags_df, by = "parameter")
  return(merged_df)
}

# Example usage
params <- c("a0", "a1", "a2", "a3", "a5", "sigma.a6", "sigma.b3")
df_compare <- extract_true_and_jags(sim$true, mod, params)

# Preserve the order in the dataframe (including duplicates)
df_compare$parameter <- fct_rev(fct_inorder(df_compare$parameter))

# Plot our estimates
ggplot(df_compare, aes(x = parameter)) +
  geom_linerange(aes(ymin = lower, ymax = upper), color = "skyblue", size = 1.5) +
  geom_point(aes(y = median), color = "blue", size = 2) +
  geom_point(aes(y = true), color = "red", shape = 18, size = 3) +
  coord_flip() +
  labs(
    y = "Parameter value",
    x = "",
    title = "JAGS posterior intervals vs true values"
  ) +
  theme_minimal(base_size = 14) + theme_bw()





# Bundled all previous work into a function. We construct count histories here
simulate_coabundance_matrix <- function(n_sites = 160, # Number of sampling units
                                        n_landscapes = 5, # Number of landscapes
                                        sites_per_landscape = c(40, 40, 20, 30, 30), # Allocation of SUs in landscapes
                                        range_land = matrix(c(1, 1, 0, 1, 1,
                                                              1, 1, 0, 1, 1), 
                                                            byrow = T,
                                                            nrow=2), # Species ranges per landscape
                                        n_species = 2, # Number of species - always 2 for co-abundance models
                                        n_rep = 10, # Number of sampling occasions or 'visits'
                                        psi = c(1, 1), # Base probability for a species to inhabit a site. If 1, defaults to the full range
                                        b0 = c(log(5), log(10)), # State formula intercept
                                        bFLLI = c(0.4, 0.4), # Forest Integrity Beta; good for both species
                                        bHFP = c(-0.2, -0.1), # HFP Beta; worse for predators, but bad for both
                                        bELEV = c(0.3, -0.2), # Elevation Beta; good for predators, a bit bad for prey usually
                                        bSIV = -2, # Species Interaction Value. Negative means Dom suppresses sub (top down) and vice versa is bottom up.
                                        sd_landscapes = 0.2,
                                        sd_source = 0.2,
                                        unmeasured_SD = 0.2,
                                        model_count_overdispersion = T){
  
  # bring in some controlling logic
  stopifnot(n_sites == sum(sites_per_landscape))
  stopifnot(n_landscapes == length(sites_per_landscape))
  # Create the landscape attribute based on # sites
  landscape <- rep(1:n_landscapes, times = sites_per_landscape)
  
  # Introduce landscapes where species are extirpated (for IZIP modelling)
  is_in_range <- matrix(NA, nrow = n_sites, ncol = n_species)
  is_in_range[,1] <- range_land[1,landscape]
  is_in_range[,2] <- range_land[2,landscape]
  
  # Simulate presence and absence of species in each landscape
  z <- matrix(0, nrow = n_sites, ncol = n_species)
  for (s in 1:n_species) {
    z[, s] <- rbinom(n_sites, 1, psi[s]) * is_in_range[,s]
  }
  
  ##### Generate covariates
  
  # Landscapes random effects
  b_landscapes <- rnorm(n_landscapes, mean = 0, sd=sd_landscapes) 
  
  # A forest integrity - type covariate; good for both species
  # Forest Landscape Integrity Index (0–10 scaled) - the actual covariate values
  flli <- as.numeric(scale(pmin(pmax(runif(n_sites, 2, 9) + rnorm(n_sites, 0, 0.7), 0), 10))) # site-level spread
  
  # A hunting/poaching type covariate, maybe worse for the predator (e.g. tiger)
  # Human Footprint Index (0–50 scaled, but right-skewed)
  hfp <- as.numeric(scale(100 * rbeta(n_sites, shape1 = 3, shape2 = 6)))
  
  # A elevation type covariate - mixed results. good for predator bad for prey (just an example)
  # Elevation: mostly lowland, with random hilltops/uplands
  elev <- as.numeric(scale(ifelse(runif(n_sites) < 0.8,
                                  runif(n_sites, 0, 300),                    # 80% lowland
                                  ifelse(runif(n_sites) < 0.7,
                                         runif(n_sites, 300, 800),           # 14% mid-elev
                                         runif(n_sites, 800, 2000)))))         # 6% upland
  
  ##### Generate Lambdas and Ns
  lambda <- matrix(NA, n_sites, n_species)
  N <- matrix(0, n_sites, n_species)
  
  
  # Dominant - always sp 1
  lambda[,1] <- exp(b0[1] + bFLLI[1] * flli + bHFP[1] * hfp + bELEV[1] * elev + b_landscapes[landscape] + rnorm(n_sites, 0, unmeasured_SD)) 
  for(i in 1:n_sites){
    if(model_count_overdispersion){
      N[i,1] <- rnbinom(1, mu = lambda[i,1] * z[i,1], size = 1.5)
    }else{
      N[i,1] <- rpois(1, lambda[i,1] * z[i,1])
    }
  }
  
  # Subordinate - always sp 2
  abu.dom <- log1p(N[,1])
  lambda[,2] <- exp(b0[2] + bFLLI[2] * flli + bHFP[2] * hfp + bELEV[2] * elev + bSIV * abu.dom + b_landscapes[landscape] + rnorm(n_sites, 0, unmeasured_SD))
  for(i in 1:n_sites){
    if(model_count_overdispersion){
    N[i,2] <- rnbinom(1, mu = lambda[i,2] * z[i,2], size = 1.5)
    }else{
      N[i,2] <- rpois(1, lambda[i,2] * z[i,2])
    }
  }
  
  ##### Generate detections and repeated visits
 
  a0 = c(0, 0)
  a1 = c(1,1)
  
  # Add some sources for our sites
  n_sources = 3
  source_per_site = c(50, 45, 65)
  stopifnot(sum(source_per_site) == n_sites)
  source <- rep(1:n_sources, times = source_per_site)
  a_source <- rnorm(n_sources, mean = 0, sd = sd_source)
  
  # Create the cams covariate
  cams_vec <- as.numeric(scale(pmin(rpois(n_sites, 2) + 1, 9)))
  cams <- matrix(rep(cams_vec, n_rep), nrow = n_sites, ncol = n_rep)
  p.dom <- matrix(NA, nrow = n_sites, ncol = n_rep)
  p.sub <- matrix(NA, nrow = n_sites, ncol = n_rep)
  # for each site 
  for(i in 1:n_sites){
    # and for each rep 
    for(j in 1:n_rep){
      eps_dom <- rnorm(1, 0, 0.1)
      eps_sub <- rnorm(1, 0, 0.1)
      ## Detection linear predictor is same for both species, w/ different intercepts/coefficents
      p.dom[i,j] = plogis(a0[1] + a1[1]*cams[i,j] + a_source[source[i]] + eps_dom)
      p.sub[i,j] = plogis(a0[2] + a1[2]*cams[i,j] + a_source[source[i]] + eps_sub)
    }
  }
  
  # Generate the fill detection matrix
  y.dom <- array(NA, c(n_sites, n_rep))
  y.sub <- array(NA, c(n_sites, n_rep))
    for (i in 1:n_sites) {
      for (j in 1:n_rep) {
        y.dom[i,j] <- rbinom(1, N[i,1], p.dom[i,j])
        y.sub[i,j] <- rbinom(1, N[i,2], p.sub[i,j])
      }
    }
  
  ##### Bundle up data for modelling
  sim <- list()
  
  sim$data <- list(
    n_sites = n_sites,
    nreps = n_rep,
    y.dom = y.dom,
    y.sub = y.sub,
    Z.dom = z[,1], #huh
    Z.sub = z[,2], #huh
    flii = flli,
    hfp = hfp,
    elev = elev,
    cams = cams,
    narea = n_landscapes,
    area = landscape,
    nsource = n_sources,
    source = source
    
  )
  
  sim$true <- list(
    a0 = b0,
    a1 = bFLLI,
    a2 = bHFP,
    a3 = bELEV,
    a5 = bSIV,
    b0 = a0,
    b2 = a1,
    sigma.a6 = sd_landscapes,
    sigma.b3 = sd_source,
    lambda.dom = lambda[,1],
    lambda.sub = lambda[,2],
    N.dom = N[,1],
    N.sub = N[,2],
    p.dom = p.dom,
    p.sub = p.sub
  )
  
  return(sim)
} 

sim <- simulate_coabundance_matrix()

Co_abundance_OD <- "simulations_vJMS/Co-abundance OD.txt"
writeLines("
model{

  # Priors for both species
  for (i in 1:2) {
      
    # Intercept
    a0[i] ~ dnorm(0, 0.01)
      
    # FLII
    a1[i] ~ dnorm(0, 0.01)
      
    # HFP
    a2[i] ~ dnorm(0, 0.01)
    
    # Elev
    a3[i] ~ dnorm(0, 0.01)
  
  
    # Det intercept
    b0[i] ~ dnorm(0, 0.01)
    
    # Cams
    b2[i] ~ dnorm(0, 0.01)
      
    ## OD params for observation model
    tau.p[i] <- pow(sd.p[i], -2) 
    sd.p[i] ~ dunif(0, 1)  #not sure how to define variance here... leaving it the same as K&R
      
    }
    
  ## Species interaction
  a5 ~ dnorm(0, 0.01)
    
  # Landscape RE hyper prior --> define it's variance
  sigma.a6 ~ dunif(0,5)
  var.a6 <- 1/(sigma.a6*sigma.a6) 
    
  #sigma.a7 ~ dunif(0,5)
  #var.a7 <- 1/(sigma.a7*sigma.a7)
    
  for (k in 1:narea) {
      
    a6[k] ~ dnorm(0,var.a6)
    #a7[k] ~ dnorm(0,var.a7)
      
  }
  
  # year RE hyper prior --> define it's variance
  #sigma.a8 ~ dunif(0,5)
  #var.a8 <- 1/(sigma.a8*sigma.a8)

  #sigma.a9 ~ dunif(0,5)
  #var.a9 <- 1/(sigma.a9*sigma.a9)

  #for (k in 1:nyear) {

    #a8[k] ~ dnorm(0,var.a8)
    #a9[k] ~ dnorm(0,var.a9)

  #}
  
    # source RE hyper prior --> define it's variance
  sigma.b3 ~ dunif(0,5)
  var.b3 <- 1/(sigma.b3*sigma.b3)

  #sigma.b4 ~ dunif(0,5)
  #var.b4 <- 1/(sigma.b4*sigma.b4)

  for (k in 1:nsource) {

    b3[k] ~ dnorm(0,var.b3)
    #b4[k] ~ dnorm(0,var.b4)

  }
 

  # Likelihood
  # Ecological model for true abundance per site
  for (j in 1:n_sites) {
    
    # Abundance of Subordinate Species w/ iZIP
    N.sub[j] ~ dpois(lambda.sub[j] * Z.sub[j])

      log(lambda.sub[j]) <- a0[2] + a1[2]*flii[j] + a2[2]*hfp[j] + a3[2]*elev[j] + a5 * log(1+N.dom[j]) + a6[area[j]] 
    
    # Abundance of Dominant Species w/ iZIP 
    N.dom[j] ~ dpois(lambda.dom[j] * Z.dom[j])
    
      log(lambda.dom[j]) <- a0[1] + a1[1]*flii[j] + a2[1]*hfp[j] + a3[1]*elev[j] + a6[area[j]] 
                      
                      
    # Observation model for counts per replicated observation with OD params 
    for (k in 1:nreps) {
    
      ## Subordinate species
      y.sub[j,k] ~ dbin(p.sub[j,k], N.sub[j])
      
        p.sub[j,k] <- 1 / (1 + exp(-lp.lim.sub[j,k]))
        
          lp.lim.sub[j,k]<- min(250, max(-250, lp.sub[j,k])) #stabilize logit
          
            lp.sub[j,k] <- b0[2] +  b2[2]*cams[j,k] + b3[source[j]] + eps.p.sub[j,k]
        
              eps.p.sub[j,k] ~ dnorm(0, tau.p[2])
          
      
      ## Dominant species
      y.dom[j,k] ~ dbin(p.dom[j,k], N.dom[j])
      
        p.dom[j,k] <- 1 / (1 + exp(-lp.lim.dom[j,k]))
        
          lp.lim.dom[j,k]<- min(250, max(-250, lp.dom[j,k])) #stabilize logit
          
            lp.dom[j,k] <- b0[1] + b2[1]*cams[j,k] + b3[source[j]] + eps.p.dom[j,k]
        
              eps.p.dom[j,k] ~ dnorm(0, tau.p[1])
         
         
    ### PPC- Subordinate
      
    ## Expected count at site j, occasion k
    exp.sub[j,k]<- N.sub[j] * p.sub[j,k]
      
    ## Discrepency from Real Data
    # small denominator added to prevent dividing by zero
    E.sub[j,k] <- pow((y.sub[j,k] - exp.sub[j,k]), 2) / (exp.sub[j,k] + 0.5) 
      
    ## Simulate new sub counts from model 
    y.rep.sub[j,k] ~ dbin(p.sub[j,k], N.sub[j])
      
    ## Discrepency from Simulation
    E.rep.sub[j,k] <- pow((y.rep.sub[j,k] - exp.sub[j,k]), 2) / (exp.sub[j,k] + 0.5)
      
      
      
    ### PPC- Dominant
      
    ## Expected count at site j, occasion k
    exp.dom[j,k]<- N.dom[j] * p.dom[j,k]
      
    ## Discrepency from Real Data
    # small denominator added to prevent dividing by zero
    E.dom[j,k] <- pow((y.dom[j,k] - exp.dom[j,k]), 2) / (exp.dom[j,k] + 0.5) 
      
    ## Simulate new dom counts from model 
    y.rep.dom[j,k] ~ dbin(p.dom[j,k], N.dom[j])
      
    ## Discrepency from Simulation
    E.rep.dom[j,k] <- pow((y.rep.dom[j,k] - exp.dom[j,k]), 2) / (exp.dom[j,k] + 0.5)
      

    } #k
  } #j

  ## Derived Parameters
  
  #Chi-Square Test Statistic- Subordinate
  fit.sub = sum(E.sub[,])
  fit.rep.sub = sum(E.rep.sub[,])

  #Chi-Square Test Statistic- Dominant
  fit.dom = sum(E.dom[,])
  fit.rep.dom = sum(E.rep.dom[,])

} 
", con = Co_abundance_OD) # end model 

# Initial values
# --- safe initial values so y <= N holds at start ---
make_inits <- function(data_list){
  
  # max values observed in matricies 
  max_y_dom <- apply(data_list$y.dom, 1, max)
  max_y_sub <- apply(data_list$y.sub, 1, max)
  
  # ensure N starts >= observed maxima, and >=1 so log(N+eps) is finite
  Ndom_init <- pmax(as.integer(max_y_dom + 1L), 1L)
  Nsub_init <- pmax(as.integer(max_y_sub + 1L), 1L)
  
  # Force zero where Z==0
  Ndom_init <- data_list$Z.dom * Ndom_init
  Nsub_init <- data_list$Z.sub * Nsub_init
  
  # bundle in a list 
  return(list(
    # key params
    a0 = rnorm(2, 0, 0.2),
    a5 = rnorm(1, 0, 0.2),
    # site covs 
    a1 = rnorm(2, 0, 0.2),
    a2 = rnorm(2, 0, 0.2),
    a3 = rnorm(2, 0, 0.2),
    a4 = rnorm(2, 0, 0.2),
    # det covs 
    alpha.p.sub = plogis(0),
    alpha.p.dom = plogis(0),
    b0 = rnorm(2, 0, 0.2),
    b2 = rnorm(2, 0, 0.2),
    sd.p = runif(2, 0, 0.3),
    # abundance inital values 
    N.dom = Ndom_init,
    N.sub = Nsub_init 
  ))
}

inits_list <- list(make_inits(sim$data), 
                   make_inits(sim$data),
                   make_inits(sim$data))

# 2) Run jagsUI with the file path
params <- c("a0","b0", "a5", "a1","a2","a3",#"a4","a5","b0","b1","sd.p", "N.dom", "N.sub",
            "sigma.a6", "sigma.b3", "lambda.sub", "y.dom", "y.sub", "p.dom", "p.sub", 
            "lambda.dom", "N.dom", "N.sub") # no 'deviance' unless DIC=TRUE

# MCMC settings
ni = 4000; nb = 1000; nt = 10; nc = 3

## set parallell processing power
options(mc.cores = nc)

# 3) call the model 
mod <- jagsUI::jags(data   = sim$data,
                    inits  = inits_list,
                    parameters.to.save = params,
                    model.file = Co_abundance_OD,
                    n.chains = nc,
                    n.burnin = nb,
                    n.iter   = ni,  
                    n.thin   = nt,
                    parallel = TRUE,  
                    DIC      = FALSE)



#### playing with spatial autocorrelation
gen_coords_by_landscape <- function(landscape, method = c("grid","random"),
                                    xlim = c(0,1), ylim = c(0,1),
                                    island_gap = 5){
  method <- match.arg(method)
  n_sites <- length(landscape)
  n_land <- length(unique(landscape))
  coords <- matrix(NA, nrow = n_sites, ncol = 2)
  
  for(L in unique(landscape)){
    idx <- which(landscape == L)
    nL <- length(idx)
    
    # coords within island L
    if(method == "grid"){
      nside <- ceiling(sqrt(nL))
      xs <- seq(xlim[1], xlim[2], length.out = nside)
      ys <- seq(ylim[1], ylim[2], length.out = nside)
      grid <- expand.grid(x = xs, y = ys)
      base_coords <- grid[seq_len(nL), ]
    } else {
      base_coords <- data.frame(
        x = runif(nL, xlim[1], xlim[2]),
        y = runif(nL, ylim[1], ylim[2])
      )
    }
    
    # shift island L by offset so it doesn’t overlap others
    offset_x <- (L-1) * island_gap
    coords[idx, ] <- cbind(base_coords$x + offset_x, base_coords$y)
  }
  
  colnames(coords) <- c("x","y")
  coords
}

coords <- gen_coords_by_landscape(sim$data$area,
                        method = 'grid')

# build dataframe for ggplot
df <- data.frame(
  x = coords[,1],
  y = coords[,2],
  landscape = factor(landscape)
)

# plot
ggplot(df, aes(x = x, y = y)) +
  geom_tile()+  # adjust binwidth for resolution
  facet_wrap(~landscape) +             # one panel per landscape (optional)
  coord_equal() +
  scale_fill_viridis_c() +
  theme_minimal() +
  labs(title = "Hex-binned site coordinates by landscape",
       x = "X", y = "Y", fill = "Count")

make_spatial_field_by_land <- function(coords, landscape, phi = 0.2, sigma = 1, scale = TRUE){
  # coords: data.frame/matrix with x,y for all sites
  # landscape: integer/factor vector length n_sites
  n <- nrow(coords)
  sp_effect <- numeric(n)
  for(L in unique(landscape)){
    idx <- which(landscape == L)
    if(length(idx) == 1){
      sp_effect[idx] <- 0
      next
    }
    sub_coords <- coords[idx, , drop=FALSE]
    d <- as.matrix(dist(sub_coords))
    Sigma <- (sigma^2) * exp(-d / phi)
    Sigma <- Sigma + diag(1e-6, nrow(Sigma))
    R <- chol(Sigma)
    fld <- as.numeric(t(R) %*% rnorm(length(idx)))
    if(scale) fld <- as.numeric(scale(fld))
    sp_effect[idx] <- fld
  }
  sp_effect
}
# example usage:
# coords <- gen_coords_by_landscape(landscape, "grid")
# spf <- make_spatial_field_by_land(coords, landscape, phi=0.15, sigma=0.8)
# lambda_log <- log(lambda[,1]) + uplift_strength * spf

spf <- make_spatial_field_by_land(coords, sim$data$area,
                                  phi = 0.5, sigma=1)


# Build dataframe for plotting
df_spf <- data.frame(
  x = coords[,1],
  y = coords[,2],
  landscape = factor(landscape),
  spf = spf
)

# Visualise as coloured points
ggplot(df_spf[1:80,], aes(x = x, y = y, fill = spf)) +
  geom_tile(size = 3) +
  facet_wrap(~landscape) +
  scale_fill_viridis_c(option = "plasma") +
  coord_equal() +
  theme_minimal() +
  labs(title = "Spatial autocorrelation field per landscape",
       colour = "Spatial effect")

library(ggplot2)
library(dplyr)

# Subset a single landscape
df1 <- df_spf %>% filter(landscape == 2)

# Plot as hex grid coloured by spf
ggplot(df1, aes(x = x, y = y, fill = spf)) +
  geom_tile() +  # aggregate spf per hex
  scale_fill_viridis_c(option = "plasma") +
  coord_equal() +
  theme_minimal() +
  labs(title = "Spatial field as hex grid (Landscape 1)",
       fill = "Spatial effect")


# Now we just paste this in to our sim function!

# test autocorr and then move to double counting






