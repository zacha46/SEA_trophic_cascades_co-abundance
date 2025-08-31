### Co-abundance models - simulation study 

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
set.seed(2025)

# Set # sites, i.e. camera
n_sites <- 100

# Set landscapes e.g. 'Singapore'
n_landscapes <- 3

# Allocate the cameras to our landscapes
sites_per_landscape <- c(40, 40, 20)

# bring in some controlling logic
stopifnot(n_sites == sum(sites_per_landscape))

# Create the landscape attribute based on # sites
landscape <- rep(1:3, times = sites_per_landscape)

# Introduce landscapes where species are extirpated (for IZIP modelling)
range_land <- c(1, 1, 0)
is_in_range <- range_land[landscape]

# Number of species - always 2 for now.
n_species <- 2

# Number of visits or sampling occasions. Start with 10
n_rep <- 10

# Simulate presence and absence of species in each landscape
psi <- c(0.7, 0.7)

z <- matrix(0, nrow = n_sites, ncol = n_species)
for (s in 1:n_species) {
  z[, s] <- rbinom(n_sites, 1, psi[s]) * is_in_range
}

##### Generate covariates

# not all landscapes are equal and some have higher base abundances due to unmeasured variation
b0 <- rep(log(2), n_species)

# Landscapes random effects
b_landscapes <- rnorm(n_landscapes, mean = 0, sd=0.3) 

# A forest integrity - type covariate; good for both species
bFLLI <- rep(1.5, n_species)

# Forest Landscape Integrity Index (0–10 scaled) - the actual covariate values
flli <- as.numeric(scale(pmin(pmax(runif(n_sites, 2, 9) + rnorm(n_sites, 0, 0.7), 0), 10))) # site-level spread

# A hunting/poaching type covariate, maybe worse for the predator (e.g. tiger)
bHFP <- c(-1.5, -0.3)

# Human Footprint Index (0–50 scaled, but right-skewed)
hfp <- as.numeric(scale(100 * rbeta(n_sites, shape1 = 3, shape2 = 6)))

# A elevation type covariate - mixed results. good for predator bad for prey (just an example)
bELEV <- c(1, -1)

# Elevation: mostly lowland, with random hilltops/uplands
elev <- as.numeric(scale(ifelse(runif(n_sites) < 0.8,
               runif(n_sites, 0, 300),                    # 80% lowland
               ifelse(runif(n_sites) < 0.7,
                      runif(n_sites, 300, 800),           # 14% mid-elev
                      runif(n_sites, 800, 2000)))))         # 6% upland


# Generate SIV. This should be a single, 'true' parameter that gets muddied by the rest. But we assume it is a natural 'law' of sorts.
bSIV <- -2

##### Generate Lambdas and Ns
lambda <- matrix(NA, n_sites, n_species)
N <- matrix(0, n_sites, n_species)

# Dominant - always sp 1
lambda[,1] <- exp(b0[1] + bFLLI[1] * flli + bHFP[1] * hfp + bELEV[1] * elev + b_landscapes[landscape])
for(i in 1:n_sites){
  N[i,1] <- rpois(1, lambda[i,1] * z[i,1])
}

# Subordinate - always sp 2
abu.dom <- log1p(N[,1])
lambda[,2] <- exp(b0[2] + bFLLI[2] * flli + bHFP[2] * hfp + bELEV[2] * elev + bSIV * abu.dom + b_landscapes[landscape])
for(i in 1:n_sites){
  N[i,2] <- rpois(1, lambda[i,2] * z[i,2])
}

##### Generate detections and repeated visits

# Set detection proabilities - vanilla to start, add covariates later
p <- c(0.6, 0.4)

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

##### Fit the model!

# Zach's co-abundance model altered for a simpler initial run.

# removed detection covariates and OD for now. 

Co_abundance_simple <- "Co-abundance Simple.txt"
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
            "a6", "lambda.sub", "lambda.dom", "N.dom", "N.sub") # no 'deviance' unless DIC=TRUE

# MCMC settings
ni = 4000; nb = 500; nt = 5; nc = 3

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

extract_jags_param <- function(jags_out, param_name) {
  # jags_out = object returned by jagsUI
  # param_name = string, e.g., "a1[1]"
  
  summary_table <- jags_out$summary  # contains mean, sd, 2.5%, 50%, 97.5%, etc.
  
  if (!param_name %in% rownames(summary_table)) {
    stop(paste("Parameter", param_name, "not found in JAGS output"))
  }
  
  param_summary <- summary_table[param_name, c("2.5%", "50%", "97.5%")]
  names(param_summary) <- c("lower", "median", "upper")
  
  return(param_summary)
}

extract_jags_param(mod,'a5')

plot_jags_param <- function(jags_out, param_names, true_values, main=NULL) {
  library(ggplot2)
  
  # Extract summaries
  summaries <- t(sapply(param_names, function(p) extract_jags_param(jags_out, p)))
  
  df <- data.frame(
    param = factor(param_names, levels=param_names),
    lower = summaries[,"lower"],
    median = summaries[,"median"],
    upper = summaries[,"upper"],
    true = true_values
  )
  
  ggplot(df, aes(x=param, y=median)) +
    geom_point(color="blue", size=3) +
    geom_errorbar(aes(ymin=lower, ymax=upper), width=0.2, color="blue") + ylim(3, -3)+
    geom_point(aes(y=true), color="red", shape=18, size=3) +
    ylab("Parameter value") +
    xlab("Parameter") +
    ggtitle(ifelse(is.null(main), "JAGS parameter estimates", main)) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle=45, hjust=1))
}

plot_jags_param(mod, 'a5', true_values = -2)


