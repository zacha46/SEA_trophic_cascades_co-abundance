### Co-abundance models - simulation study 

# Scenarios

# Scenario 1 - Top Down 
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
b_landscapes <- rnorm(n_landscapes, mean = 0, sd=0.3) # create i dimension

# A forest integrity - type covariate; good for both species
bFLII <- rep(1.5, n_species)

# A hunting/poaching type covariate, maybe worse for the predator (e.g. tiger)
bHFP <- c(-1.5, -0.3)

# A elevation type covariate - mixed results. good for predator bad for prey (just an example)
bELEV <- c(1, -1)

# Generate SIV. This should be a single, 'true' parameter that gets muddied by the rest. But we assume it is a natural 'law' of sorts.
bSIV <- 2

##### Generate Lambdas and Ns
lambda <- matrix(NA, n_sites, n_species)
N <- matrix(0, n_sites, n_species)

# Dominant - always sp 1
lambda[1] <- b0[1] + bFLLI[1] * flli + bHFP[1] * hfp + bELEV[1] * elev + b_landscapes[xxx] 
N[i,1] <- rpois(1, lambda[i,1] * z[i,1])

# Subordinate - always sp 2
lambda[2] <- b0[2] + bFLLI[2] * flli + bHFP[2] * hfp + bELEV[2] * elev + bSIV * N[1] + b_landscapes[xxx] 
N[i,2] <- rpois(1, lambda[i,2] * z[i,2])

##### Generate detections and repeated visits

# Set detection proabilities - vanilla to start
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


##### Fit the model!




