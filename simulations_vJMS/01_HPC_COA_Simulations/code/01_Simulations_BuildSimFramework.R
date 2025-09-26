#### Script for creating the simulation scenario framework which is an input for the HPC analysis

# Import libraries
library(tidyverse)
library(jagsUI)

# Declare all utility functions in this file to promote readability

#### Utility functions for simulating bias processes.
# Declaring here to promote readability

# Generate grid coordinates
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

# Use the newly created coordinates to sample a spatial covariate
make_spatial_field_by_land <- function(coords, landscape, phi = 0.2, 
                                       sigma = 1, scale = TRUE, threshold =.75){
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
    
    # thresholding to keep only strong positive values
    fld[fld <= threshold] <- 0
    fld <- fld / 1
    sp_effect[idx] <- fld
  }
  sp_effect
}

# Add a spatial spillover effect to counts (pre-detection!)

### messy doodles matrix
apply_spillover_y <- function(y_mat,
                              coords,
                              landscape,
                              sites_per_landscape,
                              spillover_rate = 0.05,
                              #cutoff = 0.17,
                              kernel = c("binary", "exp"),
                              decay_scale = .1,
                              seed = NULL,
                              frac_exposed = 0.1) {
  if(!is.null(seed)) set.seed(seed)
  kernel <- match.arg(kernel)
  
  n_sites  <- nrow(y_mat)
  n_visits <- ncol(y_mat)
  
  n_landscapes <- length(unique(landscape))
  
  cut.off_length_landscape <- rep(sites_per_landscape, times = sites_per_landscape)
  cutoff <- 1 / (ceiling(sqrt(cut.off_length_landscape))-1)
  
  coords <- as.matrix(coords)
  storage.mode(coords) <- "numeric"
  
  # distance & mask
  dmat <- as.matrix(dist(coords))
  same_land <- outer(landscape, landscape, "==")
  
  # kernel for neighbor definition
  if(kernel == "binary") {
    K <- (dmat <= cutoff) * 1
  } else {
    K <- exp(-dmat / decay_scale)
  }
  diag(K) <- 0L
  K[!same_land] <- 0
  
  # normalize to probabilities for allocation
  row_sums <- rowSums(K)
  has_nb <- row_sums > 0
  W <- matrix(0, n_sites, n_sites)
  nz <- which(has_nb)
  if(length(nz) > 0) W[nz, ] <- K[nz, , drop = FALSE] / row_sums[nz]
  
  # outputs
  y_eff       <- matrix(0L, n_sites, n_visits)
  n_spill_out <- matrix(0L, n_sites, n_visits)
  n_spill_in  <- matrix(0L, n_sites, n_visits)
  
  for(v in seq_len(n_visits)) {
    yv <- as.integer(y_mat[, v])
    
    # --- draw site-specific exposure and number spilled but only for sites with neighbors ---
    # q_i <- rbeta(n_sites, frac_exposed * 1, (1 - frac_exposed) * 1)
    spill_events <- integer(n_sites)
    if(any(has_nb)) {
      # spill_events[has_nb] <- rbinom(sum(has_nb), size = yv[has_nb], prob = pmin(q_i[has_nb] * spillover_rate,1))
      
      spill_events[has_nb] <- rbinom(
        sum(has_nb),
        size = yv[has_nb],
        prob = frac_exposed * spillover_rate
      )
      
      
      
    }
    n_spill_out[, v] <- spill_events
    
    # --- reallocate all spilled individuals to neighbors (multinomial) ---
    spill_matrix <- matrix(0L, n_sites, n_sites)
    for(i in nz) {
      szi <- spill_events[i]
      if(szi > 0L) {
        spill_matrix[i, ] <- rmultinom(1, size = szi, prob = W[i, ])[, 1]
      }
    }
    imports <- colSums(spill_matrix)
    n_spill_in[, v] <- imports
    
    # sanity check: out == in for this visit (should always hold)
    if(sum(spill_events) != sum(imports)) {
      stop(sprintf("Spill conservation violated at visit %d: out=%d, in=%d",
                   v, sum(spill_events), sum(imports)))
    }
    
    # update abundances
    y_eff[, v] <- yv - spill_events + imports
    y_eff[, v][y_eff[, v] < 0] <- 0L
  }
  
  total_before <- sum(y_mat)
  total_after  <- sum(y_eff)
  if(total_before != total_after) {
    stop(sprintf("Total conservation violated: before=%d, after=%d", total_before, total_after))
  }
  
  list(y_eff = y_eff,
       n_spill_out = n_spill_out,
       n_spill_in = n_spill_in,
       W = W,
       K = K,
       total_before = total_before,
       total_after = total_after)
}

# Bundled all previous work into a function. We construct count histories here
simulate_coabundance_matrix <- function(n_sites = 160, # Number of sampling units
                                        n_landscapes = 5, # Number of landscapes
                                        sites_per_landscape = c(40, 40, 20, 30, 30), # Allocation of SUs in landscapes
                                        range_land = matrix(c(1, 1, 0, 0, 1,
                                                              1, 0, 0, 1, 1), 
                                                            byrow = T,
                                                            nrow = 2), # Species ranges per landscape
                                        n_species = 2, # Number of species - always 2 for co-abundance models
                                        n_rep = 10, # Number of sampling occasions or 'visits'
                                        psi = c(1, 1), # Base probability for a species to inhabit a site. If 1, defaults to the full range. This actually helps us model Z-Inflation... another case!
                                        b0 = c(log(1.5), log(15)), # State formula intercept
                                        bFLLI = c(0.3, 0.3), # Forest Integrity Beta; good for both species
                                        bHFP = c(-0.4, -0.1), # HFP Beta; worse for predators, but bad for both
                                        bELEV = c(0.2, -0.2), # Elevation Beta; good for predators, a bit bad for prey usually
                                        bSIV = -1.5, # Species Interaction Value. Negative means Dom suppresses sub (top down) and vice versa is bottom up.
                                        sd_landscapes = 0.2, # RE
                                        sd_source = 0.2, # RE
                                        unmeasured_SD = 0.2, # Unmeaasured covariate/variation
                                        model_count_overdispersion = F, OD.size = 2,
                                        model_spatial_autocorrelation = T, mag.spatial = 1,
                                        model_double_counting = F, double_rate = 0.5,
                                        model_spatial_spillover = F, spillover_rate = .5, fraction_exposed = .15){
  
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
  
  # Generate species detections. If these are spatially linked we replace downstream. This is just filler
  species.detections <- scale(runif(n_sites, 1, 40))[,1]
  
  ##### Generate Lambdas and Ns
  lambda <- matrix(NA, n_sites, n_species)
  N <- matrix(0, n_sites, n_species)
  
  # Set up a spatial covariate effect for spatial autocorrelation
  # It will be 0 unless we trigger the setting
  spatial_effect <- rep(0, n_sites)
  
  coords <- gen_coords_by_landscape(landscape,
                                    method = 'grid')
  
  # Assign X and Y coordinates in an imaginary CRS, for use later in residual analysis
  x = coords[,1]
  y = coords[,2]
  
  if(model_spatial_autocorrelation){
    # Create the spatial covariate through Cholesky sampling from a Gaussian Random Field 
    # Going to leave these parameters set because they work well
    spatial_effect <- make_spatial_field_by_land(coords, landscape,
                                                 phi = 10, sigma=6.5) 
    
    # come back here if you eant to add a size parameter controlling the mangitude of uplift from spatial effect
    # also come back if we want to simulate something a different effect for the other species, currently they are same
    
    # We are now going to draw some species counts using our spatial factor and some environmental covariates as a basis
    # we assume the spatial effect 
    
    # Species specific responses
    S <- 40
    positive.bias <- .6
    n_pos <- round(positive.bias*S)
    n_neg <- S - n_pos
    
    # draw near-1/-1 values with small spread
    beta_s <- c(
      rnorm(n_pos, mean = 1, sd = 0.2),   # 70% positive
      rnorm(n_neg, mean = -1, sd = 0.2)   # 30% negative
    )
    beta_s <- sample(beta_s)  # shuffle
    
    species.detections <- scale(rowSums(matrix(rbinom(n_sites*40, 1, plogis(spatial_effect %*% t(beta_s))), n_sites, 40)))[,1]
    
  }
  
  
  # Dominant - always sp 1
  lambda[,1] <- exp(b0[1] + bFLLI[1] * flli + bHFP[1] * hfp + bELEV[1] * elev + b_landscapes[landscape] + rnorm(n_sites, 0, unmeasured_SD) + mag.spatial * spatial_effect) 
  for(i in 1:n_sites){
    if(model_count_overdispersion){
      N[i,1] <- rnbinom(1, mu = lambda[i,1] * z[i,1], size = OD.size)
    }else{
      N[i,1] <- rpois(1, lambda[i,1] * z[i,1])
    }
  }
  
  # Subordinate - always sp 2

  abu.dom <- log1p(N[,1])
  lambda[,2] <- exp(b0[2] + bFLLI[2] * flli + bHFP[2] * hfp + bELEV[2] * elev + bSIV * abu.dom + b_landscapes[landscape] + rnorm(n_sites, 0, unmeasured_SD) + mag.spatial * spatial_effect)
  for(i in 1:n_sites){
    if(model_count_overdispersion){
      N[i,2] <- rnbinom(1, mu = lambda[i,2] * z[i,2], size = OD.size)
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
  
  # For modelling spatial spillover. Important to put this before the double counting process,
  # as spillover is a 'real' ecological process and double counting is just human detection error!
  
  if(model_spatial_spillover){
    
    # bringing in our handy spillover function here. Need to call twice and spillover will be different
    # per species, which makes good sense ecologically.
    
    # fill out grid for dom
    N.matrix.temp.dom <- matrix(N[,1], nrow = n_sites, ncol = n_rep)
    
    # once over for Dominant species
    temp.spill.dom <- apply_spillover_y(y_mat = N.matrix.temp.dom,
                                        coords = coords,
                                        landscape = landscape,
                                        sites_per_landscape = sites_per_landscape,
                                        spillover_rate = spillover_rate,
                                        kernel = "binary", 
                                        decay_scale = 0.1,
                                        seed = NULL,
                                        frac_exposed = fraction_exposed)
    
    # feed back into count history
    temp.N.dom <- temp.spill.dom$y_eff
    
    spillout.matrix.dom <- temp.spill.dom$n_spill_out
    spillin.matrix.dom <- temp.spill.dom$n_spill_in
    
    # fill out grid for sub
    N.matrix.temp.sub <- matrix(N[,2], nrow = n_sites, ncol = n_rep)
    
    # once more for subordinate
    temp.spill.sub <- apply_spillover_y(y_mat = N.matrix.temp.sub,
                                        coords = coords,
                                        landscape = landscape,
                                        sites_per_landscape = sites_per_landscape,
                                        spillover_rate = spillover_rate,
                                        kernel = "binary", 
                                        decay_scale = 0.1,
                                        seed = NULL,
                                        frac_exposed = fraction_exposed)
    
    # feed back into count history
    temp.N.sub <- temp.spill.sub$y_eff
    spillout.matrix.sub <- temp.spill.sub$n_spill_out
    spillin.matrix.sub <- temp.spill.sub$n_spill_in
    
    # Finally, redraw y's given the new spillover N's
    for (i in 1:n_sites) {
      for (j in 1:n_rep) {
        #Will override previous y's
        y.dom[i,j] <- rbinom(1, temp.N.dom[i,1], p.dom[i,j])
        y.sub[i,j] <- rbinom(1, temp.N.sub[i,2], p.sub[i,j])
      }
    }
  }
  
  # Model double counting
  if(model_double_counting){
    for(i in 1:n_sites){
      for(j in 1:n_rep){
        y.dom[i,j] <- y.dom[i,j] + rbinom(1, y.dom[i,j], double_rate)
        y.sub[i,j] <- y.sub[i,j] + rbinom(1, y.sub[i,j], double_rate)
      }
    }
  }
  
  
  ##### Bundle up data for modelling
  sim <- list()
  
  sim$data <- list(
    n_sites = n_sites,
    nreps = n_rep,
    y.dom = y.dom,
    y.sub = y.sub,
    Z.dom = is_in_range[,1],
    Z.sub = is_in_range[,2],
    flii = flli,
    hfp = hfp,
    elev = elev,
    cams = cams,
    narea = n_landscapes,
    area = landscape,
    nsource = n_sources,
    source = source,
    comm_det = species.detections
    
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
    p.sub = p.sub,
    x = x,
    y = y,
    spatial_effect = spatial_effect,
    spillout.matrix.dom = if(exists("spillout.matrix.dom")) spillout.matrix.dom else NA ,
    spillin.matrix.dom = if(exists("spillin.matrix.dom")) spillin.matrix.dom else NA,
    spillin.matrix.sub = if(exists("spillin.matrix.sub")) spillin.matrix.sub else NA,
    spillout.matrix.sub =  if(exists("spillout.matrix.sub")) spillout.matrix.sub else NA
  )
  
  return(sim)
} 


#### Create Simulation Parameter Grid

# Scenario 1. Base (but with detection overdispersion)
## --- State (Ecological) Process Related
# Scenario 1 - Base case w/ OD in detection
# Sceanrio 2 - Base case w/ OD in detection + Unmeasured SD
# Scenario 3 - Base case w/ OD in detection + OD in state
# Scenario 4 - Base case w/ OD in detection + Unmeasured spatial variable

## --- Detection Process Related
# Scenario 5 - Base case w/ OD in detection + Double Counting
# Scenario 6 - Base case w/ OD in detection + Spatial Spillover process
scenario.df <- data.frame(scenario = c(1,"2a", "2b","3a", "3b","4a", "4b","5a", "5b","6a", "6b"),
                          scenario_name = c('Base',
                                       'Unmeasured Variation - Low',
                                       'Unmeasured Variation - High',
                                       'State Process Overdispersion - Low',
                                       'State Process Overdispersion - High',
                                       'Unmeasured Spatial Covariate - Low',
                                       'Unmeasured Spatial Covariate - High',
                                       'Double Counting - Low',
                                       'Double Counting - High',
                                       'Spatial Spillover - Low',
                                       'Spatial Spillover - High'),
                          unmeasured_SD = as.numeric(c(0,0.2, 0.4,0,0,0,0,0,0,0,0)),
                          model_count_overdispersion = c(F,F,F,T,T,F,F,F,F,F,F), 
                          OD.size = c(0,0,0,2,1.5,0,0,0,0,0,0),
                          model_spatial_autocorrelation = c(F,F,F,F,F,T,T,F,F,F,F),
                          mag.spatial = c(0,0,0,0,0,0.5,1.5,0,0,0,0),
                          model_double_counting = c(F,F,F,F,F,F,F,T,T,F,F),
                          double_rate = c(0,0,0,0,0,0,0,0.2,0.5,0,0),
                          model_spatial_spillover = c(F,F,F,F,F,F,F,F,F,T,T),
                          spillover_rate = c(0,0,0,0,0,0,0,0,0,0.2,0.9),
                          fraction_exposed = c(0,0,0,0,0,0,0,0,0,0.15,0.15)
                          )

# Create the SIV values we need to test
SIV.sweep <- data.frame("SIV" = seq(-1, 1, 0.1))

# Set base abundances with respect to the role they play
# When testing for bottom up, we expect there to be more prey
# When testing for top down, we expect there to be more prey
# But because sub and dom flips, we need to do this
SIV.sweep$b0.dom <- ifelse(SIV.sweep$SIV <= 0, log(0.5), log(10))
SIV.sweep$b0.sub <- ifelse(SIV.sweep$SIV <= 0, log(10), log(0.5))
SIV.sweep$bFLLI.dom <- ifelse(SIV.sweep$SIV <= 0, 0.3, 0.3)
SIV.sweep$bHFP.dom <- ifelse(SIV.sweep$SIV <= 0, -0.4, -0.1)
SIV.sweep$bELEV.dom <- ifelse(SIV.sweep$SIV <= 0, 0.2, -0.2)
SIV.sweep$bFLLI.sub <- ifelse(SIV.sweep$SIV <= 0, 0.3, 0.3)
SIV.sweep$bHFP.sub <- ifelse(SIV.sweep$SIV <= 0, -0.1, -0.4)
SIV.sweep$bELEV.sub<- ifelse(SIV.sweep$SIV <= 0, -0.2, 0.2)

# Create the number of replicates, given we have stochasticity
simulation.replicates <- data.frame("replicate" = 1:10)

# Expand the grid with replicates and the SIV
simulation.parameter.grid <- merge(merge(scenario.df, SIV.sweep), simulation.replicates) 

# Add an index column for downstream analyses

simulation.parameter.grid$index <- 1:nrow(simulation.parameter.grid)

# Save file for downstream analyses
saveRDS(simulation.parameter.grid, "simulations_vJMS/simulation_parameter_grid.rds")


# keep environment clean
rm(simulation.replicates, SIV.sweep, scenario.df)

##### Step 2. Create simulated datasets from our grid

simData <- apply(simulation.parameter.grid, 1, function(row) {
  
  # Visual tracker
  print(as.numeric(row['index']))
  
  # We use default values unless scenario inputs the specific params.
  # We are less interested in the standard params (sites etc) so have chosen
  # justifiable values for those, and we vary the params that make a difference
  # for the scenario
  
  simdat <- simulate_coabundance_matrix(unmeasured_SD = as.numeric(row['unmeasured_SD']),
                                        model_count_overdispersion = as.logical(row['model_count_overdispersion']), 
                                        OD.size = as.numeric(row['OD.size']),
                                        model_spatial_autocorrelation = as.logical(row['model_spatial_autocorrelation']),
                                        mag.spatial = as.numeric(row['mag.spatial']),
                                        model_double_counting = as.logical(row['model_double_counting']),
                                        double_rate = as.numeric(row['double_rate']),
                                        model_spatial_spillover = as.logical(row['model_spatial_spillover']),
                                        spillover_rate = as.numeric(row['spillover_rate']),
                                        fraction_exposed = as.numeric(row['fraction_exposed']),
                                        bSIV = as.numeric(row['SIV']),
                                        b0 = c(as.numeric(row['b0.dom']), as.numeric(row['b0.sub'])),
                                        bFLLI = c(as.numeric(row['bFLLI.dom']), as.numeric(row['bFLLI.sub'])),
                                        bHFP = c(as.numeric(row['bHFP.dom']), as.numeric(row['bHFP.sub'])),
                                        bELEV = c(as.numeric(row['bELEV.dom']), as.numeric(row['bELEV.sub']))
  )
  
  # Return our simulation data
  return(simdat)
})

names(simData) <- seq_len(nrow(simulation.parameter.grid))

# Write and Save output files for HPC analysis.

# Save file for downstream analyses
saveRDS(simData,  "simulations_vJMS/sim_data_list.rds")

# END











