### Co-abundance models - simulation study 

# Import libraries
library(tidyverse)
library(jagsUI)

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
  print(total_after)
  print(total_before)
  
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
                                        b0 = c(log(2), log(10)), # State formula intercept
                                        bFLLI = c(0.4, 0.4), # Forest Integrity Beta; good for both species
                                        bHFP = c(-0.2, -0.1), # HFP Beta; worse for predators, but bad for both
                                        bELEV = c(0.3, -0.2), # Elevation Beta; good for predators, a bit bad for prey usually
                                        bSIV = -1.5, # Species Interaction Value. Negative means Dom suppresses sub (top down) and vice versa is bottom up.
                                        sd_landscapes = 0.2, # RE
                                        sd_source = 0.2, # RE
                                        unmeasured_SD = 0.2, # Unmeaasured covariate/variation
                                        model_count_overdispersion = F,
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
      N[i,1] <- rnbinom(1, mu = lambda[i,1] * z[i,1], size = 2)
    }else{
      N[i,1] <- rpois(1, lambda[i,1] * z[i,1])
    }
  }
  
  # Subordinate - always sp 2
  abu.dom <- log1p(N[,1])
  lambda[,2] <- exp(b0[2] + bFLLI[2] * flli + bHFP[2] * hfp + bELEV[2] * elev + bSIV * abu.dom + b_landscapes[landscape] + rnorm(n_sites, 0, unmeasured_SD) + mag.spatial * spatial_effect)
  for(i in 1:n_sites){
    if(model_count_overdispersion){
      N[i,2] <- rnbinom(1, mu = lambda[i,2] * z[i,2], size = 2)
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

# Create initial values for JAGS model
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
    a4 = rnorm(1, 0, 0.2),
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

# Function for extracting and plotting summary statistics from jags objects
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


# Scenarios

# these are just ideas. We may need to structure main scenarios around common cam trap pitfalls such as double counting, overdispersion, autocorrelation, etc

## --- State (Ecological) Process Related
# Scenario 1 - Base case w/ OD in detection
# Sceanrio 2 - Base case w/ OD in detection + Unmeasured SD
# Scenario 3 - Base case w/ OD in detection + OD in state
# Scenario 4 - Base case w/ OD in detection + Unmeasured spatial variable

## --- Detection Process Related
# Scenario 5 - Base case w/ OD in detection + Double Counting
# Scenario 6 - Base case w/ OD in detection + Spatial Spillover process

# Getting started with a framework

# For later inspection, set random seed
set.seed(2026)

# Zach's model
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
    
  ## Community detections
  a4 ~ dnorm(0, 1)  
    
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

      log(lambda.sub[j]) <- a0[2] + a1[2]*flii[j] + a2[2]*hfp[j] + a3[2]*elev[j] + a4 * comm_det[j]  + a5 * log(1+N.dom[j]) + a6[area[j]] 
    
    # Abundance of Dominant Species w/ iZIP 
    N.dom[j] ~ dpois(lambda.dom[j] * Z.dom[j])
    
      log(lambda.dom[j]) <- a0[1] + a1[1]*flii[j] + a2[1]*hfp[j] + a3[1]*elev[j] + a4 * comm_det[j] + a6[area[j]] 
                      
                      
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

# Simulate the data
sim <- simulate_coabundance_matrix()

inits_list <- list(make_inits(sim$data), 
                   make_inits(sim$data),
                   make_inits(sim$data))

# 2) Run jagsUI with the file path
params.mod <- c("a0","b0", "a5", "a1","a2","a3","a4", #,"a5","b0","b1","sd.p", "N.dom", "N.sub",
            "sigma.a6", "sigma.b3", "lambda.sub", "y.dom", "y.sub", "p.dom", "p.sub", 
            "lambda.dom", "N.dom", "N.sub") # no 'deviance' unless DIC=TRUE

# MCMC settings
ni = 4000; nb = 1000; nt = 10; nc = 3

## set parallell processing power
options(mc.cores = nc)

# 3) call the model 
mod <- jagsUI::jags(data   = sim$data,
                    inits  = inits_list,
                    parameters.to.save = params.mod,
                    model.file = Co_abundance_OD,
                    n.chains = nc,
                    n.burnin = nb,
                    n.iter   = ni,  
                    n.thin   = nt,
                    parallel = TRUE,  
                    DIC      = FALSE)

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










for(i in 1:10){
  lgcol<-ifelse(sim$true$spillout.matrix.sub[,i] | sim$true$spillin.matrix.sub[,i]  ,adjustcolor("red", alpha.f = 0.5), adjustcolor("grey50", alpha.f = 0.5))
  plot(sim$true$N.sub+1 ~ sim$true$N.sub, col = lgcol, pch=19, type='n',
       xlab = 'True N before Spillover',
       ylab = "True N after Spillover")
  
  inpoints <- sim$true$spillin.matrix.sub[,i]
  outpoints <- sim$true$spillout.matrix.sub[,i]
  
  points((inpoints + sim$true$N.sub - outpoints)  ~sim$true$N.sub, col = lgcol, pch=19)
  Sys.sleep(1)
}






# Build dataframe for plotting
df_spf <- data.frame(
  x = sim$true$x,
  y = sim$true$y,
  landscape = factor(sim$data$area),
  spf = sim$true$spatial_effect,
  caps = sim$data$comm_det
)

# Visualise as coloured points
ggplot(df_spf, aes(x = x, y = y)) +
  geom_tile(size = 3) +
  facet_wrap(~landscape) +
  scale_fill_viridis_c(option = "plasma") +
  coord_equal() +
  theme_minimal() +
  labs(title = "Spatial autocorrelation field per landscape",
       colour = "Spatial effect")

# Subset a single landscape
df1 <- df_spf %>% filter(landscape == 1)

# Plot as hex grid coloured by spf
ggplot(df1, aes(x = x, y = y, fill = spf)) +
  geom_tile() +  # aggregate spf per hex
  geom_text(aes(label = round(caps,2) ), color = "white", size = 3)+
  scale_fill_viridis_c(option = "plasma") +
  coord_equal() +
  theme_minimal() +
  labs(title = "Spatial field as hex grid (Landscape 1)",
       fill = "Spatial effect")




apply_spillover_y <- function(y_mat,
                              coords,
                              landscape,
                              spillover_rate = 0.05,
                              cutoff = 0.17,
                              kernel = c("binary", "exp"),
                              decay_scale = 1,
                              seed = NULL) {
  if(!is.null(seed)) set.seed(seed)
  kernel <- match.arg(kernel)
  
  n_sites  <- nrow(y_mat)
  n_visits <- ncol(y_mat)
  
  coords <- as.matrix(coords)
  storage.mode(coords) <- "numeric"
  
  # distance matrix
  dmat <- as.matrix(dist(coords))
  
  # same-landscape mask
  same_land <- outer(landscape, landscape, "==")
  
  # kernel
  if(kernel == "binary") {
    K <- (dmat <= cutoff) * 1
  } else {
    K <- exp(-dmat / decay_scale)
  }
  diag(K) <- 0L
  K[!same_land] <- 0
  
  # normalize rows
  row_sums <- rowSums(K)
  row_sums[row_sums == 0] <- 1
  W <- K / row_sums
  W[row_sums == 0, ] <- 0
  
  # initialize outputs
  y_eff       <- y_mat
  n_spill_out <- matrix(0L, n_sites, n_visits)  # how many individuals exported per site
  n_spill_in  <- matrix(0L, n_sites, n_visits)  # how many received
  neigh_abund <- matrix(0, n_sites, n_visits)
  
  # loop over visits
  for(v in 1:n_visits) {
    # neighbor abundance (used to scale expected spill)
    neighbor_counts <- as.numeric(W %*% y_mat[, v])
    neigh_abund[, v] <- neighbor_counts
    
    # total spill export per site
    lambda_vec <- spillover_rate * neighbor_counts
    export_draw <- rpois(n_sites, lambda_vec)
    n_spill_out[, v] <- export_draw
    
    # allocate exports among neighbors
    spill_matrix <- matrix(0, nrow = n_sites, ncol = n_sites)
    for(i in 1:n_sites) {
      if(export_draw[i] > 0) {
        probs <- W[i, ]
        if(sum(probs) > 0) {
          spill_matrix[i, ] <- rmultinom(1, size = export_draw[i], prob = probs)
        }
      }
    }
    
    # imports = column sums of spill matrix
    import_vec <- colSums(spill_matrix)
    n_spill_in[, v] <- import_vec
    
    # update y: subtract exports, add imports
    y_eff[, v] <- y_mat[, v] - export_draw + import_vec
  }
  
  list(y_eff       = y_eff,
       n_spill_out = n_spill_out,
       n_spill_in  = n_spill_in,
       neighbor_counts = neigh_abund,
       W = W,
       K = K)
}

apply_spillover_y_conserve <- function(y_mat,
                                       coords,
                                       landscape,
                                       spillover_rate = 0.05,
                                       cutoff = 0.17,
                                       kernel = c("binary", "exp"),
                                       decay_scale = 1,
                                       seed = NULL) {
  if(!is.null(seed)) set.seed(seed)
  kernel <- match.arg(kernel)
  
  n_sites  <- nrow(y_mat)
  n_visits <- ncol(y_mat)
  
  coords <- as.matrix(coords)
  storage.mode(coords) <- "numeric"
  
  # distance & mask
  dmat <- as.matrix(dist(coords))
  same_land <- outer(landscape, landscape, "==")
  
  # kernel
  if(kernel == "binary") {
    K <- (dmat <= cutoff) * 1
  } else {
    K <- exp(-dmat / decay_scale)
  }
  diag(K) <- 0L
  K[!same_land] <- 0
  
  # row-normalize to W (rows sum to 1 for sites with neighbors)
  row_sums <- rowSums(K)
  W <- matrix(0, n_sites, n_sites)
  nz <- which(row_sums > 0)
  if(length(nz) > 0) W[nz, ] <- K[nz, , drop = FALSE] / row_sums[nz]
  
  # outputs
  y_eff       <- matrix(0L, n_sites, n_visits)
  n_spill_out <- matrix(0L, n_sites, n_visits)
  n_spill_in  <- matrix(0L, n_sites, n_visits)
  neigh_abund <- matrix(0, n_sites, n_visits)  # weighted neighbor counts (diagnostic)
  
  # loop visits
  for(v in seq_len(n_visits)) {
    yv <- as.integer(y_mat[, v])
    
    # neighbor abundance diagnostic (how many neighbors have)
    neigh_abund[, v] <- as.numeric(W %*% yv)
    
    # export draws (based on local abundance) and cap at available
    lambda_vec <- spillover_rate * yv
    export_draw <- rpois(n_sites, lambda_vec)
    export_draw <- pmin(export_draw, yv)     # cannot export more than present
    
    # if no neighbors, can't export (return to zero)
    no_nb <- which(rowSums(W) == 0)
    if(length(no_nb) > 0) export_draw[no_nb] <- 0L
    
    n_spill_out[, v] <- export_draw
    
    # build spill matrix (rows = source, cols = destination)
    spill_matrix <- matrix(0L, nrow = n_sites, ncol = n_sites)
    for(i in seq_len(n_sites)) {
      szi <- export_draw[i]
      if(szi > 0) {
        probs <- W[i, ]
        if(sum(probs) > 0) {
          # rmultinom returns matrix; take column 1
          spill_matrix[i, ] <- as.integer(rmultinom(1, size = szi, prob = probs)[,1])
        }
      }
    }
    
    imports <- colSums(spill_matrix)
    n_spill_in[, v] <- imports
    
    # update counts: subtract exports, add imports
    y_eff[, v] <- yv - export_draw + imports
    # safety floor (shouldn't be needed): ensure non-negative
    y_eff[, v][y_eff[, v] < 0] <- 0L
  }
  
  # conservation check
  total_before <- sum(y_mat)
  total_after  <- sum(y_eff)
  
  list(y_eff = y_eff,
       n_spill_out = n_spill_out,
       n_spill_in = n_spill_in,
       neighbor_counts = neigh_abund,
       W = W,
       K = K,
       total_before = total_before,
       total_after = total_after)
}

apply_spillover_y <- function(y_mat,
                              coords,
                              landscape,
                              spillover_rate = 0.05,   # mean fraction spilled
                              cutoff = 0.17,
                              kernel = c("binary", "exp"),
                              decay_scale = 1,
                              seed = NULL) {
  if(!is.null(seed)) set.seed(seed)
  kernel <- match.arg(kernel)
  
  n_sites  <- nrow(y_mat)
  n_visits <- ncol(y_mat)
  
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
  W <- matrix(0, n_sites, n_sites)
  nz <- which(row_sums > 0)
  if(length(nz) > 0) W[nz, ] <- K[nz, , drop = FALSE] / row_sums[nz]
  
  # outputs
  y_eff       <- matrix(0L, n_sites, n_visits)
  n_spill_out <- matrix(0L, n_sites, n_visits)
  n_spill_in  <- matrix(0L, n_sites, n_visits)
  
  for(v in seq_len(n_visits)) {
    yv <- as.integer(y_mat[, v])
    
    # --- step 1: draw exports (binomial fraction of current abundance) ---
    
    q_i <- rbeta(n_sites, 0.1*1, (1-0.1)*1)           # fraction exposed (~10% on average)
    export_draw <- rbinom(n_sites, size = rbinom(n_sites, yv, q_i), prob = spillover_rate)  # spill from subset
    
  #  export_draw <- rbinom(n_sites, size = yv, prob = spillover_rate)
    
    n_spill_out[, v] <- export_draw
    
    # --- step 2: allocate to neighbors ---
    spill_matrix <- matrix(0L, n_sites, n_sites)
    for(i in seq_len(n_sites)) {
      szi <- export_draw[i]
      if(szi > 0 && row_sums[i] > 0) {
        probs <- W[i, ]
        spill_matrix[i, ] <- as.integer(rmultinom(1, size = szi, prob = probs)[, 1])
      }
    }
    
    imports <- colSums(spill_matrix)
    n_spill_in[, v] <- imports
    
    # --- step 3: update abundances ---
    y_eff[, v] <- yv - export_draw + imports
    y_eff[, v][y_eff[, v] < 0] <- 0L
  }
  
  list(y_eff = y_eff,
       n_spill_out = n_spill_out,
       n_spill_in = n_spill_in,
       W = W,
       K = K,
       total_before = sum(y_mat),
       total_after = sum(y_eff))
}


ret<-apply_spillover_y(y_mat = N.matrix.temp.dom,
                  coords = coords,
                  landscape = landscape,
                  spillover_rate = .5,
                  cutoff = 0.17, # if we change the number of sites, this will need to change. Can use pythagoras to determine automatically - but keep simple for now
                  kernel = "binary", 
                  decay_scale = 0.1,
                  seed = NULL)

i <- 59
N.matrix.temp.dom[i,]
ret$y_eff[i,]
ret$n_spill_out[i,]
ret$n_spill_in[i,]


# tough recovering w/ spatial autocorrelation, explore why
cols <- ifelse(sim$true$spatial_effect > 0, "red", adjustcolor("grey50", alpha.f = 0.5))
plot(mod$mean$N.sub[100:130]~mod$mean$N.dom[100:130], col=cols, pch=19)

sim$true$N.dom

sim <- simData[[397]]


cols <- ifelse(apply(sim$data$y.dom, 1, max) > 0, adjustcolor("red", alpha.f = 0.5), adjustcolor("grey50", alpha.f = 0.5))

plot(sim$true$N.sub~sim$true$N.dom, col = cols, pch=19,
     ylab = "N - Subordinant",
     xlab = "N - Dominant")

legend("topright",
       legend = c("Predator present", "Predator absent", paste0("SIV = ", sim$true$a5)),
       col = c(adjustcolor("red", alpha.f = 0.5), 
               adjustcolor("grey50", alpha.f = 0.5),
               adjustcolor("grey50", alpha.f = 0)),
       pch = 19, pt.cex = 1.2, bty = "n")

plot(sim$true$N.sub~sim$data$flii, col = cols, pch=19,
     ylab = "N - Subordinant",
     xlab = "Forest Integrity")

legend("topright",
       legend = c("Predator present", "Predator absent", paste0("SIV = ", sim$true$a5)),
       col = c(adjustcolor("red", alpha.f = 0.5), 
               adjustcolor("grey50", alpha.f = 0.5),
               adjustcolor("grey50", alpha.f = 0)),
       pch = 19, pt.cex = 1.2, bty = "n")

domy.mod <- apply(mod$mean$y.dom, 1, max)
domy.true <- apply(sim$data$y.dom, 1, max)

plot(mod$mean$lambda.dom~sim$data$Z.dom, col=cols, pch=19)
plot(mod$mean$lambda.sub~sim$true$lambda.sub, col=cols, pch=19)

# Build dataframe for plotting
df_spf <- data.frame(
  x = sim$true$x,
  y = sim$true$y,
  landscape = factor(sim$data$area),
  spf = sim$true$spatial_effect
)

# Visualise as coloured points
ggplot(df_spf[41:80,], aes(x = x, y = y, fill = spf)) +
  geom_tile(size = 3) +
  facet_wrap(~landscape) +
  scale_fill_viridis_c(option = "plasma") +
  coord_equal() +
  theme_minimal() +
  labs(title = "Spatial autocorrelation field per landscape",
       colour = "Spatial effect")

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





#--------------------------------------------------
# Spillover function
#--------------------------------------------------
apply_spillover <- function(N = sim$true$N.dom, coords = coords, spillover_rate= 3.5){
  # N: matrix of site × species (true latent abundances)
  # coords: site coordinates (n_sites × 2)
  # movement_scale: distance decay parameter for kernel
  # spillover_strength: proportion of neighbor abundance contributing
  
  n_sites   <- nrow(N)
  n_species <- ncol(N)
  
  # distance matrix
  dmat <- as.matrix(dist(coords))
  
  # nearest-neighbor kernel (only sites exactly 1 unit away)
  K <- (dmat <= .2) * 1
  diag(K) <- 0
  # normalize rows (so contributions sum to 1)
  row_sums <- rowSums(K)
  row_sums[row_sums == 0] <- 1   # avoid division by zero
  W <- K / row_sums
  
  # initialize effective abundance
  N_eff <- N
  spill.tracker <- rep(1,length = nrow(N))
  
  for(site in 1:n_sites){
    # neighbors' total abundance
    neighbor_abund <- sum(N[,] * W[site,])
    
    # Poisson draw for spillover individuals
    n_spill <- rpois(1, lambda = spillover_rate * neighbor_abund)
    if(n_spill > 0){spill.tracker[site] <- 2}
    # add to focal site
    N_eff[site,] <- N_eff[site,s] + n_spill
  }
  
  
  return(N_eff)
}

apply_spillover()

t1 <- apply(N_eff, 1, max)
t2 <- apply(N, 1, max)

plot(t1~t2, col = spill.tracker, pch=19)

K <- matrix(0, n_sites, n_sites)

for(i in 1:n_sites){
  for(j in 1:n_sites){
    if(i == j) next
    dx <- abs(coords[i,1] - coords[j,1])
    dy <- abs(coords[i,2] - coords[j,2])
    if(diagonal){
      if((dx <= 1) & (dy <= 1)) K[i,j] <- 1
    } else {
      if((dx == 1 & dy == 0) | (dx == 0 & dy == 1)) K[i,j] <- 1
    }
  }
}
diagonal = F


apply_spillover_scalar <- function(N_vec = sim$true$N.sub, coords = data.frame(x = sim$true$x,
                                                                               y = sim$true$y), 
                                   landscape = sim$data$area,
                                   spillover_rate = 3.05,
                                   cutoff = 0.2,
                                   kernel = c("binary", "exp"),
                                   decay_scale = 0.1,
                                   seed = 2025) {
  # N_vec: numeric vector, length = n_sites
  # coords: n_sites x 2 matrix/data.frame of x,y (as from gen_coords_by_landscape)
  # landscape: integer/factor vector length n_sites
  # cutoff: max distance for binary adjacency
  # kernel: "binary" uses d <= cutoff; "exp" uses exp(-d/decay_scale) (still masked by landscape)
  if(!is.null(seed)) set.seed(seed)
  kernel <- match.arg(kernel)
  n_sites <- length(N_vec)
  
  # distance matrix
  dmat <- as.matrix(dist(coords))
  
  # same-landscape mask
  same_land <- outer(landscape, landscape, "==")
  
  # raw kernel
  if(kernel == "binary") {
    K <- (dmat <= cutoff) * 1
  } else { # exponential distance-decay
    K <- exp(-dmat / decay_scale)
  }
  diag(K) <- 0L
  
  # force zero between landscapes
  K[!same_land] <- 0
  
  # normalize rows -> W
  row_sums <- rowSums(K)
  no_nb <- which(row_sums == 0)
  row_sums[row_sums == 0] <- 1    # avoid div-by-zero
  W <- K / row_sums
  if(length(no_nb) > 0) W[no_nb, ] <- 0   # explicit zero rows for isolated cells
  
  # neighbor abundance (vectorized)
  neighbor_abund <- as.numeric(W %*% N_vec)   # n_sites x 1
  
  # spillover draws (vectorized)
  lambda_vec <- spillover_rate * neighbor_abund
  n_spill_vec <- rpois(n_sites, lambda = lambda_vec)
  
  # effective abundance (non-iterative: based on original N_vec)
  N_eff <- N_vec + n_spill_vec
  
  tracker <- ifelse(n_spill_vec > 0, 2L, 1L)  # 1 = no spill, 2 = got spill
  
  list(N_eff = N_eff,
       n_spill = n_spill_vec,
       neighbor_abund = neighbor_abund,
       W = W,
       K = K,
       tracker = tracker)
}

ret <- apply_spillover_scalar(kernel = "binary")

plot(ret$N_eff ~ sim$true$N.sub, pch=19, col = as.factor(ret$tracker))

# with Y's

apply_spillover_y <- function(y_mat=sim$data$y.dom,
                              coords= data.frame(x = sim$true$x,
                                                 y = sim$true$y),
                              landscape = sim$data$area,
                              spillover_rate = 0.05,
                              cutoff = 0.17,
                              kernel = c("binary", "exp"),
                              decay_scale = 0.1,
                              seed = NULL) {
  # y_mat: n_sites x n_visits integer matrix
  # coords: n_sites x 2 numeric coords
  # landscape: length n_sites, grouping factor
  # cutoff: distance cutoff for binary kernel
  # kernel: "binary" or "exp"
  
  if(!is.null(seed)) set.seed(seed)
  kernel <- match.arg(kernel)
  
  n_sites  <- nrow(y_mat)
  n_visits <- ncol(y_mat)
  
  coords <- as.matrix(coords)
  storage.mode(coords) <- "numeric"
  
  # distance matrix
  dmat <- as.matrix(dist(coords))
  
  # same-landscape mask
  same_land <- outer(landscape, landscape, "==")
  
  # kernel
  if(kernel == "binary") {
    K <- (dmat <= cutoff) * 1
  } else {
    K <- exp(-dmat / decay_scale)
  }
  diag(K) <- 0L
  K[!same_land] <- 0
  
  # normalize rows
  row_sums <- rowSums(K)
  row_sums[row_sums == 0] <- 1
  W <- K / row_sums
  W[row_sums == 0, ] <- 0
  
  # initialize
  y_eff      <- y_mat
  n_spill    <- matrix(0L, n_sites, n_visits)
  neigh_abund <- matrix(0, n_sites, n_visits)
  
  # loop over visits (spill calculated independently per replicate)
  for(v in 1:n_visits) {
    neighbor_counts <- as.numeric(W %*% y_mat[, v])
    neigh_abund[, v] <- neighbor_counts
    
    lambda_vec <- spillover_rate * neighbor_counts
    spill_draw <- rpois(n_sites, lambda_vec)
    
    n_spill[, v] <- spill_draw
    y_eff[, v]   <- y_eff[, v] + spill_draw
  }
  
  list(y_eff = y_eff,
       n_spill = n_spill,
       neighbor_counts = neigh_abund,
       W = W,
       K = K)
}

ret <- apply_spillover_y(kernel = "binary")




plot(ret$y_eff ~ sim$data$y.dom,
     pch = 19, col = cols,type='n',
     xlab = "Original counts (site 1)",
     ylab = "With spillover (site 1)", xlim = c(0,max(sim$data$y.dom)),
     ylim = c(0,max(ret$y_eff)))
for(i in 35){
 cols <- ifelse(ret$n_spill[i, ] > 0, adjustcolor("red", alpha.f = 0.5), adjustcolor("grey50", alpha.f = 0.5))
  
 points(ret$y_eff[i,] ~ sim$data$y.dom[i,],
        pch = 19, col = cols)
}
abline(0,1,col="grey", lty =2)

legend("topleft",
       legend = c("Spillage", "No Spillage"),
       col = c(adjustcolor("red", alpha.f = 0.5), 
               adjustcolor("grey50", alpha.f = 0.5)),
       pch = 19, pt.cex = 1.2, bty = "n")

df <- data.frame(
  site = 1:sim$data$n_sites,  # row index
  x = sim$true$x,
  y = sim$true$y,
  N = sim$true$N.dom  # site-level abundance
)

ggplot(df[1:40,], aes(x = x, y = y, fill = N)) +
  geom_tile(width = 0.2, height = 0.2) +   # small squares for heatmap
  geom_text(aes(label = site), color = "white", size = 3)+
  scale_fill_viridis_c(option = "plasma") +
  coord_equal() +
  theme_minimal() +
  labs(fill = "Abundance", title = "Dominant species abundance heatmap")

# maybe spillover should act on N. that is more ecologically sound. but ok for now
apply_spillover_y <- function(y_mat,
                              coords,
                              landscape,
                              spillover_rate = 0.05,
                              cutoff = 0.17,
                              kernel = c("binary", "exp"),
                              decay_scale = 1,
                              seed = NULL) {
  # y_mat: n_sites x n_visits integer matrix
  # coords: n_sites x 2 numeric coords
  # landscape: length n_sites, grouping factor
  # cutoff: distance cutoff for binary kernel
  # kernel: "binary" or "exp"
  
  if(!is.null(seed)) set.seed(seed)
  kernel <- match.arg(kernel)
  
  n_sites  <- nrow(y_mat)
  n_visits <- ncol(y_mat)
  
  coords <- as.matrix(coords)
  storage.mode(coords) <- "numeric"
  
  # distance matrix
  dmat <- as.matrix(dist(coords))
  
  # same-landscape mask
  same_land <- outer(landscape, landscape, "==")
  
  # kernel
  if(kernel == "binary") {
    K <- (dmat <= cutoff) * 1
  } else {
    K <- exp(-dmat / decay_scale)
  }
  diag(K) <- 0L
  K[!same_land] <- 0
  
  # normalize rows
  row_sums <- rowSums(K)
  row_sums[row_sums == 0] <- 1
  W <- K / row_sums
  W[row_sums == 0, ] <- 0
  
  # initialize
  y_eff      <- y_mat
  n_spill    <- matrix(0L, n_sites, n_visits)
  neigh_abund <- matrix(0, n_sites, n_visits)
  
  # loop over visits (spill calculated independently per replicate)
  for(v in 1:n_visits) {
    neighbor_counts <- as.numeric(W %*% y_mat[, v])
    neigh_abund[, v] <- neighbor_counts
    
    lambda_vec <- spillover_rate * neighbor_counts
    spill_draw <- rpois(n_sites, lambda_vec)
    
    n_spill[, v] <- spill_draw
    y_eff[, v]   <- y_eff[, v] + spill_draw
  }
  
  list(y_eff = y_eff,
       n_spill = n_spill,
       neighbor_counts = neigh_abund,
       W = W,
       K = K)
}


