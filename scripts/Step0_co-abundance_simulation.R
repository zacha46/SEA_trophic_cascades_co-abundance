### Exploring simulations and dataset generation to be applied for co-abundance models 
## Mainly just looking at some exploratory data for a simulation from ChatGPT 

## Zachary Amir, Z.Amir@uq.edu.au
## Code created: July 25th, 2025
## Last updated: 




### copy/paste code from ChatGPT, then edit ruthlessly 
simulate_coabundance_data <- function(
    # regular variables specificed in the data 
    n_sites = 100,
    n_reps = 3,
    n_area = 5,
    n_year = 4,
    n_source = 3,
    true_a5 = 0.3,
    # logical values to introduce noise/messy data
    add_unobserved_confounder = FALSE,
    add_covariate_noise = FALSE,
    add_double_counts = FALSE,
    inflate_detection_noise = FALSE,
    a5_vary_by_area = FALSE ## possibly remove, created by chat
) {
  # Observed covariates
  flii <- rnorm(n_sites, 1, 0.3)
  hfp <- rnorm(n_sites, 0.5, 0.2)
  elev <- rnorm(n_sites, 800, 100)
  comm_det <- rnorm(n_sites, 0, 1)
  
  # Site-level grouping
  area <- sample(1:n_area, n_sites, replace = TRUE)
  year <- sample(1:n_year, n_sites, replace = TRUE)
  source <- sample(1:n_source, n_sites, replace = TRUE)
  
  # Random effects
  a6 <- rnorm(n_area, 0, 0.5)
  a7 <- rnorm(n_area, 0, 0.5)
  a8 <- rnorm(n_year, 0, 0.5)
  a9 <- rnorm(n_year, 0, 0.5)
  b3 <- rnorm(n_source, 0, 0.5)
  b4 <- rnorm(n_source, 0, 0.5)
  
  # Optional unobserved confounding
  confounder <- if (add_unobserved_confounder) rnorm(n_sites, 0, 1) else rep(0, n_sites)
  
  # Measurement error in covariates
  flii_obs <- if (add_covariate_noise) flii + rnorm(n_sites, 0, 0.1) else flii
  hfp_obs <- if (add_covariate_noise) hfp + rnorm(n_sites, 0, 0.05) else hfp
  elev_obs <- if (add_covariate_noise) elev + rnorm(n_sites, 0, 10) else elev
  comm_det_obs <- if (add_covariate_noise) comm_det + rnorm(n_sites, 0, 0.2) else comm_det
  
  # True interaction strength: constant or varying by area
  if (a5_vary_by_area) {
    area_a5 <- rnorm(n_area, true_a5, 0.3)
    a5_used <- area_a5[area]
  } else {
    a5_used <- rep(true_a5, n_sites)
  }
  
  # Linear predictors
  lp_dom <- 0 + 0.4*flii + 0.3*hfp + 0.2*elev + 0.5*comm_det + a7[area] + a9[year] + confounder
  lambda_dom <- exp(lp_dom)
  N_dom <- rpois(n_sites, lambda_dom)
  
  lp_sub <- 0 + 0.5*flii + 0.2*hfp + 0.3*elev + 0.6*comm_det + a5_used * N_dom + a6[area] + a8[year] + confounder
  lambda_sub <- exp(lp_sub)
  N_sub <- rpois(n_sites, lambda_sub)
  
  # Simulate observation covariates (cams)
  cams <- matrix(sample(1:3, n_sites * n_reps, replace = TRUE), n_sites, n_reps)
  
  # Observation model with optional inflation
  y_dom <- matrix(NA, n_sites, n_reps)
  y_sub <- matrix(NA, n_sites, n_reps)
  
  for (j in 1:n_sites) {
    for (k in 1:n_reps) {
      sd_eps <- if (inflate_detection_noise) 1 else 0.2
      eps_dom <- rnorm(1, 0, sd_eps)
      eps_sub <- rnorm(1, 0, sd_eps)
      
      lp_p_dom <- -0.5 + 0.3 * cams[j,k] + b4[source[j]] + eps_dom
      p_dom <- plogis(lp_p_dom)
      y_dom[j,k] <- rbinom(1, N_dom[j], p_dom)
      
      lp_p_sub <- -0.5 + 0.3 * cams[j,k] + b3[source[j]] + eps_sub
      p_sub <- plogis(lp_p_sub)
      y_sub[j,k] <- rbinom(1, N_sub[j], p_sub)
    }
  }
  
  # Add double counting noise
  if (add_double_counts) {
    y_dom <- y_dom + matrix(rpois(n_sites * n_reps, 0.3), n_sites, n_reps)
    y_sub <- y_sub + matrix(rpois(n_sites * n_reps, 0.3), n_sites, n_reps)
  }
  
  return(list(
    y_dom = y_dom,
    y_sub = y_sub,
    N_dom = N_dom,
    N_sub = N_sub,
    flii = flii_obs,
    hfp = hfp_obs,
    elev = elev_obs,
    comm_det = comm_det_obs,
    area = area,
    year = year,
    source = source,
    cams = cams,
    true_a5 = a5_used,
    scenario = list(
      a5_constant = !a5_vary_by_area,
      unobserved_confounder = add_unobserved_confounder,
      covariate_noise = add_covariate_noise,
      double_counts = add_double_counts,
      detection_noise = inflate_detection_noise
    )
  ))
}


## make data 
set.seed(42)
data_sim <- simulate_coabundance_data(true_a5 = 0.4)

## inspect
head(data_sim$y_dom)
head(data_sim$comm_det)
head(data_sim$hfp)
head(data_sim$a5)

