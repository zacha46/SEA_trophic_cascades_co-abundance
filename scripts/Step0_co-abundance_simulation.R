### Exploring simulations and dataset generation to be applied for co-abundance models 
## Mainly just looking at some exploratory data for a simulation from ChatGPT 

## Zachary Amir, Z.Amir@uq.edu.au
## Code created: July 25th, 2025
## Last updated: August 5th, 2025 

# load library
library(tidyverse)     ## For lots of functions 

#### Provide example data to ChatGPT ##### 

## import example data
dat = readRDS("data_CoA_bundles/Bundled_data_for_Bayes_co-abundance_mods_community_500GB_33_species_pairs_5km_20250306.RDS")
dat = dat[[1]] # only need one

# quick inspection 
table(dat$Z.dom);table(dat$Z.sub)
summary(dat$y.dom);summary(dat$y.sub)
names(dat)
str(dat)

# Combine key variables into a data frame
meta_df <- data.frame(
  row_id = 1:length(dat$area),
  area = dat$area,
  source = dat$source,
  year = dat$year
)

# Keep rows from 5 different areas, 5 sources, and 5 years
set.seed(123)  # for reproducibility
selected_areas <- sample(unique(meta_df$area), 5)
selected_sources <- sample(unique(meta_df$source), 5)
selected_years <- sample(unique(meta_df$year), 5)

# Filter rows that match any of the selected values in each dimension
filtered_df <- meta_df %>%
  filter(area %in% selected_areas,
         source %in% selected_sources,
         year %in% selected_years)

# Take a maximum of 300 rows
dat_small_rows <- head(filtered_df$row_id, 300)

# Establish the rows and cols to subset
subset_rows <- dat_small_rows
subset_cols <- 1:10 # keep it thin! 

# Subset and save representative data
dat_small <- dat
dat_small$y.dom   <- dat$y.dom[subset_rows, subset_cols]
dat_small$y.sub   <- dat$y.sub[subset_rows, subset_cols]
dat_small$cams    <- dat$cams[subset_rows, subset_cols]
dat_small$Z.dom   <- dat$Z.dom[subset_rows]
dat_small$Z.sub   <- dat$Z.sub[subset_rows]
dat_small$flii    <- dat$flii[subset_rows]
dat_small$hfp     <- dat$hfp[subset_rows]
dat_small$elev    <- dat$elev[subset_rows]
dat_small$comm_det <- dat$comm_det[subset_rows]
dat_small$area    <- dat$area[subset_rows]
dat_small$source  <- dat$source[subset_rows]
dat_small$year    <- dat$year[subset_rows]
dat_small$nsites  <- length(subset_rows)
dat_small$nreps   <- length(subset_cols)

## Make sure we have variation! 
table(dat_small$year)
table(dat_small$area)
table(dat_small$source) # all good! 
summary(dat_small$y.dom);summary(dat_small$y.sub) # looks good. 

## Save the data for an upload
# write_rds(dat_small, "~/Dropbox/Zach PhD/Ch3 Trophic release project/SEA_TC_GitHub_data_storage/explore/example_coabundance_subset_stratified_for_ChatGPT_20250804.rds")

## Chat doesnt like RDS, so flatten to a csv
# Flatten for example, keeping site-level variables
df_export <- data.frame(
  y.dom.1 = dat_small$y.dom[, 1],
  y.sub.1 = dat_small$y.sub[, 1],
  Z.dom = dat_small$Z.dom,
  Z.sub = dat_small$Z.sub,
  flii = dat_small$flii,
  hfp = dat_small$hfp,
  elev = dat_small$elev,
  comm_det = dat_small$comm_det,
  cams.1 = dat_small$cams[, 1],
  area = dat_small$area,
  year = dat_small$year,
  source = dat_small$source
)

# write.csv(df_export, "~/Dropbox/Zach PhD/Ch3 Trophic release project/SEA_TC_GitHub_data_storage/explore/example_coabundance_subset_flat_20250804.csv", row.names = FALSE)
rm(dat_small, df_export, filtered_df, dat, meta_df, selected_areas, selected_sources, 
   selected_years, subset_cols, subset_rows, dat_small_rows)

##### Load ChatGPT data simulator funciton #####

simulate_coabundance_full <- function(
    nsites = 300,
    nreps = 20,
    narea = 10,
    nyear = 5,
    nsource = 5,
    a5 = 0,
    bias = c("none"),  # Accepts vector: e.g., c("spatial", "unmeasured_state")
    seed = 42
) {
  set.seed(seed)
  bias <- match.arg(bias, choices = c("none", "double_count", "spatial", "unmeasured_state", "unmeasured_detection"), several.ok = TRUE)
  
  # === 1. FIXED EFFECTS ===
  flii <- scale(rnorm(nsites))[,1]
  hfp <- scale(rnorm(nsites))[,1]
  elev <- scale(rnorm(nsites))[,1]
  comm_det <- scale(log(runif(nsites, 1, 100) + 1))[,1]
  
  # === 2. RANDOM EFFECT GROUPINGS ===
  area <- sample(1:narea, nsites, replace = TRUE)
  year <- sample(1:nyear, nsites, replace = TRUE)
  source <- sample(1:nsource, nsites, replace = TRUE)
  
  # === 3. OCCUPANCY FLAGS ===
  Z.dom <- rep(1, nsites)
  Z.sub <- rbinom(nsites, size = 1, prob = 0.75)
  
  # === 4. RANDOM EFFECTS ===
  a6 <- rnorm(narea, 0, 0.5)
  a7 <- rnorm(narea, 0, 0.5)
  a8 <- rnorm(nyear, 0, 0.5)
  a9 <- rnorm(nyear, 0, 0.5)
  b3 <- rnorm(nsource, 0, 0.4)
  b4 <- rnorm(nsource, 0, 0.4)
  
  # === 5. OPTIONAL UNMEASURED COVARIATES ===
  if ("unmeasured_state" %in% bias) {
    u_state <- scale(rnorm(nsites))[,1]
  } else {
    u_state <- rep(0, nsites)
  }
  
  if ("unmeasured_detection" %in% bias) {
    u_det <- scale(rnorm(nsites))[,1]
  } else {
    u_det <- rep(0, nsites)
  }
  
  # === 6. SPATIAL AUTOCORRELATION ===
  if ("spatial" %in% bias) {
    spatial_effect <- scale(stats::filter(rnorm(nsites + 10), rep(1/5, 5), sides = 2))[6:(nsites+5)]
  } else {
    spatial_effect <- rep(0, nsites)
  }
  
  # === 7. STATE MODEL: LAMBDA ===
  lambda_dom <- exp(
    1 + 0.6*flii + 0.4*hfp + 0.3*elev + 0.2*comm_det +
      a7[area] + a9[year] + u_state + spatial_effect
  )
  
  N.dom <- rpois(nsites, lambda_dom * Z.dom)
  
  lambda_sub <- exp(
    1 + 0.3*flii + 0.2*hfp + 0.1*elev + 0.1*comm_det +
      a5 * lambda_dom + a6[area] + a8[year] + u_state + spatial_effect
  )
  
  N.sub <- rpois(nsites, lambda_sub * Z.sub)
  
  # === 8. EFFORT MATRIX (cams) ===
  raw_effort <- matrix(runif(nsites * nreps, min = 0.5, max = 1.5), nsites, nreps)
  cams <- scale(raw_effort)
  
  # === 9. DOUBLE COUNTING (shared latent detection burst) ===
  if ("double_count" %in% bias) {
    burst_effect <- matrix(rnorm(nsites, 0, 0.8), nsites, nreps)
  } else {
    burst_effect <- matrix(0, nsites, nreps)
  }
  
  # === 10. DETECTION MODEL ===
  y.dom <- matrix(0, nsites, nreps)
  y.sub <- matrix(0, nsites, nreps)
  
  for (j in 1:nsites) {
    for (k in 1:nreps) {
      # Detection probabilities
      lp_dom <- -1 + 0.5 * cams[j,k] + b4[source[j]] + u_det[j] + spatial_effect[j] + burst_effect[j,k]
      lp_sub <- -1.5 + 0.3 * cams[j,k] + b3[source[j]] + u_det[j] + spatial_effect[j] + burst_effect[j,k]
      
      p_dom <- plogis(lp_dom)
      p_sub <- plogis(lp_sub)
      
      ## make sure matricies are skewed by very high numbers and cap values @ 10 
      y.dom[j,k] <- min(rbinom(1, N.dom[j], p_dom), 10)
      y.sub[j,k] <- min(rbinom(1, N.sub[j], p_sub), 10)
    }
  }
  
  # === 11. RETURN DATA STRUCTURE ===
  return(list(
    y.dom = y.dom,
    y.sub = y.sub,
    Z.dom = Z.dom,
    Z.sub = Z.sub,
    N.dom = N.dom,
    N.sub = N.sub,
    flii = flii,
    hfp = hfp,
    elev = elev,
    comm_det = comm_det,
    area = area,
    year = year,
    source = source,
    cams = cams,
    nsites = nsites,
    nreps = nreps,
    narea = narea,
    nyear = nyear,
    nsource = nsource,
    a5_true = a5,
    bias_applied = bias
  ))
}

# Clean simulation with a5 = 0
sim_clean <- simulate_coabundance_full(a5 = 0)

# Add unmeasured confounding in state only
sim_confounded <- simulate_coabundance_full(a5 = 1, bias = "unmeasured_state")

# Combine multiple biases
sim_messy <- simulate_coabundance_full(a5 = 1, bias = c("double_count", "unmeasured_detection", "spatial"))

# Inspect outputs
str(sim_messy)
image(sim_messy$y.sub)
rm(sim_clean, sim_confounded, sim_messy)

### Now we are going to create 18 datasets based on multiple conditions
# bias = none, double counts of individuals, spatial autocorrelation, unmeasured variables in state & det, and a combo of all
# SIV = -.5, 0, .5

# save biases as a list
bias_types <- list(
  "none",
  "double_count",
  "spatial",
  "unmeasured_state",
  "unmeasured_detection",
  c("double_count", "spatial", "unmeasured_state", "unmeasured_detection")  # all biases
)

# store true a5 values we want to test
a5_values <- c(-1, -.5, 0, .5, 1)

# Storage
simulated_datasets <- list()

# Generate datasets
for (a5 in a5_values) {
  for (bias in bias_types) {
    key <- paste0("a5_", a5, "_bias_", paste(bias, collapse = "+"))
    simulated_datasets[[key]] <- simulate_coabundance_full(a5 = a5, bias = bias)
  }
}
rm(a5, bias, key)

## check 
length(simulated_datasets) # 30 w/ extra a5 values  
names(simulated_datasets) # very good. 

## grab the date
day<-substr(Sys.Date(),9, 10)
month<-substr(Sys.Date(),6,7)
year<-substr(Sys.Date(),1,4)
date = paste(year, month, day, sep = "")
rm(day, month, year)

# make a path 
path = paste("~/Dropbox/Zach PhD/Ch3 Trophic release project/SEA_TC_GitHub_data_storage/data/step2_output_CoA_bundles/Bundled_data_for_bayes_co-abudance_mods_SIMULATION_", length(simulated_datasets), "_species_pairs_", 
             date, ".RDS", sep = "")

## save this as a RDS object
write_rds(simulated_datasets, path)
