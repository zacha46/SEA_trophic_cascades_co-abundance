### Exploring simulations and dataset generation to be applied for co-abundance models 
## Mainly just looking at some exploratory data for a simulation from ChatGPT 

## Zachary Amir, Z.Amir@uq.edu.au
## Code created: July 25th, 2025
## Last updated: August 5th, 2025 

# load library
library(tidyverse)     ## For lots of functions 
library(jagsUI)        ## for bayes model manipulations

##### Provide example data to ChatGPT ##### 

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
    narea = 12,
    nyear = 4,
    nsource = 6,
    a5 = 0,
    bias = c("none"),  # Accepts vector: e.g., c("spatial", "unmeasured_state")
    seed = 42
) {
  set.seed(seed)
  bias <- match.arg(bias, choices = c("none", "double_count", "spatial", "unmeasured_state", "unmeasured_detection"), several.ok = TRUE)
  
  # === 1. Assign groupings in structured blocks ===
  sites_per_area <- floor(nsites / narea)
  area <- rep(1:narea, each = sites_per_area)
  area <- c(area, sample(1:narea, nsites - length(area), replace = TRUE))  # pad remainder
  
  sites_per_year <- floor(nsites / nyear)
  year <- rep(1:nyear, each = sites_per_year)
  year <- c(year, sample(1:nyear, nsites - length(year), replace = TRUE))
  
  sites_per_source <- floor(nsites / nsource)
  source <- rep(1:nsource, each = sites_per_source)
  source <- c(source, sample(1:nsource, nsites - length(source), replace = TRUE))
  
  # === 2. Group-level covariate structure ===
  area_effect_flii <- rnorm(narea, 0, 0.7)
  area_effect_hfp  <- rnorm(narea, 0, 0.5)
  area_effect_elev <- rnorm(narea, 0, 0.3)
  
  flii_raw <- rnorm(nsites) + area_effect_flii[area]
  hfp_raw  <- rnorm(nsites) + area_effect_hfp[area]
  elev_raw <- rnorm(nsites) + area_effect_elev[area]
  
  flii <- scale(flii_raw)[,1]
  hfp  <- scale(hfp_raw)[,1]
  elev <- scale(elev_raw)[,1]
  
  comm_det <- scale(log(runif(nsites, 1, 100) + 1))[,1]
  
  
  # === 3. iZIP parameter ===
  # Ensure at least 2 areas for each combination (adjust as needed)
  n_each <- floor(narea / 4)
  remaining <- narea - 4 * n_each
  
  # Create vector of all combinations
  combos <- rep(list(c(1,1), c(1,0), c(0,1), c(0,0)), each = n_each)
  if (remaining > 0) {
    # Add a few more random combos to fill out remaining areas
    combos <- c(combos, replicate(remaining, sample(list(c(1,1), c(1,0), c(0,1), c(0,0)), 1), simplify = FALSE))
  }
  
  # Shuffle for randomness
  set.seed(seed)
  combos <- sample(combos)
  
  # Split into Z values
  landscape_Z_dom <- sapply(combos, function(x) x[1])
  landscape_Z_sub <- sapply(combos, function(x) x[2])
  
  # Assign to sites
  Z.dom <- landscape_Z_dom[area]
  Z.sub <- landscape_Z_sub[area]
  
  # === 4. Group-level random effects ===
  a_state_area <- rnorm(narea, 0, 0.5)
  a_state_yr <- rnorm(nyear, 0, 0.5)
  b_det_source <- rnorm(nsource, 0, 0.4)
  
  # === 5. OPTIONAL UNMEASURED COVARIATES ===
  u_state <- if ("unmeasured_state" %in% bias) scale(rnorm(nsites))[,1] else rep(0, nsites)
  u_det   <- if ("unmeasured_detection" %in% bias) scale(rnorm(nsites))[,1] else rep(0, nsites)
  
  # === 6. SPATIAL AUTOCORRELATION ===
  if ("spatial" %in% bias) {
    spatial_effect <- scale(stats::filter(rnorm(nsites + 10), rep(1/5, 5), sides = 2))[6:(nsites+5)]
  } else {
    spatial_effect <- rep(0, nsites)
  }
  
  # === 7. STATE MODEL: LAMBDA ===
  # use a common state model for both species, so create the base linear predictor 
  linpred <- 1 + 0.5*flii + 0.3*hfp + 0.2*elev + 0.2*comm_det +
    u_state + spatial_effect + a_state_area[area] + a_state_yr[year]
  # apply the predictor 
  lambda_dom <- exp(linpred)
  # extract abundance
  N.dom <- rpois(nsites, lambda_dom * Z.dom)
  
  # intalize lambda sub, 
  lambda_sub <- rep(NA, nsites)
  # but apply N_dom as a centered value per RE group 
  for (i in 1:nsites) {
    group_mean <- mean(N.dom[area == area[i]])
    # this is expanding upon the linear predictor to include centered N.dom 
    lambda_sub[i] <- exp(linpred[i] + a5 * (N.dom[i] - group_mean))
  }
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
  
  ## Shared detection parameters, assuming both species are detected equally 
  b_det_intercept <- -1
  b_det_effort <- 0.5
  b_source <- rnorm(nsource, 0, 0.4)  # One shared source effect vector
  
  for (j in 1:nsites) {
    for (k in 1:nreps) {
      # Common detection model for both species Detection probabilities
      lp_common <- b_det_intercept + b_det_effort * cams[j,k] +
        b_det_source[source[j]] + u_det[j] + spatial_effect[j] + burst_effect[j,k]
      
      p_dom <- plogis(lp_common)
      p_sub <- plogis(lp_common)
      
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
    bias_applied = bias,
    a_state_area = a_state_area
  ))
} # end function

### simulate positive and negative a5 w/ no bias
sim_pos = simulate_coabundance_full(a5 = 1, bias = "none")
sim_neg = simulate_coabundance_full(a5 = -1, bias = "none")
# assess correlations between N.dom and N.sub
cor(sim_pos$N.dom, sim_pos$N.sub) # positive 
cor(sim_neg$N.dom, sim_neg$N.sub) # barely negative... suspicous. 

# # Assess if N.dom varies by random effects
# boxplot(N.dom ~ area, data = sim_pos)
# boxplot(N.sub ~ area, data = sim_pos)
# 
# ## check for correlation between N.dom and area random effect 
# area_effects <- tapply(sim_pos$N.dom, sim_pos$area, mean)
# cor(area_effects, unique(sim_pos$a_state_area))  
# 
# plot(sim_neg$N.dom - tapply(sim_neg$N.dom, sim_neg$area, mean)[sim_neg$area],
#      sim_neg$N.sub, col = sim_neg$area, pch = 19)


### Now we are going to create 18 datasets based on multiple conditions
# bias = none, double counts of individuals, spatial autocorrelation, unmeasured variables in state & det, and a combo of all
# SIV = -1,-.5, 0, .5, 1

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
a5_values <- c(-2, -1, -.5, 0, .5, 1, 2)

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
length(simulated_datasets) # 30 w/ extra a5 values, but only 5 if looking at no biases
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
saveRDS(simulated_datasets, path)

##### Run models on the HPC ######

## The code that runs the simulations on the HPC is called scripts/HPC_code/HPC_co-abundance_model_simulations.R
## and the code that communicates with the HPC is called scripts/SLURM_code/MIDDLE/SLURM_co-abundance_array_MIDDLE_simulation.sh

##### Import results from HPC ######

## start fresh 
rm(list = ls())

# save working directory to where the results live in dropbox 
wd = "~/Dropbox/Zach PhD/Ch3 Trophic release project/SEA_TC_GitHub_data_storage/results/"

# and list all relevant files 
files = list.files(paste(wd, "MIDDLE_simulations_August_2025_varying_REs", sep = ""), recursive = T)

# #
# ##
# ### Coefficient data frame 
# 
# ## First, subset for coefficent results
# files_coeff = files[grepl("coefficent_dataframes/", files)]
# # import each one
# res = list()
# for(i in 1:length(files_coeff)){
#   # import the file 
#   d = read.csv(paste(wd, "MIDDLE_simulations_August_2025_v8/", files_coeff[i], sep = ""))
#   # save in the list 
#   res[[i]] = d
#   # save with the test name 
#   names(res)[i] = str_extract(files_coeff[i], "(?<=coefficents_).*(?=_\\d{8}\\.csv)")
# }
# rm(d, i, files_coeff)
# 
# ## Combine in to a DF 
# coeff = do.call(rbind, res)
# rownames(coeff) = NULL
# 
# # extract true a5 and bias as new columns
# coeff <- coeff %>%
#   mutate(
#     true_a5 = as.numeric(str_match(sim_test, "a5_([^_]+)_bias_")[,2]),
#     bias = str_match(sim_test, "bias_(.*)$")[,2]
#   )
# unique(coeff$true_a5)
# unique(coeff$bias) # both are good! 

#
##
### PPC data frame 

## First, subset for coefficent results
files_ppc = files[grepl("PPC_dataframes/", files)]
# import each one
res = list()
for(i in 1:length(files_ppc)){
  # import the file 
  d = read.csv(paste(wd, "MIDDLE_simulations_August_2025_varying_REs/", files_ppc[i], sep = ""))
  # grab the random effect test 
  d$RE_test = str_extract(files_ppc[i], "(?<=values_)[^_]+_REs?|[^_]+_RE(?=_only)")
  # save in the list 
  res[[i]] = d
  # save with the test name 
  names(res)[i] = str_extract(files_ppc[i], "a5.*(?=_\\d{8}\\.csv)")
}
rm(d, i, files_ppc)

## Combine in to a DF 
ppc = do.call(rbind, res)
rownames(ppc) = NULL

# extract true a5 and bias as new columns
ppc <- ppc %>%
  mutate(
    true_a5 = as.numeric(str_match(sim_test, "a5_([^_]+)_bias_")[,2]),
    bias = str_match(sim_test, "bias_(.*)$")[,2]
  )
unique(ppc$true_a5)
unique(ppc$bias) # both are good! 

#### Apply support levels here
## Bayes p-value
ppc$BPV_valid = "No"
ppc$BPV_valid[ppc$BPV.dom >= 0.15 & ppc$BPV.dom <= 0.85 &
                    ppc$BPV.sub >= 0.15 & ppc$BPV.sub <= 0.85] = "Yes"
table(ppc$BPV_valid[!is.na(ppc$Interaction_Estimate)]) 

## over dispersion, C-hat 
ppc$OD_valid = "No"
ppc$OD_valid[ppc$Chat.dom >= 0.95 & ppc$Chat.dom <= 1.3 &
                   ppc$Chat.sub >= 0.95 & ppc$Chat.sub <= 1.3] = "Yes"
table(ppc$OD_valid[!is.na(ppc$Interaction_Estimate)]) 

## Rhat for SIV
ppc$parameter_valid = "No"
ppc$parameter_valid[ppc$Rhat >= 0.99 & ppc$Rhat <= 1.1] = "Yes"
table(ppc$parameter_valid[!is.na(ppc$Interaction_Estimate)]) 

## Create support levels when combining w/ direction
ppc$support[ppc$BPV_valid == "No" |
              ppc$OD_valid == "No" |
              ppc$parameter_valid == "No"] = "unsupported_poor_fit" # lowest level of support --> bad mod.
ppc$support[ppc$BPV_valid == "Yes" & 
              ppc$OD_valid == "Yes" &
              ppc$parameter_valid == "Yes" &
              ppc$Significance == "Non-Significant"] = "unsupported_unclear_SIV" # mid-low support --> good mod, but not important
ppc$support[ppc$BPV_valid == "Yes" & 
              ppc$OD_valid == "Yes" &
              ppc$parameter_valid == "Yes" &
              ppc$Significance == "Significant" &
              ppc$true_a5 < 0  &
              ppc$Interaction_Estimate > 0] = "unsupported_wrong_direction" # almost supportive --> good model, significant result, but not in correct direction for hypothesis.  
ppc$support[ppc$BPV_valid == "Yes" & 
              ppc$OD_valid == "Yes" &
              ppc$parameter_valid == "Yes" &
              ppc$Significance == "Significant" &
              ppc$true_a5 > 0 &
              ppc$Interaction_Estimate < 0] = "unsupported_wrong_direction" # Same as above, but applied for bottom-up
## assign which models support our hypothesis
ppc$support[ppc$BPV_valid == "Yes" & 
              ppc$OD_valid == "Yes" &
              ppc$parameter_valid == "Yes" &
              ppc$Significance == "Significant" &
              ppc$true_a5 < 0  &
              ppc$Interaction_Estimate <= 0] = "Supported" # good model, significant result, in correct direction for hypothesis.  
ppc$support[ppc$BPV_valid == "Yes" & 
              ppc$OD_valid == "Yes" &
              ppc$parameter_valid == "Yes" &
              ppc$Significance == "Significant" &
              ppc$true_a5 > 0 &
              ppc$Interaction_Estimate >= 0] = "Supported" # Same as above, but applied for bottom-up direction. 
# Check results 
table(ppc$support) # most are a poor fit 

#
##
### abundance data frame 

# ## First, subset for coefficent results
# files_abund = files[grepl("prediction_dataframes/", files)]
# # import each one
# res = list()
# for(i in 1:length(files_abund)){
#   # import the file 
#   d = read.csv(paste(wd, "MIDDLE_simulations_August_2025_v8/", files_abund[i], sep = ""))
#   # save in the list 
#   res[[i]] = d
#   # save with the test name 
#   names(res)[i] = str_extract(files_abund[i], "(?<=comparison_).*(?=_\\d{8}\\.csv)")
# }
# rm(d, i, files_abund)
# 
# ## Combine in to a DF 
# abund = do.call(rbind, res)
# rownames(abund) = NULL
# 
# # extract true a5 and bias as new columns
# abund <- abund %>%
#   mutate(
#     true_a5 = as.numeric(str_match(sim_test, "a5_([^_]+)_bias_")[,2]),
#     bias = str_match(sim_test, "bias_(.*)$")[,2]
#   )
# unique(abund$true_a5)
# unique(abund$bias) # both are good! 

#
##
###
#### Inspect results 

## first inspect no bias -> can model recover true value?
dat = ppc[ppc$bias == "none", ]
dat[order(dat$true_a5), c("true_a5", "Interaction_Estimate","lower","upper", "Rhat", "Significance", "support", "RE_test")]
## ok, values are somewhat recovered, but still facing poor fits. 
dat = dat[order(dat$true_a5),] # some very large Rhats, insepct the true mod

# inspect one RE test at a time 
dat[dat$RE_test == "no_REs", c("true_a5", "Interaction_Estimate","lower","upper", "Rhat", "Significance", "support", "RE_test")]


## Inspect each bias level 
dat = ppc[ppc$bias == "unmeasured_state", ]
dat[order(dat$true_a5), c("true_a5", "Interaction_Estimate","lower","upper", "Significance", "support")]
# weird results 

dat = ppc[ppc$bias == "unmeasured_detection", ]
dat[order(dat$true_a5), c("true_a5", "Interaction_Estimate","lower","upper", "Significance", "support")]
# almost all poor fits 

dat = ppc[ppc$bias == "double_count", ]
dat[order(dat$true_a5), c("true_a5", "Interaction_Estimate","lower","upper", "Significance", "support")]
# all poor fits 

dat = ppc[ppc$bias == "spatial", ]
dat[order(dat$true_a5), c("true_a5", "Interaction_Estimate","lower","upper", "Significance", "support")]
# all poor fits 

dat = ppc[ppc$bias == "double_count+spatial+unmeasured_state+unmeasured_detection", ]
dat[order(dat$true_a5), c("true_a5", "Interaction_Estimate","lower","upper", "Significance", "support")]
# all poor fits 

#
#
#
#