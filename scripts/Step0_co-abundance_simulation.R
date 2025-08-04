### Exploring simulations and dataset generation to be applied for co-abundance models 
## Mainly just looking at some exploratory data for a simulation from ChatGPT 

## Zachary Amir, Z.Amir@uq.edu.au
## Code created: July 25th, 2025
## Last updated: August 4th, 2025 

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

##### Load ChatGPT data simulator funciton #####



