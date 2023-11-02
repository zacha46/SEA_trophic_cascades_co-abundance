###### Zachary Amir's co-abundance model, Z.Amir@uq.edu.au
##### Using 1) count matricies, 2) site covarites, & 3) observation covariates. 
#### Currently examining 279 pairwise interactions between 44 species based on trophic guilds

### Data has already been formatted to into the proper data bundle to run the co-abundance model on the HPC 
## But because models are too large to save ALL of them, 
# this script will extract only the relevant dataframes. 

### October 2023 Update-
## In an effort to improve convergence and speed, 
# I am attempting to run this model with NIMBLE

# Author:
# Zachary Amir, z.amir@uq.edu.au

library(tidyverse)
library(nimble)