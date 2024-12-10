#### Results file name reformatting and re-naming
### ZAchary Amir, Z.Amir@uq.edu.au
## Code initalized December 10th, 2024

## Purpose of this code is to re-organize files to minimize how many models need to be re-run. 
## ultimetly, want to copy/paste community_deteciton models as the new normal model 
## and save the old 'normal' model as a counter-factual test. 


## start fresh
rm(list = ls())

## load libraries
# library(tidyverse)  # For basic data wrangling
# library(fs)         # For file manipulation

# set working d
setwd("~/Dropbox/Zach PhD/Ch3 Trophic release project/SEA_TC_GitHub_data_storage/")

#### Locate relevant files #####

## start w/ coefficents, most important. 
coeff_files = list.files("results/LONG_5km_final_and_counterfactual_Oct2024/coefficent_dataframes/")

## counterfactual7 is the good stuff, grab that 
coeff_files = coeff_files[grepl("counterfactual7", coeff_files)]

## Make a loop to import each file, rename, and then save in a new location
for(i in 1:length(coeff_files)){
  # import one file 
  c = read.csv(paste("results/LONG_5km_final_and_counterfactual_Oct2024/coefficent_dataframes/", 
                     coeff_files[i], sep = ""))
  # change the file name 
  file = gsub("counterfactual7_new_variable_", "", coeff_files[i])
  # make a new saving path
  c_path = paste("results/LONG_5km_final_and_counterfactual_Dec2024_community_detections/coefficent_dataframes/",
                 file, sep = "")
  # save it! 
  # write.csv(c, c_path, row.names = F)
}
rm(i, c, file, c_path, coeff_files)

## now repeat for PPC files 
files = list.files("results/LONG_5km_final_and_counterfactual_Oct2024/PPC_dataframes/")
## counterfactual7 is the good stuff, grab that 
files = files[grepl("counterfactual7", files)]

## Make a loop to import each file, rename, and then save in a new location
for(i in 1:length(files)){
  # import one file 
  c = read.csv(paste("results/LONG_5km_final_and_counterfactual_Oct2024/PPC_dataframes/", 
                     files[i], sep = ""))
  # change the file name 
  file = gsub("counterfactual7_new_variable_", "", files[i])
  # make a new saving path
  c_path = paste("results/LONG_5km_final_and_counterfactual_Dec2024_community_detections/PPC_dataframes/",
                 file, sep = "")
  # save it! 
  write.csv(c, c_path, row.names = F)
}
rm(i, c, file, c_path, files)


#### Now grab the original results w/out counter-factual, and save them as a new coutnerfactual7

# coefficents first 
files = list.files("results/LONG_5km_final_and_counterfactual_Oct2024/coefficent_dataframes/")
## omit all counterfactuals 
files = files[! grepl("counterfactual", files)]

## Make a loop to import each file, rename, and then save in a new location
for(i in 1:length(files)){
  # import one file 
  c = read.csv(paste("results/LONG_5km_final_and_counterfactual_Oct2024/coefficent_dataframes/", 
                     files[i], sep = ""))
  # change the file name 
  file = gsub("G_c", "G_counterfactual7_no_community_detections_c", files[i])
  # make a new saving path
  c_path = paste("results/LONG_5km_final_and_counterfactual_Dec2024_community_detections/counterfactual_coefficent_dataframes/",
                 file, sep = "")
  # save it! 
  write.csv(c, c_path, row.names = F)
}
rm(i, c, file, c_path, files)

## now repeat for the PPC
files = list.files("results/LONG_5km_final_and_counterfactual_Oct2024/PPC_dataframes/")
## omit all counterfactuals 
files = files[! grepl("counterfactual", files)]
## and omit all plot data
files = files[!grepl("plotdata", files)]

## Make a loop to import each file, rename, and then save in a new location
for(i in 1:length(files)){
  # import one file 
  c = read.csv(paste("results/LONG_5km_final_and_counterfactual_Oct2024/PPC_dataframes/", 
                     files[i], sep = ""))
  # change the file name 
  file = gsub("G_c", "G_counterfactual7_no_community_detections_c", files[i])
  # make a new saving path
  c_path = paste("results/LONG_5km_final_and_counterfactual_Dec2024_community_detections/counterfactual_PPC_dataframes/",
                 file, sep = "")
  # save it! 
  write.csv(c, c_path, row.names = F)
}
rm(i, c, file, c_path, files)




#
#