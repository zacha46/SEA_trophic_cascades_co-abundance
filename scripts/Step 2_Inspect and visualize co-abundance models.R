##### Inspect and visualize completed co-abundance models
#### Data was prepared on R Script: Step 1_Generate co-abundance data bundles.Rmd
#### and ran on HPC on R Script: HPC_co-abundance_model.R
### Will convert this to Rmd when ready, but R script for exploring now
## Import dataframes, not complete models b/c theyre too heavy

# Zachary Amir, Z.Amir@uq.edu.au
# Last updated on Nov 27th, 2023

## start fresh
rm(list = ls())

## load libs 
library(tidyverse)
library(plyr)
library(cowplot) # for arranging plots 

# library(ggpubr) # for fancy histo + density plots 
# library(fs) # for moving files around 
# library(ggridges) # for ridgeline plots 
# library(reshape2) # for making checkerboard matrix 
# library(leaflet) # for creating concepts of study map
# library(colortools) # for nice colors in the map
# library(jagsUI) # for running co-abundance mods and inspecting jags output 
# library(traitdata) # for species trait data

## set WD (should be automatic b/c RProj, but being safe anyway)
setwd("/Users/zachary_amir/Dropbox/Zach PhD/Trophic release project/SEA_trophic_cascades_co-abundance")

## Import clean cam trap data here for referencing
og_resamp_captures = read.csv("data/send_to_HPC/clean_captures_to_make_UMFs_20230821.csv")
og_resamp_meta = read.csv("data/send_to_HPC/clean_metadata_to_make_UMFs_20230821.csv")

## should also import NON-resampled data for referencing too! 
og_captures = read.csv("/Users/zachary_amir/Dropbox/CT capture histories database/Asian ECL raw CT data/Step4_output_pre-resampling/Clean_independent_captures_20230610.csv")
og_meta = read.csv("/Users/zachary_amir/Dropbox/CT capture histories database/Asian ECL raw CT data/Step4_output_pre-resampling/Clean_independent_metadata_20230610.csv")


######## Import co-abundance coefficent dataframes ######

## first, I will import the bundled data to generate a vector of relevant species pairs 
bdata = readRDS("data/send_to_HPC/Bundled_data_for_Bayes_co-abundance_mods_260_species_pairs_20231106.RDS")
## NOTE, this was not the data sent to the HPC, but several unlikley mods were removed here.
## use this one for referencing completness. 


## save vector of relevant species pairs
all_combos = names(bdata)


# Initate data management-- 
##
###
#### If you made a mistake and forgot to move old data to the OLD folder when importing from HPC, 
### use this code below to automatically look for duplicates and move the older files to the OLD folder. 
## This will do nothing if there are no duplicates, so its ok to run thru. 

# List all files in the coefficents folder
file_list <- dir("results/coefficent_dataframes/", full.names = TRUE)

# Extract the filenames without the date
file_names <- basename(file_list) %>% str_remove("\\d{8}\\.csv$")

# Find the duplicated file names
duplicated_names <- file_names[duplicated(file_names)]

# Iterate over the duplicated names
for (name in duplicated_names) {
  # Get the files with the current name
  current_files <- file_list[file_names == name]
  
  # Find the earliest date among the duplicated files
  earliest_date <- as.character(min(as.Date(str_extract(current_files, "\\d{8}"), format = "%Y%m%d")))
  
  # replace - with nothing to match existing format
  earliest_date = gsub("-", "", earliest_date)
  
  # Get the files with the earliest date
  files_to_move <- current_files[str_extract(current_files, "\\d{8}") == earliest_date]
  
  # Move the files to the destination folder
  file_move(files_to_move, "results/coefficent_dataframes/OLD")
}
rm(current_files, earliest_date, files_to_move, name, duplicated_names, file_names, file_list)



### Import coefficent dataframes

# list all files
files = list.files("results/HALF_LONG_testing_Oct_2023/coefficent_dataframes/")
files = files[!grepl("OLD", files)] # remove any old data 

## Have a few long and all middle files in this directory, subset for one
files = files[grepl("HALF_LONG", files)]
# files = files[grepl("LONG", files)]
# files = files[grepl("MIDDLE", files)]
# files = files[grepl("SHORT", files)]


# store results here
coeff.res = list()

# loop thru each file to import into list
for(i in 1:length(files)){
  
  d = read.csv(paste("results/HALF_LONG_testing_Oct_2023/coefficent_dataframes/", files[i], sep = ""))

  ## if there are pesky row.names, remove em!
  d$X = NULL
  
  coeff.res[[i]] = d
  
  if(ncol(d) != 10){
    print(paste("The file", files[i], "with the index value of", i, 
                "has the wrong number of columns and wont combine well. INVESTIGTE!"))
  }
}
rm(d,i, files)

## bind together and inspect
coeff = do.call(rbind, coeff.res)

## quick inspection
head(coeff)
str(coeff)
rm(coeff.res)

## thin to include only updated mods from bdata 
coeff = coeff[coeff$Species_Pair %in% names(bdata), ]
length(unique(coeff$Species_Pair)) # Nov 27th, all 260 mods are completed w/ HALF_LONG settings! 


## inspect effective sample size
summary(coeff$n.eff[coeff$var == "Species_Interaction"]) # large range, but median and mean are good 

## how is convergence around the species interaction parameter?
summary(coeff$Rhat[coeff$var == "Species_Interaction"]) 
# HALF_LONG - median & 3rd quart are good! 

## Which models failed to converge?
sort(unique(coeff$Species_Pair[coeff$Rhat > 1.2 & coeff$var == "Species_Interaction"])) 
#65 models @ middle, 57 @ half long! 

## vis Rhat
hist(coeff$Rhat[coeff$var == "Species_Interaction"]) # thankfully most are small, but some biggies exist! 


## split apart species pair into sub vs dom species
coeff <- within(coeff, {
  sub_sp <- sapply(strsplit(Species_Pair, "~"), function(x) x[1])
  dom_sp <- sapply(strsplit(Species_Pair, "~"), function(x) x[2])
})
coeff$dom_sp = gsub("DOM-", "", coeff$dom_sp)
coeff$sub_sp = gsub("SUB-", "", coeff$sub_sp)


### Inspect which species pairs are present/missing, 20231127, HALF_LONG settings

setdiff(all_combos, coeff$Species_Pair) # All models present and accounted for! 

## create a dataframe to track our model performance
preform = data.frame("Species_Pair" = all_combos,
                     "sub_species" = as.character(NA),
                     "dom_species" = as.character(NA))

## add individual species to individual cols via loop 
for(i in 1:length(preform$Species_Pair)){
  
  ## gather the combo and split it 
  c = preform$Species_Pair[i]
  c = str_split(c, "~")
  
  ## add it to the DF
  preform$sub_species[i] = c[[1]][1]
  preform$dom_species[i] = c[[1]][2]
  
  ## clean it up
  preform$sub_species[i] = gsub("SUB-", "", preform$sub_species[i])
  preform$dom_species[i] = gsub("DOM-", "", preform$dom_species[i])
  
}
rm(c,i)

## inspect
head(preform)
anyNA(preform) # Must be F

## assess which models are completed
preform$mod_completion = "uncompleted"
preform$mod_completion[preform$Species_Pair %in% coeff$Species_Pair] = "completed"
table(preform$mod_completion) # 260 completed mods! 

######## Import co-abundance PPC dataframes ######

## Make sure all_combos is loaded here! 

# list all files
files = list.files("results/HALF_LONG_testing_Oct_2023/PPC_dataframes/")
files = files[!grepl("OLD", files)] # remove any old data 

## Subset files for one setting, especially if multiple are in the mix. 
files = files[grepl("HALF_LONG", files)]
# files = files[grepl("LONG", files)]
# files = files[grepl("MIDDLE", files)]
# files = files[grepl("SHORT", files)]

### Noticed a typo that included generated results for species that dont spatially overlap
## and these have been removed from bdata, but not worth re-running models b/c they take long.
# so thin files to match relevant combos and exclude pointless mods via loop
### This will be irrelevant w/ updated models and testing. 

# create a new files vector to save good files
matched_files_vals <- vector("character", length(all_combos))
matched_files_plot <- vector("character", length(all_combos))

for(i in 1:length(all_combos)){
  # if there are any TRUE matches for all_combos in files w/ PPC values
  if(any(grepl(all_combos[i], files[grepl("values", files)]))){
    ## save that file name in the new vector
    matched_files_vals[i] = files[grepl("values", files)][grepl(all_combos[i], files[grepl("values", files)])]
  }
  # Repeat this process also for PPC plotting
  if(any(grepl(all_combos[i], files[grepl("plotdata", files)]))){
    ## save that file name in the new vector
    matched_files_plot[i] = files[grepl("plotdata", files)][grepl(all_combos[i], files[grepl("plotdata", files)])]
  }
}
rm(i)

## remove matched_files w/ blank names (likely b/c mod never left HPC)
matched_files_vals = matched_files_vals[matched_files_vals != ""]
matched_files_plot = matched_files_plot[matched_files_plot != ""]

## update files w/ the good ones (for consistent naming)
files = c(matched_files_vals, matched_files_plot)
rm(matched_files_vals, matched_files_plot)

# store results here
ppc_res = list()

# loop thru each file to import into list
for(i in 1:length(files)){
  
  ## read the file 
  d = read.csv(paste("results/HALF_LONG_testing_Oct_2023/PPC_dataframes/", files[i], sep = ""))
  # d = read.csv(paste("results/MIDDLE_testing_Jul_2023/PPC_dataframes/", files[i], sep = ""))
  
  
  ## save plotdata vs values in nested list. 
  if(grepl("plotdata", files[i])){
    
    ppc_res$plotdata[[i]] = d
    
  }else{
    
    ppc_res$values[[i]] = d
    
  } # end plot vs value condition
  
}
rm(d,i, files)

## bind together and inspect
ppc_plotdat = do.call(rbind, ppc_res$plotdata)
head(ppc_plotdat) 
ppc_plotdat$X = NULL # damn row names 
str(ppc_plotdat) # looks good! 

ppc_values = do.call(rbind, ppc_res$values)
head(ppc_values) 
ppc_values$X = NULL # damn row names 
str(ppc_values) # looks good! 

## quicky inspect parameter convergence for species interaction parameter
summary(ppc_values$Rhat) # looks good! some big values tho! 
nrow(ppc_values[ppc_values$Rhat > 1.2,])/nrow(ppc_values) * 100 # 21% of mods have bad interaction param

# how about sample size?
summary(coeff$n.eff[coeff$var == "Species_Interaction"]) # big range
length(coeff$n.eff[coeff$var == "Species_Interaction" & 
                     coeff$n.eff < 100]) / length(coeff$n.eff[coeff$var == "Species_Interaction"]) * 100 
# 76% of mods have small (< 100) sample size on short settings
# 36% of mods have small (< 100) sample size on middle settings
# 28% of mods have small (< 100) sample size on HALF_LONG settings


## summary of SIV
summary(ppc_values$Interaction_Estimate) # who are the crazy ones?
ppc_values[ppc_values$Interaction_Estimate < -3, ] # SUB-Cuon_alpinus~DOM-Tapirus_indicus & SUB-Panthera_pardus~DOM-Paradoxurus_hermaphroditus
ppc_values[ppc_values$Interaction_Estimate > 3, ] # nothing on positive side! 

# remove list now that were all stored in dfs. 
rm(ppc_res)

### Inspect which species pairs are present/missing, 20231127
setdiff(all_combos, ppc_plotdat$Species_Pair) #all accounted for! 


######## Import co-abundance prediction dataframes ######

#### UPDATE--> In order to save ram and files, I have stopped generating these. 
#### When presenting final results, I can re-make these on the HPC later. 
### We will only want to generate very specific plots anyway, no need for all. 

## I forgot to delete old data in predictions results folder on HPC
## and now all middle/long results are mingling with the new short ones
## but all of those are saved here anyway: results/OLD_20230617/prediction_dataframes
#### SO, delete all middle and long files! 

# ## grab all files 
# files = dir("results/prediction_dataframes/", full.names = TRUE)
# # now middle
# a = files[grepl("MIDDLE", files)]
# # now long
# b = files[grepl("LONG", files)]
# # and combine 
# rm = c(a,b)
# rm(a,b,files)
# 
# ## simple loop to delete all of them 
# for(i in 1:length(rm)){
#   file.remove(rm[i])
# }
# rm(i,rm)

## ONLY UNHASH IF YOU MUST DELETE FILES, there is no control + Z here. 

# # List all files in the coefficents folder
# file_list <- dir("results/prediction_dataframes/", full.names = TRUE)
# 
# # Extract the filenames without the date
# file_names <- basename(file_list) %>% str_remove("\\d{8}\\.csv$")
# 
# # Find the duplicated file names
# duplicated_names <- file_names[duplicated(file_names)]
# 
# # Iterate over the duplicated names
# for (name in duplicated_names) {
#   # Get the files with the current name
#   current_files <- file_list[file_names == name]
#   
#   # Find the earliest date among the duplicated files
#   earliest_date <- as.character(min(as.Date(str_extract(current_files, "\\d{8}"), format = "%Y%m%d")))
#   
#   # replace - with nothing to match existing format
#   earliest_date = gsub("-", "", earliest_date)
#   
#   # Get the files with the earliest date
#   files_to_move <- current_files[str_extract(current_files, "\\d{8}") == earliest_date]
#   
#   # Move the files to the destination folder
#   file_move(files_to_move, "results/coefficent_dataframes/OLD")
# }
# rm(current_files, earliest_date, files_to_move, name, duplicated_names, file_names, file_list)



## Make sure all_combos is loaded here! 

# list all files
files = list.files("results/MIDDLE_testing_Jul_2023/prediction_dataframes/")
files = files[!grepl("OLD", files)] # remove any old data 

## Subset files for one setting, especially if multiple are in the mix. 
# files = files[grepl("LONG", files)]
files = files[grepl("MIDDLE", files)]
# files = files[grepl("SHORT", files)]

### Noticed a typo that included generated results for species that dont spatially overlap
## and these have been removed from bdata, but not worth re-running models b/c they take long.
# so thin files to match relevant combos and exclude pointless mods via loop
### This will be irrelevant w/ updated models and testing. 

# create a new files vector to save good files
matched_files_est <- vector("character", length(all_combos))
matched_files_pred <- vector("character", length(all_combos))


for(i in 1:length(all_combos)){
  # if there are any TRUE matches for all_combos in files w/ estimated abundance
  if(any(grepl(all_combos[i], files[grepl("estimated", files)]))){
    ## save that file name in the new vector
    matched_files_est[i] = files[grepl("estimated", files)][grepl(all_combos[i], files[grepl("estimated", files)])]
  }
  # Repeat this process also for predicted abundance
  if(any(grepl(all_combos[i], files[grepl("predicted", files)]))){
    ## save that file name in the new vector
    matched_files_pred[i] = files[grepl("predicted", files)][grepl(all_combos[i], files[grepl("predicted", files)])]
  }
}
rm(i)

## remove matched_files w/ blank names (likely b/c mod never left HPC)
matched_files_est = matched_files_est[matched_files_est != ""]
matched_files_pred = matched_files_pred[matched_files_pred != ""]

## update files w/ the good ones (for consistent naming)
files = c(matched_files_est, matched_files_pred)
rm(matched_files_est, matched_files_pred)

# store results here
predict_res = list()

# loop thru each file to import into list
for(i in 1:length(files)){
  
  ## read the file 
  d = read.csv(paste("results/MIDDLE_testing_Jul_2023/prediction_dataframes/", files[i], sep = ""))
  
  ## save estimated vs predicted abundance in nested list. 
  if(grepl("estimated", files[i])){
    
    predict_res$estimated[[i]] = d
    
  }else{
    
    predict_res$predicted[[i]] = d
    
  } # end estimated vs predicted condition
  
}
rm(d,i, files)

## bind together and inspect
est_abund = do.call(rbind, predict_res$estimated)
head(est_abund) 
est_abund$X = NULL # damn row names 
str(est_abund) # looks good! 

pred_abund = do.call(rbind, predict_res$predicted)
head(pred_abund) 
pred_abund$X = NULL # damn row names 
str(pred_abund) # looks good! 

## Remove list now that data is good to go
rm(predict_res)


### Inspect which species pairs are present/missing,  20230724
setdiff(all_combos, pred_abund$Species_Pair) #14 missing models based off name spelling differences from old to new bundled data
## this wont matter anyway, dont worry about them. 


######## Import guild data and dietary preferences of large carnivores #####

#load guild data
guilds = read.csv(("data/species_traits/clean_44_species_trait_data_20230821.csv"))
head(guilds)
str(guilds)
## looks good, dietary preferences are ready to go! 

## Address two NA typo made in step 1
guilds$tiger_pref[grepl("Panthera", guilds$scientificNameStd)] = "No" # tigers dont preferentially prey upon other big cats

#
##
###
#### Combine species trait/guild data w/ model performance

## make sure all species are present
setdiff(preform$sub_species, guilds$scientificNameStd) # these match, add subordinate guild first

add = select(guilds, scientificNameStd, TrophicGuild)
names(add) = c("sub_species", "SUB_guild")

preform = merge(preform, add, by = "sub_species")
head(preform)

## now do the same for dominant speices
setdiff(preform$dom_species, guilds$scientificNameStd) # these match
add = select(guilds, scientificNameStd, TrophicGuild)
names(add) = c("dom_species", "DOM_guild")

preform = merge(preform, add, by = "dom_species")
head(preform)

## combine into a guild pair 
preform$guild_pair = paste("SUB-", preform$SUB_guild, "~",
                           "DOM-", preform$DOM_guild, sep = "")
table(preform$guild_pair)
## what can we actually work with rn
table(preform$guild_pair[preform$mod_completion == "completed"]) # a lot!

## clean up 
rm(add)

######## Track model performance and validity #######

## Will use the preform data.frame (generated in coefficents import), 
## to track which models are robust, working, and converging. 

## inspect data
head(preform)
str(preform) 
# need to add sub + dom BPV/Chat & species interaction, its significance, and the Rhat value

### inspect available data 
head(ppc_values) # its all here! 

## merge it! 
preform = merge(preform, ppc_values, by = "Species_Pair", all.x = T) # make sure to keep all values from preform b/c some could be missing from PPC! 


## inspect
head(preform) # looks good 
head(preform[preform$mod_completion == "uncompleted",]) # all NA values here for BPV,Chat,Rhat,etc for the models that havent finished 
str(preform) # looks good! 

## Add a column denoting if the BPV values converged based on OG values
preform$BPV_valid = "No"
preform$BPV_valid[preform$BPV.dom >= 0.25 & preform$BPV.dom <= 0.75 &
                    preform$BPV.sub >= 0.25 & preform$BPV.sub <= 0.75] = "Yes"
table(preform$BPV_valid[!is.na(preform$Interaction_Estimate)]) 
# mostly invalid rn, but thats w/ short settings. 
## Still mostly invalid with middle settings. 76% are not a good fit. 
## Still mostly invalid with HALF_LONG settings. 
table(preform$BPV_valid[!is.na(preform$Interaction_Estimate)])[1]/
  length(preform$Species_Pair[!is.na(preform$Interaction_Estimate)]) * 100
#69% are not a good fit. 


## Add another column that is a more liberal BPV intepreation --> explore these values!
preform$BPV_valid_liberal = "No"
preform$BPV_valid_liberal[preform$BPV.dom >= 0.15 & preform$BPV.dom <= 0.85 &
                            preform$BPV.sub >= 0.15 & preform$BPV.sub <= 0.85] = "Yes"
table(preform$BPV_valid_liberal[!is.na(preform$Interaction_Estimate)]) # yes is the majority! 
table(preform$BPV_valid_liberal[!is.na(preform$Interaction_Estimate)])[1]/
  length(preform$Species_Pair[!is.na(preform$Interaction_Estimate)]) * 100
#21.5% are not a good fit. 

## add a column denoting if overdispersion remains
preform$OD_valid = "No"
preform$OD_valid[preform$Chat.dom >= 0.98 & preform$Chat.dom <= 1.1 &
                     preform$Chat.sub >= 0.98 & preform$Chat.sub <= 1.1] = "Yes"
table(preform$OD_valid[!is.na(preform$Interaction_Estimate)]) 
## Short settings have majority with OD. TBH, better preformance than expected w/ low settings. 
## I bet this is because active_cams + RE for source are actually good for detection models... less reliance on ODRE
## Middle settings are very good, 98% have no OD
## HALF_LONG only have 3 invalid mods --> great! 

## Add another column for a more liberal OD parameter
preform$OD_valid_liberal = "No"
preform$OD_valid_liberal[preform$Chat.dom >= 0.95 & preform$Chat.dom <= 1.3 &
                           preform$Chat.sub >= 0.95 & preform$Chat.sub <= 1.3] = "Yes"
table(preform$OD_valid_liberal[!is.na(preform$Interaction_Estimate)]) # only YES!


## add a column denoting if the species interaction parameter converged 
preform$parameter_valid = "No"
preform$parameter_valid[preform$Rhat >= 0.99 & preform$Rhat <= 1.2] = "Yes"
table(preform$parameter_valid[!is.na(preform$Interaction_Estimate)]) # majority is valid!
## HALF_LONG settings are good, 79% have good interaction parameters 
table(preform$parameter_valid[!is.na(preform$Interaction_Estimate)])[1]/
  length(preform$Species_Pair[!is.na(preform$Interaction_Estimate)]) * 100
## only 21% have invalid parameters

## Which species pairs involved large carnivores significantly regulating the abundance of other critters?

## conservative settings-
preform$Species_Pair[preform$Interaction_Estimate < 0 & 
                       preform$Significance == "Significant" &
                       preform$parameter_valid == "Yes" & 
                       preform$BPV_valid == "Yes" & 
                       grepl("DOM-Large_Carnivore", preform$guild_pair)] 
## liberal settings
preform$Species_Pair[preform$Interaction_Estimate < 0 & 
                       preform$Significance == "Significant" &
                       preform$parameter_valid == "Yes" & 
                       preform$BPV_valid_liberal == "Yes" & 
                       grepl("DOM-Large_Carnivore", preform$guild_pair)] 
## much more to work with here! 


## Are there models with a good fit, but invalid SIV?
preform$Species_Pair[preform$BPV_valid == "Yes" &
                       preform$parameter_valid == "No"] # 14! This is odd

## Add the direction of the relationship
preform$direction[preform$sub_species %in% c("Panthera_tigris", "Panthera_pardus", 
                                             "Neofelis_genus", "Cuon_alpinus")] = "bottom-up"
preform$direction[preform$dom_species %in% c("Panthera_tigris", "Panthera_pardus", 
                                             "Neofelis_genus", "Cuon_alpinus")] = "top-down" # this should allow large carnivore pairings to be top-down
## address the large carnivore pairings as top-down
preform$Species_Pair[preform$guild_pair == "SUB-Large_Carnivore~DOM-Large_Carnivore"] # they are all top-down already... 


## Assign 3 levels of non-support --> unsupported_1, unsupported_2, unsupported_3
preform$support[preform$BPV_valid == "No" |
                  preform$OD_valid == "No" |
                  preform$parameter_valid == "No"] = "Unsupported_3" # lowest level of support --> bad mod.
preform$support[preform$BPV_valid == "Yes" & 
                  preform$OD_valid == "Yes" &
                  preform$parameter_valid == "Yes" &
                  preform$Significance == "Non-Significant"] = "Unsupported_2" # mid-low support --> good mod, but not important
preform$support[preform$BPV_valid == "Yes" & 
                  preform$OD_valid == "Yes" &
                  preform$parameter_valid == "Yes" &
                  preform$Significance == "Significant" &
                  preform$direction == "top-down" &
                  preform$Interaction_Estimate > 0] = "Unsupported_1" # almost supportive --> good model, significant result, but not in correct direction for hypothesis.  
preform$support[preform$BPV_valid == "Yes" & 
                  preform$OD_valid == "Yes" &
                  preform$parameter_valid == "Yes" &
                  preform$Significance == "Significant" &
                  preform$direction == "bottom-up" &
                  preform$Interaction_Estimate < 0] = "Unsupported_1" # Same as above, but applied for bottom-up direction. 
## assign which models support our hypothesis
preform$support[preform$BPV_valid == "Yes" & 
                  preform$OD_valid == "Yes" &
                  preform$parameter_valid == "Yes" &
                  preform$Significance == "Significant" &
                  preform$direction == "top-down" &
                  preform$Interaction_Estimate <= 0] = "Supported" # good model, significant result, in correct direction for hypothesis.  
preform$support[preform$BPV_valid == "Yes" & 
                  preform$OD_valid == "Yes" &
                  preform$parameter_valid == "Yes" &
                  preform$Significance == "Significant" &
                  preform$direction == "bottom-up" &
                  preform$Interaction_Estimate >= 0] = "Supported" # Same as above, but applied for bottom-up direction. 
## inspect
table(preform$support);anyNA(preform$support) # must be F for NA! 
# who is supported?
preform[preform$support == "Supported", ] # mostly bottom up! 


#### Do the same thing, but with the liberal values
## Assign 3 levels of non-support --> unsupported_1, unsupported_2, unsupported_3
preform$support_liberal[preform$BPV_valid_liberal == "No" |
                          preform$OD_valid_liberal == "No" |
                          preform$parameter_valid == "No"] = "Unsupported_3" # lowest level of support --> bad mod.
preform$support_liberal[preform$BPV_valid_liberal == "Yes" & 
                          preform$OD_valid_liberal == "Yes" &
                          preform$parameter_valid == "Yes" &
                          preform$Significance == "Non-Significant"] = "Unsupported_2" # mid-low support --> good mod, but not important
preform$support_liberal[preform$BPV_valid_liberal == "Yes" & 
                          preform$OD_valid_liberal == "Yes" &
                          preform$parameter_valid == "Yes" &
                          preform$Significance == "Significant" &
                          preform$direction == "top-down" &
                          preform$Interaction_Estimate > 0] = "Unsupported_1" # almost supportive --> good model, significant result, but not in correct direction for hypothesis.  
preform$support_liberal[preform$BPV_valid_liberal == "Yes" & 
                          preform$OD_valid_liberal == "Yes" &
                          preform$parameter_valid == "Yes" &
                          preform$Significance == "Significant" &
                          preform$direction == "bottom-up" &
                          preform$Interaction_Estimate < 0] = "Unsupported_1" # Same as above, but applied for bottom-up direction. 
## assign which models support our hypothesis
preform$support_liberal[preform$BPV_valid_liberal == "Yes" & 
                          preform$OD_valid_liberal == "Yes" &
                          preform$parameter_valid == "Yes" &
                          preform$Significance == "Significant" &
                          preform$direction == "top-down" &
                          preform$Interaction_Estimate <= 0] = "Supported" # good model, significant result, in correct direction for hypothesis.  
preform$support_liberal[preform$BPV_valid_liberal == "Yes" & 
                          preform$OD_valid_liberal == "Yes" &
                          preform$parameter_valid == "Yes" &
                          preform$Significance == "Significant" &
                          preform$direction == "bottom-up" &
                          preform$Interaction_Estimate >= 0] = "Supported" # Same as above, but applied for bottom-up direction. 
## inspect
table(preform$support_liberal[!is.na(preform$Interaction_Estimate)]);anyNA(preform$support[!is.na(preform$Interaction_Estimate)]) 
# must be F for NA! 
## MANY more fall into the supported category, and unsupported _3 and _2 are very similar 
# who is supported?
preform[preform$support_liberal == "Supported", ] # mostly bottom up! Also lots of random species pairings here... 


### This is a super useful dataframe, will need to share this in the supplmentary tables section
## save it! 
# grab the date
day<-str_sub(Sys.Date(),-2)
month<-str_sub(Sys.Date(),-5,-4)
year<-str_sub(Sys.Date(),-10,-7)
date = paste(year,month,day, sep = "")
rm(year,month,day)

# write.csv(preform, paste("results/summarized_results/model_performance_and_SIVs_", date, ".csv", sep =""))


######## Prepare plotting data to visualize all SIVs ######

### Prepare the data for plotting-

## grab top-down mods
td_plot_dat = preform[grepl("DOM-Large_Carnivore", preform$guild_pair),]

## create a new col that is the relevant species
td_plot_dat$rel_species = td_plot_dat$dom_species

## double down for all_carns 
a = td_plot_dat
a$rel_species = "All_large_carnivores"
td_plot_dat = rbind(td_plot_dat, a)

## add direction
td_plot_dat$direction = "top-down"


## grab bottom-up mods
bu_plot_dat = preform[grepl("SUB-Large_Carnivore", preform$guild_pair),]

## create a new col that is the relevant species
bu_plot_dat$rel_species = bu_plot_dat$sub_species

## double down for all_canrs 
a = bu_plot_dat
a$rel_species = "All_large_carnivores"
bu_plot_dat = rbind(bu_plot_dat, a)

## add direction
bu_plot_dat$direction = "bottom-up"

## remove any models where there is also a dominant large carnivore
bu_plot_dat = bu_plot_dat[!grepl("DOM-Large_Carnivore", bu_plot_dat$guild_pair),]

## and combine them
plot_dat = rbind(td_plot_dat, bu_plot_dat)
rm(td_plot_dat, bu_plot_dat)

### verify all large-carnivore pairings are top-down only
table(plot_dat$direction[plot_dat$guild_pair == "SUB-Large_Carnivore~DOM-Large_Carnivore"]) # good! 

## set order manually b/c order is silly on graph
plot_dat$rel_species = factor(plot_dat$rel_species, levels = c("Neofelis_genus","Cuon_alpinus",
                                                               "Panthera_pardus","Panthera_tigris",
                                                                 "All_large_carnivores"))

## create a new column for simple groupings: top-down, bottom-up and unsupported
plot_dat$support_simple = plot_dat$direction
plot_dat$support_simple_liberal = plot_dat$direction

### edit support_simple to allow for unsupported mods
plot_dat$support_simple[grepl("Unsupported", plot_dat$support)] = "unsupported"
table(plot_dat$support_simple)

### and do the same for liberal values
plot_dat$support_simple_liberal[grepl("Unsupported", plot_dat$support_liberal)] = "unsupported"
table(plot_dat$support_simple_liberal)


## remove NA values from plot_dat if present --> from non-completed mods
if(any(is.na(plot_dat$Interaction_Estimate))){
  
  plot_dat = plot_dat[! is.na(plot_dat$Interaction_Estimate), ]
  
}

## Check if there are extreme species interaction values (> |2|)
summary(plot_dat$Interaction_Estimate) 
plot_dat[plot_dat$Interaction_Estimate < -3,] 
## 2 Problematic species pairs: SUB-Panthera_pardus~DOM-Paradoxurus_hermaphroditus & SUB-Cuon_alpinus~DOM-Tapirus_indicus 
## both are a good model fit, but unsupported_1 (SIV in wrong direction). Not preferred prey anyway

## check for even less extreme values
plot_dat[plot_dat$Interaction_Estimate < -2,] 
## SUB-Macaca_fascicularis~DOM-Cuon_alpinus, SUB-Trichys_fasciculata~DOM-Panthera_tigris
## but these are supported top-down interactions!! They are large effects, but not in preferred prey 

## for plotting sake, remove anything < -3
# plot_dat = plot_dat[plot_dat$Interaction_Estimate > -3,]
## NOT REMOVING THESE FROM THE DATA --> these will only be supplementary figures anyway b/c it wont make it into preferred prey figs. 

## set factors so unsupported is in the back 
plot_dat$support_simple = factor(plot_dat$support_simple, levels = c("top-down", "bottom-up", "unsupported"))
plot_dat$support_simple_liberal = factor(plot_dat$support_simple_liberal, levels = c("top-down", "bottom-up", "unsupported"))

## catch a tricky typo! 
plot_dat$Species_Pair[plot_dat$support_simple_liberal == "bottom-up" & plot_dat$Interaction_Estimate < 0] # these should be top-down! 
plot_dat$support_simple_liberal[plot_dat$Species_Pair == "SUB-Neofelis_genus~DOM-Cuon_alpinus"] = "top-down"
plot_dat$support_liberal[plot_dat$Species_Pair == "SUB-Neofelis_genus~DOM-Cuon_alpinus"] = "top-down"

## make the support column more informative
plot_dat$support[plot_dat$support == "Supported" & plot_dat$direction == "top-down"] = "top-down"
plot_dat$support_liberal[plot_dat$support_liberal == "Supported" & plot_dat$direction == "top-down"] = "top-down"
plot_dat$support[plot_dat$support == "Supported" & plot_dat$direction == "bottom-up"] = "bottom-up"
plot_dat$support_liberal[plot_dat$support_liberal == "Supported" & plot_dat$direction == "bottom-up"] = "bottom-up"

## make letters all lowercase
plot_dat$support= tolower(plot_dat$support)
plot_dat$support_liberal = tolower(plot_dat$support_liberal)

## and set factor level 
plot_dat$support = factor(plot_dat$support, levels = c("top-down", "bottom-up", 
                                                       "unsupported_1", "unsupported_2",
                                                       "unsupported_3"))
plot_dat$support_liberal = factor(plot_dat$support_liberal, levels = c("top-down", "bottom-up", 
                                                                       "unsupported_1", "unsupported_2",
                                                                       "unsupported_3"))

### set up colors for all values 
new_colors = c("unsupported_1" = "gray87",
               "unsupported_2" = "gray55",
               "unsupported_3" = "gray32",
               "unsupported" = "gray32",
               "top-down" = "darkorange4",
               "bottom-up" = "springgreen4")


###### Explore the different kinds of plots that can be generated, using all large carn as the example
{
  # ## first grab today's date for saving plots
  # day<-str_sub(Sys.Date(),-2)
  # month<-str_sub(Sys.Date(),-5,-4)
  # year<-str_sub(Sys.Date(),-10,-7)
  # date = paste(year,month,day, sep = "")
  # rm(year,month,day)
  
  ### load ggridges library for these to work. 
  
  # # The classic ridge line I am used to 
  # p = 
  #   ggplot(ridge_dat[ridge_dat$rel_species ==  "All_large_carnivores", ], 
  #          aes(x = Interaction_Estimate, y = rel_species, fill = direction)) +
  #   geom_density_ridges(scale=1, alpha = .8, jittered_points = TRUE, point_alpha=.5)+
  #   geom_vline(aes(xintercept = 0), color = "black") +
  #   theme_test() +
  #   scale_fill_manual(values = new_colors) +
  #   labs(x = "Species Interaction Value", y = NULL, fill = NULL, color = NULL, title = "all mods")
  # # ggsave(paste("explore/SIV_distribution_ridgeline_non_scaled_",date, ".png", sep = ""),
  # #        p, width = 7, height = 3, units = "in")
  
  # # the 'scaled' version w/ ridgeline
  # p =
  #   ggplot(ridge_dat2[ridge_dat2$rel_species ==  "All_large_carnivores", ], 
  #          aes(x = Interaction_Estimate, y = rel_species, fill = direction, height = stat(density))) +
  #   geom_density_ridges(stat = "density", scale=1, alpha = .8)+ # default bandwidth is 1 (adjust = 1 in code), high makes distribution more round. 
  #   geom_vline(aes(xintercept = 0), color = "black") +
  #   theme_test() +
  #   scale_fill_manual(values = new_colors) +
  #   labs(x = "Species Interaction Value", y = NULL, fill = NULL, color = NULL, title = "all mods")
  # # ggsave(paste("explore/SIV_distribution_ridgeline_scaled_",date, ".png", sep = ""),
  # #        p, width = 7, height = 3, units = "in")
  
  ## A frequency distribution, colored per category
  # p= 
  ggplot(plot_dat[plot_dat$rel_species == "All_large_carnivores", ], aes(x = Interaction_Estimate, fill = direction_liberal))+
    geom_histogram(aes(y=after_stat(count)), colour="black", bins = 20, alpha = 1)+#, position = "dodge")+
    # geom_density(alpha=.6, bw = .1)+
    theme_test()+
    scale_fill_manual(values = new_colors) +
    labs(x = "Species Interaction Value", y = NULL, fill = NULL, color = NULL, title = "all mods")
  # ggsave(paste("explore/SIV_distribution_frequency_plot_",date, ".png", sep = ""),
  #        p, width = 7, height = 3, units = "in")
  
  ## A frequency distribution with density overlaid
  # p=
  ggplot(plot_dat[plot_dat$rel_species == "All_large_carnivores", ], aes(x = Interaction_Estimate, fill = direction))+
    geom_histogram(aes(y=after_stat(count)), colour="black", bins = 80, alpha = 1)+
    geom_density(alpha=.8, bw = .1)+
    theme_test()+
    scale_fill_manual(values = new_colors) +
    labs(x = "Species Interaction Value", y = "Frequency", fill = NULL, color = NULL, title = "all mods")
  # ggsave(paste("explore/SIV_distribution_frequency_plot_with_density_",date, ".png", sep = ""),
  #        p, width = 7, height = 3, units = "in")
  
  ## a regular gg density plot
  # p=
  ggplot(plot_dat[plot_dat$rel_species == "All_large_carnivores", ], aes(x = Interaction_Estimate, fill = direction))+
    geom_density(bw = .1, alpha=.8, )+
    theme_test()+
    scale_fill_manual(values = new_colors) +
    labs(x = "Species Interaction Value", y = "Density", fill = NULL, color = NULL, title = "all mods")
  # ggsave(paste("explore/SIV_distribution_density_plot_",date, ".png", sep = ""),
  #        p, width = 7, height = 3, units = "in")
  
  ## density plot with only one category, so the area under the curve MUST be 1
  # p=
  ggplot(plot_dat[plot_dat$rel_species == "All_large_carnivores", ], aes(x = Interaction_Estimate))+#, fill = direction))+
    geom_density(bw = .1, alpha=.8, )+
    # geom_histogram(aes(y=after_stat(count)), colour="black", bins = 80, alpha = 1)+
    theme_test()+
    # scale_fill_manual(values = new_colors) +
    labs(x = "Species Interaction Value", y = "Density", fill = NULL, color = NULL, title = "all mods")
  # ggsave(paste("explore/SIV_distribution_density_plot_one_level_",date, ".png", sep = ""),
  #        p, width = 7, height = 3, units = "in")
  
  
  ### Try to overlay the one-category density plot with the frequency plot 
  
  # Frequency plot 
  # f_plot = 
  ggplot(plot_dat[plot_dat$rel_species == "All_large_carnivores", ], 
         aes(x = Interaction_Estimate, fill = direction))+
    geom_histogram(aes(y=after_stat(count)), colour="black", bins = 80, alpha = 1)+ #, position = "dodge")+
    scale_fill_manual(values = new_colors) +
    geom_vline(aes(xintercept = 0), linetype = "dashed", alpha = 0.6)+
    theme_classic()+
    labs(x = "Species Interaction Value", y = "Frequency") 
  
  ### LOAD ggpubr and cowplot libraries to make this work. 
  # # Density plot
  # d_plot = ggdensity(plot_dat[plot_dat$rel_species == "All_large_carnivores", ], 
  #                    
  #                    x = "Interaction_Estimate")+
  #   geom_vline(aes(xintercept = 0), linetype = "dashed", alpha = 0.6)+
  #   theme_classic()+
  #   scale_y_continuous(expand = expansion(mult = c(0, 0.05)), position = "right")+
  #   rremove("x.axis")+
  #   rremove("xlab") +
  #   rremove("x.text") +
  #   rremove("x.ticks")
  # 
  # ## put them together! 
  # # f_plot + d_plot
  # # plot_grid(f_plot, d_plot, align = "hv", axis = "tblr")
  # # plot_grid(f_plot, d_plot, align = "v", axis = "l")
  # ggarrange(f_plot, d_plot, ncol = 1, nrow = 2)
  # 
  # ## This code is supposed to overlay the two plots, 
  # ## but will not work --> 2nd plot always overwrite the first. 
  # aligned_plots = align_plots(f_plot,d_plot, align="hv", axis="tblr")
  # ggdraw(aligned_plots[[1]]) + draw_plot(aligned_plots[[2]])
  # 
  # rm(f_plot, d_plot, aligned_plots, p, a, date)
  # rm(date,p,f_plot)
}

### histogram is the final move! 
ggplot(plot_dat, aes(x = Interaction_Estimate, fill = support_simple_liberal))+
  geom_histogram(aes(y=after_stat(count)), colour="black", binwidth = .08, alpha = 1)+#, position = "dodge")+
  # theme_minimal()+
  # theme_test()+
  theme_classic()+
  geom_vline(aes(xintercept = 0), linetype = "dashed", color = "firebrick4", alpha = .5)+
  scale_fill_manual(values = new_colors) +
  # coord_cartesian(xlim = c(min_val, max_val))+
  labs(x = "Species Interaction Value", y = NULL, fill = NULL, color = NULL)



############# Visualize unsupported SIVs using frequency plots ######

#
##
###
#### Top-down non-support first 

## subset rel_species for top-down models 
table(plot_dat$direction)
td_dat = plot_dat[plot_dat$direction == "top-down",]
table(td_dat$support); table(td_dat$support_liberal)

## frequency plot version --> make one for each species 
plot_list = list()
plot_list_lib = list()
order = c("All_large_carnivores","Panthera_tigris","Panthera_pardus", 
          "Cuon_alpinus", "Neofelis_genus")
# grab min and max vals 
min_val = min(td_dat$Interaction_Estimate[grepl("unsupport", td_dat$support) ])
max_val = max(td_dat$Interaction_Estimate[grepl("unsupport", td_dat$support) ])
min_val_lib = min(td_dat$Interaction_Estimate[grepl("unsupport", td_dat$support_liberal) ])
max_val_lib = max(td_dat$Interaction_Estimate[grepl("unsupport", td_dat$support_liberal) ])

for(i in 1:length(order)){
  
  # subset the data
  d = td_dat[grepl("unsupport", td_dat$support) &
               td_dat$rel_species == order[i], ]
  
  # make the plot 
  p = ggplot(d, aes(x = Interaction_Estimate, fill = support))+
    geom_histogram(aes(y=after_stat(count)), colour="black", binwidth = .1, alpha = 1)+#, position = "dodge")+
    # theme_minimal()+
    # theme_test()+
    theme_classic()+
    geom_vline(aes(xintercept = 0), linetype = "dashed", color = "firebrick4", alpha = .5)+
    scale_fill_manual(values = new_colors) +
    coord_cartesian(xlim = c(min_val, max_val))+
    labs(x = "Species Interaction Value", y = NULL, fill = NULL, color = NULL, title = order[i])
  
  # save it! 
  plot_list[[i]] = p
  
  ### Repeat for liberal values
  
  # subset the data
  d = td_dat[grepl("unsupport", td_dat$support_liberal) &
               td_dat$rel_species == order[i], ]
  
  # make the plot 
  p = ggplot(d, aes(x = Interaction_Estimate, fill = support_liberal))+
    geom_histogram(aes(y=after_stat(count)), colour="black", binwidth = .1, alpha = 1)+#, position = "dodge")+
    # theme_minimal()+
    # theme_test()+
    theme_classic()+
    geom_vline(aes(xintercept = 0), linetype = "dashed", color = "firebrick4", alpha = .5)+
    scale_fill_manual(values = new_colors) +
    coord_cartesian(xlim = c(min_val_lib, max_val_lib))+
    labs(x = "Species Interaction Value", y = NULL, fill = NULL, color = NULL, title = order[i])
  
  # save it! 
  plot_list_lib[[i]] = p
  
}
rm(d,p,i, min_val, max_val, min_val_lib, max_val_lib)


## inspect
plot_list[[1]]

## arrange them vertically
p = plot_grid(plotlist = plot_list, align = "v",
              ncol = 1, nrow = length(plot_list))
p

## grab today's date for saving plots
day<-str_sub(Sys.Date(),-2)
month<-str_sub(Sys.Date(),-5,-4)
year<-str_sub(Sys.Date(),-10,-7)
date = paste(year,month,day, sep = "")
rm(year,month,day)

# and save it!
# ggsave(paste("figures/Grouped Histograms/Histogram_UNSUPPORTED_&_CONSERVATIVE_top-down_SIV_HALF_LONG_",date, ".png", sep = ""),
#                p, width = 10, height = 10, units = "in")

## repeat for the liberal plots 
p = plot_grid(plotlist = plot_list_lib, align = "v",
              ncol = 1, nrow = length(plot_list))
p

# and save it!
# ggsave(paste("figures/Grouped Histograms/Histogram_UNSUPPORTED_&_LIBERAL_top-down_SIV_HALF_LONG_",date, ".png", sep = ""),
#                p, width = 10, height = 10, units = "in")


## summarize percentages in each non-supported category, which will be added to graphs manually in ppt
# first how many mods per species per support
top_down_nonsupport = ddply(td_dat[grepl("unsupported", td_dat$support), ], .(rel_species, support), summarize,
                            num_lvl = length(Species_Pair))
# then total number of mods per species
top_down_nonsupport$num_total[top_down_nonsupport$rel_species == "Neofelis_genus"] = 
  sum(top_down_nonsupport$num_lvl[top_down_nonsupport$rel_species == "Neofelis_genus"])
top_down_nonsupport$num_total[top_down_nonsupport$rel_species == "Cuon_alpinus"] = 
  sum(top_down_nonsupport$num_lvl[top_down_nonsupport$rel_species == "Cuon_alpinus"])
top_down_nonsupport$num_total[top_down_nonsupport$rel_species == "Panthera_pardus"] = 
  sum(top_down_nonsupport$num_lvl[top_down_nonsupport$rel_species == "Panthera_pardus"])
top_down_nonsupport$num_total[top_down_nonsupport$rel_species == "Panthera_tigris"] = 
  sum(top_down_nonsupport$num_lvl[top_down_nonsupport$rel_species == "Panthera_tigris"])
top_down_nonsupport$num_total[top_down_nonsupport$rel_species == "All_large_carnivores"] = 
  sum(top_down_nonsupport$num_lvl[top_down_nonsupport$rel_species == "All_large_carnivores"])
# finally calulate the percent-
top_down_nonsupport$percent = top_down_nonsupport$num_lvl / top_down_nonsupport$num_total * 100
## add the PPC settings
top_down_nonsupport$PPC_setting = "conservative"
## and the direction
top_down_nonsupport$direction = "top-down"


### remake w/ liberal values
top_down_nonsupport_lib = ddply(td_dat[grepl("unsupported", td_dat$support_liberal), ], .(rel_species, support_liberal), summarize,
                            num_lvl = length(Species_Pair))
# then total number of mods per species
top_down_nonsupport_lib$num_total[top_down_nonsupport_lib$rel_species == "Neofelis_genus"] = 
  sum(top_down_nonsupport_lib$num_lvl[top_down_nonsupport_lib$rel_species == "Neofelis_genus"])
top_down_nonsupport_lib$num_total[top_down_nonsupport_lib$rel_species == "Cuon_alpinus"] = 
  sum(top_down_nonsupport_lib$num_lvl[top_down_nonsupport_lib$rel_species == "Cuon_alpinus"])
top_down_nonsupport_lib$num_total[top_down_nonsupport_lib$rel_species == "Panthera_pardus"] = 
  sum(top_down_nonsupport_lib$num_lvl[top_down_nonsupport_lib$rel_species == "Panthera_pardus"])
top_down_nonsupport_lib$num_total[top_down_nonsupport_lib$rel_species == "Panthera_tigris"] = 
  sum(top_down_nonsupport_lib$num_lvl[top_down_nonsupport_lib$rel_species == "Panthera_tigris"])
top_down_nonsupport_lib$num_total[top_down_nonsupport_lib$rel_species == "All_large_carnivores"] = 
  sum(top_down_nonsupport_lib$num_lvl[top_down_nonsupport_lib$rel_species == "All_large_carnivores"])
# finally calulate the percent-
top_down_nonsupport_lib$percent = top_down_nonsupport_lib$num_lvl / top_down_nonsupport_lib$num_total * 100
## add the PPC settings
top_down_nonsupport_lib$PPC_setting = "liberal"
## and the direction
top_down_nonsupport_lib$direction = "top-down"

## standardize support column for easy rbinding
names(top_down_nonsupport)[2] = "support_levels"
names(top_down_nonsupport_lib)[2] = "support_levels"

## and combine 
top_down_nonsupport = rbind(top_down_nonsupport, top_down_nonsupport_lib)
## will combien with bottom-up later and save it! 



## keep it clean 
rm(a, p, plot_list, plot_list_lib, top_down_nonsupport_lib)


#
##
###
#### Bottom-up non-support next 

## subset plot_dat for bottom-up models 
table(plot_dat$direction)
bu_dat = plot_dat[plot_dat$direction == "bottom-up",]
table(bu_dat$support); table(bu_dat$support_liberal)

## frequency plot version --> make one for each species 
plot_list = list()
plot_list_lib = list()
order = c("All_large_carnivores","Panthera_tigris","Panthera_pardus", 
          "Cuon_alpinus", "Neofelis_genus")
# grab min and max vals 
min_val = min(bu_dat$Interaction_Estimate[grepl("unsupport", bu_dat$support)])
max_val = max(bu_dat$Interaction_Estimate[grepl("unsupport", bu_dat$support)])
min_val_lib = min(bu_dat$Interaction_Estimate[grepl("unsupport", bu_dat$support_liberal)])
max_val_lib = max(bu_dat$Interaction_Estimate[grepl("unsupport", bu_dat$support_liberal)])

for(i in 1:length(order)){
  
  # subset the data
  d = bu_dat[grepl("unsupport", bu_dat$support) &
               bu_dat$rel_species == order[i], ]
  
  # make the plot 
  p = ggplot(d, aes(x = Interaction_Estimate, fill = support))+
    geom_histogram(aes(y=after_stat(count)), colour="black", binwidth = .1, alpha = 1)+#, position = "dodge")+
    # theme_minimal()+
    # theme_test()+
    theme_classic()+
    geom_vline(aes(xintercept = 0), linetype = "dashed", color = "firebrick4", alpha = .5)+
    scale_fill_manual(values = new_colors) +
    coord_cartesian(xlim = c(min_val, max_val))+
    labs(x = "Species Interaction Value", y = NULL, fill = NULL, color = NULL, title = order[i])
  
  # save it! 
  plot_list[[i]] = p
  
  ### Repeat for liberal values
  
  # subset the data
  d = bu_dat[grepl("unsupport", bu_dat$support_liberal) &
               bu_dat$rel_species == order[i], ]
  
  # make the plot 
  p = ggplot(d, aes(x = Interaction_Estimate, fill = support_liberal))+
    geom_histogram(aes(y=after_stat(count)), colour="black", binwidth = .1, alpha = 1)+#, position = "dodge")+
    # theme_minimal()+
    # theme_test()+
    theme_classic()+
    geom_vline(aes(xintercept = 0), linetype = "dashed", color = "firebrick4", alpha = .5)+
    scale_fill_manual(values = new_colors) +
    coord_cartesian(xlim = c(min_val_lib, max_val_lib))+
    labs(x = "Species Interaction Value", y = NULL, fill = NULL, color = NULL, title = order[i])
  
  # save it! 
  plot_list_lib[[i]] = p
  
}
rm(d,p,i, min_val, max_val, min_val_lib, max_val_lib)

## arrange them vertically 
p = plot_grid(plotlist = plot_list, align = "v", 
              ncol = 1, nrow = length(plot_list))
# and save it! 
# ggsave(paste("figures/Grouped Histograms/Histogram_UNSUPPORTED_&_CONSERVATIVE_bottom-up_SIV_HALF_LONG_",date, ".png", sep = ""),
#                p, width = 10, height = 10, units = "in")

## repeat for the liberal plots 
p = plot_grid(plotlist = plot_list_lib, align = "v",
              ncol = 1, nrow = length(plot_list))
p

# and save it!
# ggsave(paste("figures/Grouped Histograms/Histogram_UNSUPPORTED_&_LIBERAL_bottom-up_SIV_HALF_LONG_",date, ".png", sep = ""),
#                p, width = 10, height = 10, units = "in")


## summarize percentages in each non-supported category, which will be added to graphs manually in ppt
# first how many mods per species per support
bottom_up_nonsupport = ddply(bu_dat[grepl("unsupport", bu_dat$support), ], .(rel_species, support), summarize,
                            num_lvl = length(Species_Pair))
# then total number of mods per species
bottom_up_nonsupport$num_total[bottom_up_nonsupport$rel_species == "Neofelis_genus"] = 
  sum(bottom_up_nonsupport$num_lvl[bottom_up_nonsupport$rel_species == "Neofelis_genus"])
bottom_up_nonsupport$num_total[bottom_up_nonsupport$rel_species == "Cuon_alpinus"] = 
  sum(bottom_up_nonsupport$num_lvl[bottom_up_nonsupport$rel_species == "Cuon_alpinus"])
bottom_up_nonsupport$num_total[bottom_up_nonsupport$rel_species == "Panthera_pardus"] = 
  sum(bottom_up_nonsupport$num_lvl[bottom_up_nonsupport$rel_species == "Panthera_pardus"])
bottom_up_nonsupport$num_total[bottom_up_nonsupport$rel_species == "Panthera_tigris"] = 
  sum(bottom_up_nonsupport$num_lvl[bottom_up_nonsupport$rel_species == "Panthera_tigris"])
bottom_up_nonsupport$num_total[bottom_up_nonsupport$rel_species == "All_large_carnivores"] = 
  sum(bottom_up_nonsupport$num_lvl[bottom_up_nonsupport$rel_species == "All_large_carnivores"])
# calulate the percent-
bottom_up_nonsupport$percent = bottom_up_nonsupport$num_lvl / bottom_up_nonsupport$num_total * 100
## add the PPC settings
bottom_up_nonsupport$PPC_setting = "conservative"
## and the direction
bottom_up_nonsupport$direction = "bottom-up"


#### Repeat w/ liberal groupings 
## summarize percentages in each non-supported category, which will be added to graphs manually in ppt
# first how many mods per species per support
bottom_up_nonsupport_lib = ddply(bu_dat[grepl("unsupport", bu_dat$support_liberal), ], .(rel_species, support_liberal), summarize,
                             num_lvl = length(Species_Pair))
# then total number of mods per species
bottom_up_nonsupport_lib$num_total[bottom_up_nonsupport_lib$rel_species == "Neofelis_genus"] = 
  sum(bottom_up_nonsupport_lib$num_lvl[bottom_up_nonsupport_lib$rel_species == "Neofelis_genus"])
bottom_up_nonsupport_lib$num_total[bottom_up_nonsupport_lib$rel_species == "Cuon_alpinus"] = 
  sum(bottom_up_nonsupport_lib$num_lvl[bottom_up_nonsupport_lib$rel_species == "Cuon_alpinus"])
bottom_up_nonsupport_lib$num_total[bottom_up_nonsupport_lib$rel_species == "Panthera_pardus"] = 
  sum(bottom_up_nonsupport_lib$num_lvl[bottom_up_nonsupport_lib$rel_species == "Panthera_pardus"])
bottom_up_nonsupport_lib$num_total[bottom_up_nonsupport_lib$rel_species == "Panthera_tigris"] = 
  sum(bottom_up_nonsupport_lib$num_lvl[bottom_up_nonsupport_lib$rel_species == "Panthera_tigris"])
bottom_up_nonsupport_lib$num_total[bottom_up_nonsupport_lib$rel_species == "All_large_carnivores"] = 
  sum(bottom_up_nonsupport_lib$num_lvl[bottom_up_nonsupport_lib$rel_species == "All_large_carnivores"])
# finally calulate the percent-
bottom_up_nonsupport_lib$percent = bottom_up_nonsupport_lib$num_lvl / bottom_up_nonsupport_lib$num_total * 100
## add the PPC settings
bottom_up_nonsupport_lib$PPC_setting = "liberal"
## and the direction
bottom_up_nonsupport_lib$direction = "bottom-up"

## standardize col names for easier rbinding
names(bottom_up_nonsupport)[2] = "support_levels"
names(bottom_up_nonsupport_lib)[2] = "support_levels"

## rbind them together
bottom_up_nonsupport = rbind(bottom_up_nonsupport, bottom_up_nonsupport_lib)

## combine with top-down 
nonsupport = rbind(top_down_nonsupport, bottom_up_nonsupport)
table(nonsupport$PPC_setting)
table(nonsupport$direction)

## order by PPC setting before saving
nonsupport = nonsupport[order(nonsupport$PPC_setting),]

# save it for fig making 
# write.csv(nonsupport, paste("results/summarized_results/percentages_for_all_unsupported_histograms_", date, ".csv", sep = ""), row.names = F)


## keep it clean 
rm(p, plot_list, plot_list_lib, bottom_up_nonsupport_lib, bottom_up_nonsupport, 
   top_down_nonsupport, bu_dat, td_dat)


############# Visualize ALL SIVs using frequency plots ######

#
##
###
#### Plots that show ALL models (i.e., top-down, bottom-up and all levels of unsupported)

### Verify plot_data has the correct info 
table(plot_dat$rel_species[plot_dat$direction == "top-down"])
table(plot_dat$rel_species[plot_dat$direction == "bottom-up"]) # all critters accounted for 
## double check support
table(plot_dat$support[plot_dat$direction == "top-down"])
table(plot_dat$support_liberal[plot_dat$direction == "top-down"])
table(plot_dat$support[plot_dat$direction == "bottom-up"])
table(plot_dat$support_liberal[plot_dat$direction == "bottom-up"])


## Frequency plot version
plot_list = list()
plot_list_lib = list()

order = c("All_large_carnivores","Panthera_tigris","Panthera_pardus", 
          "Cuon_alpinus", "Neofelis_genus")
# grab min and max vals 
min_val = min(plot_dat$Interaction_Estimate)
max_val = max(plot_dat$Interaction_Estimate)

for(i in 1:length(order)){
  
  # subset the data
  d = plot_dat[plot_dat$rel_species == order[i], ]
  
  # make the plot 
  p = ggplot(d, aes(x = Interaction_Estimate, fill = support_simple))+
    geom_histogram(aes(y=after_stat(count)), colour="black", binwidth = .1, alpha = 1)+#, position = "dodge")+
    # theme_minimal()+
    # theme_test()+
    theme_classic()+
    geom_vline(aes(xintercept = 0), linetype = "dashed", color = "firebrick4", alpha = .5)+
    scale_fill_manual(values = new_colors) +
    coord_cartesian(xlim = c(min_val, max_val))+
    labs(x = "Species Interaction Value", y = NULL, fill = NULL, color = NULL, title = order[i])
  
  # save it! 
  plot_list[[i]] = p
  
  ### Repeat for liberal values

  # make the plot 
  p = ggplot(d, aes(x = Interaction_Estimate, fill = support_simple_liberal))+
    geom_histogram(aes(y=after_stat(count)), colour="black", binwidth = .1, alpha = 1)+#, position = "dodge")+
    # theme_minimal()+
    # theme_test()+
    theme_classic()+
    geom_vline(aes(xintercept = 0), linetype = "dashed", color = "firebrick4", alpha = .5)+
    scale_fill_manual(values = new_colors) +
    coord_cartesian(xlim = c(min_val, max_val))+
    labs(x = "Species Interaction Value", y = NULL, fill = NULL, color = NULL, title = order[i])
  
  # save it! 
  plot_list_lib[[i]] = p
}
rm(d,p,i, min_val, max_val)

## arrange them vertically 
p = plot_grid(plotlist = plot_list, align = "v", 
              ncol = 1, nrow = length(plot_list)) 
# and save it! 
# ggsave(paste("figures/Grouped Histograms/Histogram_SUPP_&_NONSUPP_&_CONSERVATIVE_SIVs_all_mods_HALF_LONG_",date, ".png", sep = ""),
#                p, width = 8, height = 10, units = "in")

## arrange them vertically 
p = plot_grid(plotlist = plot_list_lib, align = "v", 
              ncol = 1, nrow = length(plot_list_lib)) 
p

# and save it! 
# ggsave(paste("figures/Grouped Histograms/Histogram_SUPP_&_NONSUPP_&_LIBERAL_SIVs_all_mods_HALF_LONG_",date, ".png", sep = ""),
# p, width = 8, height = 10, units = "in")

### MAKE CUSTOM SIZED PLOTS
{
  
  # ## so all_large carnivore is long, and other species are small underneath. 
  # p = ggplot(plot_dat[plot_dat$rel_species == "All_large_carnivores",], aes(x = Interaction_Estimate, fill = support))+
  #   geom_histogram(aes(y=after_stat(count)), colour="black", binwidth = .1, alpha = 1)+#, position = "dodge")+
  #   theme_classic()+
  #   geom_vline(aes(xintercept = 0), linetype = "dashed", color = "firebrick4", alpha = .5)+
  #   scale_fill_manual(values = new_colors) +
  #   labs(x = "Species Interaction Value", y = NULL, fill = NULL, color = NULL, title = "All_large_carnivores all models")+
  #   theme(axis.text.x = element_text(size = 18),
  #         axis.text.y = element_text(size = 18),
  #         axis.text = element_text(color = "black"),
  #         text = element_text(family = "Helvetica"))
  # ## save it! 
  # # ggsave(paste("figures/Grouped Histograms/custom_size_histo/Histogram_all_large_carnivores_multi_unsupported_and_supported_species_interaction_all_mods_HALF_LONG_",date, ".png", sep = ""),
  # #        p, width = 10, height = 4, units = "in")
  # 
  # ## and same with liberal!! 
  # p = ggplot(plot_dat[plot_dat$rel_species == "All_large_carnivores",], aes(x = Interaction_Estimate, fill = support_liberal))+
  #   geom_histogram(aes(y=after_stat(count)), colour="black", binwidth = .1, alpha = 1)+#, position = "dodge")+
  #   theme_classic()+
  #   geom_vline(aes(xintercept = 0), linetype = "dashed", color = "firebrick4", alpha = .5)+
  #   scale_fill_manual(values = new_colors) +
  #   labs(x = "Species Interaction Value", y = NULL, fill = NULL, color = NULL, title = "All_large_carnivores all models liberal")+
  #   theme(axis.text.x = element_text(size = 18),
  #         axis.text.y = element_text(size = 18),
  #         axis.text = element_text(color = "black"),
  #         text = element_text(family = "Helvetica"))
  # ## save it! 
  # ggsave(paste("figures/Grouped Histograms/custom_size_histo/Histogram_all_large_carnivores_LIBERAL_multi_unsupported_and_supported_species_interaction_all_mods_HALF_LONG_",date, ".png", sep = ""),
  #        p, width = 10, height = 4, units = "in")
  # 
  # for(i in 1:length(unique(plot_dat$rel_species))){
  #   
  #   ## select one species 
  #   r =unique(plot_dat$rel_species)[i]
  #   if(r == "All_large_carnivores"){next} # skip all large carns 
  #   
  #   # make the plot 
  #   p = ggplot(plot_dat[plot_dat$rel_species == r,], aes(x = Interaction_Estimate, fill = support))+
  #     geom_histogram(aes(y=after_stat(count)), colour="black", binwidth = .1, alpha = 1)+#, position = "dodge")+
  #     theme_classic()+
  #     geom_vline(aes(xintercept = 0), linetype = "dashed", color = "firebrick4", alpha = .5)+
  #     scale_fill_manual(values = new_colors) +
  #     labs(x = "Species Interaction Value", y = NULL, color = NULL, title = paste(r, "all models"))+
  #     guides(fill = "none") + 
  #     theme(axis.text.x = element_text(size = 16),
  #           axis.text.y = element_text(size = 16),
  #           axis.text = element_text(color = "black"),
  #           text = element_text(family = "Helvetica"))
  # 
  #   # and save it!
  #   # ggsave(paste("figures/Grouped Histograms/custom_size_histo/Histogram_", r, "_multi_unsupported_and_supported_species_interaction_all_mods_HALF_LONG_",date, ".png", sep = ""),
  #   #               p, width = 5, height = 2, units = "in")
  #   
  #   ## Repeat for the libs! 
  #   p = ggplot(plot_dat[plot_dat$rel_species == r,], aes(x = Interaction_Estimate, fill = support_liberal))+
  #     geom_histogram(aes(y=after_stat(count)), colour="black", binwidth = .1, alpha = 1)+#, position = "dodge")+
  #     theme_classic()+
  #     geom_vline(aes(xintercept = 0), linetype = "dashed", color = "firebrick4", alpha = .5)+
  #     scale_fill_manual(values = new_colors) +
  #     labs(x = "Species Interaction Value", y = NULL, color = NULL, title = paste(r, "all models liberal"))+
  #     guides(fill = "none") + 
  #     theme(axis.text.x = element_text(size = 16),
  #           axis.text.y = element_text(size = 16),
  #           axis.text = element_text(color = "black"),
  #           text = element_text(family = "Helvetica"))
  #   
  #   # and save it!
  #   # ggsave(paste("figures/Grouped Histograms/custom_size_histo/Histogram_", r, "_LIBERAL_multi_unsupported_and_supported_species_interaction_all_mods_HALF_LONG_",date, ".png", sep = ""),
  #   #        p, width = 5, height = 2, units = "in")
  #   
  # }
  # # clean it up
  # rm(i,r,p)
  # 
}


### determine percentages in each category
perc_table = ddply(plot_dat, .(rel_species, support_simple), summarize,
                   num_lvl = length(Species_Pair))
# then total number of mods per species
perc_table$num_total[perc_table$rel_species == "Neofelis_genus"] = 
  sum(perc_table$num_lvl[perc_table$rel_species == "Neofelis_genus"])
perc_table$num_total[perc_table$rel_species == "Cuon_alpinus"] = 
  sum(perc_table$num_lvl[perc_table$rel_species == "Cuon_alpinus"])
perc_table$num_total[perc_table$rel_species == "Panthera_pardus"] = 
  sum(perc_table$num_lvl[perc_table$rel_species == "Panthera_pardus"])
perc_table$num_total[perc_table$rel_species == "Panthera_tigris"] = 
  sum(perc_table$num_lvl[perc_table$rel_species == "Panthera_tigris"])
perc_table$num_total[perc_table$rel_species == "All_large_carnivores"] = 
  sum(perc_table$num_lvl[perc_table$rel_species == "All_large_carnivores"])
# finally calulate the percent-
perc_table$percent = perc_table$num_lvl / perc_table$num_total * 100
## add the PPC settings
perc_table$PPC_setting = "conservative"


### repeat for the liberal grouping
perc_table_lib = ddply(plot_dat, .(rel_species, support_simple_liberal), summarize,
                   num_lvl = length(Species_Pair))
# then total number of mods per species
perc_table_lib$num_total[perc_table_lib$rel_species == "Neofelis_genus"] = 
  sum(perc_table_lib$num_lvl[perc_table_lib$rel_species == "Neofelis_genus"])
perc_table_lib$num_total[perc_table_lib$rel_species == "Cuon_alpinus"] = 
  sum(perc_table_lib$num_lvl[perc_table_lib$rel_species == "Cuon_alpinus"])
perc_table_lib$num_total[perc_table_lib$rel_species == "Panthera_pardus"] = 
  sum(perc_table_lib$num_lvl[perc_table_lib$rel_species == "Panthera_pardus"])
perc_table_lib$num_total[perc_table_lib$rel_species == "Panthera_tigris"] = 
  sum(perc_table_lib$num_lvl[perc_table_lib$rel_species == "Panthera_tigris"])
perc_table_lib$num_total[perc_table_lib$rel_species == "All_large_carnivores"] = 
  sum(perc_table_lib$num_lvl[perc_table_lib$rel_species == "All_large_carnivores"])
# finally calulate the percent-
perc_table_lib$percent = perc_table_lib$num_lvl / perc_table_lib$num_total * 100
## add the PPC settings
perc_table_lib$PPC_setting = "liberal"

## match col names for easier binding
names(perc_table)[2] = "support_levels"
names(perc_table_lib)[2] = "support_levels"

## combine
perc_table = rbind(perc_table, perc_table_lib)

# save it for fig making 
# write.csv(perc_table, paste("results/summarized_results/percentages_for_3-level_support_histograms_all_models_", date, ".csv", sep = ""), row.names = F)
rm(perc_table, perc_table_lib)

## make the graph, omititng the very negative values 
p = 
  ggplot(plot_dat[plot_dat$rel_species == "All_large_carnivores" &
                    plot_dat$Interaction_Estimate > -2,], aes(x = Interaction_Estimate, fill = support_simple_liberal))+
  geom_histogram(aes(y=after_stat(count)), colour="black", binwidth = .1, alpha = 1)+#, position = "dodge")+
  theme_classic()+
  geom_vline(aes(xintercept = 0), linetype = "dashed", color = "firebrick4", alpha = .5)+
  scale_fill_manual(values = new_colors) +
  labs(x = "Species Interaction Value", y = "Count", color = NULL, title = paste("all_large_carnivores all prey liberal"))+
  guides(fill = "none") + 
  theme(axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        axis.text = element_text(color = "black"),
        text = element_text(family = "Helvetica"))
## save it! 
# ggsave(paste("figures/Grouped Histograms/Histogram_all_large_carnivores_LIBERAL_unsupported_and_supported_all_mods_", date, ".png", sep = ""), p,
#        width = 12, height = 4, units = "in")


#
##
###
####
##### Repeat for preferred prey species too! 
td_dat = plot_dat[plot_dat$direction == "top-down",]

## subset for each pref prey and then combine all
a = td_dat[td_dat$rel_species == "Neofelis_genus" &
                   td_dat$sub_species %in% guilds$scientificNameStd[guilds$CL_pref == "Yes"],]
b = td_dat[td_dat$rel_species == "Cuon_alpinus" &
                   td_dat$sub_species %in% guilds$scientificNameStd[guilds$dhole_pref == "Yes"],]
c = td_dat[td_dat$rel_species == "Panthera_pardus" &
                   td_dat$sub_species %in% guilds$scientificNameStd[guilds$leopard_pref == "Yes"],]
d = td_dat[td_dat$rel_species == "Panthera_tigris" &
                   td_dat$sub_species %in% guilds$scientificNameStd[guilds$tiger_pref == "Yes"],]
e = rbind(a,b,c,d)

## double down for all_carns 
f = e
f$rel_species = "All_large_carnivores"
td_dat_pref = rbind(f,e)
rm(a,b,c,d,e,f)

## inspect
table(td_dat_pref$support);table(td_dat_pref$support_liberal) # all good! 


#
##
### Repeat the process for bottom-up 
bu_dat = plot_dat[plot_dat$direction == "bottom-up",]

## subset for each pref prey and then combine all
a = bu_dat[bu_dat$rel_species == "Neofelis_genus" &
                   bu_dat$dom_species %in% guilds$scientificNameStd[guilds$CL_pref == "Yes"],]
b = bu_dat[bu_dat$rel_species == "Cuon_alpinus" &
                   bu_dat$dom_species %in% guilds$scientificNameStd[guilds$dhole_pref == "Yes"],]
c = bu_dat[bu_dat$rel_species == "Panthera_pardus" &
                   bu_dat$dom_species %in% guilds$scientificNameStd[guilds$leopard_pref == "Yes"],]
d = bu_dat[bu_dat$rel_species == "Panthera_tigris" &
                   bu_dat$dom_species %in% guilds$scientificNameStd[guilds$tiger_pref == "Yes"],]
e = rbind(a,b,c,d)

## double down for all_carns 
f = e
f$rel_species = "All_large_carnivores"
bu_dat_pref = rbind(f,e)
rm(a,b,c,d,e,f)

## inspect
table(bu_dat_pref$support);table(bu_dat_pref$support_liberal) # all good! 

## and combine them
pref_dat = rbind(td_dat_pref, bu_dat_pref)
rm(td_dat_pref, bu_dat_pref)

## inspect
summary(pref_dat$Interaction_Estimate) # small range
anyNA(pref_dat) # should be F
table(pref_dat$rel_species) #looks good! 

## set factor levels for order on graph
pref_dat$support_simple = factor(pref_dat$support_simple, 
                           levels = c("top-down", "bottom-up","unsupported"))
pref_dat$support_simple_liberal = factor(pref_dat$support_simple_liberal, 
                                         levels = c("top-down", "bottom-up","unsupported"))

table(pref_dat$support_simple);table(pref_dat$support_simple_liberal)
anyNA(pref_dat$support_simple);anyNA(pref_dat$support_simple_liberal) # must be F! 

## Frequency plots 
plot_list = list()
plot_list_lib = list()

order = c("All_large_carnivores","Panthera_tigris","Panthera_pardus", 
          "Cuon_alpinus", "Neofelis_genus")

# # grab min and max vals 
min_val = min(pref_dat$Interaction_Estimate)
max_val = max(pref_dat$Interaction_Estimate)
for(i in 1:length(order)){
  
  # subset the data
  d = pref_dat[pref_dat$rel_species == order[i], ]
  
  # make the plot 
  p = ggplot(d, aes(x = Interaction_Estimate, fill = support_simple))+
    geom_histogram(aes(y=after_stat(count)), colour="black", binwidth = .1, alpha = 1)+#, position = "dodge")+
    # theme_minimal()+
    # theme_test()+
    theme_classic()+
    geom_vline(aes(xintercept = 0), linetype = "dashed", color = "firebrick4", alpha = .5)+
    scale_fill_manual(values = new_colors) +
    coord_cartesian(xlim = c(min_val, max_val))+
    labs(x = "Species Interaction Value", y = NULL, fill = NULL, color = NULL, title = order[i])
  
  # save it! 
  plot_list[[i]] = p
  
  ### Repeat for liberal values
  
  # make the plot 
  p = ggplot(d, aes(x = Interaction_Estimate, fill = support_simple_liberal))+
    geom_histogram(aes(y=after_stat(count)), colour="black", binwidth = .1, alpha = 1)+#, position = "dodge")+
    # theme_minimal()+
    # theme_test()+
    theme_classic()+
    geom_vline(aes(xintercept = 0), linetype = "dashed", color = "firebrick4", alpha = .5)+
    scale_fill_manual(values = new_colors) +
    coord_cartesian(xlim = c(min_val, max_val))+
    labs(x = "Species Interaction Value", y = NULL, fill = NULL, color = NULL, title = order[i])+
    theme(axis.text.x = element_text(size = 16),
          axis.text.y = element_text(size = 16),
          axis.text = element_text(color = "black"),
          text = element_text(family = "Helvetica"))
  
  # save it! 
  plot_list_lib[[i]] = p
}
rm(d,p,i, min_val, max_val)

## arrange them vertically 
p = plot_grid(plotlist = plot_list, align = "v",
              ncol = 1, nrow = length(plot_list)) #
# and save it! 
ggsave(paste("figures/Grouped Histograms/Histogram_SUPP_&_NONSUPP_&_CONSERVATIVE_SIVs_preferred_prey_HALF_LONG_",date, ".png", sep = ""),
               p, width = 8, height = 10, units = "in")

## repeat for liberal grouping
p = plot_grid(plotlist = plot_list_lib, align = "v", 
              ncol = 1, nrow = length(plot_list)) # 
# and save it! 
# ggsave(paste("figures/Grouped Histograms/Histogram_SUPP_&_NONSUPP_&_LIBERAL_SIVs_preferred_prey_HALF_LONG_",date, ".png", sep = ""),
#                p, width = 8, height = 10, units = "in")


## custom-make all large carnivore graph for sizing
p = 
  ggplot(pref_dat[pref_dat$rel_species == "All_large_carnivores" ,], aes(x = Interaction_Estimate, fill = support_simple_liberal))+
  geom_histogram(aes(y=after_stat(count)), colour="black", binwidth = .1, alpha = 1)+#, position = "dodge")+
  theme_classic()+
  geom_vline(aes(xintercept = 0), linetype = "dashed", color = "firebrick4", alpha = .5)+
  scale_fill_manual(values = new_colors) +
  labs(x = "Species Interaction Value", y = "Count", color = NULL, title = paste("all_large_carnivores preferred prey liberal"))+
  guides(fill = "none") + 
  theme(axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        axis.text = element_text(color = "black"),
        text = element_text(family = "Helvetica"))
## save it! 
# ggsave(paste("figures/Grouped Histograms/Histogram_all_large_carnivores_LIBERAL_unsupported_and_supported_preferred_prey_", date, ".png", sep = ""), p,
#        width = 12, height = 4, units = "in")


### determine percentages in each category
pref_perc_table = ddply(pref_dat, .(rel_species, support_simple), summarize,
                   num_lvl = length(Species_Pair))
# then total number of mods per species
pref_perc_table$num_total[pref_perc_table$rel_species == "Neofelis_genus"] = 
  sum(pref_perc_table$num_lvl[pref_perc_table$rel_species == "Neofelis_genus"])
pref_perc_table$num_total[pref_perc_table$rel_species == "Cuon_alpinus"] = 
  sum(pref_perc_table$num_lvl[pref_perc_table$rel_species == "Cuon_alpinus"])
pref_perc_table$num_total[pref_perc_table$rel_species == "Panthera_pardus"] = 
  sum(pref_perc_table$num_lvl[pref_perc_table$rel_species == "Panthera_pardus"])
pref_perc_table$num_total[pref_perc_table$rel_species == "Panthera_tigris"] = 
  sum(pref_perc_table$num_lvl[pref_perc_table$rel_species == "Panthera_tigris"])
pref_perc_table$num_total[pref_perc_table$rel_species == "All_large_carnivores"] = 
  sum(pref_perc_table$num_lvl[pref_perc_table$rel_species == "All_large_carnivores"])
# finally calulate the percent-
pref_perc_table$percent = pref_perc_table$num_lvl / pref_perc_table$num_total * 100
# add the PPC setting
pref_perc_table$PPC_setting = "conservative"

### Repeat for liberal groupings
pref_perc_table_lib = ddply(pref_dat, .(rel_species, support_simple_liberal), summarize,
                   num_lvl = length(Species_Pair))
# then total number of mods per species
pref_perc_table_lib$num_total[pref_perc_table_lib$rel_species == "Neofelis_genus"] = 
  sum(pref_perc_table_lib$num_lvl[pref_perc_table_lib$rel_species == "Neofelis_genus"])
pref_perc_table_lib$num_total[pref_perc_table_lib$rel_species == "Cuon_alpinus"] = 
  sum(pref_perc_table_lib$num_lvl[pref_perc_table_lib$rel_species == "Cuon_alpinus"])
pref_perc_table_lib$num_total[pref_perc_table_lib$rel_species == "Panthera_pardus"] = 
  sum(pref_perc_table_lib$num_lvl[pref_perc_table_lib$rel_species == "Panthera_pardus"])
pref_perc_table_lib$num_total[pref_perc_table_lib$rel_species == "Panthera_tigris"] = 
  sum(pref_perc_table_lib$num_lvl[pref_perc_table_lib$rel_species == "Panthera_tigris"])
pref_perc_table_lib$num_total[pref_perc_table_lib$rel_species == "All_large_carnivores"] = 
  sum(pref_perc_table_lib$num_lvl[pref_perc_table_lib$rel_species == "All_large_carnivores"])
# finally calulate the percent-
pref_perc_table_lib$percent = pref_perc_table_lib$num_lvl / pref_perc_table_lib$num_total * 100
## add the PPC setting
pref_perc_table_lib$PPC_setting = "liberal"

## re-name cols to save
names(pref_perc_table_lib)[2] = "support_levels"
names(pref_perc_table)[2] = "support_levels"

## rbind them 
pref_perc_table = rbind(pref_perc_table, pref_perc_table_lib)

# save it for fig making 
# write.csv(pref_perc_table, paste("results/summarized_results/percentages_for_3-level_support_histogram_preferred_prey_", date, ".csv", sep = ""), row.names = F)


## Clean up enviro after finishing plots! 
rm(p, pref_perc_table_lib,pref_perc_table, perc_table,
   plot_dat, bu_dat, td_dat,pref_dat,new_colors, order,
   date, plot_list, plot_list_lib)


############# Visualize coefficient effect sizes ###### 

#### ONLY RUN THRU THIS ONCE ALL MODELS ARE COMPLETED! 


### Will be making a graph comparing the effect sizes for each state and detection variable
## Repeated for each species pair

### Will need to make directories in the figures folder for relevant guild pairs
sort(unique(preform$guild_pair)) # 10 guilds is much more reasonable than 610 species pairs. 

# create each directory
for(i in 1:length(unique(preform$guild_pair))){
  
  #construct the path
  path = paste("figures/", unique(preform$guild_pair)[i], sep = "")
  
  # and make the directory
  dir.create(path)
  
} # end per guild pair 
rm(i, path)
## will get warnings if the directory is already created, no dramas tho. 

## create a directory for each species pair in the relevant guild pair b/c each species pair will have multiple plots 
for(i in 1:length(unique(preform$Species_Pair))){
  
  #construct the path
  path = paste("figures/", preform$guild_pair[preform$Species_Pair == unique(preform$Species_Pair)[i]],
               "/", unique(preform$Species_Pair)[i], sep = "")
  
  # and make the directory
  dir.create(path)
}
rm(i, path)

#
##
###
#### Loop through all species pairs to make a plot for each!

sp_pair_coeff_plots = list() # save em here. 

for(i in 1:length(unique(coeff$Species_Pair))){
  
  ## subset for a single species 
  dat = coeff[coeff$Species_Pair == unique(coeff$Species_Pair)[i],]
  
  ## subset for relevant variables 
  # dat = dat[dat$var %in% c("Abundance_intercept","FLII","HFP","Elevation","Hunting",
  #                          "Species_Interaction","Detection_intercept","Active_cams"),]
  dat = dat[dat$var %in% c("FLII","HFP","Elevation","Hunting",
                           "Species_Interaction","Active_cams"),] # removing intercepts b/c they skew graphs
  
  ## replace _ with space for vars for clean x-vars
  dat$var = gsub("_", " ", dat$var)
  
  ## Change active cams to effort for easier understanding
  dat$var[dat$var == "Active cams"] = "Effort"
  
  ## order the factor levels for each variable
  dat$var = factor(dat$var, levels = c("Species Interaction","FLII","HFP", "Hunting", # state vars
                                       "Elevation", # more state vars
                                       "Effort")) # then det vars. 
  
  # Add a column for significance marker
  dat$marker <- ifelse(dat$sig == "Significant", "*", "")
  
  # Find the position of "Abundance intercept" on the x-axis to add the dashed line. 
  # split_position <- which(levels(dat$var) == "Abundance intercept")
  split_position <- which(levels(dat$var) == "Elevation")
  
  
  # specify the distance were spacing error bars and stars
  dodge_width = 0.8
  
  # Extract the prefix from species to determine sub or dom
  dat$prefix <- sub("^(SUB|DOM).*", "\\1", dat$species)
  
  # specify the colors for dominant and subordinate species based on prefix 
  dat$color[dat$prefix == "SUB"] = "wheat3"
  dat$color[dat$prefix == "DOM"] = "tan3"
  
  # Clean up species name for clarity 
  dat$species_name = sub("^(SUB|DOM)-", "", dat$species)
  
  # Reorder factor levels of species_name
  dat$species_name <- factor(dat$species_name, 
                             levels = c(unique(dat$species_name[dat$prefix == "DOM"]), 
                                        unique(dat$species_name[dat$prefix == "SUB"])))
  
  # make a title 
  title = paste("Dominant species:", unique(dat$species_name[dat$prefix == "DOM"]), 
                "affects subordinate species:", unique(dat$species_name[dat$prefix == "SUB"]))
  
  
  # Plot the grouped bar chart
  p =
    ggplot(dat, aes(x = var, y = mean, fill = species_name)) +
    geom_bar(stat = "identity", position = position_dodge(width = dodge_width), color = "black", width = dodge_width) +
    geom_errorbar(aes(ymin = lower, ymax = upper), position = position_dodge(width = dodge_width), width = 0, color = "black") +
    geom_text(aes(y = mean+.3*sign(mean), label = marker), position = position_dodge(width = dodge_width), size = 8) +
    geom_vline(xintercept = split_position + 0.5, color = "red", linetype = "dashed") +
    geom_hline(yintercept = 0)+
    scale_fill_manual(values = dat$color) +
    labs(x = NULL, y = "Mean Effect Size", fill = NULL, title = title) +
    theme_test()+
    theme(
      axis.text.x = element_text(size = 12, angle = 45, hjust = 1),   # Increase x-axis label font size and tilt it
      axis.text.y = element_text(size = 12)                           # Increase y-axis tick label font size
    )
  ## looks good! 
  
  ## save it 
  sp_pair_coeff_plots[[i]] = p
  names(sp_pair_coeff_plots)[i] = unique(coeff$Species_Pair)[i]
  
}
# keep it clean! 
rm(dodge_width, i, split_position, title, p, dat) 

## Save em to the relevant file directories (based off guilds)
for(i in 1:length(sp_pair_coeff_plots)){
  
  ## grab one plot 
  p = sp_pair_coeff_plots[[i]]
  #and the name
  n = names(sp_pair_coeff_plots)[i]
  
  # find the matching guild pair
  gp = preform$guild_pair[preform$Species_Pair == n]
  
  # grab today's date
  day<-str_sub(Sys.Date(),-2)
  month<-str_sub(Sys.Date(),-5,-4)
  year<-str_sub(Sys.Date(),-10,-7)
  date = paste(year,month,day, sep = "")
  
  # construct a saving path 
  path = paste("figures/", paste(gp, n, sep = "/"), "/model_coefficents_", n, "_", date, ".png", sep = "")
  
  ## save it! 
  ggsave(path, p, width = 10.5, height = 6.5, units = "in")
  
} 
rm(p,n,gp,i,path,day,month,year,date)

# dont need these plots anymore as they are saved
rm(sp_pair_coeff_plots)



############# Visualize prediction plots ########

## Goal is to visualze the co-abundance relationship between every pairwise interaction
## With so much data now, we will not use est_abund b/c it makes the graphs too messy. 

## Inspect data
head(est_abund) #estimated abundance per sampling unit --> no longer using 
str(est_abund) # already numeric

head(pred_abund) #Predicted change in subordinate abundance as a function of dominant species abundance --> this is what we want! 
str(pred_abund) # already numeric

## set colors 
three_colors = c("unsupported" = "gray32",
                 "top-down" = "darkorange4",
                 "bottom-up" = "springgreen4")


## repeat for each species pair 
for(i in 1:length(unique(pred_abund$Species_Pair))){
  
  # save one pair name
  n = unique(pred_abund$Species_Pair)[i]
  
  ## subset est and pred
  est = est_abund[est_abund$Species_Pair == n,]
  pred = pred_abund[pred_abund$Species_Pair == n,]
  
  ## grab the relevant guild pair for saving later
  gp = preform$guild_pair[preform$Species_Pair == n]
  
  ## and the date while were at it
  day<-str_sub(Sys.Date(),-2)
  month<-str_sub(Sys.Date(),-5,-4)
  year<-str_sub(Sys.Date(),-10,-7)
  date = paste(year,month,day, sep = "")
  
  ## split species apart for nicer names 
  n_split = str_split(n, "~")[[1]]
  n_split = gsub("SUB-", "", n_split)
  n_split = gsub("DOM-", "", n_split)
  
  
  ## make a nice title 
  title = paste("Dominant species:", n_split[2], 
                "affects subordinate species:", n_split[1])
  
  ## make plots according to support level color 
  if(preform$support_liberal[preform$Species_Pair == n] == "Supported"){
    
    ## make the top-down plot 
    if(preform$direction[preform$Species_Pair == n] == "top-down"){
      
      p = 
        ggplot(pred, aes(y = Predicted.N, x = var))+
        geom_line(color = "darkorange4", linewidth = 2)+
        geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.1, color = "black", linetype = "dashed")+
        labs(y = paste0(n_split[1], " Abundance"), 
             x = paste0(n_split[2], " Abundance"),
             title = title, color = NULL)+ 
        theme_test()+
        theme(axis.text.x = element_text(size = 22),
              axis.text.y = element_text(size = 22),
              axis.text = element_text(color = "black"),
              text = element_text(family = "Helvetica"))
      
    } # end top-down
    
    # make the bottom-up plot
    if(preform$direction[preform$Species_Pair == n] == "bottom-up"){
      
      
      p = 
        ggplot(pred, aes(y = Predicted.N, x = var))+
        geom_line(color = "springgreen4", linewidth = 2)+
        geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.1, color = "black", linetype = "dashed")+
        labs(y = paste0(n_split[1], " Abundance"), 
             x = paste0(n_split[2], " Abundance"),
             title = title, color = NULL)+ 
        theme_test()+
        theme(axis.text.x = element_text(size = 22),
              axis.text.y = element_text(size = 22),
              axis.text = element_text(color = "black"),
              text = element_text(family = "Helvetica"))
      
    } # end bottom-up
    
  }else{
    
    ## make the unsupported plot 
    p = 
      ggplot(pred, aes(y = Predicted.N, x = var))+
      geom_line(color = "black", linewidth = 2, linetype = "dashed")+
      geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.1, color = "black")+
      labs(y = paste0(n_split[1], " Abundance"), 
           x = paste0(n_split[2], " Abundance"),
           title = title, color = NULL)+ 
      guides(color = "none") + # way too many landscapes to list them now. 
      theme_test()+
      theme(axis.text.x = element_text(size = 22),
            axis.text.y = element_text(size = 22),
            axis.text = element_text(color = "black"),
            text = element_text(family = "Helvetica"))
    
  }
  
  ## construct the saving path
  path = paste("figures/", paste(gp, n, sep = "/"), "/Prediction_plot_", n, "_", date, ".png", sep = "")
  
  ## save the dang thang 
  ggsave(path, p, width = 5, height = 4, units = "in")
  
}
rm(path, p,i, title, n, n_split, gp, est, pred, day, month, year, date, 
   three_colors)





############# Visualize PPC plots ########

## Will generate two plots, 2 w/ the scatter (one for dom and one for sub) and 1 w/ the bars 

## scatter first 
# add guild_pair to the data
add = select(preform, Species_Pair, guild_pair)
ppc_plotdat = merge(ppc_plotdat, add, by = "Species_Pair")

## loop thru each species pair to make the scatter plots per species and save them directly. 
for(i in 1:length(unique(ppc_plotdat$Species_Pair))){
  
  ## grab species pair name 
  n = unique(ppc_plotdat$Species_Pair)[i]
  
  ## subset for a single species 
  dat = ppc_plotdat[ppc_plotdat$Species_Pair == n,]
  
  # grab today's date
  day<-str_sub(Sys.Date(),-2)
  month<-str_sub(Sys.Date(),-5,-4)
  year<-str_sub(Sys.Date(),-10,-7)
  date = paste(year,month,day, sep = "")
  
  ## Create titles based on values
  val = ppc_values[ppc_values$Species_Pair == n,]
  sub.bpv = paste("Bayesian P-value =", round(val$BPV.sub, 3),
                  "& C-hat =", round(val$Chat.sub, 3))
  dom.bpv = paste("Bayesian P-value =", round(val$BPV.dom, 3),
                  "& C-hat =",round(val$Chat.dom, 3))
  
  ## Subordinate plot
  plot=
    ggplot(dat, aes( y = sub.sim, x = sub.real))+
    geom_point(aes(color = sub.label))+
    geom_abline(intercept = 0, slope = 1, color = "black", size = 2)+
    theme_classic()+
    scale_color_manual(values = c("red", "blue"))+
    labs(title = sub.bpv, subtitle = str_split(n, "~")[[1]][1]  , 
         y = "Fit Statistic Simulated Data", x = "Fit Statistic Real Data")+
    theme(legend.position = "None")
  
  ## make a saving path
  path = paste("figures/", unique(dat$guild_pair), "/", n, "/Subordinate_PPC_scatter_plot_", date, ".png", sep = "")
  
  # save it! 
  ggsave(path, plot, width = 7.5, height = 7, units = "in")
  
  
  ## Dominant plot
  plot = 
    ggplot(dat, aes( y = dom.sim, x = dom.real))+
    geom_point(aes(color = dom.label))+
    geom_abline(intercept = 0, slope = 1, color = "black", size = 2)+
    theme_classic()+
    scale_color_manual(values = c("red", "blue"))+
    labs(title = dom.bpv, subtitle = str_split(n, "~")[[1]][2], 
         y = "Fit Statistic Simulated Data", x = "Fit Statistic Real Data")+
    theme(legend.position = "None")
  
  ## make a saving path
  path = paste("figures/", unique(dat$guild_pair), "/", n, "/Dominant_PPC_scatter_plot_", date, ".png", sep = "")
  
  # save it! 
  ggsave(path, plot, width = 7.5, height = 7, units = "in")
  
  
}
rm(path, day,month,year,date, plot, n, sub.bpv, dom.bpv, val, dat)

# add guild_pair to the values
add = select(preform, Species_Pair, guild_pair)
ppc_values = merge(ppc_values, add, by = "Species_Pair")

## Loop thru each species pair in the values DF to make bar graphs for BPV and Chat values
for(i in 1:length(unique(ppc_values$Species_Pair))){
  
  ## grab species pair name 
  n = unique(ppc_values$Species_Pair)[i]
  
  ## subset for a single species 
  dat = ppc_values[ppc_values$Species_Pair == n,]
  
  ## slightly transform the data to plot better 
  plot_dat = data.frame("BPV" = c(dat$BPV.dom, dat$BPV.sub),
                        "Chat" = c(dat$Chat.dom, dat$Chat.sub),
                        "position" = c("Dominant", "Subordinate"),
                        "species" = str_split(n, "~")[[1]])
  
  # Define labels and positions for text
  text_data <- data.frame(
    label = c("liberal","liberal", "conservative", "conservative", "perfect fit"),
    y = c(0.85, 0.15, 0.75, 0.25, 0.5)  # Adjust these values based on your needs
  )
  
  ## BPV plot
  bpv =
  ggplot(plot_dat, aes(x = species, y = BPV, fill = position))+
    geom_col(color = "black", position = "dodge", color = 'black')+
    geom_hline(yintercept = 0.5)+
    geom_hline(yintercept = 0.75, linetype = "dashed", color = "gray32", alpha = 0.5)+
    geom_hline(yintercept = 0.25, linetype = "dashed", color = "gray32", alpha = 0.5)+
    geom_hline(yintercept = 0.85, linetype = "dashed", color = "red")+
    geom_hline(yintercept = 0.15, linetype = "dashed", color = "red")+
    geom_text(data = text_data, aes(x = Inf, y = y, label = label), hjust = 1.2, vjust = -1, size = 5, inherit.aes = F)+
    theme_test()+
    labs(y = "Bayesian P-Value", x = NULL, fill = NULL)+
    coord_cartesian(ylim = c(0,1))+
    scale_fill_manual(values = c("tan3","wheat3"))+
    theme(axis.text.x = element_text(size = 12, face = "bold"))
  
  # grab today's date
  day<-str_sub(Sys.Date(),-2)
  month<-str_sub(Sys.Date(),-5,-4)
  year<-str_sub(Sys.Date(),-10,-7)
  date = paste(year,month,day, sep = "")
  
  ## make a saving path
  path = paste("figures/", paste(dat$guild_pair, dat$Species_Pair, sep = "/"), "/", "Bayes_P-value_bar_plot_", date, ".png", sep = "")
  
  # save it! 
  ggsave(path, bpv, width = 7.5, height = 7, units = "in")
  
  # Define labels and positions for text
  text_data <- data.frame(
    label = c("liberal", "conservative", "perfect fit"),
    y = c(1.3, 1.1, 1)  # Adjust these values based on your needs
  )
  
  ##Chat plot
  chat = 
    ggplot(plot_dat, aes(x = species, y = Chat, fill = position))+
    geom_col(color = "black", position = "dodge", color = 'black')+
    geom_hline(yintercept = 1)+
    geom_hline(yintercept = 1.1, linetype = "dashed", color = "gray32", alpha = 0.5)+
    geom_hline(yintercept = 1.3, linetype = "dashed", color = "red")+
    geom_text(data = text_data, aes(x = Inf, y = y, label = label), hjust = 1.2, vjust = -1, size = 5, inherit.aes = F)+
    theme_test()+
    labs(y = "Magnitude of Over-Dispersion (C-hat)", x = NULL, fill = NULL)+ 
    coord_cartesian(ylim = c(0.9,1.35))+
    scale_fill_manual(values = c("tan3","wheat3"))+
    theme(axis.text.x = element_text(size = 12, face = "bold"))
  
  # grab today's date
  day<-str_sub(Sys.Date(),-2)
  month<-str_sub(Sys.Date(),-5,-4)
  year<-str_sub(Sys.Date(),-10,-7)
  date = paste(year,month,day, sep = "")
  
  ## make a saving path
  path = paste("figures/", paste(dat$guild_pair, dat$Species_Pair, sep = "/"), "/", "C-hat_overdispersion_bar_plot_", date, ".png", sep = "")
  
  # save it! 
  ggsave(path, chat, width = 7.5, height = 7, units = "in")
  
}
rm(bpv, chat, i, path, day,month,year,date,plot_dat, 
   n, dat,add, text_data)



######## Summary statistics to describe overall dataset ######

## how many detections did we get for our large carnivores?
og_captures$Species = gsub(" ", "_", og_captures$Species)
ddply(og_captures[og_captures$Species %in% c("Panthera_pardus", "Panthera_tigris", 
                                             "Cuon_alpinus", "Neofelis_nebulosa", 
                                             "Neofelis_diardi" ),], .(Species), summarize,
      num_cap = length(Individuals))

### effort per survey?
og_meta$camera_end.date = as.Date(og_meta$camera_end.date, format = "%Y-%m-%d")
og_meta$camera_start.date = as.Date(og_meta$camera_start.date, format = "%Y-%m-%d")

eff = ddply(og_meta, .(survey_id), summarize,
            total_eff = sum(effort),
            dur = as.numeric(difftime(max(camera_end.date), min(camera_start.date), units = "days")))

summary(eff$dur) # mean = 80.25
sd(eff$dur) # sd = 30.16

rm(eff)

## wait, why do these number differ?!
summary(og_meta$effort)
summary(og_resamp_meta$Cell_Effort)
summary(og_resamp_meta$Avg_effort) # this is probably the most faithful data to include based on what was analyzed

### how many detections of our study species 
og_captures$Species = gsub(" ", "_", og_captures$Species)
nrow(og_captures[og_captures$Species %in% preform$sub_species,]) # 140,334

## how many surveys before removal?
length(unique(og_meta$survey_id)) #361

## how many surveys after removal 
length(unique(og_resamp_meta$survey_id)) #284

## how many landscapes made the cut in the end?
length(unique(og_resamp_meta$Landscape)) #46

## how many cameras did we deploy at each landscape?
info = ddply(og_resamp_meta, .(Landscape), summarize,
             num_cam = length(unique(cell_id_3km)))
summary(info$num_cam)
mean(info$num_cam) #103
sd(info$num_cam) #141 --> lots of variation! 


############# Visualize conceptual study map to send to Jon ###### 

#### We can make a study map (conceptually) w/ the original metadata
### Where we can color our points based off which (or how many) large carnivores are present

## create a back up of og_captures for manipulation
caps = distinct(select(og_captures, deployment_id, survey_id, Landscape, Species))

## Make each column manually by hand
caps$tig_pres = "No"
caps$tig_pres[caps$Landscape %in% unique(caps$Landscape[caps$Species == "Panthera tigris"])] = "Yes"
caps$leo_pres = "No"
caps$leo_pres[caps$Landscape %in% unique(caps$Landscape[caps$Species == "Panthera pardus"])] = "Yes"
caps$CL_pres = "No"
caps$CL_pres[caps$Landscape %in% unique(caps$Landscape[grepl("Neofelis", caps$Species )])] = "Yes"
caps$dho_pres = "No"
caps$dho_pres[caps$Landscape %in% unique(caps$Landscape[caps$Species == "Cuon alpinus"])] = "Yes"

## reduce to minimal info 
caps = distinct(select(caps, Landscape, tig_pres, leo_pres, CL_pres, dho_pres))

# Create a new column to denote which species are present, thx Chat GPT
species_present <- apply(caps[, -1], 1, function(row) {
  species_names <- names(caps)[-1][row == "Yes"]
  if (length(species_names) > 0) {
    return(paste(species_names, collapse = ", "))
  } else {
    return("None")
  }
})

# Add the new column to the dataframe
caps$LC_pres = species_present
rm(species_present)

## clean up names to be nicer
caps$LC_pres = gsub("tig_pres", "Panthera tigris", caps$LC_pres)
caps$LC_pres = gsub("leo_pres", "Panthera pardus", caps$LC_pres)
caps$LC_pres = gsub("CL_pres", "Neofelis genus", caps$LC_pres)
caps$LC_pres = gsub("dho_pres", "Cuon alpinus", caps$LC_pres)

## create a new col to count how many LC present
caps$LC_count = str_count(caps$LC_pres, ",")
## add 1 to those greater than zero to get the trailing species name 
caps$LC_count[caps$LC_count > 0] = caps$LC_count[caps$LC_count > 0] + 1
## and for those with only one species
caps$LC_count[caps$LC_pres != "None" & caps$LC_count == 0] = 1

## now thin to relevant info 
LC_pres = distinct(select(caps, Landscape, LC_pres, LC_count))
rm(caps)

## calculate average coords per landsacpe to visualize 
avg_land = ddply(og_meta, .(Landscape), summarize,
                 avg_lat = mean(Latitude),
                 avg_long = mean(Longitude))

# verify all landscapes match
setdiff(avg_land$Landscape, LC_pres$Landscape)
setdiff(LC_pres$Landscape, avg_land$Landscape) # full match on both sides! 
# now merge w/ LC pres
avg_land = merge(avg_land, LC_pres, by = "Landscape")


## first make a map for the whole study area

#Genereate colour scheme as factor
category <- "LC_pres" #What should we use as our category to generate colours
colour <- "lightgreen"# Define a colour from the R options to base the colourscheme

## Manually set factor level so it displays nicely it a factor now 
avg_land$LC_pres <- factor(avg_land$LC_pres, levels = c("Panthera tigris, Panthera pardus, Neofelis genus, Cuon alpinus", 
                                                        "Panthera tigris, Panthera pardus, Cuon alpinus",
                                                        "Panthera tigris, Neofelis genus, Cuon alpinus",
                                                        "Panthera tigris, Neofelis genus",
                                                        "Panthera pardus, Neofelis genus",
                                                        "Neofelis genus, Cuon alpinus",
                                                        "Cuon alpinus", "Neofelis genus", "None"))
## Set the colors manually based off color brewer website 
col.cat = c("Panthera tigris, Panthera pardus, Neofelis genus, Cuon alpinus" = "#7f2704", 
            "Panthera tigris, Panthera pardus, Cuon alpinus" = "#a63603",
            "Panthera tigris, Neofelis genus, Cuon alpinus"= "#d94801",
            "Panthera tigris, Neofelis genus" = "#f16913",
            "Panthera pardus, Neofelis genus" = "#fd8d3c",
            "Neofelis genus, Cuon alpinus" = "#fdae6b",
            "Cuon alpinus" = "#fdd0a2", "Neofelis genus" = "#fee6ce", "None" = "#fff5eb")
avg_land$Cols <- col.cat[avg_land$LC_pres]

## exclude data from E_Indo
avg_land = avg_land[!grepl("E_Indo", avg_land$Landscape),]

## make the map! 
leaflet() %>%
  # Add base map data
  addProviderTiles(providers$Esri.WorldImagery, group="Satellite") %>%  
  addCircleMarkers(lng= avg_land$avg_long, #Add Longitude
                   lat= avg_land$avg_lat, #Add latitude
                   color= avg_land$Cols
                   ) %>%
  # Add a legend
  addLegend(position = "topright",
            colors = col.cat,
            labels = levels(avg_land$LC_pres),
            title = "Large carnivores present",
            opacity = 1,
            group = "Legend") %>%
  # Add a scale bar
  addScaleBar(
    position = c("bottomleft")) %>%
  # Layers control
  addLayersControl(
    overlayGroups = "Legend",
    options = layersControlOptions(collapsed = FALSE)
  )
## Screen-shot'ed! 
rm(col.cat, colour, category)

### Save avg land for jon to do some map making
# write.csv(avg_land, "results/summarized_results/average_landscapes_with_predators_present_20230702.csv", row.names = F)


### Now subset the data for just the Thai E forest complex
thai = og_meta[grepl("Thai_E_FC_", og_meta$Landscape),]

### But goal is to visualize the hexagons for resampling, so summarize data by cell_id
thai_cell = ddply(thai, .(cell_id_3km), summarize,
                  num_cam = length(unique(deployment_id)),
                  avg_lat = mean(Latitude),
                  avg_long = mean(Longitude))

# library(hexbin) # for hex bins
library(sf) # make grids

## Convert metadata to a spatial object using the coordinates
sf_data = st_as_sf(thai_cell, coords = c("avg_long", "avg_lat"), 
                 crs = "+proj=longlat +datum=WGS84")

# Ensure EPSG is 4087 via transformation
sf_data = st_transform(sf_data, 4087)

## 3km hex
hex3 = st_make_grid(sf_data, cellsize = 1861.2, square = FALSE) %>%
  st_sf() 

## This will suffice for concept now 
ggplot()+
  geom_sf(data = hex3)+
  geom_sf(data = sf_data)#, aes(color = num_cam))#[shape$deployment_id == p,])
## screen shott'ed

rm(hex3, sf_data, thai, thai_cell)

### now make a histogram of the captures for Thai E forest complex to show how we brought it all together
thai_caps = og_captures[grepl("Thai_E_FC", og_captures$Landscape),]
## format dates
thai_caps$Photo.Date = as.Date(thai_caps$Photo.Date, format = "%Y-%m-%d")
summary(thai_caps$Photo.Date) # all g

## set factor of landscape so it shows up nicely
thai_caps$Landscape = factor(thai_caps$Landscape, levels = c("Thailand_Thai_E_FC_West", 
                                                             "Thailand_Thai_E_FC_Center",
                                                             "Thailand_Thai_E_FC_East"))

## want to merge data source here too
add = distinct(select(og_meta[og_meta$deployment_id %in% thai_caps$deployment_id,], deployment_id, source))
thai_caps = merge(add, thai_caps, by = "deployment_id")

length(unique(thai_caps$survey_id)) ## 84 surveys is not feasible to graph! use source instead

# grab surv start and end dates tho 
surv_dates = ddply(thai_caps, .(survey_id), summarize,
                   start_date = min(Photo.Date), end_date = max(Photo.Date))

# make the plot 
p = 
ggplot(thai_caps, aes(x = Photo.Date))+
  geom_histogram(aes(y=after_stat(count), fill = source), colour="black", bins = 50)+
  geom_density(alpha=.2, fill="#FF6666")+
  theme_test()+
  facet_wrap(~Landscape)+
  labs(x = "Duration of Deployment", y="Number of Independent Detections", fill = "Data provider") 
# ggsave("results/summarized_results/Thai_E_FC_histogram_20230802.png", p, width = 12, height = 3, units = "in")
rm(thai_caps, p, surv_dates)



####### Do the same, but for BBS since there are less surveys but from multiple sources

## Screeshot from leaflet first
bbs = og_meta[grepl("Bukit_Barisan_Sel", og_meta$Landscape),]

## thin to relevant info
bbs_map = distinct(select(bbs, placename, source, Longitude, Latitude))

## add colors for source
bbs_map$cols = "blue"
bbs_map$cols[bbs_map$source == "Luskin"] = "orange"

## make the map! 
leaflet() %>%
  # Add base map data
  addProviderTiles(providers$Esri.WorldImagery, group="Satellite") %>%  
  addCircleMarkers(lng= bbs_map$Longitude, #Add Longitude
                   lat= bbs_map$Latitude, #Add latitude
                   color= (bbs_map$cols)
  ) %>%
  # Add a scale bar
  addScaleBar(
    position = c("bottomleft"))
## screenshotted! 

### Now subset the data for just the Thai E forest complex
bbs = og_meta[grepl("Bukit_Barisan_Sel", og_meta$Landscape),]

### But goal is to visualize the hexagons for resampling, so summarize data by cell_id
bbs_cell = ddply(bbs, .(cell_id_3km), summarize,
                  num_cam = length(unique(deployment_id)),
                  avg_lat = mean(Latitude),
                  avg_long = mean(Longitude))

# library(hexbin) # for hex bins
library(sf) # make grids

## Convert metadata to a spatial object using the coordinates
sf_data = st_as_sf(bbs_cell, coords = c("avg_long", "avg_lat"), 
                   crs = "+proj=longlat +datum=WGS84")

# Ensure EPSG is 4087 via transformation
sf_data = st_transform(sf_data, 4087)

## 3km hex
hex3 = st_make_grid(sf_data, cellsize = 1861.2, square = FALSE) %>%
  st_sf() 

## This will suffice for concept now 
ggplot()+
  geom_sf(data = hex3)+
  geom_sf(data = sf_data)#, aes(color = num_cam))#[shape$deployment_id == p,])
## screen shott'ed

rm(hex3, sf_data, bbs, bbs_cell)

### now make a histogram of the captures for Thai E forest complex to show how we brought it all together
bbs_caps = og_captures[grepl("Bukit_Barisan_Sel", og_captures$Landscape),]
## format dates
bbs_caps$Photo.Date = as.Date(bbs_caps$Photo.Date, format = "%Y-%m-%d")
summary(bbs_caps$Photo.Date) # all g

## want to merge data source here too
add = distinct(select(og_meta[og_meta$deployment_id %in% bbs_caps$deployment_id,], deployment_id, source))
bbs_caps = merge(add, bbs_caps, by = "deployment_id")

length(unique(bbs_caps$survey_id)) ## 15 surveys is more feasible to graph! 

## clean up caps tho
bbs_caps$survey_id = gsub('Bukit_Barisan_Selatan_NP', "BBSNP", bbs_caps$survey_id)

## make the graph 
p = 
  ggplot(bbs_caps, aes(x = Photo.Date))+
  geom_histogram(aes(y=after_stat(count), fill = survey_id), colour="black", bins = 50)+
  geom_density(alpha=.2, fill="#FF6666")+
  theme_test()+
  facet_wrap(~Landscape)+
  labs(x = "Duration of Deployment", y="Number of Independent Detections", fill = "survey_id") +
  theme(legend.position="bottom")
# ggsave("results/summarized_results/BBS_histogram_20230802.png", p, width = 10, height = 4, units = "in")

rm(p, category, col.cat, colour, LC_pres, avg_land, add, bbs_map, bbs_caps)



############# Visualize conceptual figures ######

### Essentially the goal here is to create new conceptual figures 
## that describe our hypothesis. 

### Will generate two histograms, and both will have two distributions:
# H1 --> negative distribution + norm distribution
# H2 --> positive distribution + norm distribution

## Use the sn pacakge to generate specified distributions
library(sn)

# Set a seed for reproducibility
set.seed(123)

# Central Normal distribution with mean=0
normal_data <- rnorm(100, mean = 0, sd = 6) # equal
# normal_data <- rnorm(150, mean = 0, sd = 6) # realistic


# Negative-skewed Normal distribution with mean= -2
# Using shape parameter to achieve negative skewness
negative_dat <- rsn(100, xi = -2, omega = 5, alpha = -3) # equal
# negative_dat <- rsn(50, xi = -2, omega = 5, alpha = -3)  # realistic

# Positive-skewed Normal distribution with mean=2
# Using shape parameter to achieve positive skewness
positive_dat <- rsn(100, xi = 2, omega = 5, alpha = 3) # equal
# positive_dat <- rsn(50, xi = 2, omega = 5, alpha = 3)  # realistic


# Combine data into one data frame
df <- data.frame(
  value = c(positive_dat, negative_dat, normal_data),
  distribution = factor(c(rep("positive_dat", length(positive_dat)),
                          rep("negative_dat", length(negative_dat)),
                          rep("Normal", length(normal_data)))),
  baseline = 0
)
rm(normal_data, positive_dat, negative_dat)

## make sure normal distribution goes behind the non-normal one!
df$distribution = factor(df$distribution, levels = c("positive_dat", "negative_dat","Normal"))

## and set the color scheme approprietly 
concept_col = c("Normal" = "gray32",
                "negative_dat" = "darkorange4",
                "positive_dat" = "springgreen4")

## make the x-axis closer to what the results are
df$value = df$value / 8 # smaller range

## make sure no negative values are positive!
df$value[df$distribution == "negative_dat" & df$value > 0] = (df$value[df$distribution == "negative_dat" & df$value > 0]) * -1


## Plotting the distributions
#H1
p =
  ggplot(df[df$distribution %in% c("negative_dat", "Normal"),], aes(x = value, fill = distribution))+
  geom_histogram(aes(y=after_stat(count)), colour="black")+#, binwidth = .1, alpha = 1)+#, position = "dodge")+
  theme_classic()+
  geom_vline(aes(xintercept = 0), linetype = "dashed", color = "firebrick4", alpha = .5)+
  scale_fill_manual(values = concept_col) +
  labs(x = "Species Interaction Value", y = NULL, color = NULL, title = "top-down hypothesis")+
  guides(fill = "none") + 
  coord_cartesian((xlim = c(-2,2)))+
  theme(axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        axis.text = element_text(color = "black"),
        text = element_text(family = "Helvetica"))
# ## save it! 
# ggsave("explore/H1_negative_concept_many_unsupported_histogram_20231116.png", p,
#        width = 12, height = 4, units = "in")

## remade plot, but with equal weighting between groups
# ggsave("explore/H1_negative_concept_equal_unsupported_histogram_20231116.png", p,
#        width = 6, height = 4, units = "in")

## make sure no positive values are negative!
df$value[df$distribution == "positive_dat" & df$value < 0] = (df$value[df$distribution == "positive_dat" & df$value < 0]) * -1

## repeat but for the positive ones
p = 
  ggplot(df[df$distribution %in% c("positive_dat", "Normal"),], aes(x = value, fill = distribution))+
  geom_histogram(aes(y=after_stat(count)), colour="black")+#, binwidth = .1, alpha = 1)+#, position = "dodge")+
  theme_classic()+
  geom_vline(aes(xintercept = 0), linetype = "dashed", color = "firebrick4", alpha = .5)+
  scale_fill_manual(values = concept_col) +
  labs(x = "Species Interaction Value", y = NULL, color = NULL, title = "bottom-up hypothesis")+
  guides(fill = "none") + 
  coord_cartesian((xlim = c(-2,2)))+
  theme(axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        axis.text = element_text(color = "black"),
        text = element_text(family = "Helvetica"))
# save it! 
# ggsave("explore/H1_positive_concept_many_unsupported_histogram_20231116.png", p,
#        width = 12, height = 4, units = "in")

## remade plot, but with equal weighting between groups
# ggsave("explore/H2_positive_concept_equal_unsupported_histogram_20231116.png", p,
#        width = 6, height = 4, units = "in")


## run a plot with both top down and bottom up 
p = 
ggplot(df, aes(x = value, fill = distribution))+
  geom_histogram(aes(y=after_stat(count)), colour="black")+#, binwidth = .1, alpha = 1)+#, position = "dodge")+
  theme_classic()+
  geom_vline(aes(xintercept = 0), linetype = "dashed", color = "firebrick4", alpha = .5)+
  scale_fill_manual(values = concept_col) +
  labs(x = "Species Interaction Value", y = "Count", color = NULL, title = "both hypothesis")+
  guides(fill = "none") + 
  coord_cartesian((xlim = c(-2,2)))+
  theme(axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        axis.text = element_text(color = "black"),
        text = element_text(family = "Helvetica"))
# ## save it!
# ggsave("explore/H1&2_neg_and_pos_concept_histogram_wide_20231117.png", p,
#        width = 12, height = 4, units = "in")
# # and again but narrower
# ggsave("explore/H1&2_neg_and_pos_concept_histogram_narrow_20231117.png", p,
#        width = 6, height = 4, units = "in")



##### Now make a conceptual trend line to compliment the histograms

# 1. Simulate data with a logarithmic trend
set.seed(123)
n <- 100
x <- seq(0.1, 5, length.out = n)  # x should be > 0 for logarithmic relationship
y <- 10 - 3 * log(x) + rnorm(n, 0, 5)  # Logarithmic relationship with some noise

df <- data.frame(x = x, y = y)
rm(n,x,y)

# 2. Fit a logarithmic regression model
model <- lm(y ~ log(x), data = df)

# ## thin df to not include the positive tilt at the end 
# df = df[df$x < 1,]

# 3. Plot the fitted line with confidence interval
p =
ggplot(df, aes(x = x, y = y)) +
  geom_smooth(method = "lm", formula = y ~ log(x), se = TRUE, col = "darkorange4") +
  labs(title = "H1: Negative predator-prey relationship", 
       x = "Predator abundance", y = "Prey abundance") +
  guides(color = "none") +
  theme_test() +  
  # coord_cartesian(xlim = c(0.1, 3))+
  theme(
    axis.text = element_blank(),  # Removes axis values
    axis.ticks = element_blank(), # Removes axis ticks
    axis.title.x = element_text(size = 20, family = "Helvetica"),  # Adjusts x-axis title properties
    axis.title.y = element_text(size = 20, family = "Helvetica")   # Adjusts y-axis title properties
  )

# 4. save it!
ggsave("explore/H1_negative_concept_pred-prey_relationship_20231116.png", p,
       width = 6, height = 5, units = "in")


#### Now reverse it for the positive trend! 

# 1. Simulate data with an exponential trend
set.seed(123)
n <- 100
x <- seq(0.1, 5, length.out = n)
y <- 2 + 1.5 * exp(0.5 * x) + rnorm(n, 0, 7)  # Exponential relationship with some noise

df <- data.frame(x = x, y = y)
rm(n,x,y)

# 2. Fit an exponential regression model using a transformation
model <- lm(log(y - min(y) + 1) ~ x, data = df)

# 3. Plot the fitted line with confidence interval
p=
ggplot(df, aes(x = x, y = y)) +
  geom_smooth(method = "lm", formula = y ~ exp(x), se = TRUE, col = "springgreen4") +
  labs(title = "H2: Positive prey-predator relationship", 
       x = "Prey abundance", y = "Predator abundance") +
  guides(color = "none") +
  theme_test() +  
  theme(
    axis.text = element_blank(),  # Removes axis values
    axis.ticks = element_blank(), # Removes axis ticks
    axis.title.x = element_text(size = 20, family = "Times New Roman"),  # Adjusts x-axis title properties
    axis.title.y = element_text(size = 20, family = "Times New Roman")   # Adjusts y-axis title properties
  )

## 4. save it! 
# ggsave("explore/H2_positive_concept_pred-prey_relationship_20231116.png", p,
#        width = 6, height = 5, units = "in")

## keep it clean!
rm(df, concept_col, model, p)

############# Visualize pairwise interaction comparison from opposing directions ######

#
##
### THESE ARE NOT USED ANY MORE, FULLY HASHED OUT FOR NOW
##
#

# ### The goal here is to assess both directions of a species pair, e.g.:
# ## SUB-Tapirus_indicus~DOM-Panthera_tigris vs SUB-Panthera_tigris~DOM-Tapirus_indicus
# # To determine the overall direciton of the relationship. 
# 
# preform[preform$Species_Pair == "SUB-Tapirus_indicus~DOM-Panthera_tigris", 
#         c("Interaction_Estimate","lower","upper", "Significance")] # clear positive effect of tigers on taps 
# 
# preform[preform$Species_Pair == "SUB-Panthera_tigris~DOM-Tapirus_indicus", 
#         c("Interaction_Estimate","lower","upper", "Significance")] # unclear negative effect of taps on tigers. 
# 
# ## We would move in favor of the larger and more significant response --> 
# ## tigers possibly bottom-up limited by taps, but both could be sharing responses to unmeasured covarites.  
# 
# ## what is the difference in the interaction?
# preform$Interaction_Estimate[preform$Species_Pair == "SUB-Tapirus_indicus~DOM-Panthera_tigris"] -
#   preform$Interaction_Estimate[preform$Species_Pair == "SUB-Panthera_tigris~DOM-Tapirus_indicus"] 
# # Again, this is a larger positive difference in interactions, suggesting a positive (bottom-up relaitonship)
# # but make sure to calculate difference when large carnivores are dominant FIRST! (maybe?)
# 
# 
# ## Make a plot comparing both interactions (i.e. when large carnivores are dominant AND subordinate)
# ## to determine direction of relationship for each species pairs 
# plots = list() #store plots here
# 
# 
# for(i in 1:length(unique(coeff$Species_Pair))){
#   
#   ## grab the name of the species pair 
#   n = unique(coeff$Species_Pair)[i]
#   
#   ##subset coeff for it 
#   a = coeff[coeff$Species_Pair == n, ]
#   
#   # and for the the complementary species pair data.
#   b = subset(coeff, sub_sp == unique(a$dom_sp) & dom_sp == unique(a$sub_sp))
#   
#   ## then combine
#   dat = rbind(a[a$var == "Species_Interaction",],
#               b[b$var == "Species_Interaction",])
#   
#   ## bypass models that are missing their complement 
#   if(nrow(dat)<2){
#     next
#   }
#   
#   ## add a col for which position by looking for all dominant large carnivores
#   dat$dom_position = ifelse(dat$dom_sp %in% c("Canis_lupus_familiaris","Panthera_tigris",
#                                               "Panthera_pardus","Neofelis_genus",
#                                               "Cuon_alpinus"),
#                             "Large_carnivore", "Prey")
#   
#   ## set the dom_pos as a factor 
#   # dat$dom_position = factor(dat$dom_position, levels = c( "Large_carnivore", "Prey"))
#   
#   # Add a column for significance marker
#   dat$marker <- ifelse(dat$sig == "Significant", "*", "")
#   
#   ## if there is a large carnivore pairwise comparison,
#   if(length(unique(dat$dom_position)) == 1){
#     
#     # Calculate the difference in effect sizes just by 1-2
#     dat$diff = dat$mean[1] - dat$mean[2]
#     
#     # and change dom_guild to the large carnivore species names
#     dat$dom_position = dat$dom_sp
#     
#     # but if its not two large carnivores, 
#   }else{
#     
#     # Calculate the difference in effect sizes following a standard order
#     dat$diff = dat$mean[dat$dom_position == "Large_carnivore"] - dat$mean[dat$dom_position == "Prey"]
#     
#   }
#   
#   ## combine both species pairs for the ID
#   id = paste(dat$Species_Pair[1], dat$Species_Pair[2], sep = " & ")
#   
#   ## combine dom and sub species for a title 
#   title = paste(dat$dom_sp[1], dat$sub_sp[1], sep = " & ")
#   
#   ## set bar dodge width 
#   dodge_width = 0.8
#   
#   # Assign colors based on unique dominant positions
#   colors <- c("Prey" = "mediumslateblue", "Large_carnivore" = "mediumseagreen", 
#               "Cuon_alpinus" = "red1", "Panthera_pardus" = "purple1", 
#               "Canis_lupus_familiaris" = "grey", "Panthera_tigris" = "orange1", 
#               "Neofelis_genus" = "blue1")
#   colors <- colors[unique(dat$dom_position)]
#   
#   ## make the plot
#   p =
#     ggplot(dat, aes(x = var, y = mean, fill = dom_position)) +
#     geom_bar(stat = "identity", position = position_dodge(width = dodge_width), color = "black", width = dodge_width) +
#     geom_errorbar(aes(ymin = lower, ymax = upper), position = position_dodge(width = dodge_width), width = 0, color = "black") +
#     geom_text(aes(y = mean+.2*sign(mean), label = marker), position = position_dodge(width = dodge_width), size = 8) +
#     geom_text(aes(y = .04, label = round(Rhat, 2)), position = position_dodge(width = dodge_width), size = 3,  hjust = -.3) +
#     geom_text(aes(y = max(upper), x = 1, label = paste("Difference =", round(unique(diff),2))), size = 4, vjust = 1, hjust = 1) +  
#     geom_hline(yintercept = 0)+
#     scale_fill_manual(values = colors) +  # Assign dynamically generated color scale
#     labs(x = "Species Interaction Estimate", y = "Mean Effect Size", fill = "Dominant Guild", title = title) +
#     theme_test()+
#     theme(
#       axis.text.y = element_text(size = 12), # Increase y-axis tick label font size
#       axis.text.x = element_blank() 
#     )
#   
#   ## save it
#   plots[[i]] = p
#   names(plots)[i] = id
#   
# }
# rm(p,i,dodge_width, a,b, title, n, id, dat, colors)
# 
# names(plots)
# plots[[387]] # "" in title means there was no complementary model 
# 
# ## remove plots w/ NA
# plots = plots[names(plots) != ""]
# ## need to remove redundant plots tho 
# 
# ## Chat GPT to the rescue?
# # Function to sort and concatenate species pairs in canonical form
# getCanonicalForm <- function(label) {
#   sorted_pairs <- sort(strsplit(label, " & ")[[1]])
#   return(paste(sorted_pairs, collapse = " & "))
# }
# 
# # Get the canonical forms of the plot labels
# canonical_labels <- sapply(names(plots), getCanonicalForm)
# 
# # Identify duplicate plots based on canonical labels
# duplicate_indices <- duplicated(canonical_labels)
# 
# # Remove duplicate plots
# plots <- plots[!duplicate_indices]
# 
# names(plots) # seems better! 
# plots$`SUB-Canis_lupus_familiaris~DOM-Presbytis_rubicunda & SUB-Presbytis_rubicunda~DOM-Canis_lupus_familiaris` # good + interesting result
# plots$`SUB-Neofelis_genus~DOM-Cuon_alpinus & SUB-Cuon_alpinus~DOM-Neofelis_genus` # its working! 
# 
# ## keep it clean
# rm(getCanonicalForm, canonical_labels, duplicate_indices)
# 
# 
# ## save each of these graphs! 
# for(i in 1:length(plots)){
#   
#   # select one plot 
#   p = plots[[i]]
#   
#   # Isolate the relevant species pair
#   id = str_split(names(plots)[i], " & ")[[1]][1]
#   id = paste((str_split(id, "~")[[1]]), collapse=" & ")
#   id = gsub("SUB-", "", id)
#   id = gsub("DOM-", "", id)
#   
#   # grab today's date
#   day<-str_sub(Sys.Date(),-2)
#   month<-str_sub(Sys.Date(),-5,-4)
#   year<-str_sub(Sys.Date(),-10,-7)
#   date = paste(year,month,day, sep = "")
#   
#   # construct a saving path 
#   path = paste("figures/Double pairwise species interaction comparison/", id, "_", date, ".png", sep = "")
#   
#   # save it! 
#   ggsave(path, p, width = 6, height = 5, units = "in")
#   
#   
# }
# rm(p, path, day, month, year, date, id, i)



##### Quick trouble shooting pt1, 20230614 ####

## Based on coefficent dataframe import, ~20230614

### Inspect which species pairs are present/missing, 20230615
setdiff(all_combos, coeff$Species_Pair) # only 20 missing models!! This is an improvement, but still check. 
which(all_combos == "SUB-Paradoxurus_hermaphroditus~DOM-Cuon_alpinus") # currently re-running on HPC
which(all_combos == "SUB-Arctictis_binturong~DOM-Canis_lupus_familiaris") # 2nd fix made this work! 


### use the 3 examples listed above to debug 

# pt1
check = bdata$`SUB-Prionailurus_bengalensis~DOM-Panthera_tigris` # Error file: one node produced an error: Error in node y.sub[691,3]

# inspect problematic cell of matrix
as.numeric(check$y.sub[691,3]) # zero
as.numeric(check$y.dom[691,3]) # also zero

# inspect the rest of the row 
check$y.dom[691,]
check$y.sub[691,] # seems ok... 

# which SU is this?
rownames(check$y.dom)[691] # "Sabah_DVCA_163_2019_c"...

# chekc SiteCovs
check$elev[691]
check$flii[691]
check$hfp[691] # all seem fine

# check obsCov
check$cams[691,3] # woah, crazy value...
check$cams[691,] # but not necessarily more/less crazy than the others... 

# check iZIP
check$Z.dom[691]
check$Z.sub[691] # these make sense. 

# check RE
check$area[691] # 5
check$year[691] # 15
check$source[691] # 6, all seem fine 

# make sure sites are fine. 
check$nsites #4658
dim(check$y.dom)
dim(check$y.sub) # both 4658

#
##
### pt2
check = bdata$'SUB-Rusa_unicolor~DOM-Cuon_alpinus' #one node produced an error: Error in node y.sub[675,3]

# inspect problematic cell of matrix
as.numeric(check$y.dom[675,3]) # zero
as.numeric(check$y.sub[675,3]) # 1

# inspect the rest of the row 
check$y.dom[675,]
check$y.sub[675,] # lots of numbers here

# which SU is this?
rownames(check$y.dom)[675] # SAME SU ==> "Sabah_DVCA_163_2019_c"...

# chekc SiteCovs
check$elev[675]
check$flii[675]
check$hfp[675] # all seem fine

# check obsCov
check$cams[675,3] # woah, crazy value...
check$cams[675,] # but not necessarily more/less crazy than the others... 

# check iZIP
check$Z.dom[675]
check$Z.sub[675] # these make sense. 

# check RE
check$area[675] # 4
check$year[675] # 15
check$source[675] # 6, all seem fine 

# make sure sites are fine. 
check$nsites #4516
dim(check$y.dom)
dim(check$y.sub) # both 4658

#
##
###
#### pt 3
check = bdata$"SUB-Neofelis_genus~DOM-Paguma_larvata" # one node produced an error: Error in node y.dom[675,4]

# inspect problematic cell of matrix
as.numeric(check$y.dom[675,3]) # 2
as.numeric(check$y.sub[675,3]) # 1

# inspect the rest of the row 
check$y.dom[675,]
check$y.sub[675,] # looks fine

# which SU is this?
rownames(check$y.dom)[675] # SAME SU ==> "Sabah_DVCA_163_2019_c"... AGAIN! 
### Three times is NOT random! 
## there is an error w/ this sampling unit! 
# Investigate! 

meta = read.csv("data/send_to_HPC/clean_metadata_to_make_UMFs_20230615.csv")
prob_meta = meta[meta$cell_id_3km == "Sabah_DVCA_163_2019_c",]
prob_meta$cameras_included # this is alot, figure out how many!
prob_surv = meta[meta$survey_id == "Danum_Valley_CA_2019_c_Luskin",]
prob_surv$cameras_included # problematic SU is an outlier! Probably exclude SUs w/ this many active cams --> skews models. 

## need to find out how many cameras are stored per SU
prob_surv$cams_included_count = str_count(prob_surv$cameras_included, " - ") + 1
head(prob_surv[, c("cameras_included", "cams_included_count")])

## how many cams in the problem SU?
prob_surv$cams_included_count[prob_surv$cell_id_3km == "Sabah_DVCA_163_2019_c"] # 14

## check how this applys across all data
meta$cams_included_count = str_count(meta$cameras_included, " - ") + 1
summary(meta$cams_included_count) # 3rd quarter is 2, so 14 is big!!
hist(meta$cams_included_count) # classic Poisson 

## which are the problematic SUs?
bad_cams = meta$cell_id_3km[meta$cams_included_count >= 14]

## how many SUs do we lose by exclusing them?
nrow(meta[meta$cell_id_3km %in% bad_cams,]) / nrow(meta) * 100
## lose only 1% of all SUs --> good enough! 


### Double-check a good model 
check_good = bdata[[preform$Species_Pair[1]]]

rownames(check_good$y.dom)[grepl("DVCA", rownames(check_good$y.dom))] # No DVCA sites here! 
check_good$y.dom["Sabah_DVCA_163_2019_c"]

# also char
str(check_good$y.sub) # also char

##### Quick trouble shooting pt2, 20230615 ####### 

## Based off PPC data import on 20230615 w/ 44 missing pairs. 

which(all_combos == "SUB-Panthera_tigris~DOM-Rusa_unicolor") # Error file: one node produced an error: Error in node y.dom[683,1], Node inconsistent with parents
which(all_combos == "SUB-Tragulus_genus~DOM-Cuon_alpinus") # Error file: one node produced an error: Error in node y.sub[682,1], Node inconsistent with parents
which(all_combos == "SUB-Echinosorex_gymnura~DOM-Panthera_pardus") # Error file: one node produced an error: Error in node y.sub[596,1]


# pt1
check = bdata$"SUB-Panthera_tigris~DOM-Rusa_unicolor" # Error file: one node produced an error: Error in node y.dom[683,1]

# inspect problematic cell of matrix
as.numeric(check$y.sub[683,1]) # zero
as.numeric(check$y.dom[683,1]) # 3

# inspect the rest of the row 
check$y.dom[683,] # lots of numbers 
check$y.sub[683,] # all zero 

# which SU is this?
rownames(check$y.dom)[683] # "Sabah_DVCA_163_2020_a" --> same site as last time, but different time! 

# chekc SiteCovs
check$elev[683]
check$flii[683]
check$hfp[683] # all seem fine

# check obsCov
check$cams[683,3] # woah, crazy value...
check$cams[683,] #biiiiig values here! 




# pt2
check = bdata$"SUB-Tragulus_genus~DOM-Cuon_alpinus" # Error file: one node produced an error: Error in node y.sub[682,1], Node inconsistent with parents

# which SU is this?
rownames(check$y.dom)[682] # "Sabah_DVCA_163_2020_a" --> same site as above

# check obsCov
check$cams[682,] #biiiiig values here! 




# pt3
check = bdata$"SUB-Echinosorex_gymnura~DOM-Panthera_pardus" # Error file: one node produced an error: Error in node y.sub[596,1]

# which SU is this?
rownames(check$y.dom)[596] # "Sabah_DVCA_163_2020_a" --> same site as above

# check obsCov
check$cams[596,] #biiiiig values here! 


#### inspect problem cam 
meta = read.csv("data/send_to_HPC/clean_metadata_to_make_UMFs_20230616.csv")
prob_meta = meta[meta$cell_id_3km == "Sabah_DVCA_163_2020_a",]
prob_meta$cams_included_count #9 cams active here. 
## which are the problematic SUs?
bad_cams = meta$cell_id_3km[meta$cams_included_count >= 9]
## how many SUs do we lose by exclusing them?
nrow(meta[meta$cell_id_3km %in% bad_cams,]) / nrow(meta) * 100
## lose 3% of all SUs --> notideal but ok. 

## compare active cams below and above 9 threshold 
summary(meta$cams_included_count)
summary(meta$cams_included_count[meta$cams_included_count < 9])

## check how much data we lose if were even more conservative and go to 7 cams
bad_cams = meta$cell_id_3km[meta$cams_included_count >= 7]
## how many SUs do we lose by exclusing them?
nrow(meta[meta$cell_id_3km %in% bad_cams,]) / nrow(meta) * 100
# lose 5.6% of all SUs --> go w/ this as an acceptable cut off. 

## are all landscapes still represented?
unique(meta$Landscape[! meta$cell_id_3km %in% bad_cams]) # yes! 


#### One more trouble shoot on 20230724

check = bdata$"SUB-Cuon_alpinus~DOM-Paradoxurus_hermaphroditus"
check$y.dom[1688,] # seems ok... 
check$y.dom[1689,]
rownames(check$y.dom)[1688:1690] # Singapore_NAS_139
unique(og_meta$Landscape[og_meta$land == "NAS"])
nrow(og_captures[og_captures$Landscape == "Singapore" & og_captures$Species == "Paradoxurus hermaphroditus",]) # plenty of detections! 

## check sampling effort b/c its been problematic in the past
check$cams[1688,] # looooots of active cams here! 
## how many?
summary(og_resamp_captures$num_cams_active_at_date[og_resamp_captures$cell_id_3km == rownames(check$y.dom)[1688]])
## maxes out at 6.... odd. 
og_resamp_captures$num_cams_active_at_date[og_resamp_captures$cell_id_3km == rownames(check$y.dom)[1688] &
                                             og_resamp_captures$Species == "Paradoxurus_hermaphroditus"]
## ahhh, but its ONLY 5s and 6s....
## just accept this as a bad mod! 

## double check covariates
check$flii[1688];summary(check$flii) # very low FLII
check$elev[1688];summary(check$elev) # very low ELEV
check$hfp[1688];summary(check$hfp) # very high HFP
check$hunt[1688];summary(check$hunt) # higher range HUNT

## Generally consider this an outlier along multiple fronts, accept model failure. 

