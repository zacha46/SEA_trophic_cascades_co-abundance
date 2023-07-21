##### Inspect and visualize completed co-abundance models
#### Data was prepared on R Script: Step 1_Generate co-abundance data bundles.Rmd
#### and ran on HPC on R Script: HPC_co-abundance_model.R
### Will convert this to Rmd when ready, but R script for exploring now
## Import dataframes, not complete models b/c theyre too heavy

# Zachary Amir, Z.Amir@uq.edu.au
# Last updated on July 20th, 2023

## start fresh
rm(list = ls())

## load libs 
library(tidyverse)
library(plyr)
library(reshape2) # for making checkerboard matrix 
library(fs) # for moving files around 
library(ggridges) # for ridgeline plots 


# library(leaflet) # for creating concepts of study map
# library(colortools) # for nice colors in the map 
# library(jagsUI) 
# library(traitdata) # for species trait data

## set WD (should be automatic b/c RProj, but being safe anyway)
setwd("/Users/zachary_amir/Dropbox/Zach PhD/Trophic release project/SEA_trophic_cascades_co-abundance")

## Import clean cam trap data here for referencing
og_resamp_captures = read.csv("data/send_to_HPC/clean_captures_to_make_UMFs_20230718.csv")
og_resamp_meta = read.csv("data/send_to_HPC/clean_metadata_to_make_UMFs_20230718.csv")

## should also import NON-resampled data for referencing too! 
og_captures = read.csv("/Users/zachary_amir/Dropbox/CT capture histories database/Asian ECL raw CT data/Step4_output_pre-resampling/Clean_independent_captures_20230610.csv")
og_meta = read.csv("/Users/zachary_amir/Dropbox/CT capture histories database/Asian ECL raw CT data/Step4_output_pre-resampling/Clean_independent_metadata_20230610.csv")


######## Import co-abundance coefficent dataframes ######

## first, I will import the bundled data to generate a vector of relevant species pairs 
bdata = readRDS("data/send_to_HPC/Bundled_data_for_Bayes_co-abundance_mods_394_species_pairs_20230720.RDS")
all_combos = names(bdata)


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
files = list.files("results/coefficent_dataframes/")
files = files[!grepl("OLD", files)] # remove any old data 

## Have a few long and all middle files in this directory, subset for one
# files = files[grepl("LONG", files)]
# files = files[grepl("MIDDLE", files)]
files = files[grepl("SHORT", files)]


### Noticed a typo that included generated results for species that dont spatially overlap
## and these have been removed from bdata, but not worth re-running models b/c they take long.
# so thin files to match relevant combos and exclude pointless mods via loop
### This will be irrelevant w/ updated models and testing. 

# create a new files vector to save good files
matched_files <- vector("character", length(all_combos))

for(i in 1:length(all_combos)){
  # if there are any TRUE matches for all_combos in files
  if(any(grepl(all_combos[i], files))){
    ## save that file name in the new vector
    matched_files[i] = files[grepl(all_combos[i], files)]
  }
}
rm(i)

## remove matched_files w/ blank names (likely b/c mod never left HPC)
matched_files = matched_files[matched_files != ""]

## update files w/ the good ones (for consistent naming)
files = matched_files
rm(matched_files)


# store results here
coeff.res = list()

# loop thru each file to import into list
for(i in 1:length(files)){

    d = read.csv(paste("results/coefficent_dataframes/", files[i], sep = ""))
    
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
head(coeff)
str(coeff)
rm(coeff.res)

summary(coeff$n.eff[coeff$var == "Species_Interaction"]) # not good! short settings tho
summary(coeff$Rhat[coeff$var == "Species_Interaction"]) # The Meidan is acceptable! but some high values exist. 


### Inspect which species pairs are present/missing, 20230719, SHORT
setdiff(all_combos, coeff$Species_Pair) # only 6 missing models!! This is a huge improvement, but still check. 
which(all_combos == "SUB-Panthera_tigris~DOM-Helarctos_malayanus")       # exceeded 4 hour wall time, cancelled at 4.5 hours
which(all_combos == "SUB-Canis_lupus_familiaris~DOM-Lophura_diardi")     # exceeded 4 hour wall time, cancelled at 4.5 hours
which(all_combos == "SUB-Neofelis_genus~DOM-Arctonyx_collaris")          # exceeded 4 hour wall time, cancelled at 4.5 hours
which(all_combos == "SUB-Canis_lupus_familiaris~DOM-Argusianus_argus")   # exceeded 4 hour wall time, cancelled at 4.5 hours
which(all_combos == "SUB-Herpestes_genus~DOM-Neofelis_genus")            # exceeded 4 hour wall time, cancelled at 4.5 hours
which(all_combos == "SUB-Neofelis_genus~DOM-Bos_taurus")                 # exceeded 4 hour wall time, cancelled at 4.5 hours

### Long-story-short, I found out they all had problem SUs from DVCA that contained 14+ active cams 
## --> removed all SUs w/ > 14 cameras and re-run analysis. 
## 2nd update --> removed all SUs w/ > 7 cameras and re-run analysis with MUCH better performance. 
## 3rd update --> Seems to be going well! only probelmatic models are excluded (i.e. took too long to run)

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
table(preform$mod_completion) # 1.5% of models have failed! Epic improvement. 

#
##
###
#### Combine species trait/guild data w/ model performance

#load guild data
guilds = read.csv(("data/species_traits/clean_51_species_trait_data_20230718.csv"))
head(guilds)
str(guilds)

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

# ## Set factor levels for guilds for later re-ordering
# preform$SUB_guild = factor(preform$SUB_guild, levels = c("Large_Carnivore", "Small_Carnivore",
#                                                          "Large_Omnivore", "Small_Omnivore",
#                                                          "Large_Herbivore", "Small_Herbivore"))
# preform$DOM_guild = factor(preform$DOM_guild, levels = c("Large_Carnivore", "Small_Carnivore",
#                                                          "Large_Omnivore", "Small_Omnivore",
#                                                          "Large_Herbivore", "Small_Herbivore"))
# ## set the order as the dom species
# preform = preform[order(preform$DOM_guild, preform$dom_species),]

## clean up 
rm(add)


######## Import co-abundance PPC dataframes ######

## Make sure all_combos is loaded here! 

# list all files
files = list.files("results/PPC_dataframes/")
files = files[!grepl("OLD", files)] # remove any old data 

## Subset files for one setting, especially if multiple are in the mix. 
# files = files[grepl("LONG", files)]
# files = files[grepl("MIDDLE", files)]
files = files[grepl("SHORT", files)]

### Noticed a typo that included generated results for species that dont spatially overlap
## and these have been removed from bdata, but not worth re-running models b/c they take long.
# so thin files to match relevant combos and exclude pointless mods via loop
### This will be irrelevant w/ updated models and testing. 

# create a new files vector to save good files
matched_files <- vector("character", length(all_combos))

for(i in 1:length(all_combos)){
  # if there are any TRUE matches for all_combos in files
  if(any(grepl(all_combos[i], files))){
    ## save that file name in the new vector
    matched_files[i] = files[grepl(all_combos[i], files)]
  }
}
rm(i)

## remove matched_files w/ blank names (likely b/c mod never left HPC)
matched_files = matched_files[matched_files != ""]

## update files w/ the good ones (for consistent naming)
files = matched_files
rm(matched_files)

# store results here
ppc_res = list()

# loop thru each file to import into list
for(i in 1:length(files)){
  
  ## read the file 
  d = read.csv(paste("results/PPC_dataframes/", files[i], sep = ""))
  
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

### UPDATE 2023-07-19, typo in updated co-abundance model did not let this one generate! 
## It has been fixed for next round tho. 
ppc_values = do.call(rbind, ppc_res$values)
head(ppc_values) 
ppc_values$X = NULL # damn row names 
str(ppc_values) # looks good! 

## quicky inspect parameter convergence for species interaction parameter
summary(ppc_values$Rhat) # looks decent! 
nrow(ppc_values[ppc_values$Rhat > 1.2,])/nrow(ppc_values) * 100 # 3/4ths of mods have good interaction param. 
# how about sample size?
summary(coeff$n.eff[coeff$var == "Species_Interaction"]) # this could certainly be improved! 
length(coeff$n.eff[coeff$var == "Species_Interaction" & 
                     coeff$n.eff < 100]) / length(coeff$n.eff[coeff$var == "Species_Interaction"]) * 100 # 76% of mods have small (< 100) sample size. 

### Seems like middle settings are decent, but definitely not up to snuff! 

# remove list now that were all stored in dfs. 
rm(ppc_res)

### Inspect which species pairs are present/missing, 20230719
setdiff(all_combos, ppc_plotdat$Species_Pair) #6 missing models, same as above


######## Import co-abundance prediction dataframes ######


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

# List all files in the coefficents folder
file_list <- dir("results/prediction_dataframes/", full.names = TRUE)

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



## Make sure all_combos is loaded here! 

# list all files
files = list.files("results/prediction_dataframes/")
files = files[!grepl("OLD", files)] # remove any old data 

## Subset files for one setting, especially if multiple are in the mix. 
# files = files[grepl("LONG", files)]
# files = files[grepl("MIDDLE", files)]
files = files[grepl("SHORT", files)]

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
  d = read.csv(paste("results/prediction_dataframes/", files[i], sep = ""))
  
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


### Inspect which species pairs are present/missing,  20230719
setdiff(all_combos, pred_abund$Species_Pair) #6 missing models, same as above


######### Track model performance and validity #######

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
head(preform[preform$mod_completion == "uncompleted",]) # all NA values here for BPV,Chat,Rhat,etc b/c the models never finished! 
str(preform) # looks good! 

## Add a column denoting if the BPV values converged
preform$BPV_valid = "No"
preform$BPV_valid[preform$BPV.dom >= 0.25 & preform$BPV.dom <= 0.75 &
                    preform$BPV.sub >= 0.25 & preform$BPV.sub <= 0.75] = "Yes"
table(preform$BPV_valid) # mostly invalid rn, but thats w/ short settings. 
## Still mostly invalid with middle settings. 70% are not a good fit. 

## add a column denoting if overdispersion remains
preform$OD_valid = "No"
preform$OD_valid[preform$Chat.dom >= 0.98 & preform$Chat.dom <= 1.1 &
                     preform$Chat.sub >= 0.98 & preform$Chat.sub <= 1.1] = "Yes"
table(preform$OD_valid) # mostly OD, but thats w/ short settings. TBH, better preformance than expected w/ low settings. 
## Middle settings are very good, 91% have no OD
## I bet this is because active_cams + RE for source are actually good for detection models... less reliance on ODRE

## add a column denoting if the species interaction parameter converged 
preform$parameter_valid = "No"
preform$parameter_valid[preform$Rhat >= 0.99 & preform$Rhat <= 1.2] = "Yes"
table(preform$parameter_valid) # majority is valid! Thats exciting! Should be better w/ more reps. 
## Middle settings are good, 76% have good interaction parameters 

## Finally, use all of this information to denote if the overall model is valid and should be trusted
preform$mod_valid = "No"
preform$mod_valid[preform$BPV_valid == "Yes" &
                    preform$OD_valid == "Yes" &
                    preform$parameter_valid == "Yes"] = "Yes"
table(preform$mod_valid) # Vast majority are NOT valid, but this checks out b/c of short settings. 
## middle settings have 25% of all mods as totally valid --> improvement! 

## Save which species pairs are valid to examine data later
valid_mods = preform$Species_Pair[preform$mod_valid == "Yes"]

## Which species pairs involved large carnivores significantly regulating the abundance of other critters?
preform$Species_Pair[preform$Interaction_Estimate < 0 & 
                       preform$Significance == "Significant" &
                       preform$mod_valid == "Yes" &
                       grepl("DOM-Large_Carnivore", preform$guild_pair)] 

# only 9 with middle settings 
# c("SUB-Lophura_nycthemera~DOM-Panthera_tigris", "SUB-Sus_barbatus~DOM-Cuon_alpinus", "SUB-Trichys_fasciculata~DOM-Panthera_tigris")
## but a lot more when we exclude the mod_valid part! 

## Curious to see what happens w/ mods w/ higher MCMC settings. 



######### Compare pairwise interactions in different directions ######

### The goal here is to assess both directions of a species pair, e.g.:
## SUB-Tapirus_indicus~DOM-Panthera_tigris vs SUB-Panthera_tigris~DOM-Tapirus_indicus
# To determine the overal direciton of the relationship. 

preform[preform$Species_Pair == "SUB-Tapirus_indicus~DOM-Panthera_tigris", 
        c("Interaction_Estimate","lower","upper", "Significance")] # clear positive effect of tigers on taps 

preform[preform$Species_Pair == "SUB-Panthera_tigris~DOM-Tapirus_indicus", 
        c("Interaction_Estimate","lower","upper", "Significance")] # unclear negative effect of taps on tigers. 

## We would move in favor of the larger and more significant response --> 
## tigers possibly bottom-up limited by taps, but both could be sharing responses to unmeasured covarites.  

## what is the difference in the interaction?
preform$Interaction_Estimate[preform$Species_Pair == "SUB-Tapirus_indicus~DOM-Panthera_tigris"] -
  preform$Interaction_Estimate[preform$Species_Pair == "SUB-Panthera_tigris~DOM-Tapirus_indicus"] 
# Again, this is a larger positive difference in interactions, suggesting a positive (bottom-up relaitonship)
# but make sure to calculate difference when large carnivores are dominant FIRST! (maybe?)


## Make a plot comparing both interactions (i.e. when large carnivores are dominant AND subordinate)
## to determine direction of relationship for each species pairs 
plots = list() #store plots here

## First, create a dom + sub species column 
coeff <- within(coeff, {
  sub_sp <- sapply(strsplit(Species_Pair, "~"), function(x) x[1])
  dom_sp <- sapply(strsplit(Species_Pair, "~"), function(x) x[2])
})
coeff$dom_sp = gsub("DOM-", "", coeff$dom_sp)
coeff$sub_sp = gsub("SUB-", "", coeff$sub_sp)


for(i in 1:length(unique(coeff$Species_Pair))){
  
  ## grab the name of the species pair 
  n = unique(coeff$Species_Pair)[i]
  
  ##subset coeff for it 
  a = coeff[coeff$Species_Pair == n, ]
  
  # and for the the complementary species pair data.
  b = subset(coeff, sub_sp == unique(a$dom_sp) & dom_sp == unique(a$sub_sp))
 
  ## then combine
  dat = rbind(a[a$var == "Species_Interaction",],
              b[b$var == "Species_Interaction",])
  
  ## bypass models that are missing their complement 
  if(nrow(dat)<2){
    next
  }
  
  ## add a col for which position by looking for all dominant large carnivores
  dat$dom_position = ifelse(dat$dom_sp %in% c("Canis_lupus_familiaris","Panthera_tigris",
                                              "Panthera_pardus","Neofelis_genus",
                                              "Cuon_alpinus"),
                            "Large_carnivore", "Prey")
  
  ## set the dom_pos as a factor 
  # dat$dom_position = factor(dat$dom_position, levels = c( "Large_carnivore", "Prey"))
  
  # Add a column for significance marker
  dat$marker <- ifelse(dat$sig == "Significant", "*", "")
  
  ## if there is a large carnivore pairwise comparison,
  if(length(unique(dat$dom_position)) == 1){
    
    # Calculate the difference in effect sizes just by 1-2
    dat$diff = dat$mean[1] - dat$mean[2]
    
    # and change dom_guild to the large carnivore species names
    dat$dom_position = dat$dom_sp
    
    # but if its not two large carnivores, 
  }else{
    
    # Calculate the difference in effect sizes following a standard order
    dat$diff = dat$mean[dat$dom_position == "Large_carnivore"] - dat$mean[dat$dom_position == "Prey"]
    
  }
  
  ## combine both species pairs for the ID
  id = paste(dat$Species_Pair[1], dat$Species_Pair[2], sep = " & ")
  
  ## combine dom and sub species for a title 
  title = paste(dat$dom_sp[1], dat$sub_sp[1], sep = " & ")
  
  ## set bar dodge width 
  dodge_width = 0.8
  
  # Assign colors based on unique dominant positions
  colors <- c("Prey" = "mediumslateblue", "Large_carnivore" = "mediumseagreen", 
              "Cuon_alpinus" = "red1", "Panthera_pardus" = "purple1", 
              "Canis_lupus_familiaris" = "grey", "Panthera_tigris" = "orange1", 
              "Neofelis_genus" = "blue1")
  colors <- colors[unique(dat$dom_position)]
  
  ## make the plot
  p =
    ggplot(dat, aes(x = var, y = mean, fill = dom_position)) +
    geom_bar(stat = "identity", position = position_dodge(width = dodge_width), color = "black", width = dodge_width) +
    geom_errorbar(aes(ymin = lower, ymax = upper), position = position_dodge(width = dodge_width), width = 0, color = "black") +
    geom_text(aes(y = mean+.2*sign(mean), label = marker), position = position_dodge(width = dodge_width), size = 8) +
    geom_text(aes(y = .04, label = round(Rhat, 2)), position = position_dodge(width = dodge_width), size = 3,  hjust = -.3) +
    geom_text(aes(y = max(upper), x = 1, label = paste("Difference =", round(unique(diff),2))), size = 4, vjust = 1, hjust = 1) +  
    geom_hline(yintercept = 0)+
    scale_fill_manual(values = colors) +  # Assign dynamically generated color scale
    labs(x = "Species Interaction Estimate", y = "Mean Effect Size", fill = "Dominant Guild", title = title) +
    theme_test()+
    theme(
      axis.text.y = element_text(size = 12), # Increase y-axis tick label font size
      axis.text.x = element_blank() 
    )
  
  ## save it
  plots[[i]] = p
  names(plots)[i] = id
  
}
rm(p,i,dodge_width, a,b, title, n, id, dat, colors)

names(plots)
plots[[383]] # "" in title means there was no complementary model 

## remove plots w/ NA
plots = plots[names(plots) != ""]
## need to remove redundant plots tho 
# plots$`SUB-Geokichla_citrina~DOM-Cuon_alpinus & SUB-Cuon_alpinus~DOM-Geokichla_citrina` # its working! 
# plots$`SUB-Cuon_alpinus~DOM-Geokichla_citrina & SUB-Geokichla_citrina~DOM-Cuon_alpinus` # also working, but need to remove redundant graphs! 

## Chat GPT to the rescue?
# Function to sort and concatenate species pairs in canonical form
getCanonicalForm <- function(label) {
  sorted_pairs <- sort(strsplit(label, " & ")[[1]])
  return(paste(sorted_pairs, collapse = " & "))
}

# Get the canonical forms of the plot labels
canonical_labels <- sapply(names(plots), getCanonicalForm)

# Identify duplicate plots based on canonical labels
duplicate_indices <- duplicated(canonical_labels)

# Remove duplicate plots
plots <- plots[!duplicate_indices]

names(plots) # seems better! 
plots$`SUB-Canis_lupus_familiaris~DOM-Presbytis_rubicunda & SUB-Presbytis_rubicunda~DOM-Canis_lupus_familiaris` # good + interesting result
plots$`SUB-Neofelis_genus~DOM-Cuon_alpinus & SUB-Cuon_alpinus~DOM-Neofelis_genus` # its working! 

## keep it clean
rm(getCanonicalForm, canonical_labels, duplicate_indices)


## save each of these graphs! 
for(i in 1:length(plots)){
  
  # select one plot 
  p = plots[[i]]
  
  # Isolate the relevant species pair
  id = str_split(names(plots)[i], " & ")[[1]][1]
  id = paste((str_split(id, "~")[[1]]), collapse=" & ")
  id = gsub("SUB-", "", id)
  id = gsub("DOM-", "", id)
  
  # grab today's date
  day<-str_sub(Sys.Date(),-2)
  month<-str_sub(Sys.Date(),-5,-4)
  year<-str_sub(Sys.Date(),-10,-7)
  date = paste(year,month,day, sep = "")
  
  # construct a saving path 
  path = paste("figures/Double pairwise species interaction comparison/", id, "_", date, ".png", sep = "")
  
  # save it! 
  ggsave(path, p, width = 6, height = 5, units = "in")

  
}
rm(p, path, day, month, year, date, id, i)

#
##

######### Add data on dietary preferences for large carnivores ######

#### TBH, this should probably get added to step 1 when I generate the guilds DF anyway. 

### The goal here is to determine if large carnivores have different effects 
## when looking at 'preferred' prey species vs non-preferred species. 
# Ideally, just add a yes-no column to the guilds dataframe for each large carnivore. 

# Diet changes so much based on habitat, interactions, and time of year... 
# All species seem to prefer to some degree Sus scrofa based on readings. 

#### Tigers
# https://zslpublications.onlinelibrary.wiley.com/doi/full/10.1111/j.1469-7998.2011.00871.x
# suggests preferred weight range of 60-250 kg, and wil eat the largest available prey.... So anything >= 60 kg to accomodate tapir, gaur, and cattle
guilds$tiger_pref = "No"
guilds$tiger_pref[guilds$AdultBodyMass_g > 60000] = "Yes"
guilds$tiger_pref[guilds$scientificNameStd == "Panthera_tigris"] = "NA" # no cannibals allowed. 
# who is preferred?
guilds$scientificNameStd[guilds$tiger_pref == "Yes"] # 8 species, ok
# black bear made it in, and this MS says they are preyed upon by tigers in Laos: https://onlinelibrary.wiley.com/doi/full/10.1002/ece3.9067
#who is missing?
sort(guilds$scientificNameStd[guilds$tiger_pref == "No"]) # seems ok, 42 species

#### Leopards 
#  https://doi.org/10.1111/j.1469-7998.2006.00139.x
# suggest preferred weight range of 10-40 kg. 
# https://www.sciencedirect.com/science/article/pii/S0006320712005149 # langurs/primates, pigs when tigers are absent. 
guilds$leopard_pref = "No"
guilds$leopard_pref[guilds$AdultBodyMass_g > 10000 & guilds$AdultBodyMass_g <= 40000] = "Yes"
guilds$tiger_pref[guilds$scientificNameStd == "Panthera_pardus"] = "NA" # no cannibals allowed. 
# who is preferred?
guilds$scientificNameStd[guilds$leopard_pref == "Yes"] # 5 species, mostly carnivores
## exclude clouded leopards and dholes, leave dogs tho (evidence from india)
guilds$leopard_pref[guilds$scientificNameStd %in% c("Cuon_alpinus","Neofelis_genus")] = "No"
## also add primates and pigs based off readings
guilds$leopard_pref[grepl("Macac", guilds$scientificNameStd)] = "Yes" # leopards dont overlap w/ langurs in our dataset, so only include macaques. 
guilds$leopard_pref[guilds$scientificNameStd == "Sus_scrofa"] = "Yes"
# who is preferred now?
guilds$scientificNameStd[guilds$leopard_pref == "Yes"] # 8 species, much better! 


#### Dholes
# https://zslpublications.onlinelibrary.wiley.com/doi/10.1111/jzo.12171
## Hayward style pref wieght range of 130–190 kg, and highlights sambar as important. 
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9388674/
## scat analysis that emphasizes large deer (sambar + chital), pigs and muntjac were moderetly avoided but still consumed.
# https://doi.org/10.1016/j.mambio.2013.08.007
## lit review of scats and abundance: 40 and 60 kg weight range, sambar preferred most. 
guilds$dhole_pref = "No"
guilds$dhole_pref[guilds$AdultBodyMass_g > 40000 & guilds$AdultBodyMass_g < 190000 &
                    guilds$TrophicLevel != "Carnivore"] = "Yes"
# who is preferred?
guilds$scientificNameStd[guilds$dhole_pref == "Yes"] # Some species here dont make sense: bears --> remove! Also some dont spatially overlap (orangutan)
# remove bears + orangutans, no evidence for this! 
guilds$dhole_pref[guilds$scientificNameStd %in% c("Helarctos_malayanus", "Ursus_thibetanus", "Pongo_pygmaeus")] = "No"
## unsure if muntjac should be added here... adding for now because they are consumed
guilds$dhole_pref[guilds$scientificNameStd %in% c("Muntiacus_genus")] = "Yes"
# who is preferred?
guilds$scientificNameStd[guilds$dhole_pref == "Yes"] # 5 species now... just ok. 


#### Clouded leopards
# https://onlinelibrary.wiley.com/doi/full/10.1002/ece3.9067 
# Serow! Ungulates in particular. No evidence for primates in this citation,
# but evidence from elsewhere for primates, and generally evidence for a varied and wide diet. 
guilds$CL_pref = "No"
guilds$CL_pref[guilds$AdultBodyMass_g > 7000 & guilds$AdultBodyMass_g < 190000 &
                           guilds$TrophicGuild != "Large_Carnivore"] = "Yes"
# who is preferred?
guilds$scientificNameStd[guilds$CL_pref == "Yes"] # Remove some of the larger carnivores (e.g. bears, cats), but leaving binturong b/c evidence exists! 
guilds$CL_pref[guilds$scientificNameStd %in% c("Catopuma_temminckii", "Helarctos_malayanus", "Ursus_thibetanus")] = "No"
## this has 14 species, many more than other predators, but little info exists for them! 


######### Visualize coefficient effect sizes ###### 

### Will be making two kinds of graphs:
# 1 comparing the effect sizes for each state variable within each species pair
# 1 histogram comparing the species interaction parameter for many species pairs per guild pair 

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
  dat = dat[dat$var %in% c("Abundance_intercept","FLII","HFP","Elevation","Hunting",
                           "Species_Interaction","Detection_intercept","Active_cams"),]
  
  ## replace _ with space for vars for clean x-vars
  dat$var = gsub("_", " ", dat$var)
  
  ## Change active cams to effort for easier understanding
  dat$var[dat$var == "Active cams"] = "Effort"
  
  ## order the factor levels for each variable
  dat$var = factor(dat$var, levels = c("Species Interaction","FLII","HFP", "Hunting", # state vars
                                       "Elevation", "Abundance intercept", # more state vars
                                       "Detection intercept","Effort")) # then det vars. 
  
  # Add a column for significance marker
  dat$marker <- ifelse(dat$sig == "Significant", "*", "")
  
  # Find the position of "Abundance intercept" on the x-axis to add the dashed line. 
  split_position <- which(levels(dat$var) == "Abundance intercept")
  
  # specify the distance were spacing error bars and stars
  dodge_width = 0.8
  
  # Extract the prefix from species to determine sub or dom
  dat$prefix <- sub("^(SUB|DOM).*", "\\1", dat$species)
  
  # specify the colors for dominant and subordinate species based on prefix 
  dat$color[dat$prefix == "SUB"] = "mediumslateblue"
  dat$color[dat$prefix == "DOM"] = "mediumseagreen"
  
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


######### Visualize species interaction coefficient via histograms #######

##### First, examine at the guild-level 

# first we need to add guild pairs to coeff (probs should have been done earlier)
add = distinct(select(preform, Species_Pair, guild_pair))
coeff = merge(coeff, add, by = "Species_Pair")
rm(add)

## add a column if interaction is significant or not and negative or positve
coeff$neg_vs_pos = "Positive & Non-Significant"
coeff$neg_vs_pos[coeff$mean < 0 & coeff$sig == "Significant"] = "Negative & Significant"
coeff$neg_vs_pos[coeff$mean < 0 & coeff$sig == "Non-Significant"] = "Negative & Non-Significant"
coeff$neg_vs_pos[coeff$mean > 0 & coeff$sig == "Significant"] = "Positive & Significant"

# specify colors manually
color_values_full <- c("Positive & Non-Significant" = "darkolivegreen1",
                       "Positive & Significant" = "darkgreen",
                       "Negative & Non-Significant" = "goldenrod1",
                       "Negative & Significant" = "firebrick4")

## Create a new sig level that combines non-sig into one
coeff$neg_vs_pos2 = coeff$neg_vs_pos
coeff$neg_vs_pos2[grepl("Non-Sig", coeff$neg_vs_pos2)] = "Non-Significant"
coeff$neg_vs_pos2[grepl("Positive", coeff$neg_vs_pos2)] = "Positive"
coeff$neg_vs_pos2[grepl("Negative", coeff$neg_vs_pos2)] = "Negative" # clean up names

# specify colors manually
color_values_short <- c("Positive" = "darkgreen",
                        "Negative" = "firebrick4",
                        "Non-Significant" = "gray75")

#### Make a test plot 

# subset the data
dat = coeff[coeff$var ==  "Species_Interaction",]

# Calculate the median value
median_value <- median(dat$mean)


# test the plot 
p =
ggplot(dat, aes(x = mean))+
  geom_histogram(aes(y = after_stat(density), fill = neg_vs_pos), bins = 100)+
  geom_density(aes(fill = neg_vs_pos), alpha = .75) +  # Add density plot with opacity
  scale_fill_manual(values = color_values_full)+
  geom_vline(aes(xintercept = 0), linetype = "dashed")+
  geom_vline(aes(xintercept = median(mean)), color = "purple") +
  annotate("text", x = Inf, y = Inf, label = paste("Median =", round(median(dat$mean), 2)),
           vjust = 10, hjust = 1.5, color = "black", size = 6) +
  theme_test()+
  labs(x = "Mean Species Interaction Value", y = "Density", fill = NULL, title = "All species pairs")

## save it!
# ggsave("figures/Grouped Histograms/Histogram_species_interaction_value_ALL_species_pairs_20230720.png", p,
#        width = 12, height = 8, units = "in")
rm(median_value, y_position, p, pos)

### Loop through all guild pairs, and save directly 
for(i in 1:length(unique(coeff$guild_pair))){
  
  ## save the name of one guild pair
  gp = unique(coeff$guild_pair)[i]
  
  ## subset data for 1 pair and the relevant data 
  dat = coeff[coeff$guild_pair == gp & coeff$var == "Species_Interaction",]
  
  # Calculate the median value
  median_value <- median(dat$mean)
  
  # # count how many values we have that will be graphed
  # pos = ddply(dat, .(neg_vs_pos), summarize,
  #             count = length(neg_vs_pos))
  # # make the label just below the largest count. 
  # y_position <- 0.75 * max(pos$count)
  ### THIS^^ is only needed if looking at raw counts, but using density now. 
  
  # make the plot 
  p =
    ggplot(dat, aes(x = mean))+
    geom_histogram(aes(y = after_stat(density), fill = neg_vs_pos), bins = 100)+
    geom_density(aes(fill = neg_vs_pos), alpha = .75) +  # Add density plot with opacity
    scale_fill_manual(values = color_values_full)+
    geom_vline(aes(xintercept = 0), linetype = "dashed")+
    geom_vline(aes(xintercept = median(mean)), color = "purple") +
    annotate("text", x = Inf, y = Inf, label = paste("Median =", round(median(dat$mean), 2)),
             vjust = 10, hjust = 1.5, color = "black", size = 6) +
    theme_test()+
    labs(x = "Mean Species Interaction Value", y = "Density", fill = NULL, title = gp)
  
  ## create a relevant file directory to save it
  # grab today's date
  day<-str_sub(Sys.Date(),-2)
  month<-str_sub(Sys.Date(),-5,-4)
  year<-str_sub(Sys.Date(),-10,-7)
  date = paste(year,month,day, sep = "")
  
  # construct a saving path 
  path = paste("figures/Grouped Histograms/Histogram_species_interaction_values_", gp, "_", date, ".png", sep = "")
  
  ## save it 
  ggsave(path, p, width = 12, height = 8, units = "in")
  
}
rm(p,i,gp,pos, y_position, dat, median_value, 
   day, month, year, date,path)


### Now make two more histos --> one w/ sub large carn and dom large carn
## Pt1 --> dominant
# subset the data
dat = coeff[coeff$var ==  "Species_Interaction" & grepl("DOM-Large_Carnivore", coeff$guild_pair),]

# # count how many values we have that will be graphed
# pos = ddply(dat, .(neg_vs_pos), summarize,
#             count = length(neg_vs_pos))
# # make the label just below the largest count. 
# y_position <- 0.5 * max(pos$count)

# make the plot 
p =
  ggplot(dat, aes(x = mean))+
  geom_histogram(aes(y = after_stat(density), fill = neg_vs_pos), bins = 100)+
  geom_density(aes(fill = neg_vs_pos), alpha = .75) +  # Add density plot with opacity
  scale_fill_manual(values = color_values_full)+
  geom_vline(aes(xintercept = 0), linetype = "dashed")+
  geom_vline(aes(xintercept = median(mean)), color = "purple") +
  annotate("text", x = Inf, y = Inf, label = paste("Median =", round(median(dat$mean), 2)),
           vjust = 10, hjust = 1.5, color = "black", size = 6) +
  theme_test()+
  labs(x = "Mean Species Interaction Value", y = "Density", fill = NULL, title = "Large carnivores are dominant")
# # save it!
# ggsave("figures/Grouped Histograms/Histogram_species_interaction_value_large_carnivores_dominant_20230720.png", p,
#        width = 12, height = 8, units = "in")


## Do the same, but exclude dogs for exploratory purposes. 
dat = dat[!grepl("Canis_lupus", dat$dom_sp),]
p =
  ggplot(dat, aes(x = mean))+
  geom_histogram(aes(y = after_stat(density), fill = neg_vs_pos), bins = 100)+
  geom_density(aes(fill = neg_vs_pos), alpha = .75) +  # Add density plot with opacity
  scale_fill_manual(values = color_values_full)+
  geom_vline(aes(xintercept = 0), linetype = "dashed")+
  geom_vline(aes(xintercept = median(mean)), color = "purple") +
  annotate("text", x = Inf, y = Inf, label = paste("Median =", round(median(dat$mean), 2)),
           vjust = 10, hjust = 1.5, color = "black", size = 6) +
  theme_test()+
  labs(x = "Mean Species Interaction Value", y = "Density", fill = NULL, title = "Large carnivores are dominant but no dogs")
# # save it!
# ggsave("figures/Grouped Histograms/Histogram_species_interaction_value_large_carnivores_dominant_NO_DOGS_20230720.png", p,
#        width = 12, height = 8, units = "in")


## Pt2 --> subordinate
# subset the data
dat = coeff[coeff$var ==  "Species_Interaction" & grepl("SUB-Large_Carnivore", coeff$guild_pair),]

# # Calculate the median value
# median_value <- median(dat$mean)
# 
# # count how many values we have that will be graphed
# pos = ddply(dat, .(neg_vs_pos), summarize,
#             count = length(neg_vs_pos))
# # make the label just below the largest count. 
# y_position <- 0.75 * max(pos$count)

# test the plot 
p = 
  ggplot(dat, aes(x = mean))+
  geom_histogram(aes(y = after_stat(density), fill = neg_vs_pos), bins = 100)+
  geom_density(aes(fill = neg_vs_pos), alpha = .75) +  # Add density plot with opacity
  scale_fill_manual(values = color_values_full)+
  geom_vline(aes(xintercept = 0), linetype = "dashed")+
  geom_vline(aes(xintercept = median(mean)), color = "purple") +
  annotate("text", x = Inf, y = Inf, label = paste("Median =", round(median(dat$mean), 2)),
           vjust = 10, hjust = 1.5, color = "black", size = 6) +
  theme_test()+
  labs(x = "Mean Species Interaction Value", y = "Density", fill = NULL, title = "Large carnivores are subordinate")
# # save it!
# ggsave("figures/Grouped Histograms/Histogram_species_interaction_value_large_carnivores_subordinate_20230720.png", p,
#        width = 12, height = 8, units = "in")

## Do the same, but exclude dogs for exploratory purposes. 
dat = dat[!grepl("Canis_lupus", dat$sub_sp),]
p = 
  ggplot(dat, aes(x = mean))+
  geom_histogram(aes(y = after_stat(density), fill = neg_vs_pos), bins = 100)+
  geom_density(aes(fill = neg_vs_pos), alpha = .75) +  # Add density plot with opacity
  scale_fill_manual(values = color_values_full)+
  geom_vline(aes(xintercept = 0), linetype = "dashed")+
  geom_vline(aes(xintercept = median(mean)), color = "purple") +
  annotate("text", x = Inf, y = Inf, label = paste("Median =", round(median(dat$mean), 2)),
           vjust = 10, hjust = 1.5, color = "black", size = 6) +
  theme_test()+
  labs(x = "Mean Species Interaction Value", y = "Density", fill = NULL, title = "Large carnivores are subordinate but no dogs")
# ggsave("figures/Grouped Histograms/Histogram_species_interaction_value_large_carnivores_subordinate_NO_DOGS_20230720.png", p,
#        width = 12, height = 8, units = "in")


## keep it clean
rm(p, y_position, pos, median_value, dat, i, day, month, year, date, path,gp)

#
##
###
#### Next, examine at the species-level 
##### Generate Ridgeline histograms per large carnivore

## The goal here is to condense several of those histograms into a single figure. 

## first subset coeff to just look at the species interaction value 
dat = coeff[coeff$var == "Species_Interaction",]

# ## remove dom from species name for a cleaner look 
dat$dom_sp = gsub("DOM-", "", dat$dom_sp)

## create two data subsets, 
a = dat[grepl("DOM-Large_Carn", dat$guild_pair),] # one with all large carnivores as dominant
b = a
unique(b$dom_sp) # and a second where all of these are combined into one 
b$dom_sp = "All_large_carnivores"

## and combine
dat = rbind(a,b)

## Set factor level based off weights from guild dataframe 
order = arrange(guilds[guilds$TrophicGuild == "Large_Carnivore",], desc(AdultBodyMass_g))
order$scientificNameStd
## do it manually b/c order is silly on graph
dat$dom_sp = factor(dat$dom_sp, levels = c("Neofelis_genus","Cuon_alpinus","Canis_lupus_familiaris", 
                                           "Panthera_pardus","Panthera_tigris","All_large_carnivores"))
rm(order)

## dominant large carnivore histogram
p = 
  ggplot(dat, aes(x = mean, y = dom_sp, fill = neg_vs_pos2)) +
  geom_density_ridges(scale = 1, alpha = .8) +
  geom_vline(aes(xintercept = 0), color = "black") +
  theme_ridges() +
  scale_fill_manual(values = color_values_short) +
  labs(x = "Mean Species Interaction Value", y = NULL, fill = NULL, title = "Large carnivores are dominant")
## save it! 
# ggsave("figures/Grouped Histograms/Histogram_Ridgeline_species_interaction_value_large_carnivores_dominant_20230720.png",p,
#        width = 10, height = 8, units = "in")

## make the same plot but without dogs 
## first subset coeff to just look at the species interaction value 
dat = coeff[coeff$var == "Species_Interaction",]

## create two data subsets, 
a = dat[grepl("DOM-Large_Carn", dat$guild_pair),] # one with all large carnivores as dominant
## but remove dogs here !
a = a[a$dom_sp != "Canis_lupus_familiaris",]
b = a
unique(b$dom_sp) # and a second where all of these are combined into one 
b$dom_sp = "All_large_carnivores"
# combine 
dat = rbind(a,b)
# and set factor levels
dat$dom_sp = factor(dat$dom_sp, levels = c("Neofelis_genus","Cuon_alpinus",
                                           "Panthera_pardus","Panthera_tigris","All_large_carnivores"))

## make the plot
p = 
  ggplot(dat, aes(x = mean, y = dom_sp, fill = neg_vs_pos2)) +
  geom_density_ridges(scale = 1, alpha = .8) +
  geom_vline(aes(xintercept = 0), color = "black") +
  theme_ridges() +
  scale_fill_manual(values = color_values_short) +
  labs(x = "Mean Species Interaction Value", y = NULL, fill = NULL, title = "Large carnivores are dominant but no dogs")
# ggsave("figures/Grouped Histograms/Histogram_Ridgeline_species_interaction_value_large_carnivores_dominant_NO_DOGS_20230720.png",p,
#        width = 10, height = 8, units = "in")


#### Remake this graph, but ONLY w/ preferred prey species
dat = coeff[coeff$var == "Species_Interaction",]
# First, gather and combine relevant data 
a = dat[dat$dom_sp == "Neofelis_genus" & 
          dat$sub_sp %in% guilds$scientificNameStd[guilds$CL_pref == "Yes"],]
b = dat[dat$dom_sp == "Panthera_tigris" & 
          dat$sub_sp %in% guilds$scientificNameStd[guilds$tiger_pref == "Yes"],]
c = dat[dat$dom_sp == "Panthera_pardus" & 
          dat$sub_sp %in% guilds$scientificNameStd[guilds$leopard_pref == "Yes"],]
d = dat[dat$dom_sp == "Canis_lupus_familiaris" & 
          dat$sub_sp %in% guilds$scientificNameStd[guilds$CL_pref == "Yes"],] # no accurate dog preferences, going w/ most broad group
e = dat[dat$dom_sp == "Cuon_alpinus" & 
          dat$sub_sp %in% guilds$scientificNameStd[guilds$dhole_pref == "Yes"],]
## also adding all large carnivores to keep sizing consistent
f = dat[grepl("DOM-Large_Carn", dat$guild_pair),]
f$dom_sp = "All_large_carnivores"

ridge_dat = rbind(a,b,c,d,e,f)
rm(a,b,c,d,e,f)

## set factor levels 
ridge_dat$dom_sp = factor(ridge_dat$dom_sp, levels = c("Neofelis_genus","Cuon_alpinus","Canis_lupus_familiaris", 
                                                       "Panthera_pardus","Panthera_tigris","All_large_carnivores"))
unique(ridge_dat$dom_sp) # good, want 6 spots 

# make the plot! 
p =
  ggplot(ridge_dat, aes(x = mean, y = dom_sp, fill = neg_vs_pos2)) +
  geom_density_ridges(scale = 1, alpha = .8) +
  geom_vline(aes(xintercept = 0), color = "black") +
  theme_ridges() +
  coord_cartesian(xlim = c(-7,2))+ # to remove all large carn fluff
  scale_fill_manual(values = color_values_short) +
  labs(x = "Mean Species Interaction Value", y = NULL, fill = NULL, title = "Large carnivores are dominant upon preferred prey")
# ggsave("figures/Grouped Histograms/Histogram_Ridgeline_species_interaction_value_large_carnivores_dominant_w_preferred_prey_20230721.png",p,
#        width = 10, height = 8, units = "in")

### now make a second one w/out dogs
p =
  ggplot(ridge_dat[ridge_dat$dom_sp != "Canis_lupus_familiaris",], aes(x = mean, y = dom_sp, fill = neg_vs_pos2)) +
  geom_density_ridges(scale = 1, alpha = .8) +
  geom_vline(aes(xintercept = 0), color = "black") +
  theme_ridges() +
  coord_cartesian(xlim = c(-7,2))+ # to remove all large carn fluff
  scale_fill_manual(values = color_values_short) +
  labs(x = "Mean Species Interaction Value", y = NULL, fill = NULL, title = "Large carnivores are dominant upon preferred prey")
# ggsave("figures/Grouped Histograms/Histogram_Ridgeline_species_interaction_value_large_carnivores_dominant_w_preferred_prey_and_NO_DOGS_20230721.png",p,
#        width = 10, height = 8, units = "in")


### Which preferred prey species are significant suppressed by tigers?
dat = coeff[coeff$var == "Species_Interaction",]
dat[dat$mean < 0 &
      dat$dom_sp == "Panthera_tigris" &
      dat$sig == "Significant" &
      dat$sub_sp %in% guilds$scientificNameStd[guilds$tiger_pref == "Yes"],]
## now what about leopards 
dat[dat$mean < 0 &
      dat$dom_sp == "Panthera_pardus" &
      dat$sig == "Significant" &
      dat$sub_sp %in% guilds$scientificNameStd[guilds$leopard_pref == "Yes"],] # macaques! 



#
##
### The graphs looks WAY better when we exclude models > -4, but were excluding 25 mods... 11 mods after removing small critters. 
length(unique(dat$Species_Pair[grepl("DOM-Large_Carn", dat$guild_pair) & dat$mean < -4]))

## inspect 
check = dat[dat$mean < -4 & dat$sig == "Significant",]
summary(check$sd) # very large SD! Not good! 
summary(check$Rhat) # not bad rhats tho 
table(preform$mod_valid[preform$Species_Pair %in% check$Species_Pair]) # 6 of these mods w/ huge estimates are valid. 

preform[preform$Species_Pair %in% preform$Species_Pair[preform$Species_Pair %in% check$Species_Pair &
                                                         preform$mod_valid == "Yes"],]
## These almost all involve the emerald dove! but there are other critters here too. 
rm(check)



####
###
##
# Do the same but reverse it! 

## first subset coeff to just look at the species interaction value 
dat = coeff[coeff$var == "Species_Interaction",]

## create two data subsets, 
a = dat[grepl("SUB-Large_Carn", dat$guild_pair),] # one with all large carnivores as dominant
b = a
unique(b$sub_sp) # and a second where all of these are combined into one 
b$sub_sp = "All_large_carnivores"

## and combine
dat = rbind(a,b)

## Set factor level to match other graphs 
dat$sub_sp = factor(dat$sub_sp, levels = c("Neofelis_genus","Cuon_alpinus","Canis_lupus_familiaris", 
                                           "Panthera_pardus","Panthera_tigris","All_large_carnivores"))

p =
  ggplot(dat, aes(x = mean, y = sub_sp, fill = neg_vs_pos2)) +
    geom_density_ridges(scale = 1, alpha = .8) +
    theme_ridges() +
    geom_vline(aes(xintercept = 0), color = "black") +
    scale_fill_manual(values = color_values_short) +
    labs(x = "Mean Species Interaction Value", y = NULL, fill = NULL, title = "Large carnivores are subordinate")
# ggsave("figures/Grouped Histograms/Histogram_Ridgeline_species_interaction_value_large_carnivores_subordinate_20230720.png",p,
#        width = 10, height = 8, units = "in")

##### do the same w/out dogs
## first subset coeff to just look at the species interaction value 
dat = coeff[coeff$var == "Species_Interaction",]

## create two data subsets, 
a = dat[grepl("SUB-Large_Carn", dat$guild_pair),] # one with all large carnivores as dominant
# but remove dogs!
a = a[a$sub_sp != "Canis_lupus_familiaris",]
b = a
unique(b$sub_sp) # and a second where all of these are combined into one 
b$sub_sp = "All_large_carnivores"

## and combine
dat = rbind(a,b)

## Set factor level to match other graphs 
dat$sub_sp = factor(dat$sub_sp, levels = c("Neofelis_genus","Cuon_alpinus", 
                                           "Panthera_pardus","Panthera_tigris","All_large_carnivores"))
## make the plot!
p =
  ggplot(dat, aes(x = mean, y = sub_sp, fill = neg_vs_pos2)) +
  geom_density_ridges(scale = 1, alpha = .8) +
  theme_ridges() +
  geom_vline(aes(xintercept = 0), color = "black") +
  scale_fill_manual(values = color_values_short) +
  labs(x = "Mean Species Interaction Value", y = NULL, fill = NULL, title = "Large carnivores are subordinate")
# ggsave("figures/Grouped Histograms/Histogram_Ridgeline_species_interaction_value_large_carnivores_subordinate_NO_DOGS_20230720.png",p,
#        width = 10, height = 8, units = "in")

####
###
##
# now subset for preferred prey species

# First, start fresh
dat = coeff[coeff$var == "Species_Interaction",]
# then, gather and combine relevant data 
a = dat[dat$sub_sp == "Neofelis_genus" & 
          dat$dom_sp %in% guilds$scientificNameStd[guilds$CL_pref == "Yes"],]
b = dat[dat$sub_sp == "Panthera_tigris" & 
          dat$dom_sp %in% guilds$scientificNameStd[guilds$tiger_pref == "Yes"],]
c = dat[dat$sub_sp == "Panthera_pardus" & 
          dat$dom_sp %in% guilds$scientificNameStd[guilds$leopard_pref == "Yes"],]
d = dat[dat$sub_sp == "Canis_lupus_familiaris" & 
          dat$dom_sp %in% guilds$scientificNameStd[guilds$CL_pref == "Yes"],]# no accurat dog preferences, going w/ most broad group
e = dat[dat$sub_sp == "Cuon_alpinus" & 
          dat$dom_sp %in% guilds$scientificNameStd[guilds$dhole_pref == "Yes"],]
f = dat[grepl("SUB-Large_Carn", dat$guild_pair),]
f$sub_sp = "All_large_carnivores"

## combine
ridge_dat = rbind(a,b,c,d,e,f)
rm(a,b,c,d,e,f)

## Set factor level to match other graphs 
ridge_dat$sub_sp = factor(ridge_dat$sub_sp, levels = c("Neofelis_genus","Cuon_alpinus","Canis_lupus_familiaris", 
                                                       "Panthera_pardus","Panthera_tigris","All_large_carnivores"))


## Make the plots 
p=
  ggplot(ridge_dat, aes(x = mean, y = sub_sp, fill = neg_vs_pos2)) +
  geom_density_ridges(scale = 1, alpha = .8) +
  theme_ridges() +
  coord_cartesian(xlim = c(-4,2))+ # to remove all large carn fluff
  geom_vline(aes(xintercept = 0), color = "black") + 
  scale_fill_manual(values = color_values_short) +
  labs(x = "Mean Species Interaction Value", y = NULL, fill = NULL, title = "Large carnivores are subordinate to preferred prey")
# ggsave("figures/Grouped Histograms/Histogram_Ridgeline_species_interaction_value_large_carnivores_subordinate_with_preferred_prey_20230721.png",p,
#        width = 10, height = 8, units = "in")

## do the same w/out dogs
p =
  ggplot(ridge_dat[ridge_dat$sub_sp != "Canis_lupus_familiaris",], aes(x = mean, y = sub_sp, fill = neg_vs_pos2)) +
  geom_density_ridges(scale = 1, alpha = .8) +
  theme_ridges() +
  coord_cartesian(xlim = c(-4,2))+ # to remove all large carn fluff
  geom_vline(aes(xintercept = 0), color = "black") +
  scale_fill_manual(values = color_values_short) +
  labs(x = "Mean Species Interaction Value", y = NULL, fill = NULL, title = "Large carnivores are subordinate to preferred prey")
# ggsave("figures/Grouped Histograms/Histogram_Ridgeline_species_interaction_value_large_carnivores_subordinate_with_preferred_prey_NO_DOGS_20230721.png",p,
#        width = 10, height = 8, units = "in")

### Which preferred prey species are significant supporting tigers?
dat = coeff[coeff$var == "Species_Interaction",]
dat[dat$mean > 0 &
      dat$sub_sp == "Panthera_tigris" &
      dat$sig == "Significant" &
      dat$dom_sp %in% guilds$scientificNameStd[guilds$tiger_pref == "Yes"],]
## now what about leopards 
dat[dat$mean > 0 &
      dat$sub_sp == "Panthera_pardus" &
      dat$sig == "Significant" &
      dat$dom_sp %in% guilds$scientificNameStd[guilds$leopard_pref == "Yes"],] 
## now what about dholes 
dat[dat$mean > 0 &
      dat$sub_sp == "Cuon_alpinus" &
      dat$sig == "Significant" &
      dat$dom_sp %in% guilds$scientificNameStd[guilds$dhole_pref == "Yes"],] # macaques! 




rm(dat, ridge_dat,order, p, color_values_short, color_values_full, a,b)

#
######### Visualize PPC plots ########

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
                        "Chat" = c(dat$Chat.dom, dat$BPV.sub),
                        "position" = c("Dominant", "Subordinate"),
                        "species" = str_split(n, "~")[[1]])
  
  ## BPV plot
  bpv = 
  ggplot(plot_dat, aes(x = species, y = BPV, fill = position))+
    geom_col(color = "black", position = "dodge", color = 'black')+
    geom_hline(yintercept = 0.5)+
    geom_hline(yintercept = 0.75, linetype = "dashed", color = "red")+
    geom_hline(yintercept = 0.25, linetype = "dashed", color = "red")+
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
  
  ##Chat plot
  chat = 
    ggplot(plot_dat, aes(x = species, y = Chat, fill = position))+
    geom_col(color = "black", position = "dodge", color = 'black')+
    geom_hline(yintercept = 1)+
    geom_hline(yintercept = 1.1, linetype = "dashed", color = "red")+
    theme_test()+
    labs(y = "Magnitude of Over-Dispersion (C-hat)", x = NULL, fill = NULL)+ 
    coord_cartesian(ylim = c(0.9,1.2))+
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
rm(bpv, chat, i, path, day,month,year,date,plot_dat, n, dat)



######### Visualize prediction plots ########

## Should ideally subset availible mods to the valid_mods species pairs before running this! 

## Inspect data
head(est_abund) #estimated abundance per sampling unit
str(est_abund) # already numeric

head(pred_abund) #Predicted change in subordinate abundance as a function of dominant species abundance
str(pred_abund) # already numeric


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
  
  
  ## make the plot! 
  # Hashed out lines are different styles you can test out. 
  p = 
  ggplot(pred, aes(y = Predicted.N, x = var))+
    # geom_point(data = est, aes(y = Sub_abundance, x = Dom_abundance, color = Landscape), alpha = 0.5, inherit.aes = F)+
    geom_jitter(data = est, aes(y = Sub_abundance, x = Dom_abundance, color = Landscape), width = .5)+
    geom_line(color = "black", size = 2)+
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.1, color = "black", linetype = "dashed")+
    # geom_errorbar(data = est, aes(y = Sub_abundance, x = Dom_abundance,
    #                                  ymin = lower_sub, ymax = upper_sub), alpha = 0.2, inherit.aes = F)+
    # geom_errorbarh(data = est, aes(y = Sub_abundance, x = Dom_abundance,
    #                                   xmin = lower_dom, xmax = upper_dom), alpha = 0.2, inherit.aes = F)+
    labs(y = paste0(n_split[1], " Abundance"), 
         x = paste0(n_split[2], " Abundance"),
         title = title, color = NULL)+ 
    guides(color = "none") + # way too many landscapes to list them now. 
    theme_test()+
    coord_cartesian(xlim = c(0, max(est$Dom_abundance)),
                    ylim = c(0, max(est$Sub_abundance)))
  
  ## construct the saving path
  path = paste("figures/", paste(gp, n, sep = "/"), "/Prediction_plot_", n, "_", date, ".png", sep = "")
  
  ## save the dang thang 
  ggsave(path, p, width = 9, height = 6, units = "in")
 
 }
rm(path, p,i, title, n, n_split, gp, est, pred, day, month, year, date)

######## Summary statistics to describe overall interactions #####

#### Summarize results by guild to clearly present them-
sum = ddply(preform, .(guild_pair), summarize,
            num_sig_neg_int = sum(Interaction_Estimate < 0 & Significance == "Significant", na.rm = T), # number of significant negative interactions
            num_sig_pos_int = sum(Interaction_Estimate > 0 & Significance == "Significant", na.rm = T), # number of significant positive interactions
            num_nonsig_neg_int = sum(Interaction_Estimate < 0 & Significance == "Non-Significant", na.rm = T), # number of NON-significant negative interactions
            num_nonsig_pos_int = sum(Interaction_Estimate > 0 & Significance == "Non-Significant", na.rm = T), # number of NON-significant positive interactions
            num_combos = length(unique(Species_Pair)))     # number of different species pairs per guild 

## what is the percent of negative or positive interactions per guild?
sum$percent_sig_neg = sum$num_sig_neg_int / sum$num_combos * 100
sum$percent_sig_pos = sum$num_sig_pos_int / sum$num_combos * 100
sum$percent_nonsig_neg = sum$num_nonsig_neg_int / sum$num_combos * 100
sum$percent_nonsig_pos = sum$num_nonsig_pos_int / sum$num_combos * 100
sum
## Key results! for every guild pairing and significant interactions, 
## there is a higher percentage of positive interactions than negative interactions! 


## Subset preform for only large carnivores as dominant
# dat = preform[preform$DOM_guild == "Large_Carnivore",]
# unique(dat$dom_species) # good. 

## Preform failed in recent test, using coeff rn
dat = coeff[grepl("DOM-Large_Carnivore", coeff$guild_pair) & coeff$var == "Species_Interaction",]
# rename cols to match 
names(dat)[c(2, 10, 11)] = c("Interaction_Estimate", "Significance", "dom_species")


#### Summarize results by large carnivore to clearly present them
sum_dom_large_carn = ddply(dat, .(dom_species), summarize,
                           num_sig_neg_int = sum(Interaction_Estimate < 0 & Significance == "Significant", na.rm = T), # number of significant negative interactions
                           num_sig_pos_int = sum(Interaction_Estimate > 0 & Significance == "Significant", na.rm = T), # number of significant positive interactions
                           num_nonsig_neg_int = sum(Interaction_Estimate < 0 & Significance == "Non-Significant", na.rm = T), # number of NON-significant negative interactions
                           num_nonsig_pos_int = sum(Interaction_Estimate > 0 & Significance == "Non-Significant", na.rm = T), # number of NON-significant positive interactions
                           num_combos = length(unique(Species_Pair)))     # number of different species pairs per dom large carn

## what is the percent of negative or positive interactions per dominant large carnivore?
sum_dom_large_carn$percent_sig_neg = sum_dom_large_carn$num_sig_neg_int / sum_dom_large_carn$num_combos * 100
sum_dom_large_carn$percent_sig_pos = sum_dom_large_carn$num_sig_pos_int / sum_dom_large_carn$num_combos * 100
sum_dom_large_carn$percent_nonsig_neg = sum_dom_large_carn$num_nonsig_neg_int / sum_dom_large_carn$num_combos * 100
sum_dom_large_carn$percent_nonsig_pos = sum_dom_large_carn$num_nonsig_pos_int / sum_dom_large_carn$num_combos * 100
sum_dom_large_carn
## Very interesting! Canis lupus is the only species with more significant negative interactions than any other! 
## But now leopards have an equal number of significant negative and positive interactions.  

## which species are feral dogs negetivly affecting?
preform$sub_species[preform$dom_species == "Canis_lupus_familiaris" &
                      preform$Interaction_Estimate < 0 &
                      # preform$mod_valid == "Yes" & # There are no sub_species if checking for valid mods! 
                      preform$Significance == "Significant"]
## very interesting results here! Tigers are suppressed by feral dogs... 

## which guilds are negativly affected by feral dogs?
table(preform$SUB_guild[preform$dom_species == "Canis_lupus_familiaris" &
                          preform$Interaction_Estimate < 0 &
                          preform$Significance == "Significant"])
## All guilds! But particularly small omnivores. 

## inspect feral dog records
dog = bdata$`SUB-Amaurornis_genus~DOM-Canis_lupus_familiaris`
# convert matrix to a df 
dog = (dog$y.dom)
dog = apply(dog, 2, as.numeric)
dog = as.data.frame(dog)
# grab row names
dog$SU = rownames(bdata$`SUB-Amaurornis_genus~DOM-Canis_lupus_familiaris`$y.dom)
# get total counts
dog$sum_count = rowSums(dog[,1:20], na.rm = T)
dog$landscape = sub("_\\d+.*", "", dog$SU)

##summarize it all
dog_sum = ddply(dog, .(landscape), summarize,
                total_dog_count = sum(sum_count))
dog_sum
## There is something here! 

######## Summary statistics to describe overall dataset ######

### effort per survey?
og_meta$camera_end.date = as.Date(og_meta$camera_end.date, format = "%Y-%m-%d")
og_meta$camera_start.date = as.Date(og_meta$camera_start.date, format = "%Y-%m-%d")

eff = ddply(og_meta, .(survey_id), summarize,
            total_eff = sum(effort),
            dur = as.numeric(difftime(max(camera_end.date), min(camera_start.date), units = "days")))

summary(eff$dur) # mean = 80.25
sd(eff$dur) # sd = 30.16

rm(eff)

### how many detections of our study species 
nrow(og_captures[og_captures$Species %in% preform$sub_species,])


#### We can also make a study map (coneptually) w/ the original metadata

### Import NON-spatially re-sampled metadata for the map 
non_resamp_meta = read_csv("/Users/zachary_amir/Dropbox/CT capture histories database/Asian ECL raw CT data/Step4_output_pre-resampling/Clean_independent_metadata_20230610.csv")

## thin this to be a bit easier to work with 
non_resamp_meta = select(non_resamp_meta, deployment_id, cell_id_3km, survey_id, Landscape, country, Latitude, Longitude)

## calculate average coords per landsacpe to visualize 
avg_land = ddply(non_resamp_meta, .(Landscape), summarize,
                 avg_lat = mean(Latitude),
                 avg_long = mean(Longitude))

## first make a map for the whole study area

#Genereate colour scheme as factor
category <- "Landscape" #What should we use as our category to generate colours for now survey ID
colour <- "lightgreen"# Define a colour from the R options to base the colourscheme

## make it a factor now 
avg_land$Landscape <- as.factor(avg_land$Landscape)
col.cat <- wheel(colour, num = length(levels(avg_land$Landscape)))
avg_land$Cols <- col.cat[avg_land$Landscape]


## make the map! 
leaflet() %>%
  # Add base map data
  addProviderTiles(providers$Esri.WorldImagery, group="Satellite") %>%  
  addCircleMarkers(lng= avg_land$avg_long, #Add Longitude
                   lat= avg_land$avg_lat, #Add latitude
                   color= avg_land$Cols
                   ) %>%
  # Add a legend
  addLegend("topleft", colors = col.cat,
            labels = levels(avg_land$Landscape),
            title = category,
            labFormat = labelFormat(prefix = "$"),
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

### Now subset the data for just the Thai E forest complex
thai = non_resamp_meta[grepl("Thai_E_FC_", non_resamp_meta$Landscape),]

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




##### Visualize interaction across community via matrix ######

### THIS PLOT SUCKS, not imformative and hard to make
## leaving code here for now, but will probably delete later. 

### will be working with coeff
### want a square matrix w/ dominant as cols and subordinat as rows
### fill in squares based on interaction strength, significance, and direction. 

# isolate the relevant info
dat = coeff[coeff$var == "Species_Interaction",]

# create a new col for dominant species
dat$dom_species = gsub("DOM-", "", sapply(str_split(dat$Species_Pair, "~"), "[", 2))
sort(unique(dat$dom_species)) # 55 species!

## are we missing any species?
setdiff(guilds$scientificNameStd[guilds$TrophicGuild != "Small_Carnivore"], dat$dom_species) # only 1, a lil bird. 

## rename other species col to sub and remove SUB from name
names(dat)[8] = "sub_species"
dat$sub_species = gsub("SUB-", "", dat$sub_species)


## need to create a matrix where rownames are sub and col names are dom
int_mat = matrix(NA, nrow = length(unique(dat$sub_species)),
                 ncol = length(unique(dat$dom_species)),
                 dimnames = list(unique(dat$sub_species), 
                                 unique(dat$dom_species)))

## and create a second matrix to store significance info
sig_mat =  matrix(NA, nrow = length(unique(dat$sub_species)),
                  ncol = length(unique(dat$dom_species)),
                  dimnames = list(unique(dat$sub_species), unique(dat$dom_species)))

## fill it in, row by row
for(i in 1:nrow(int_mat)){
  
  # select data for a single subordinate
  d = dat[dat$sub_species == rownames(int_mat)[i],]
  sub = unique(d$sub_species)
  
  for(l in 1:length(d$dom_species)){
    
    # select a single dom species
    dom = d$dom_species[l]
    
    # fill in the matricies
    int_mat[sub,dom] = dat$mean[dat$sub_species == sub & dat$dom_species == dom]
    sig_mat[sub,dom] = dat$sig[dat$sub_species == sub & dat$dom_species == dom]
    
  } # end per dom
  
} # end per sub
rm(dom, sub, d, i,l)

##### This is code that attempts to order the matrix, but I cant figure it out now
#### Hash it out! 
{
# ## susbset preform to extract relevant row order
# ord = distinct(select(preform, dom_species, DOM_guild))
# 
# ## order the rows of the matrix based on guild order
# col_order = data.frame(
#   species = rownames(int_mat),
#   species_order = match(rownames(int_mat), ord$dom_species), 
#   stringsAsFactors = F)
# 
# # Sort the guild order data frame
# col_order <- col_order[order(col_order$species_order), ]
# # Remove any thing not in order --> dont want it on graph. 
# # col_order = col_order[!is.na(col_order$species_order),]
# 
# ## subset preform to extract relevant col order
# ord = distinct(select(preform, sub_species, SUB_guild))
# 
# ## order the cols of the matrix based on guild order
# row_order = data.frame(
#   species = colnames(int_mat),
#   species_order = match(colnames(int_mat), ord$sub_species), 
#   stringsAsFactors = F)
# 
# # Sort the guild order data frame
# row_order <- row_order[order(row_order$species_order), ]
# # Remove any thing not in order --> dont want it on graph. 
# # row_order = row_order[!is.na(row_order$species_order),]
# 
# 
# ## create a new mat that matches a logical order
# new_int_mat =  matrix(NA, nrow = length(col_order$species_order),
#                       ncol = length(row_order$species_order),
#                       dimnames = list(col_order$species,row_order$species))
# 
# ## fill in the matrix with the right values 
# for (i in 1:nrow(int_mat)) {
#   for (j in 1:ncol(int_mat)) {
#     ## select the proper cell in the new matrix 
#     new_int_mat[row_order$species[i], col_order$species[j]] <- 
#       ## and fill in based on species names in old matrix 
#       int_mat[row_order$species[i], col_order$species[j]]
#   }
# }
# 
# ## do the same thing with sig_mat
# new_sig_mat =  matrix(NA, nrow = nrow(sig_mat),
#                       ncol = ncol(sig_mat),
#                       dimnames = list(guild_order$species, 
#                                       guild_order$species))
# ## fill in the matrix with the right values 
# for (i in 1:nrow(sig_mat)) {
#   for (j in 1:ncol(sig_mat)) {
#     ## select the proper cell in the new matrix 
#     new_sig_mat[guild_order$species[i], guild_order$species[j]] <- 
#       ## and fill in based on species names in old matrix 
#       sig_mat[guild_order$species[i], guild_order$species[j]]
#   }
# }
}

### Melt 'em!
## I think this turns a matrix into a dataframe (?)
int = melt(int_mat)
names(int)[1:2] = c("sub_species", "dom_species")

sig = melt(sig_mat)
names(sig)[1:2] = c("sub_species", "dom_species")
# rm(int_mat, new_int_mat, sig_mat, new_sig_mat)

## Add species_pair to both DFs 
sig$Species_Pair = paste("SUB-", sig$sub_species[], "~",
                         "DOM-", sig$dom_species[], sep = "")
int$Species_Pair = paste("SUB-", int$sub_species[], "~",
                         "DOM-", int$dom_species[], sep = "")

## inspect
setdiff(preform$Species_Pair, int$Species_Pair) # preform has 10 pairs not included here --> probs faile dmods. 
setdiff(int$Species_Pair, preform$Species_Pair) # 2000 + missing options! 

int = merge(int, select(preform, Species_Pair, guild_pair, mod_completion,
                        SUB_guild, DOM_guild), by = "Species_Pair")
table(int$guild_pair)

## combine them  combine them 
names(sig)[3] = "Significance"
dat = merge(int, sig, by = c("sub_species", "dom_species", "Species_Pair"))


## thin to completed mods 
# dat = dat[dat$mod_completion == "completed",] # doesnt matter, uncompleted mods still show up 

## Create a new column for simple positive vs negative values
dat$direction = "insignificant"
dat$direction[dat$Significance == "Significant" & dat$value > 0] = "Positive"
dat$direction[dat$Significance == "Significant" & dat$value < 0] = "Negative"




## make a plot object
p=
  ggplot(dat, aes(x = sub_species,
                  y = dom_species, 
                  fill = direction)) +
  geom_tile() +
  # scale_y_discrete(limits=rev)+ # reverse y-axis if wanted 
  # coord_equal() + # to make it into squares instead of rectangles 
  theme_classic() + 
  geom_text(aes(x = sub_species, y= dom_species, 
                label = ifelse(Significance == "Significant", "*", "")), inherit.aes = F) +
  labs(y = "Subordinate species", x = "Dominant Species")+
  scale_fill_manual(values = c("insignificant" = "grey",
                               "Negative" = "red",
                               "Positive" = "green" ))+  # use this when looking at parameter direction
  # scale_fill_continuous(type = "viridis")+ # use this when looking at parameter value
  theme(axis.text.x = element_text(angle = 40, vjust = 1, hjust= .95),
        text = element_text(size=16))
### TBH, this looks like crap rn... 



## looks like a great start! 
# setwd("~/Dropbox/Zach PhD/Trophic release project/Trophic release FIGURES/")
# ggsave("Pre-liminary species interaction result matrix plot_20230221.png",p,
#        height = 10, width = 12, units = "in")


## examine all mods where large carnivores are dominant
tmp = dat[grepl("DOM-Large_Carnivore", dat$guild_pair),]

## add sub and dom guilds
for(i in 1:length(tmp$guild_pair)){
  
  ## gather the combo and split it 
  c = tmp$guild_pair[i]
  c = str_split(c, "~")
  
  ## add it to the DF
  tmp$sub_guild[i] = c[[1]][1]
  tmp$dom_guild[i] = c[[1]][2]
  
  ## clean it up
  tmp$sub_guild[i] = gsub("SUB-", "", tmp$sub_guild[i])
  tmp$dom_guild[i] = gsub("DOM-", "", tmp$dom_guild[i])
  
}
rm(i,c)

## change the order so the Y axis starts w/ herbivores
# tmp = tmp[order(tmp$sub_guild),]
tmp$sub_species = factor(tmp$sub_species, levels = c(unique(tmp$sub_species[tmp$sub_guild == "Large_Herbivore"]),
                                                     unique(tmp$sub_species[tmp$sub_guild == "Small_Herbivore"]),
                                                     unique(tmp$sub_species[tmp$sub_guild == "Large_Omnivore"]),
                                                     unique(tmp$sub_species[tmp$sub_guild == "Small_Omnivore"]),
                                                     unique(tmp$sub_species[tmp$sub_guild == "Large_Carnivore"]),
                                                     unique(tmp$sub_species[tmp$sub_guild == "Small_Carnivore"])))

ggplot(tmp, aes(y = sub_species, x = dom_species, fill = direction)) +
  geom_tile(height = 5) +
  geom_text(aes(y = sub_species, x= dom_species, 
                label = ifelse(Significance == "Significant", "*", "")), inherit.aes = F) +
  labs(y = "Subordinate species", x = "Dominant Species")+
  scale_fill_manual(values = c("insignificant" = "grey",
                               "Negative" = "red",
                               "Positive" = "green" ))+  # use this when looking at parameter direction
  theme(axis.text.x = element_text(angle = 40, vjust = 1, hjust= .95),
        text = element_text(size=12))
### I dont think this is accurate!! 
## Nor do I think its useful... 

rm(p, sig,sig_mat, int,int_mat, dat, tmp)







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
