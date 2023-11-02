##### Inspect and visualize completed co-abundance models
#### Data was prepared on R Script: Step 1_Generate co-abundance data bundles.Rmd
#### and ran on HPC on R Script: HPC_co-abundance_model.R
### Will convert this to Rmd when ready, but R script for exploring now
## Import dataframes, not complete models b/c theyre too heavy

# Zachary Amir, Z.Amir@uq.edu.au
# Last updated on Sep 4th, 2023

## start fresh
rm(list = ls())

## load libs 
library(tidyverse)
library(plyr)
library(fs) # for moving files around 
library(ggridges) # for ridgeline plots 


# library(reshape2) # for making checkerboard matrix 
# library(leaflet) # for creating concepts of study map
# library(colortools) # for nice colors in the map
# library(jagsUI) 
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
bdata = readRDS("data/send_to_HPC/Bundled_data_for_Bayes_co-abundance_mods_279_species_pairs_20230821.RDS")

# ## I renamed a few species to move from broad genera to specific species --> fix here to match old results
# names(bdata) = gsub("Hystrix_brachyura", "Hystrix_genus", names(bdata))
# names(bdata) = gsub("Herpestes_urva", "Herpestes_genus", names(bdata)) ## herpestes will be tricky because its just the genera level in old code, but split in two species in new code. 

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


#
##
###
#### If you have imported files from the HPC but they are in the wrong folder,
#### use this code to move them around 

## key word to look for-
key = "HALF_LONG_3Chain"

# List all files in the relevant folders
coef_files <- dir("results/LONG_testing_Sep_2023/coefficent_dataframes/", full.names = TRUE)
ppc_files <- dir("results/LONG_testing_Sep_2023/PPC_dataframes/", full.names = TRUE)

# extract relevant files
coef_files = coef_files[grepl(key, coef_files)]
ppc_files = ppc_files[grepl(key, ppc_files)]

# combine into a list for easier looping
files = list("coef" = coef_files, "ppc" = ppc_files)

## for both directories
for(i in 1:length(files)){
  
  # chose one directory
  dir = files[[i]]
  
  for(l in 1:length(dir)){
    
    # select one file 
    f = dir[[l]]
    
    # and move it to the relevant folder
    if(names(files)[i] == "coef"){
      
      # and move it! 
      file_move(f, "results/HALF_LONG_testing_Oct_2023/coefficent_dataframes/")
      
    }
    
    if(names(files)[i] == "ppc"){

      # move it! 
      file_move(f, "results/HALF_LONG_testing_Oct_2023/PPC_dataframes//")
    }
    
    
  } # end l 
  
} # end i 
rm(i,l,f,dir,files,key,coef_files,ppc_files)



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

## inspect effective sample size
summary(coeff$n.eff[coeff$var == "Species_Interaction"]) # not good! Getting better w/ mid settings tho 
## which models lack sufficent sample size for the species interaction parameter?
sort(unique(coeff$Species_Pair[coeff$n.eff < 100 & coeff$var == "Species_Interaction"])) #163 models! 
## vis n.eff
hist(coeff$n.eff[coeff$var == "Species_Interaction"]) # a lot of small ones! Not good! 

summary(coeff$Rhat[coeff$var == "Species_Interaction"]) # The Meidan is good! but some high values exist, 3rd quart is still too high. 
sort(unique(coeff$Species_Pair[coeff$Rhat > 1.2 & coeff$var == "Species_Interaction"])) #120 models! 

## vis Rhat
hist(coeff$Rhat[coeff$var == "Species_Interaction"]) # thankfully most are small, but some biggies exist! 


## split apart species pair into sub vs dom species
coeff <- within(coeff, {
  sub_sp <- sapply(strsplit(Species_Pair, "~"), function(x) x[1])
  dom_sp <- sapply(strsplit(Species_Pair, "~"), function(x) x[2])
})
coeff$dom_sp = gsub("DOM-", "", coeff$dom_sp)
coeff$sub_sp = gsub("SUB-", "", coeff$sub_sp)


### Inspect which species pairs are present/missing, 20230829, MIDDLE settings

setdiff(all_combos, coeff$Species_Pair) # only 2 missing models!! This is a huge improvement, but still check. 
which(all_combos == "SUB-Varanus_salvator~DOM-Neofelis_genus")    # Check OE file number 86
which(all_combos == "SUB-Paradoxurus_hermaphroditus~DOM-Panthera_tigris")  # Check OE file number 157

### Long-story-short, I found out they all had problem SUs from DVCA that contained 14+ active cams 
## --> removed all SUs w/ > 14 cameras and re-run analysis. 
## 2nd update --> removed all SUs w/ > 7 cameras and re-run analysis with MUCH better performance. 
## 3rd update --> Seems to be going well! only probelmatic models are excluded (i.e. took too long to run)
## 4th update --> 2 failed models w/ palm civet might be because too many active cams and outlier covariates --> accept failure! 


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
table(preform$mod_completion) # 0.5% of models have failed! Epic improvement. 

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
summary(ppc_values$Rhat) # looks decent! some big values tho! 
nrow(ppc_values[ppc_values$Rhat > 1.2,])/nrow(ppc_values) * 100 # 23% of mods have bad interaction param --> better after removing non-relevant species! 
# how about sample size?
summary(coeff$n.eff[coeff$var == "Species_Interaction"]) # this could certainly be improved! 
length(coeff$n.eff[coeff$var == "Species_Interaction" & 
                     coeff$n.eff < 100]) / length(coeff$n.eff[coeff$var == "Species_Interaction"]) * 100 
# 76% of mods have small (< 100) sample size on short settings
# 36% of mods have small (< 100) sample size on middle settings

# remove list now that were all stored in dfs. 
rm(ppc_res)

### Inspect which species pairs are present/missing, 20230724
setdiff(all_combos, ppc_plotdat$Species_Pair) #2 missing models, same as above


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


### Inspect which species pairs are present/missing,  20230724
setdiff(all_combos, pred_abund$Species_Pair) #2 missing models, same as above


######## Import guild data and dietary preferences of large carnivores #####

#load guild data
guilds = read.csv(("data/species_traits/clean_44_species_trait_data_20230821.csv"))
head(guilds)
str(guilds)
## looks good, dietary preferences are ready to go! 

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
head(preform[preform$mod_completion == "uncompleted",]) # all NA values here for BPV,Chat,Rhat,etc b/c the models never finished! 
str(preform) # looks good! 

## Add a column denoting if the BPV values converged
preform$BPV_valid = "No"
preform$BPV_valid[preform$BPV.dom >= 0.25 & preform$BPV.dom <= 0.75 &
                    preform$BPV.sub >= 0.25 & preform$BPV.sub <= 0.75] = "Yes"
table(preform$BPV_valid) # mostly invalid rn, but thats w/ short settings. 
## Still mostly invalid with middle settings. 76% are not a good fit. 

## add a column denoting if overdispersion remains
preform$OD_valid = "No"
preform$OD_valid[preform$Chat.dom >= 0.98 & preform$Chat.dom <= 1.1 &
                     preform$Chat.sub >= 0.98 & preform$Chat.sub <= 1.1] = "Yes"
table(preform$OD_valid) # mostly OD, but thats w/ short settings. TBH, better preformance than expected w/ low settings. 
## Middle settings are very good, 93% have no OD
## I bet this is because active_cams + RE for source are actually good for detection models... less reliance on ODRE

## add a column denoting if the species interaction parameter converged 
preform$parameter_valid = "No"
preform$parameter_valid[preform$Rhat >= 0.99 & preform$Rhat <= 1.2] = "Yes"
table(preform$parameter_valid) # majority is valid! Should be better w/ more reps. 
## Middle settings are good, 69% have good interaction parameters 

## Finally, use all of this information to denote if the overall model is valid and should be trusted
preform$mod_valid = "No"
preform$mod_valid[preform$BPV_valid == "Yes" &
                    preform$OD_valid == "Yes" &
                    preform$parameter_valid == "Yes"] = "Yes"
table(preform$mod_valid) # Vast majority are NOT valid, but this checks out b/c of short/middle settings. 
## middle settings have 20% of all mods as totally valid

## Save which species pairs are valid to examine data later
valid_mods = preform$Species_Pair[preform$mod_valid == "Yes"]

## Which species pairs involved large carnivores significantly regulating the abundance of other critters?
preform$Species_Pair[preform$Interaction_Estimate < 0 & 
                       preform$Significance == "Significant" &
                       preform$mod_valid == "Yes" &
                       grepl("DOM-Large_Carnivore", preform$guild_pair)] 
## Curious to see what happens w/ mods w/ higher MCMC settings. 



######### Visualize coefficient effect sizes ###### 

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

# first we need to add guild pairs to coeff
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

#### Make a plot with ALL models in it  

# subset the data
dat = coeff[coeff$var ==  "Species_Interaction",]

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

## Grab today's date
day<-str_sub(Sys.Date(),-2)
month<-str_sub(Sys.Date(),-5,-4)
year<-str_sub(Sys.Date(),-10,-7)
date = paste(year,month,day, sep = "")
rm(year,month,day)

## and save it!
# ggsave(paste("figures/Grouped Histograms/Histogram_species_interaction_value_ALL_species_pairs_",date,".png", sep = ""), p,
#        width = 12, height = 8, units = "in")
rm(p, date)

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
rm(p,i,gp, dat, median_value, 
   day, month, year, date,path)


### Now make two more histos --> one w/ sub large carn and dom large carn
## Pt1 --> dominant
# subset the data
dat = coeff[coeff$var ==  "Species_Interaction" & grepl("DOM-Large_Carnivore", coeff$guild_pair),]

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

## Grab today's date
day<-str_sub(Sys.Date(),-2)
month<-str_sub(Sys.Date(),-5,-4)
year<-str_sub(Sys.Date(),-10,-7)
date = paste(year,month,day, sep = "")
rm(year,month,day)

## and save it!
# ggsave(paste("figures/Grouped Histograms/Histogram_species_interaction_value_large_carnivores_dominant_", date, ".png", sep = ""),
#        p, width = 12, height = 8, units = "in")

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
# save it!
# ggsave(paste("figures/Grouped Histograms/Histogram_species_interaction_value_large_carnivores_dominant_NO_DOGS_", date, ".png", sep = ""),
#        p, width = 12, height = 8, units = "in")


## Pt2 --> subordinate
# subset the data
dat = coeff[coeff$var ==  "Species_Interaction" & grepl("SUB-Large_Carnivore", coeff$guild_pair),]

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
  labs(x = "Mean Species Interaction Value", y = "Density", fill = NULL, title = "Large carnivores are subordinate")
# # save it!
# ggsave(paste("figures/Grouped Histograms/Histogram_species_interaction_value_large_carnivores_subordinate_", date, ".png", sep = ""),
#        p, width = 12, height = 8, units = "in")

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
# Save it! 
# ggsave(paste("figures/Grouped Histograms/Histogram_species_interaction_value_large_carnivores_subordinate_NO_DOGS_", date, ".png"),
#        p, width = 12, height = 8, units = "in")


## keep it clean
rm(p, dat, date)

#
##
###
#### Next, examine at the species-level 
##### Generate Ridgeline histograms per large carnivore
### The goal here is to condense several of those histograms into a single figure. 

## first subset coeff to just look at the species interaction value 
ridge_dat = coeff[coeff$var == "Species_Interaction",]

## create two data subsets, 
a = ridge_dat[grepl("DOM-Large_Carn", ridge_dat$guild_pair),] # one with all large carnivores as dominant
b = a
unique(b$dom_sp) # and a second where all of these are combined into one 
b$dom_sp = "All_large_carnivores"

## and combine
ridge_dat = rbind(a,b)
rm(a,b)

## Set factor level based off weights from guild dataframe 
order = arrange(guilds[guilds$TrophicGuild == "Large_Carnivore",], desc(AdultBodyMass_g))
order$scientificNameStd
## do it manually b/c order is silly on graph
ridge_dat$dom_sp = factor(ridge_dat$dom_sp, levels = c("Neofelis_genus","Cuon_alpinus", 
                                                       "Panthera_pardus","Panthera_tigris",
                                                       "All_large_carnivores"))
rm(order)

## try to group into two levels to reduce few data points
ridge_dat$hyp_support = "Unsupported"
ridge_dat$hyp_support[ridge_dat$mean < 0 & ridge_dat$sig == "Significant"] = "Supported"

## make sure supported comes first! 
ridge_dat$hyp_support = factor(ridge_dat$hyp_support, levels = c("Unsupported", "Supported"))

## and specify the colors for those levels 
hyp_colors = c("Unsupported" = "gray68",
               "Supported" = "firebrick4")


###### UPDATE --> geom_density_ridges() will EXCLUDE data if there are <= 2 data points for it! 
#### Potential solution is to graph those points WITH the regular ridgeline, 
## but include them just as points, NOT in the full density distribution. 

##### Subset records that wont fit in density estimation (i.e. <= 2 values in neg_vs_pos2 OR hyp_support)

## first grab all dom_species and take note of which need subsetting. 
# small_dat = ddply(ridge_dat, .(dom_sp, neg_vs_pos2), summarize,
#                   num_levls = length((neg_vs_pos2)))
small_dat = ddply(ridge_dat, .(dom_sp, hyp_support), summarize,
                  num_levls = length((hyp_support)))
## then remove all that dont need subsetting
small_dat = small_dat[small_dat$num_levls <= 2,]

## merge small_dat and ridge_dat based on common columns dom_sp and neg_vs_pos2 to ensure no repeats are added
# ridge_points <- merge(ridge_dat[ridge_dat$dom_sp %in% small_dat$dom_sp,], small_dat, by = c("dom_sp", "neg_vs_pos2"))
ridge_points <- merge(ridge_dat[ridge_dat$dom_sp %in% small_dat$dom_sp,], small_dat, by = c("dom_sp", "hyp_support"))

## verify we have the proper number of records 
nrow(ridge_points) == sum(small_dat$num_levls) ## MUST BE TRUE 
rm(small_dat)

## Thin ridge_dat to include only valid mods! 
# ridge_dat = ridge_dat[ridge_dat$Species_Pair %in% valid_mods,]

## dominant large carnivore histogram
p = 
  ggplot(ridge_dat, aes(x = mean, y = dom_sp, fill = hyp_support)) +
  geom_density_ridges(scale=.9, alpha = .8, jittered_points = TRUE, point_alpha=.5)+
  geom_point(data = ridge_points, aes(x = mean, y = dom_sp, color = hyp_support), 
             size = 5, shape = "circle", alpha = .8, position = position_nudge(y = .08), inherit.aes = F)+
  geom_vline(aes(xintercept = 0), color = "black") +
  theme_ridges() +
  # coord_cartesian(xlim = c(min,max))+ # to remove all large carn fluff
  # scale_fill_manual(values = color_values_short) +
  # scale_color_manual(values = color_values_short) +
  scale_fill_manual(values = hyp_colors) +
  scale_color_manual(values = hyp_colors) + 
  labs(x = "Mean Species Interaction Value", y = NULL, fill = NULL, color = NULL, 
       title = "Large carnivores are dominant")

## Grab today's date
day<-str_sub(Sys.Date(),-2)
month<-str_sub(Sys.Date(),-5,-4)
year<-str_sub(Sys.Date(),-10,-7)
date = paste(year,month,day, sep = "")
rm(year,month,day)

## and save it!
# ggsave(paste("figures/Grouped Histograms/Histogram_Ridgeline_species_interaction_value_large_carnivores_dominant_2_levels_",date, ".png", sep = ""),
#        p, width = 10, height = 8, units = "in")

## Make a second plot scaled per number of observations-
p = 
  ggplot(ridge_dat, aes(x = mean, y = dom_sp, fill = hyp_support, height = after_stat(density)))+
  geom_density_ridges(stat = "density", scale=.9, alpha = .8, jittered_points = TRUE, point_alpha=.5)+
  geom_point(data = ridge_points, aes(x = mean, y = dom_sp, color = hyp_support), 
             size = 5, shape = "circle", alpha = .8, position = position_nudge(y = .08), inherit.aes = F)+
  geom_vline(aes(xintercept = 0), color = "black") +
  theme_ridges() +
  scale_fill_manual(values = hyp_colors) +
  scale_color_manual(values = hyp_colors) +
  labs(x = "Mean Species Interaction Value", y = NULL, fill = NULL, color = NULL, 
       title = "Large carnivores are dominant")
# Save it! 
# ggsave(paste("figures/Grouped Histograms/Histogram_Ridgeline_species_interaction_value_large_carnivores_dominant_2_levels_DENSITY_SCALED_", date, ".png", sep = ""),p,
#        width = 10, height = 8, units = "in")



### Hashed out code below to make same graph, but in three levels --> positive, negative, non-significant

# ##### Subset records that wont fit in density estimation (i.e. <= 2 values in neg_vs_pos2 OR hyp_support)
# 
# ## first grab all dom_species and take note of which need subsetting. 
# small_dat = ddply(ridge_dat, .(dom_sp, neg_vs_pos2), summarize,
#                   num_levls = length((neg_vs_pos2)))
# 
# ## then remove all that dont need subsetting
# small_dat = small_dat[small_dat$num_levls <= 2,]
# 
# ## merge small_dat and ridge_dat based on common columns dom_sp and neg_vs_pos2 to ensure no repeats are added
# ridge_points <- merge(ridge_dat[ridge_dat$dom_sp %in% small_dat$dom_sp,], small_dat, by = c("dom_sp", "neg_vs_pos2"))
# 
# ## verify we have the proper number of records 
# nrow(ridge_points) == sum(small_dat$num_levls) ## MUST BE TRUE 
# rm(small_dat)
# 
# ##Thin ridge_dat to include only valid mods!
# ridge_dat = ridge_dat[ridge_dat$Species_Pair %in% valid_mods,]
# 
# ## make the plot
# p = 
#   ggplot(ridge_dat, aes(x = mean, y = dom_sp, fill = neg_vs_pos2)) +
#   geom_density_ridges(scale=.9, alpha = .8, jittered_points = TRUE, point_alpha=.5)+
#   geom_point(data = ridge_points, aes(x = mean, y = dom_sp, color = neg_vs_pos2), 
#              size = 5, shape = "circle", alpha = .8, position = position_nudge(y = .08), inherit.aes = F)+
#   geom_vline(aes(xintercept = 0), color = "black") +
#   theme_ridges() +
#   scale_fill_manual(values = color_values_short) +
#   scale_color_manual(values = color_values_short) +
#   labs(x = "Mean Species Interaction Value", y = NULL, fill = NULL, color = NULL, 
#        title = "Large carnivores are dominant")
# ## Save it!
# # ggsave(paste("figures/Grouped Histograms/Histogram_Ridgeline_species_interaction_value_large_carnivores_dominant_3_levels_", date, ".png", sep = ""),
# #        p, width = 10, height = 8, units = "in")



#
##
###
#### Remake this graph, but ONLY w/ preferred prey species
dat = coeff[coeff$var == "Species_Interaction",]
# First, gather and combine relevant data 
a = dat[dat$dom_sp == "Neofelis_genus" & 
          dat$sub_sp %in% guilds$scientificNameStd[guilds$CL_pref == "Yes"],]
b = dat[dat$dom_sp == "Panthera_tigris" & 
          dat$sub_sp %in% guilds$scientificNameStd[guilds$tiger_pref == "Yes"],]
c = dat[dat$dom_sp == "Panthera_pardus" & 
          dat$sub_sp %in% guilds$scientificNameStd[guilds$leopard_pref == "Yes"],]
d = dat[dat$dom_sp == "Cuon_alpinus" & 
          dat$sub_sp %in% guilds$scientificNameStd[guilds$dhole_pref == "Yes"],]
## also adding all large carnivores w/ preferred prey
e = rbind(a,b,c,d)
e$dom_sp = "All_large_carnivores"

## combine
ridge_dat = rbind(a,b,c,d,e)
rm(a,b,c,d,e)

## set factor levels
ridge_dat$dom_sp = factor(ridge_dat$dom_sp, levels = c("Neofelis_genus","Cuon_alpinus",
                                                       "Panthera_pardus","Panthera_tigris",
                                                       "All_large_carnivores"))
unique(ridge_dat$dom_sp) # good, want 5 spots

## try to group into two levels to reduce few data points
ridge_dat$hyp_support = "Unsupported"
ridge_dat$hyp_support[ridge_dat$mean < 0 & ridge_dat$sig == "Significant"] = "Supported"

## make sure supported comes first! 
ridge_dat$hyp_support = factor(ridge_dat$hyp_support, levels = c("Unsupported", "Supported"))

## and specify the colors for those levels 
hyp_colors = c("Unsupported" = "gray68",
               "Supported" = "firebrick4")



### subset records that wont fit in density estimation (i.e. <= 2 values in hyp_support)

## first grab all dom_species and take note of which need subsetting. 
small_dat = ddply(ridge_dat, .(dom_sp, hyp_support), summarize,
                  num_levls = length((hyp_support)))
## then remove all that dont need subsetting
small_dat = small_dat[small_dat$num_levls <= 2,]
## merge small_dat and ridge_dat based on common columns dom_sp and hyp_support to ensure no repeats are added
ridge_points <- merge(ridge_dat[ridge_dat$dom_sp %in% small_dat$dom_sp,], small_dat, by = c("dom_sp", "hyp_support"))
## verify we have the proper number of records 
nrow(ridge_points) == sum(small_dat$num_levls) ## MUST BE TRUE 
rm(small_dat)

# ## Thin ridge_dat to include only valid mods! 
# ridge_dat = ridge_dat[ridge_dat$Species_Pair %in% valid_mods,]

# make the plot! 
p =
  ggplot(ridge_dat, aes(x = mean, y = dom_sp, fill = hyp_support)) +
  geom_density_ridges(scale=.9, alpha = .8, jittered_points = TRUE, point_alpha=.5)+
  geom_point(data = ridge_points, aes(x = mean, y = dom_sp, color = hyp_support), 
             size = 5, shape = "circle", alpha = .8, position = position_nudge(y = .08), inherit.aes = F)+
  geom_vline(aes(xintercept = 0), color = "black") +
  theme_ridges() +
  scale_fill_manual(values = hyp_colors) +
  scale_color_manual(values = hyp_colors) +
  labs(x = "Mean Species Interaction Value", y = NULL, fill = NULL, color = NULL, 
       title = "Large carnivores are dominant upon preferred prey")
# Save it! 
# ggsave(paste("figures/Grouped Histograms/Histogram_Ridgeline_species_interaction_value_large_carnivores_dominant_w_preferred_prey_2_levels_", date, ".png", sep = ""),
#        p, width = 10, height = 8, units = "in")

## Make a second plot scaled per number of observations-
p = 
  ggplot(ridge_dat, aes(x = mean, y = dom_sp, fill = hyp_support, height = stat(density)))+
  geom_density_ridges(stat = "density", scale=.9, alpha = .8, jittered_points = TRUE, point_alpha=.5)+
  geom_point(data = ridge_points, aes(x = mean, y = dom_sp, color = hyp_support), 
             size = 5, shape = "circle", alpha = .8, position = position_nudge(y = .08), inherit.aes = F)+
  geom_vline(aes(xintercept = 0), color = "black") +
  theme_ridges() +
  scale_fill_manual(values = hyp_colors) +
  scale_color_manual(values = hyp_colors) +
  labs(x = "Mean Species Interaction Value", y = NULL, fill = NULL, color = NULL, 
       title = "Large carnivores are dominant upon preferred prey")
# Save it! 
# ggsave(paste("figures/Grouped Histograms/Histogram_Ridgeline_species_interaction_value_large_carnivores_dominant_w_preferred_prey_2_levels_DENSITY_SCALED_", date, ".png", sep = ""),p,
#        width = 10, height = 8, units = "in")



#
##
### now make a second one w/ 3 levels

# ##### Subset records that wont fit in density estimation (i.e. <= 2 values in neg_vs_pos2 OR hyp_support)
# 
# ## first grab all dom_species and take note of which need subsetting. 
# small_dat = ddply(ridge_dat, .(dom_sp, neg_vs_pos2), summarize,
#                   num_levls = length((neg_vs_pos2)))
# 
# ## then remove all that dont need subsetting
# small_dat = small_dat[small_dat$num_levls <= 2,]
# 
# ## merge small_dat and ridge_dat based on common columns dom_sp and neg_vs_pos2 to ensure no repeats are added
# ridge_points <- merge(ridge_dat[ridge_dat$dom_sp %in% small_dat$dom_sp,], small_dat, by = c("dom_sp", "neg_vs_pos2"))
# 
# ## verify we have the proper number of records 
# nrow(ridge_points) == sum(small_dat$num_levls) ## MUST BE TRUE 
# rm(small_dat)
# 
## Thin ridge_dat to include only valid mods! 
# ridge_dat = ridge_dat[ridge_dat$Species_Pair %in% valid_mods,]
# 
# ## make the plot
# p = 
#   ggplot(ridge_dat, aes(x = mean, y = dom_sp, fill = neg_vs_pos2)) +
#   geom_density_ridges(scale=.9, alpha = .8, jittered_points = TRUE, point_alpha=.5)+
#   geom_point(data = ridge_points, aes(x = mean, y = dom_sp, color = neg_vs_pos2), 
#              size = 5, shape = "circle", alpha = .8, position = position_nudge(y = .08), inherit.aes = F)+
#   geom_vline(aes(xintercept = 0), color = "black") +
#   theme_ridges() +
#   scale_fill_manual(values = color_values_short) +
#   scale_color_manual(values = color_values_short) +
#   labs(x = "Mean Species Interaction Value", y = NULL, fill = NULL, color = NULL, 
#        title = "Large carnivores are dominant upon preferred prey")
# ## Save it!
# # ggsave(paste("figures/Grouped Histograms/Histogram_Ridgeline_species_interaction_value_large_carnivores_dominant_w_preferred_prey_3_levels_", date, ".png", sep = ""),
# #        p, width = 10, height = 8, units = "in")

### Which preferred prey species are significant suppressed by tigers?
dat = coeff[coeff$var == "Species_Interaction",]
dat[dat$mean < 0 &
      dat$dom_sp == "Panthera_tigris" &
      dat$sig == "Significant" &
      dat$sub_sp %in% guilds$scientificNameStd[guilds$tiger_pref == "Yes"],]
## Serow and Muntjac, classic! 

## now what about leopards 
dat[dat$mean < 0 &
      dat$dom_sp == "Panthera_pardus" &
      dat$sig == "Significant" &
      dat$sub_sp %in% guilds$scientificNameStd[guilds$leopard_pref == "Yes"],] # macaques! both pig and long tailed 

## and one should be suppresed by dholes, who?
dat[dat$mean < 0 &
      dat$dom_sp == "Cuon_alpinus" &
      dat$sig == "Significant" &
      dat$sub_sp %in% guilds$scientificNameStd[guilds$dhole_pref == "Yes"],] # serow! 
## UPDATE, 2023-08-29 --> this model failed to converge on middle settings, much lower n.eff


#
##
### The graphs looks WAY better when we exclude models > -4, but were excluding 25 mods... 11 mods after removing small critters. 
length(unique(dat$Species_Pair[grepl("DOM-Large_Carn", dat$guild_pair) & dat$mean < -4])) # only 6 mods now (2023-08-29)

## inspect 
check = dat[dat$mean < -4 & dat$sig == "Significant",]
summary(check$sd) # very large SD! Not good! 
summary(check$Rhat) # not bad rhats tho 
table(preform$mod_valid[preform$Species_Pair %in% check$Species_Pair]) # 2 of these mods w/ huge estimates are valid. 

preform[preform$Species_Pair %in% preform$Species_Pair[preform$Species_Pair %in% check$Species_Pair &
                                                         preform$mod_valid == "Yes"],]
## small carnivores subordinate to large carns... 
rm(check, min,max)



####
###
##
# Do the same but reverse it! 

## Grab a fresh date just to be safe
day<-str_sub(Sys.Date(),-2)
month<-str_sub(Sys.Date(),-5,-4)
year<-str_sub(Sys.Date(),-10,-7)
date = paste(year,month,day, sep = "")
rm(year,month,day)

## first subset coeff to just look at the species interaction value 
dat = coeff[coeff$var == "Species_Interaction",]

## create two data subsets, 
a = dat[grepl("SUB-Large_Carn", dat$guild_pair),] # one with all large carnivores as dominant
b = a
unique(b$sub_sp) # and a second where all of these are combined into one 
b$sub_sp = "All_large_carnivores"

## and combine
ridge_dat = rbind(a,b)

## Set factor level to match other graphs 
ridge_dat$sub_sp = factor(ridge_dat$sub_sp, levels = c("Neofelis_genus","Cuon_alpinus", 
                                                       "Panthera_pardus","Panthera_tigris",
                                                       "All_large_carnivores"))

#### try to group into two levels to reduce few data points
ridge_dat$hyp_support = "Unsupported"
ridge_dat$hyp_support[ridge_dat$mean > 0 & ridge_dat$sig == "Significant"] = "Supported"

## make sure supported comes first! 
ridge_dat$hyp_support = factor(ridge_dat$hyp_support, levels = c("Unsupported", "Supported"))

## and specify the colors for those levels 
hyp_colors = c("Unsupported" = "gray68",
               "Supported" = "darkgreen")


##### Subset records that wont fit in density estimation (i.e. <= 2 values in hyp_support)
## first grab all sub_species and take note of which need subsetting. 
small_dat = ddply(ridge_dat, .(sub_sp, hyp_support), summarize,
                  num_levls = length((hyp_support)))
## then remove all that dont need subsetting
small_dat = small_dat[small_dat$num_levls <= 2,]
## merge small_dat and ridge_dat based on common columns sub_sp and hyp_support to ensure no repeats are added
ridge_points <- merge(ridge_dat[ridge_dat$sub_sp %in% small_dat$sub_sp,], small_dat, by = c("sub_sp", "hyp_support"))
## verify we have the proper number of records 
nrow(ridge_points) == sum(small_dat$num_levls) ## MUST BE TRUE 
rm(small_dat)

## Thin ridge_dat to include only valid mods! 
# ridge_dat = ridge_dat[ridge_dat$Species_Pair %in% valid_mods,]

## make the plot 
p =
  ggplot(ridge_dat, aes(x = mean, y = sub_sp, fill = hyp_support)) +
  geom_density_ridges(scale=.9, alpha = .8, jittered_points = TRUE, point_alpha=.5)+
  geom_point(data = ridge_points, aes(x = mean, y = sub_sp, color = hyp_support), 
             size = 5, shape = "circle", alpha = .8, position = position_nudge(y = .08), inherit.aes = F)+
  geom_vline(aes(xintercept = 0), color = "black") +
  theme_ridges() +
  scale_fill_manual(values = hyp_colors) +
  scale_color_manual(values = hyp_colors) +
  labs(x = "Mean Species Interaction Value", y = NULL, fill = NULL, color = NULL, 
       title = "Large carnivores are subordinate")
# Save it! 
# ggsave(paste("figures/Grouped Histograms/Histogram_Ridgeline_species_interaction_value_large_carnivores_subordinate_2_levels_", date, ".png", sep = ""),p,
#        width = 10, height = 8, units = "in")

## Make a second plot scaled per number of observations-
p = 
  ggplot(ridge_dat, aes(x = mean, y = sub_sp, fill = hyp_support, height = stat(density)))+
  geom_density_ridges(stat = "density", scale=.9, alpha = .8, jittered_points = TRUE, point_alpha=.5)+
  geom_point(data = ridge_points, aes(x = mean, y = sub_sp, color = hyp_support), 
             size = 5, shape = "circle", alpha = .8, position = position_nudge(y = .08), inherit.aes = F)+
  geom_vline(aes(xintercept = 0), color = "black") +
  theme_ridges() +
  scale_fill_manual(values = hyp_colors) +
  scale_color_manual(values = hyp_colors) +
  labs(x = "Mean Species Interaction Value", y = NULL, fill = NULL, color = NULL, 
       title = "Large carnivores are subordinate")
# Save it! 
# ggsave(paste("figures/Grouped Histograms/Histogram_Ridgeline_species_interaction_value_large_carnivores_subordinate_2_levels_DENSITY_SCALED_", date, ".png", sep = ""),p,
#        width = 10, height = 8, units = "in")

##### do the same w/ 3 levels
# 
# ##### Subset records that wont fit in density estimation (i.e. <= 2 values in neg_vs_pos2 OR hyp_support)
# 
# ## first grab all sub_species and take note of which need subsetting. 
# small_dat = ddply(ridge_dat, .(sub_sp, neg_vs_pos2), summarize,
#                   num_levls = length((neg_vs_pos2)))
# ## then remove all that dont need subsetting
# small_dat = small_dat[small_dat$num_levls <= 2,]
# ## merge small_dat and ridge_dat based on common columns sub_sp and neg_vs_pos2 to ensure no repeats are added
# ridge_points <- merge(ridge_dat[ridge_dat$sub_sp %in% small_dat$sub_sp,], small_dat, by = c("sub_sp", "neg_vs_pos2"))
# ## verify we have the proper number of records 
# nrow(ridge_points) == sum(small_dat$num_levls) ## MUST BE TRUE 
# rm(small_dat)
# 
#
## Thin ridge_dat to include only valid mods! 
# ridge_dat = ridge_dat[ridge_dat$Species_Pair %in% valid_mods,]
#
# ## make the plot
# p = 
#   ggplot(ridge_dat, aes(x = mean, y = sub_sp, fill = neg_vs_pos2)) +
#   geom_density_ridges(scale=.9, alpha = .8, jittered_points = TRUE, point_alpha=.5)+
#   geom_point(data = ridge_points, aes(x = mean, y = sub_sp, color = neg_vs_pos2), 
#              size = 5, shape = "circle", alpha = .8, position = position_nudge(y = .08), inherit.aes = F)+
#   geom_vline(aes(xintercept = 0), color = "black") +
#   theme_ridges() +
#   scale_fill_manual(values = color_values_short) +
#   scale_color_manual(values = color_values_short) +
#   labs(x = "Mean Species Interaction Value", y = NULL, fill = NULL, color = NULL, 
#        title = "Large carnivores are subordinate")
# # Save it! 
# # ggsave(paste("figures/Grouped Histograms/Histogram_Ridgeline_species_interaction_value_large_carnivores_subordinate_3_levels_", date, ".png", sep = ""),p,
# #        width = 10, height = 8, units = "in")

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
d = dat[dat$sub_sp == "Cuon_alpinus" & 
          dat$dom_sp %in% guilds$scientificNameStd[guilds$dhole_pref == "Yes"],]
## combine all for all_large_carn
e = rbind(a,b,c,d)
e$sub_sp = "All_large_carnivores"

## combine
ridge_dat = rbind(a,b,c,d,e)
rm(a,b,c,d,e)

## Set factor level to match other graphs 
ridge_dat$sub_sp = factor(ridge_dat$sub_sp, levels = c("Neofelis_genus","Cuon_alpinus",
                                                       "Panthera_pardus","Panthera_tigris",
                                                       "All_large_carnivores"))
## try to group into two levels to reduce few data points
ridge_dat$hyp_support = "Unsupported"
ridge_dat$hyp_support[ridge_dat$mean > 0 & ridge_dat$sig == "Significant"] = "Supported"

## make sure supported comes first! 
ridge_dat$hyp_support = factor(ridge_dat$hyp_support, levels = c("Unsupported", "Supported"))


##### Subset records that wont fit in density estimation (i.e. <= 2 values in neg_vs_pos2 OR hyp_supported)
## first grab all dom_species and take note of which need subsetting. 
small_dat = ddply(ridge_dat, .(sub_sp, hyp_support), summarize,
                  num_levls = length((hyp_support)))
## then remove all that dont need subsetting
small_dat = small_dat[small_dat$num_levls <= 2,]
## merge small_dat and ridge_dat based on common columns sub_sp and hyp_support to ensure no repeats are added
ridge_points <- merge(ridge_dat[ridge_dat$sub_sp %in% small_dat$sub_sp,], small_dat, by = c("sub_sp", "hyp_support"))
## verify we have the proper number of records 
nrow(ridge_points) == sum(small_dat$num_levls) ## MUST BE TRUE 
rm(small_dat)

## Thin ridge_dat to include only valid mods! 
# ridge_dat = ridge_dat[ridge_dat$Species_Pair %in% valid_mods,]

## Make the plots 
p=
  ggplot(ridge_dat, aes(x = mean, y = sub_sp, fill = hyp_support)) +
  geom_density_ridges(scale=.9, alpha = .8, jittered_points = TRUE, point_alpha=.5)+
  geom_point(data = ridge_points, aes(x = mean, y = sub_sp, color = hyp_support), 
             size = 5, shape = "circle", alpha = .8, position = position_nudge(y = .08), inherit.aes = F)+
  geom_vline(aes(xintercept = 0), color = "black") +
  theme_ridges() +
  scale_fill_manual(values = hyp_colors) +
  scale_color_manual(values = hyp_colors) +
  labs(x = "Mean Species Interaction Value", y = NULL, fill = NULL, color = NULL, 
       title = "Large carnivores are subordinate to preferred prey")
# Save it! 
# ggsave(paste("figures/Grouped Histograms/Histogram_Ridgeline_species_interaction_value_large_carnivores_subordinate_with_preferred_prey_2_levels_", date, ".png", sep = ""),p,
#        width = 10, height = 8, units = "in")

## Make a second plot scaled per number of observations-
p = 
  ggplot(ridge_dat, aes(x = mean, y = sub_sp, fill = hyp_support, height = stat(density)))+
  geom_density_ridges(stat = "density", scale=.9, alpha = .8, jittered_points = TRUE, point_alpha=.5)+
  geom_point(data = ridge_points, aes(x = mean, y = sub_sp, color = hyp_support), 
             size = 5, shape = "circle", alpha = .8, position = position_nudge(y = .08), inherit.aes = F)+
  geom_vline(aes(xintercept = 0), color = "black") +
  theme_ridges() +
  scale_fill_manual(values = hyp_colors) +
  scale_color_manual(values = hyp_colors) +
  labs(x = "Mean Species Interaction Value", y = NULL, fill = NULL, color = NULL, 
       title = "Large carnivores are subordinate to preferred prey")
# Save it! 
# ggsave(paste("figures/Grouped Histograms/Histogram_Ridgeline_species_interaction_value_large_carnivores_subordinate_with_preferred_prey_2_levels_DENSITY_SCALED_", date, ".png", sep = ""),p,
#        width = 10, height = 8, units = "in")




## remake the plot w/ three levels 

# ##### Subset records that wont fit in density estimation (i.e. <= 2 values in neg_vs_pos2 OR hyp_support)
# 
# ## first grab all sub_species and take note of which need subsetting. 
# small_dat = ddply(ridge_dat, .(sub_sp, neg_vs_pos2), summarize,
#                   num_levls = length((neg_vs_pos2)))
# ## then remove all that dont need subsetting
# small_dat = small_dat[small_dat$num_levls <= 2,]
# ## merge small_dat and ridge_dat based on common columns sub_sp and neg_vs_pos2 to ensure no repeats are added
# ridge_points <- merge(ridge_dat[ridge_dat$sub_sp %in% small_dat$sub_sp,], small_dat, by = c("sub_sp", "neg_vs_pos2"))
# ## verify we have the proper number of records 
# nrow(ridge_points) == sum(small_dat$num_levls) ## MUST BE TRUE 
# rm(small_dat)
# 
# ## make the plot
# p = 
#   ggplot(ridge_dat, aes(x = mean, y = sub_sp, fill = neg_vs_pos2)) +
#   geom_density_ridges(scale=.9, alpha = .8, jittered_points = TRUE, point_alpha=.5)+
#   geom_point(data = ridge_points, aes(x = mean, y = sub_sp, color = neg_vs_pos2), 
#              size = 5, shape = "circle", alpha = .8, position = position_nudge(y = .08), inherit.aes = F)+
#   geom_vline(aes(xintercept = 0), color = "black") +
#   theme_ridges() +
#   scale_fill_manual(values = color_values_short) +
#   scale_color_manual(values = color_values_short) +
#   labs(x = "Mean Species Interaction Value", y = NULL, fill = NULL, color = NULL, 
#        title = "Large carnivores are subordinate to preferred prey")
# # Save it! 
# # ggsave(paste("figures/Grouped Histograms/Histogram_Ridgeline_species_interaction_value_large_carnivores_subordinate_with_preferred_prey_3_levels_", date, ".png", sep = ""),p,
# #        width = 10, height = 8, units = "in")


### Which preferred prey species are significant supporting tigers?
dat = coeff[coeff$var == "Species_Interaction",]
dat[dat$mean > 0 &
      dat$sub_sp == "Panthera_tigris" &
      dat$sig == "Significant" &
      dat$dom_sp %in% guilds$scientificNameStd[guilds$tiger_pref == "Yes"],] 
# 3 species --> Rusa, Sus scrofa, and Ursus thibetanus

## now what about leopards 
dat[dat$mean > 0 &
      dat$sub_sp == "Panthera_pardus" &
      dat$sig == "Significant" &
      dat$dom_sp %in% guilds$scientificNameStd[guilds$leopard_pref == "Yes"],]  
# 3 species --> Muntjac, Macaca_arctoides, & Sus scrofa

## now what about dholes 
dat[dat$mean > 0 &
      dat$sub_sp == "Cuon_alpinus" &
      dat$sig == "Significant" &
      dat$dom_sp %in% guilds$scientificNameStd[guilds$dhole_pref == "Yes"],] 
# 3 species --> Muntjac, Rusa, and Sus scrofa

## and what about SCL?
dat[dat$mean > 0 &
      dat$sub_sp == "Neofelis_genus" &
      dat$sig == "Significant" &
      dat$dom_sp %in% guilds$scientificNameStd[guilds$CL_pref == "Yes"],] 
# 9 species! --> Arctictis_binturong, Arctonyx_collaris, Capricornis_genus, Hystrix_brachyura,
# Macaca_nemestrina, Muntiacus_genus, Sus_scrofa, Viverra_tangalunga, Viverra_zibetha

### interesting, Sus scrofa is a key prey species for ALL large carnivores...
# They are the keystone species that matters in SEA! 

## who is the one negative for SCL?
dat[dat$mean < 0 &
      dat$sub_sp == "Neofelis_genus" &
      dat$sig == "Significant" &
      dat$dom_sp %in% guilds$scientificNameStd[guilds$CL_pref == "Yes"],] # Sambar! interesting! 


#
##
###
####
##### Make a unified plot w/ top-down, bottom-up, and unsupported distributions

## first grab all species interaction values
dat = coeff[coeff$var == "Species_Interaction",]

## create top-down data next-
# grab two data subsets, 
a = dat[grepl("DOM-Large_Carn", dat$guild_pair),] # one with all large carnivores as dominant
b = a
unique(b$dom_sp) # and a second where all of these are combined into one 
b$dom_sp = "All_large_carnivores"

# and combine
td = rbind(a,b)
rm(a,b)

## create a column for hypothesis support or not
td$hyp_support = "Unsupported"
td$hyp_support[td$mean < 0 & td$sig == "Significant"] = "Top-Down"

## seperate out the large carns into their own col
td$large_carnivores = td$dom_sp



## create bottom-up data
# grab two data subsets, 
a = dat[grepl("SUB-Large_Carn", dat$guild_pair),] # one with all large carnivores as dominant
b = a
unique(b$sub_sp) # and a second where all of these are combined into one 
b$sub_sp = "All_large_carnivores"

# and combine
bu = rbind(a,b)
rm(a,b)

## create a column for hypothesis support or not
bu$hyp_support = "Unsupported"
bu$hyp_support[bu$mean > 0 & bu$sig == "Significant"] = "Bottom-Up"

## seperate out the large carns into their own col
bu$large_carnivores = bu$sub_sp

## Combine bottom-up w/ top-down
tc = rbind(td, bu)
rm(td,bu)

## do it manually b/c order is silly on graph
tc$large_carnivores = factor(tc$large_carnivores, levels = c("Neofelis_genus","Cuon_alpinus", 
                                                             "Panthera_pardus","Panthera_tigris",
                                                             "All_large_carnivores"))

## and specify the colors for those levels 
hyp_colors = c("Unsupported" = "gray68",
               "Top-Down" = "firebrick4", 
               "Bottom-Up" = "darkgreen")

## make sure supported comes first! 
tc$hyp_support = factor(tc$hyp_support, levels = c("Unsupported", "Bottom-Up", "Top-Down"))


##### Subset records that wont fit in density estimation (i.e. <= 2 values in neg_vs_pos2 OR hyp_supported)
## first grab all dom_species and take note of which need subsetting. 
small_dat = ddply(tc, .(large_carnivores, hyp_support), summarize,
                  num_levls = length((hyp_support)))
## then remove all that dont need subsetting
small_dat = small_dat[small_dat$num_levls <= 2,]
## merge small_dat and ridge_dat based on common columns large_carnivores and hyp_support to ensure no repeats are added
ridge_points <- merge(tc[tc$large_carnivores %in% small_dat$large_carnivores,], 
                      small_dat, 
                      by = c("large_carnivores", "hyp_support"))
## verify we have the proper number of records 
nrow(ridge_points) == sum(small_dat$num_levls) ## MUST BE TRUE 
rm(small_dat)

## Thin ridge_dat to include only valid mods! 
# ridge_dat = ridge_dat[ridge_dat$Species_Pair %in% valid_mods,]

## Make the plots 
p=
  ggplot(tc, aes(x = mean, y = large_carnivores, fill = hyp_support)) +
  geom_density_ridges(scale=.9, alpha = .8, jittered_points = TRUE, point_alpha=.5)+
  geom_point(data = ridge_points, aes(x = mean, y = large_carnivores, color = hyp_support), 
             size = 5, shape = "circle", alpha = .8, position = position_nudge(y = .08), inherit.aes = F)+
  geom_vline(aes(xintercept = 0), color = "black") +
  theme_ridges() +
  scale_fill_manual(values = hyp_colors) +
  scale_color_manual(values = hyp_colors) +
  labs(x = "Mean Species Interaction Value", y = NULL, fill = NULL, color = NULL, 
       title = "Large carnivores dominant (red) and subordinate (green)")
# Save it! 
ggsave(paste("figures/Grouped Histograms/Histogram_Ridgeline_species_interaction_value_large_carnivores_dominant_and_subordinate_3_levels_", date, ".png", sep = ""),p,
       width = 10, height = 8, units = "in")

## Make a second plot scaled per number of observations-
p =
  ggplot(tc, aes(x = mean, y = large_carnivores, fill = hyp_support, height = stat(density)))+
  geom_density_ridges(stat = "density", scale=.9, alpha = .8, jittered_points = TRUE, point_alpha=.5)+
  geom_point(data = ridge_points, aes(x = mean, y = large_carnivores, color = hyp_support), 
             size = 5, shape = "circle", alpha = .8, position = position_nudge(y = .08), inherit.aes = F)+
  geom_vline(aes(xintercept = 0), color = "black") +
  theme_ridges() +
  scale_fill_manual(values = hyp_colors) +
  scale_color_manual(values = hyp_colors) +
  labs(x = "Mean Species Interaction Value", y = NULL, fill = NULL, color = NULL, 
       title = "Large carnivores dominant (red) and subordinate (green)")
# Save it! 
# ggsave(paste("figures/Grouped Histograms/Histogram_Ridgeline_species_interaction_value_large_carnivores_dominant_and_subordinate_3_levels_DENSITY_SCALED_", date, ".png", sep = ""),p,
#        width = 10, height = 8, units = "in")


#
##
###
####
##### Re-make unified plot w/ preferred prey species 

## first grab all species interaction values
dat = coeff[coeff$var == "Species_Interaction",]

# then, gather and combine relevant data 
a = dat[dat$sub_sp == "Neofelis_genus" & 
          dat$dom_sp %in% guilds$scientificNameStd[guilds$CL_pref == "Yes"],]
b = dat[dat$sub_sp == "Panthera_tigris" & 
          dat$dom_sp %in% guilds$scientificNameStd[guilds$tiger_pref == "Yes"],]
c = dat[dat$sub_sp == "Panthera_pardus" & 
          dat$dom_sp %in% guilds$scientificNameStd[guilds$leopard_pref == "Yes"],]
d = dat[dat$sub_sp == "Cuon_alpinus" & 
          dat$dom_sp %in% guilds$scientificNameStd[guilds$dhole_pref == "Yes"],]
## combine all for all_large_carn
e = rbind(a,b,c,d)
e$sub_sp = "All_large_carnivores"

## combine
bu = rbind(a,b,c,d,e)
rm(a,b,c,d,e)

## create a column for hypothesis support or not
bu$hyp_support = "Unsupported"
bu$hyp_support[bu$mean > 0 & bu$sig == "Significant"] = "Bottom-Up"

## seperate out the large carns into their own col
bu$large_carnivores = bu$sub_sp



#### Remake this graph, but ONLY w/ preferred prey species
dat = coeff[coeff$var == "Species_Interaction",]
# First, gather and combine relevant data 
a = dat[dat$dom_sp == "Neofelis_genus" & 
          dat$sub_sp %in% guilds$scientificNameStd[guilds$CL_pref == "Yes"],]
b = dat[dat$dom_sp == "Panthera_tigris" & 
          dat$sub_sp %in% guilds$scientificNameStd[guilds$tiger_pref == "Yes"],]
c = dat[dat$dom_sp == "Panthera_pardus" & 
          dat$sub_sp %in% guilds$scientificNameStd[guilds$leopard_pref == "Yes"],]
d = dat[dat$dom_sp == "Cuon_alpinus" & 
          dat$sub_sp %in% guilds$scientificNameStd[guilds$dhole_pref == "Yes"],]
## also adding all large carnivores w/ preferred prey
e = rbind(a,b,c,d)
e$dom_sp = "All_large_carnivores"

## combine
td = rbind(a,b,c,d,e)
rm(a,b,c,d,e)

## create a column for hypothesis support or not
td$hyp_support = "Unsupported"
td$hyp_support[td$mean < 0 & td$sig == "Significant"] = "Top-Down"

## seperate out the large carns into their own col
td$large_carnivores = td$dom_sp


## Combine bottom-up w/ top-down
tc = rbind(td, bu)
rm(td,bu)

## do it manually b/c order is silly on graph
tc$large_carnivores = factor(tc$large_carnivores, levels = c("Neofelis_genus","Cuon_alpinus", 
                                                             "Panthera_pardus","Panthera_tigris",
                                                             "All_large_carnivores"))

## and specify the colors for those levels 
hyp_colors = c("Unsupported" = "gray68",
               "Top-Down" = "firebrick4", 
               "Bottom-Up" = "darkgreen")

## make sure supported comes first! 
tc$hyp_support = factor(tc$hyp_support, levels = c("Unsupported", "Bottom-Up", "Top-Down"))


##### Subset records that wont fit in density estimation (i.e. <= 2 values in neg_vs_pos2 OR hyp_supported)
## first grab all dom_species and take note of which need subsetting. 
small_dat = ddply(tc, .(large_carnivores, hyp_support), summarize,
                  num_levls = length((hyp_support)))
## then remove all that dont need subsetting
small_dat = small_dat[small_dat$num_levls <= 2,]
## merge small_dat and ridge_dat based on common columns large_carnivores and hyp_support to ensure no repeats are added
ridge_points <- merge(tc[tc$large_carnivores %in% small_dat$large_carnivores,], 
                      small_dat, 
                      by = c("large_carnivores", "hyp_support"))
## verify we have the proper number of records 
nrow(ridge_points) == sum(small_dat$num_levls) ## MUST BE TRUE 
rm(small_dat)

## Thin ridge_dat to include only valid mods! 
# tc = tc[tc$Species_Pair %in% valid_mods,]

## Make the plots 
p=
  ggplot(tc, aes(x = mean, y = large_carnivores, fill = hyp_support)) +
  geom_density_ridges(scale=.9, alpha = .8, jittered_points = TRUE, point_alpha=.5)+
  geom_point(data = ridge_points, aes(x = mean, y = large_carnivores, color = hyp_support), 
             size = 5, shape = "circle", alpha = .8, position = position_nudge(y = .08), inherit.aes = F)+
  geom_vline(aes(xintercept = 0), color = "black") +
  theme_ridges() +
  scale_fill_manual(values = hyp_colors) +
  scale_color_manual(values = hyp_colors) +
  labs(x = "Mean Species Interaction Value", y = NULL, fill = NULL, color = NULL, 
       title = "Large carnivores dominant (red) and subordinate (green) to preferred prey")
# Save it! 
# ggsave(paste("figures/Grouped Histograms/Histogram_Ridgeline_species_interaction_value_large_carnivores_dominant_and_subordinate_with_preferred_prey_3_levels_", date, ".png", sep = ""),p,
#        width = 10, height = 8, units = "in")

## Make a second plot scaled per number of observations-
p =
  ggplot(tc, aes(x = mean, y = large_carnivores, fill = hyp_support, height = stat(density)))+
  geom_density_ridges(stat = "density", scale=.9, alpha = .8, jittered_points = TRUE, point_alpha=.5)+
  geom_point(data = ridge_points, aes(x = mean, y = large_carnivores, color = hyp_support), 
             size = 5, shape = "circle", alpha = .8, position = position_nudge(y = .08), inherit.aes = F)+
  geom_vline(aes(xintercept = 0), color = "black") +
  theme_ridges() +
  scale_fill_manual(values = hyp_colors) +
  scale_color_manual(values = hyp_colors) +
  labs(x = "Mean Species Interaction Value", y = NULL, fill = NULL, color = NULL, 
       title = "Large carnivores dominant (red) and subordinate (green) to preferred prey")
# Save it! 
# ggsave(paste("figures/Grouped Histograms/Histogram_Ridgeline_species_interaction_value_large_carnivores_dominant_and_subordinate_with_preferred_prey_3_levels_DENSITY_SCALED_", date, ".png", sep = ""),p,
#        width = 10, height = 8, units = "in")




## keep it clean
rm(dat, ridge_dat, ridge_points, order, p, color_values_short,
   color_values_full, a,b, min,max,date, hyp_colors, tc, bu, td)


######### Visually compare pairwise interactions in different directions ######

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
plots[[387]] # "" in title means there was no complementary model 

## remove plots w/ NA
plots = plots[names(plots) != ""]
## need to remove redundant plots tho 

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
                        "Chat" = c(dat$Chat.dom, dat$Chat.sub),
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
rm(bpv, chat, i, path, day,month,year,date,plot_dat, n, dat,add)



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

### firstly, how many top-down models did we run?
length(preform$Species_Pair[grepl("DOM-Large_Carn",preform$guild_pair)]) #152

### and how many bottom-up models?
length(preform$Species_Pair[grepl("SUB-Large_Carn",preform$guild_pair)]) #127

## how many subordinate species in top-down models per guild?
ddply(preform[grepl("DOM-Large_Carn",preform$guild_pair), ], .(SUB_guild), summarise,
      num_sp = length(unique(sub_species)))

## how many detections did we get for our large carnivores?
og_captures$Species = gsub(" ", "_", og_captures$Species)
ddply(og_captures[og_captures$Species %in% c("Panthera_pardus", "Panthera_tigris", 
                                             "Cuon_alpinus", "Neofelis_nebulosa", 
                                             "Neofelis_diardi" ),], .(Species), summarize,
      num_cap = length(Individuals))



#### Summarize results by guild to clearly present them-
sum = ddply(preform, .(guild_pair), summarize,
            num_sig_neg_int = sum(Interaction_Estimate < 0 & Significance == "Significant", na.rm = T), # number of significant negative interactions
            num_sig_pos_int = sum(Interaction_Estimate > 0 & Significance == "Significant", na.rm = T), # number of significant positive interactions
            num_nonsig_neg_int = sum(Interaction_Estimate < 0 & Significance == "Non-Significant", na.rm = T), # number of NON-significant negative interactions
            num_nonsig_pos_int = sum(Interaction_Estimate > 0 & Significance == "Non-Significant", na.rm = T), # number of NON-significant positive interactions
            num_combos = length(unique(Species_Pair)))     # number of different species pairs per guild 

## what is the percent of negative or positive interactions per guild?
sum$percent_sig_neg = round(sum$num_sig_neg_int / sum$num_combos * 100, 3)
sum$percent_sig_pos = round(sum$num_sig_pos_int / sum$num_combos * 100, 3)
sum$percent_nonsig_neg = round(sum$num_nonsig_neg_int / sum$num_combos * 100, 3)
sum$percent_nonsig_pos = round(sum$num_nonsig_pos_int / sum$num_combos * 100, 3)
sum
## Key results! for every guild pairing and significant interactions, 
## there is a higher percentage of positive interactions than negative interactions! 

### Middle teir testing on 2023-07-25 has ONE exception 
## with slightly more negative interactions: SUB-Small_Herbivore~DOM-Large_Carnivore
### And this could be contested w/ some species for removal eg:
#### Hystrix_genus, orangutan + langur, + pheasents, tufted ground squirell, etc. 

# grab the date
day<-str_sub(Sys.Date(),-2)
month<-str_sub(Sys.Date(),-5,-4)
year<-str_sub(Sys.Date(),-10,-7)
date = paste(year,month,day, sep = "")
rm(year,month,day)

## and save these results! 
# write.csv(sum, paste("results/summarized_results/guild_pair_interactions_", date, ".csv", sep = ""), row.names = F)
rm(sum)



#
##
###
#### investigate per large carnivore 

## Subset preform for only large carnivores as dominant
dat = preform[preform$DOM_guild == "Large_Carnivore",]
unique(dat$dom_species) # good.

#### Summarize results by large carnivore to clearly present them
sum_dom_large_carn = ddply(dat, .(dom_species), summarize,
                           num_sig_neg_int = sum(Interaction_Estimate < 0 & Significance == "Significant", na.rm = T), # number of significant negative interactions
                           num_sig_pos_int = sum(Interaction_Estimate > 0 & Significance == "Significant", na.rm = T), # number of significant positive interactions
                           num_nonsig_neg_int = sum(Interaction_Estimate < 0 & Significance == "Non-Significant", na.rm = T), # number of NON-significant negative interactions
                           num_nonsig_pos_int = sum(Interaction_Estimate > 0 & Significance == "Non-Significant", na.rm = T), # number of NON-significant positive interactions
                           num_combos = length(unique(Species_Pair)))     # number of different species pairs per dom large carn

## what is the percent of negative or positive interactions per dominant large carnivore?
sum_dom_large_carn$percent_sig_neg = round(sum_dom_large_carn$num_sig_neg_int / sum_dom_large_carn$num_combos * 100, 3)
sum_dom_large_carn$percent_sig_pos = round(sum_dom_large_carn$num_sig_pos_int / sum_dom_large_carn$num_combos * 100, 3)
sum_dom_large_carn$percent_nonsig_neg = round(sum_dom_large_carn$num_nonsig_neg_int / sum_dom_large_carn$num_combos * 100, 3)
sum_dom_large_carn$percent_nonsig_pos = round(sum_dom_large_carn$num_nonsig_pos_int / sum_dom_large_carn$num_combos * 100, 3)
sum_dom_large_carn
## Very interesting! Canis lupus is the only species with more significant negative interactions than any other! 
## But now leopards have an equal number of significant negative and positive interactions.  

## which species are feral dogs negetivly affecting?
preform$sub_species[preform$dom_species == "Canis_lupus_familiaris" &
                      preform$Interaction_Estimate < 0 &
                      # preform$mod_valid == "Yes" & # There are no sub_species if checking for valid mods!
                      preform$Significance == "Significant"]
## very interesting results here! Tigers are suppressed by feral dogs... 

# ## which guilds are negativly affected by feral dogs?
# table(preform$SUB_guild[preform$dom_species == "Canis_lupus_familiaris" &
#                           preform$Interaction_Estimate < 0 &
#                           preform$Significance == "Significant"])
# ## All guilds! But particularly small omnivores/herbivores. 
# 
# ## inspect feral dog records
# dog = bdata$`SUB-Pardofelis_marmorata~DOM-Canis_lupus_familiaris`
# # convert matrix to a df 
# dog = (dog$y.dom)
# dog = apply(dog, 2, as.numeric)
# dog = as.data.frame(dog)
# # grab row names
# dog$SU = rownames(bdata$`SUB-Pardofelis_marmorata~DOM-Canis_lupus_familiaris`$y.dom)
# # get total counts
# dog$sum_count = rowSums(dog[,1:20], na.rm = T)
# dog$landscape = sub("_\\d+.*", "", dog$SU)
# 
# ##summarize it all
# dog_sum = ddply(dog, .(landscape), summarize,
#                 total_dog_count = sum(sum_count))
# dog_sum[order(dog_sum$total_dog_count),]
# ## There is something here! 
# rm(dog,dog_sum)

#
##
###
#### remake this for only preferred prey species 
dat = preform[preform$DOM_guild == "Large_Carnivore",]
unique(dat$dom_species) # good.

## create a new dataframe
sum_dom_large_carn_pref = data.frame("dom_species" = unique(dat$dom_species))
### and fill in row by row, but via loop
prefs = c("tiger_pref","leopard_pref","dhole_pref","CL_pref") 
LC = c("Panthera_tigris", "Panthera_pardus","Cuon_alpinus","Neofelis_genus")

for(i in 1:length(LC)){
  
  #sig negative
  sum_dom_large_carn_pref$num_sig_neg_int_pref[sum_dom_large_carn_pref$dom_species == LC[i]] = 
    nrow(dat[dat$Interaction_Estimate < 0 & 
               dat$Significance == "Significant" &
               dat$dom_species == LC[i] &
               dat$sub_species %in% guilds$scientificNameStd[guilds[[prefs[i]]] == "Yes"],])
  #sig positive
  sum_dom_large_carn_pref$num_sig_pos_int_pref[sum_dom_large_carn_pref$dom_species == LC[i]] = 
    nrow(dat[dat$Interaction_Estimate > 0 & 
               dat$Significance == "Significant" &
               dat$dom_species == LC[i] &
               dat$sub_species %in% guilds$scientificNameStd[guilds[[prefs[i]]] == "Yes"],])
  #nonsig negative
  sum_dom_large_carn_pref$num_nonsig_neg_int_pref[sum_dom_large_carn_pref$dom_species == LC[i]] = 
    nrow(dat[dat$Interaction_Estimate < 0 & 
               dat$Significance == "Non-Significant" &
               dat$dom_species == LC[i] &
               dat$sub_species %in% guilds$scientificNameStd[guilds[[prefs[i]]] == "Yes"],])
  #nonsig positive
  sum_dom_large_carn_pref$num_nonsig_pos_int_pref[sum_dom_large_carn_pref$dom_species == LC[i]] = 
    nrow(dat[dat$Interaction_Estimate > 0 & 
               dat$Significance == "Non-Significant" &
               dat$dom_species == LC[i] &
               dat$sub_species %in% guilds$scientificNameStd[guilds[[prefs[i]]] == "Yes"],])
  #total num pref species
  sum_dom_large_carn_pref$num_pref_sp[sum_dom_large_carn_pref$dom_species == LC[i]] = 
    length(guilds$scientificNameStd[guilds[[prefs[i]]] == "Yes"])
  
}
rm(i,LC, prefs)

## what is the percent of negative or positive interactions per subordinate large carnivore?
sum_dom_large_carn_pref$percent_nonsig_neg_pref = round(sum_dom_large_carn_pref$num_nonsig_neg_int_pref / sum_dom_large_carn_pref$num_pref_sp * 100, 3)
sum_dom_large_carn_pref$percent_nonsig_pos_pref = round(sum_dom_large_carn_pref$num_nonsig_pos_int_pref / sum_dom_large_carn_pref$num_pref_sp * 100, 3)
sum_dom_large_carn_pref$percent_sig_neg_pref = round(sum_dom_large_carn_pref$num_sig_neg_int_pref / sum_dom_large_carn_pref$num_pref_sp * 100, 3)
sum_dom_large_carn_pref$percent_sig_pos_pref = round(sum_dom_large_carn_pref$num_sig_pos_int_pref / sum_dom_large_carn_pref$num_pref_sp * 100, 3)
sum_dom_large_carn_pref
## perfect! 

## now merge w/ all results 
sum_dom_LC = merge(sum_dom_large_carn, sum_dom_large_carn_pref, by = "dom_species")

## but calculate percentages for ALL large carnivores too,
sum_dom_LC[5, "dom_species"] = "All_large_carnivores"
sum_dom_LC[5, c(2:6, 11:15)] = colSums(sum_dom_LC[1:4, c(2:6, 11:15)])
# percentages
sum_dom_LC$percent_sig_neg[sum_dom_LC$dom_species == "All_large_carnivores"] = 
  round(sum_dom_LC$num_sig_neg_int[sum_dom_LC$dom_species == "All_large_carnivores"] / 
          sum_dom_LC$num_combos[sum_dom_LC$dom_species == "All_large_carnivores"] * 100, 3)
sum_dom_LC$percent_sig_pos[sum_dom_LC$dom_species == "All_large_carnivores"] = 
  round(sum_dom_LC$num_sig_pos_int[sum_dom_LC$dom_species == "All_large_carnivores"] / 
          sum_dom_LC$num_combos[sum_dom_LC$dom_species == "All_large_carnivores"] * 100, 3)
sum_dom_LC$percent_nonsig_neg[sum_dom_LC$dom_species == "All_large_carnivores"] = 
  round(sum_dom_LC$num_nonsig_neg_int[sum_dom_LC$dom_species == "All_large_carnivores"] / 
          sum_dom_LC$num_combos[sum_dom_LC$dom_species == "All_large_carnivores"] * 100, 3)
sum_dom_LC$percent_nonsig_pos[sum_dom_LC$dom_species == "All_large_carnivores"] = 
  round(sum_dom_LC$num_nonsig_pos_int[sum_dom_LC$dom_species == "All_large_carnivores"] / 
          sum_dom_LC$num_combos[sum_dom_LC$dom_species == "All_large_carnivores"] * 100, 3)
## preferred now
sum_dom_LC$percent_nonsig_pos_pref[sum_dom_LC$dom_species == "All_large_carnivores"] = 
  round(sum_dom_LC$num_nonsig_pos_int_pref[sum_dom_LC$dom_species == "All_large_carnivores"] / 
          sum_dom_LC$num_pref_sp[sum_dom_LC$dom_species == "All_large_carnivores"] * 100, 3)
sum_dom_LC$percent_nonsig_neg_pref[sum_dom_LC$dom_species == "All_large_carnivores"] = 
  round(sum_dom_LC$num_nonsig_neg_int_pref[sum_dom_LC$dom_species == "All_large_carnivores"] / 
          sum_dom_LC$num_pref_sp[sum_dom_LC$dom_species == "All_large_carnivores"] * 100, 3)
sum_dom_LC$percent_sig_pos_pref[sum_dom_LC$dom_species == "All_large_carnivores"] = 
  round(sum_dom_LC$num_sig_pos_int_pref[sum_dom_LC$dom_species == "All_large_carnivores"] / 
          sum_dom_LC$num_pref_sp[sum_dom_LC$dom_species == "All_large_carnivores"] * 100, 3)
sum_dom_LC$percent_sig_neg_pref[sum_dom_LC$dom_species == "All_large_carnivores"] = 
  round(sum_dom_LC$num_sig_neg_int_pref[sum_dom_LC$dom_species == "All_large_carnivores"] / 
          sum_dom_LC$num_pref_sp[sum_dom_LC$dom_species == "All_large_carnivores"] * 100, 3)
sum_dom_LC # good! 

## grab today's date
day<-str_sub(Sys.Date(),-2)
month<-str_sub(Sys.Date(),-5,-4)
year<-str_sub(Sys.Date(),-10,-7)
date = paste(year,month,day, sep = "")
rm(year,month,day)

## Save it! 
# write.csv(sum_dom_LC, paste("results/summarized_results/Dominant_large_carnivore_interactions_", date, ".csv", sep = ""), row.names = F)
rm(sum_dom_large_carn, sum_dom_large_carn_pref, date)


#
##
###
#### Repeat the process, but examine large carnivores as subordinate

## Subset preform for only large carnivores as subordinate
dat = preform[preform$SUB_guild == "Large_Carnivore",]
unique(dat$sub_species) # good.

#### Summarize results by large carnivore to clearly present them
sum_sub_large_carn = ddply(dat, .(sub_species), summarize,
                           num_sig_neg_int = sum(Interaction_Estimate < 0 & Significance == "Significant", na.rm = T), # number of significant negative interactions
                           num_sig_pos_int = sum(Interaction_Estimate > 0 & Significance == "Significant", na.rm = T), # number of significant positive interactions
                           num_nonsig_neg_int = sum(Interaction_Estimate < 0 & Significance == "Non-Significant", na.rm = T), # number of NON-significant negative interactions
                           num_nonsig_pos_int = sum(Interaction_Estimate > 0 & Significance == "Non-Significant", na.rm = T), # number of NON-significant positive interactions
                           num_combos = length(unique(Species_Pair)))  # number of different species pairs per dom large carn
                           
## what is the percent of negative or positive interactions per subordinate large carnivore?
sum_sub_large_carn$percent_sig_neg = round(sum_sub_large_carn$num_sig_neg_int / sum_sub_large_carn$num_combos * 100, 3)
sum_sub_large_carn$percent_sig_pos = round(sum_sub_large_carn$num_sig_pos_int / sum_sub_large_carn$num_combos * 100, 3)
sum_sub_large_carn$percent_nonsig_neg = round(sum_sub_large_carn$num_nonsig_neg_int / sum_sub_large_carn$num_combos * 100, 3)
sum_sub_large_carn$percent_nonsig_pos = round(sum_sub_large_carn$num_nonsig_pos_int / sum_sub_large_carn$num_combos * 100, 3)
sum_sub_large_carn
# Very interesting! Canis lupus is the only species with more significant negative interactions than any other! 


#
##
###
#### remake this for only preferred prey species 
dat = preform[preform$SUB_guild == "Large_Carnivore",]
unique(dat$sub_species) # good.

## create a new dataframe
sum_sub_large_carn_pref = data.frame("sub_species" = unique(dat$sub_species))
### and fill in row by row, but via loop
prefs = c("tiger_pref","leopard_pref","dhole_pref","CL_pref") 
LC = c("Panthera_tigris", "Panthera_pardus","Cuon_alpinus","Neofelis_genus")

for(i in 1:length(LC)){
  
  #sig negative
  sum_sub_large_carn_pref$num_sig_neg_int_pref[sum_sub_large_carn_pref$sub_species == LC[i]] = 
    nrow(dat[dat$Interaction_Estimate < 0 & 
               dat$Significance == "Significant" &
               dat$sub_species == LC[i] &
               dat$dom_species %in% guilds$scientificNameStd[guilds[[prefs[i]]] == "Yes"],])
  #sig positive
  sum_sub_large_carn_pref$num_sig_pos_int_pref[sum_sub_large_carn_pref$sub_species == LC[i]] = 
    nrow(dat[dat$Interaction_Estimate > 0 & 
               dat$Significance == "Significant" &
               dat$sub_species == LC[i] &
               dat$dom_species %in% guilds$scientificNameStd[guilds[[prefs[i]]] == "Yes"],])
  #nonsig negative
  sum_sub_large_carn_pref$num_nonsig_neg_int_pref[sum_sub_large_carn_pref$sub_species == LC[i]] = 
    nrow(dat[dat$Interaction_Estimate < 0 & 
               dat$Significance == "Non-Significant" &
               dat$sub_species == LC[i] &
               dat$dom_species %in% guilds$scientificNameStd[guilds[[prefs[i]]] == "Yes"],])
  #nonsig positive
  sum_sub_large_carn_pref$num_nonsig_pos_int_pref[sum_sub_large_carn_pref$sub_species == LC[i]] = 
    nrow(dat[dat$Interaction_Estimate > 0 & 
               dat$Significance == "Non-Significant" &
               dat$sub_species == LC[i] &
               dat$dom_species %in% guilds$scientificNameStd[guilds[[prefs[i]]] == "Yes"],])
  #total num pref species
  sum_sub_large_carn_pref$num_pref_sp[sum_sub_large_carn_pref$sub_species == LC[i]] = 
    length(guilds$scientificNameStd[guilds[[prefs[i]]] == "Yes"])
  
}
rm(i,LC, prefs)

## what is the percent of negative or positive interactions per subordinate large carnivore?
sum_sub_large_carn_pref$percent_nonsig_neg_pref = sum_sub_large_carn_pref$num_nonsig_neg_int_pref / sum_sub_large_carn_pref$num_pref_sp * 100
sum_sub_large_carn_pref$percent_nonsig_pos_pref = sum_sub_large_carn_pref$num_nonsig_pos_int_pref / sum_sub_large_carn_pref$num_pref_sp * 100
sum_sub_large_carn_pref$percent_sig_neg_pref = sum_sub_large_carn_pref$num_sig_neg_int_pref / sum_sub_large_carn_pref$num_pref_sp * 100
sum_sub_large_carn_pref$percent_sig_pos_pref = sum_sub_large_carn_pref$num_sig_pos_int_pref / sum_sub_large_carn_pref$num_pref_sp * 100
sum_sub_large_carn_pref
## perfect! 

## now merge w/ all results 
sum_sub_LC = merge(sum_sub_large_carn, sum_sub_large_carn_pref, by = "sub_species")

## but calculate percentages for ALL large carnivores too, just not preferred
sum_sub_LC[5, "sub_species"] = "All_large_carnivores"
sum_sub_LC[5, c(2:6, 11:15)] = colSums(sum_sub_LC[1:4, c(2:6, 11:15)])
# percentages
sum_sub_LC$percent_sig_neg[sum_sub_LC$sub_species == "All_large_carnivores"] = 
  round(sum_sub_LC$num_sig_neg_int[sum_sub_LC$sub_species == "All_large_carnivores"] / 
          sum_sub_LC$num_combos[sum_sub_LC$sub_species == "All_large_carnivores"] * 100, 3)
sum_sub_LC$percent_sig_pos[sum_sub_LC$sub_species == "All_large_carnivores"] = 
  round(sum_sub_LC$num_sig_pos_int[sum_sub_LC$sub_species == "All_large_carnivores"] / 
          sum_sub_LC$num_combos[sum_sub_LC$sub_species == "All_large_carnivores"] * 100, 3)
sum_sub_LC$percent_nonsig_neg[sum_sub_LC$sub_species == "All_large_carnivores"] = 
  round(sum_sub_LC$num_nonsig_neg_int[sum_sub_LC$sub_species == "All_large_carnivores"] / 
          sum_sub_LC$num_combos[sum_sub_LC$sub_species == "All_large_carnivores"] * 100, 3)
sum_sub_LC$percent_nonsig_pos[sum_sub_LC$sub_species == "All_large_carnivores"] = 
  round(sum_sub_LC$num_nonsig_pos_int[sum_sub_LC$sub_species == "All_large_carnivores"] / 
          sum_sub_LC$num_combos[sum_sub_LC$sub_species == "All_large_carnivores"] * 100, 3)
## preferred now
sum_sub_LC$percent_nonsig_pos_pref[sum_sub_LC$sub_species == "All_large_carnivores"] = 
  round(sum_sub_LC$num_nonsig_pos_int_pref[sum_sub_LC$sub_species == "All_large_carnivores"] / 
          sum_sub_LC$num_pref_sp[sum_sub_LC$sub_species == "All_large_carnivores"] * 100, 3)
sum_sub_LC$percent_nonsig_neg_pref[sum_sub_LC$sub_species == "All_large_carnivores"] = 
  round(sum_sub_LC$num_nonsig_neg_int_pref[sum_sub_LC$sub_species == "All_large_carnivores"] / 
          sum_sub_LC$num_pref_sp[sum_sub_LC$sub_species == "All_large_carnivores"] * 100, 3)
sum_sub_LC$percent_sig_pos_pref[sum_sub_LC$sub_species == "All_large_carnivores"] = 
  round(sum_sub_LC$num_sig_pos_int_pref[sum_sub_LC$sub_species == "All_large_carnivores"] / 
          sum_sub_LC$num_pref_sp[sum_sub_LC$sub_species == "All_large_carnivores"] * 100, 3)
sum_sub_LC$percent_sig_neg_pref[sum_sub_LC$sub_species == "All_large_carnivores"] = 
  round(sum_sub_LC$num_sig_neg_int_pref[sum_sub_LC$sub_species == "All_large_carnivores"] / 
          sum_sub_LC$num_pref_sp[sum_sub_LC$sub_species == "All_large_carnivores"] * 100, 3)
sum_sub_LC # good! 

## grab today's date
day<-str_sub(Sys.Date(),-2)
month<-str_sub(Sys.Date(),-5,-4)
year<-str_sub(Sys.Date(),-10,-7)
date = paste(year,month,day, sep = "")
rm(year,month,day)

## Save it! 
# write.csv(sum_sub_LC, paste("results/summarized_results/Subordinate_large_carnivore_interactions_", date, ".csv", sep = ""), row.names = F)

rm(sum_sub_large_carn, sum_sub_large_carn_pref, date)


### Try to combine both sum_sub_LC and sum_dom_LC for total interacitons
sum_LC = data.frame("Large_carnivore" = sum_dom_LC$dom_species, # grab all carnivore species
                    "num_total_combos" = sum_dom_LC$num_combos + sum_sub_LC$num_combos, # all combos
                    "num_bottom_up" = sum_sub_LC$num_sig_pos_int, # number of significant positive interactions where preds are subordinate
                    "num_top_down" = sum_dom_LC$num_sig_neg_int, # number of significant negative interactions where preds are dominant
                    "num_unsupported" = (sum_dom_LC$num_sig_pos_int + sum_sub_LC$num_sig_neg_int +
                                            sum_dom_LC$num_nonsig_neg_int + sum_sub_LC$num_nonsig_neg_int +
                                            sum_dom_LC$num_nonsig_pos_int + sum_sub_LC$num_nonsig_pos_int), # ALL interactions that are irrelevant to hypotheses. 
                    "num_pref_combos" = sum_dom_LC$num_pref_sp + sum_sub_LC$num_pref_sp, # all preferred prey species combos 
                    "num_pref_bottom_up" = sum_sub_LC$num_sig_pos_int_pref, # number of significant positive interactions where preds are subordinate to preferred prey species! 
                    "num_pref_top_down" = sum_dom_LC$num_sig_neg_int_pref,  # number of significant positive interactions where preds are subordinate to preferred prey species! 
                    "num_pref_unsupported" = (sum_dom_LC$num_sig_pos_int_pref + sum_sub_LC$num_sig_neg_int_pref +
                                                sum_dom_LC$num_nonsig_neg_int_pref + sum_sub_LC$num_nonsig_neg_int_pref +
                                                sum_dom_LC$num_nonsig_pos_int_pref + sum_sub_LC$num_nonsig_pos_int_pref)# ALL interactions that are irrelevant to hypotheses thinned to preferred prey!
                    )

## calc percentages
sum_LC$percent_unsupported = round(sum_LC$num_unsupported / sum_LC$num_total_combos * 100, 3)
sum_LC$percent_bottom_up = round(sum_LC$num_bottom_up / sum_LC$num_total_combos * 100, 3)
sum_LC$percent_top_down = round(sum_LC$num_top_down / sum_LC$num_total_combos * 100, 3)
sum_LC$percent_unsupported_pref = round(sum_LC$num_pref_unsupported / sum_LC$num_pref_combos * 100, 3)
sum_LC$percent_bottom_up_pref = round(sum_LC$num_pref_bottom_up / sum_LC$num_pref_combos * 100, 3)
sum_LC$percent_top_down_pref = round(sum_LC$num_pref_top_down/ sum_LC$num_pref_combos * 100, 3)

## grab today's date
day<-str_sub(Sys.Date(),-2)
month<-str_sub(Sys.Date(),-5,-4)
year<-str_sub(Sys.Date(),-10,-7)
date = paste(year,month,day, sep = "")
rm(year,month,day)

## Save it! 
# write.csv(sum_LC, paste("results/summarized_results/COMBO_Dom_and_Sub_large_carnivore_interactions_", date, ".csv", sep = ""), row.names = F)

rm(sum_sub_LC, sum_dom_LC, sum_LC)

#
##
###
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

## wait, why do these number differ?!
summary(og_meta$effort)
summary(og_resamp_meta$effort)
ummary(og_resamp_meta$Cell_Effort)
summary(og_resamp_meta$Avg_effort) # this is probably the most faithful data to include based on what was analyzed

### how many detections of our study species 
og_captures$Species = gsub(" ", "_", og_captures$Species)
nrow(og_captures[og_captures$Species %in% preform$sub_species,])

## how many surveys before removal?
length(unique(og_meta$survey_id)) #361

## how many surveys after removal 
length(unique(og_resamp_meta$survey_id)) #284

## how many landscapes made the cut in the end?
length(unique(og_resamp_meta$Landscape))

## how mant cameras did we deploy at each landscape?
info = ddply(og_resamp_meta, .(Landscape), summarize,
             num_cam = length(unique(cell_id_3km)))
summary(info$num_cam)
mean(info$num_cam)
sd(info$num_cam)

#### We can also make a study map (coneptually) w/ the original metadata
### But we should color our points based off which (or how many) large carnivores are present

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

######## Visualize conceptual figures ######

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
normal_data <- rnorm(100, mean = 0, sd = 6)

# Negative-skewed Normal distribution with mean= -2
# Using shape parameter to achieve negative skewness
negative_dat <- rsn(100, xi = -2, omega = 5, alpha = -3)

# Positive-skewed Normal distribution with mean=2
# Using shape parameter to achieve positive skewness
positive_dat <- rsn(100, xi = 2, omega = 5, alpha = 3)


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
df$distribution = factor(df$distribution, levels = c("Normal","positive_dat", "negative_dat"))

## and set the color scheme approprietly 
concept_col = c("Normal" = "gray68",
                "negative_dat" = "firebrick4",
                "positive_dat" = "darkgreen")

## Plotting the distributions
#H1
p = 
ggplot(df[df$distribution %in% c("negative_dat", "Normal"),], aes(x = value, fill = distribution)) +
  geom_density(aes(y = ..density..), alpha = 0.9) +
  labs(title = "H1 concept", x = "Species Interaction Value", y = "Density") +
  guides(fill = "none") +
  scale_fill_manual(values = concept_col)+
  theme_ridges()

## save it! 
# ggsave("explore/H1_negative_concept_histogram_20230904.png", p, 
#        width = 12, height = 4, units = "in")

## repeat but for the positive ones
p = 
  ggplot(df[df$distribution %in% c("positive_dat", "Normal"),], aes(x = value)) +
  geom_density(aes(y = ..density.., fill = distribution), alpha = 0.9) +
  labs(title = "H2 concept", x = "Species Interaction Value", y = "Density") +
  guides(fill = "none") +
  scale_fill_manual(values = concept_col)+
  theme_ridges()

# ggsave("explore/H2_positive_concept_histogram_20230904.png", p,
#        width = 12, height = 4, units = "in")



##### Now make a conceptual trend line to compliment the histograms

# 1. Simulate data with a logarithmic trend
set.seed(123)
n <- 100
x <- seq(0.1, 5, length.out = n)  # x should be > 0 for logarithmic relationship
y <- 10 - 3 * log(x) + rnorm(n, 0, 3)  # Logarithmic relationship with some noise

df <- data.frame(x = x, y = y)
rm(n,x,y)

# 2. Fit a logarithmic regression model
model <- lm(y ~ log(x), data = df)

# 3. Plot the fitted line with confidence interval
p = 
ggplot(df, aes(x = x, y = y)) +
  geom_smooth(method = "lm", formula = y ~ x + I(x^2), se = TRUE, col = "firebrick4") +
  labs(title = "H1: Negative predator-prey relationship", 
       x = "Predator abundance", y = "Prey abundance") +
  guides(color = "none") +
  theme_test() +  
  theme(
    axis.text = element_blank(),  # Removes axis values
    axis.ticks = element_blank(), # Removes axis ticks
    axis.title.x = element_text(size = 20, family = "Times New Roman"),  # Adjusts x-axis title properties
    axis.title.y = element_text(size = 20, family = "Times New Roman")   # Adjusts y-axis title properties
  )

# 4. save it!
# ggsave("explore/H1_negative_concept_pred-prey_relationship_20230904.png", p,
#        width = 6, height = 5, units = "in")


#### Now reverse it for the positive trend! 

# 1. Simulate data with an exponential trend
set.seed(123)
n <- 100
x <- seq(0.1, 5, length.out = n)
y <- 2 + 1.5 * exp(0.5 * x) + rnorm(n, 0, 5)  # Exponential relationship with some noise

df <- data.frame(x = x, y = y)
rm(n,x,y)

# 2. Fit an exponential regression model using a transformation
model <- lm(log(y - min(y) + 1) ~ x, data = df)

# 3. Plot the fitted line with confidence interval
p=
ggplot(df, aes(x = x, y = y)) +
  geom_smooth(method = "lm", formula = y ~ exp(x), se = TRUE, col = "darkgreen") +
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
# ggsave("explore/H2_positive_concept_pred-prey_relationship_20230904.png", p,
#        width = 6, height = 5, units = "in")

## keep it clean!
rm(df, concept_col, model)

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

