##### Inspect and visualize completed co-abundance models
#### Data was prepared on R Script: Step 1_Generate co-abundance data bundles.Rmd
#### and ran on HPC on R Script: HPC_co-abundance_model.R
### Will convert this to Rmd when ready, but R script for exploring now
## Import dataframes, not complete models b/c theyre too heavy

# Zachary Amir, Z.Amir@uq.edu.au
# Last updated on June 16th, 2023

## start fresh
rm(list = ls())

## load libs 
library(jagsUI)
library(tidyverse)
library(plyr)
library(reshape2)
library(fs) # for moving files around 
# library(traitdata) # for species trait data

## set WD (should be automatic b/c RProj, but being safe anyway)
setwd("/Users/zachary_amir/Dropbox/Zach PhD/Trophic release project/SEA_trophic_cascades_co-abundance")


######## Import co-abundance coefficent dataframes ######

## first, I will import the bundled data to generate a vector of relevant species pairs 
bdata = readRDS("data/send_to_HPC/Bundled_data_for_Bayes_co-abundance_mods_610_species_pairs_20230617.RDS")
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


### Inspect which species pairs are present/missing, 20230615
setdiff(all_combos, coeff$Species_Pair) # only 20 missing models!! This is an improvement, but still check. 
which(all_combos == "SUB-Paradoxurus_hermaphroditus~DOM-Cuon_alpinus") # currently re-running on HPC
which(all_combos == "SUB-Arctictis_binturong~DOM-Canis_lupus_familiaris") # 2nd fix made this work! 

### Long-story-short, I found out they all had problem SUs from DVCA that contained 14+ active cams 
## --> removed all SUs w/ > 14 cameras and re-run analysis. 
## 2nd update --> removed all SUs w/ > 7 cameras and re-run analysis with MUCH better performance. 

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
table(preform$mod_completion) # only 3.4% of models have failed! Epic improvement. 

#
##
###
#### Combine species trait/guild data w/ model performance

#load guild data
guilds = read.csv(("data/species_traits/clean_72_species_trait_data_20230615.csv"))
head(guilds)
str(guilds)

## make sure all speceis are present
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

ppc_values = do.call(rbind, ppc_res$values)
head(ppc_values) 
ppc_values$X = NULL # damn row names 
str(ppc_values) # looks good! 

# remove list now that were all stored in dfs. 
rm(ppc_res)

### Inspect which species pairs are present/missing, 20230615
setdiff(all_combos, ppc_values$Species_Pair) #44 missing models, certainly an improvement! 7% of mods are missing. 
## See Quick trouble shooting pt2, 20230615 for trouble shooting missing combos
## Long story short, I needed to reduce num_active_cams to < 7 for models to converge --> on the way to better results! 


######## Import co-abundance prediction dataframes ######

## Make sure all_combos is loaded here! 

# list all files
files = list.files("results/prediction_dataframes/")
files = files[!grepl("OLD", files)] # remove any old data 

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


### Inspect which species pairs are present/missing, 20230615
setdiff(all_combos, pred_abund$Species_Pair) #496 missing mods! But this is because many are still running on HPC. 


###### Visualize coefficient effect sizes ###### 

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
  
}
rm(i, path)
## will get warnings if the directory is already created, no dramas tho. 


#
##
###
#### Loop through all species pairs to make a plot for each!

sp_pair_coeff_plots = list() # save em here. 

for(i in 1:length(unique(coeff$Species_Pair))){
  
  ## subset for a single species 
  dat = coeff[coeff$Species_Pair == unique(coeff$Species_Pair)[i],]
  
  ## subset for relevant variables 
  dat = dat[dat$var %in% c("Abundance_intercept","FLII","HFP","Elevation",
                           "Species_Interaction","Detection_intercept","Active_cams"),]
  
  ## replace _ with space for vars for clean x-vars
  dat$var = gsub("_", " ", dat$var)
  
  ## Change active cams to effort for easier understanding
  dat$var[dat$var == "Active cams"] = "Effort"
  
  ## order the factor levels for each variable
  dat$var = factor(dat$var, levels = c("Species Interaction","FLII","HFP","Elevation", "Abundance intercept", # state vars
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
  path = paste("figures/", gp, "/model_coefficents_", n, "_", date, ".png", sep = "")
  
  ## save it! 
  ggsave(path, p, width = 10.5, height = 6.5, units = "in")
  
}
rm(p,n,gp,i,path,day,month,year,date)

# dont need these plots anymore as they are saved
rm(sp_pair_coeff_plots)

#
##
###
#### Generate histograms per guild pair + overall 

# first we need to add guild pairs to coeff (probs should have been done earlier)
add = distinct(select(preform, Species_Pair, guild_pair))
coeff = merge(coeff, add, by = "Species_Pair")
rm(add)

## add a column if interaction is significant or not and negative or positve
coeff$neg_vs_pos = "Positive & Non-Significant"
coeff$neg_vs_pos[coeff$mean < 0 & coeff$sig == "Significant"] = "Negative & Significant"
coeff$neg_vs_pos[coeff$mean < 0 & coeff$sig == "Non-Significant"] = "Negative & Non-Significant"
coeff$neg_vs_pos[coeff$mean > 0 & coeff$sig == "Significant"] = "Positive & Significant"

# ## add colors 
# coeff$neg_pos_color[coeff$neg_vs_pos == "Positive_Non-Significant"] = 'darkolivegreen1'
# coeff$neg_pos_color[coeff$neg_vs_pos == "Positive_Significant"] = "darkgreen"
# coeff$neg_pos_color[coeff$neg_vs_pos == "Negative_Non-Significant"] = "goldenrod1"
# coeff$neg_pos_color[coeff$neg_vs_pos == "Negative_Significant"] = "firebrick4"
# specify manually
color_values <- c("Positive & Non-Significant" = "darkolivegreen1",
                  "Positive & Significant" = "darkgreen",
                  "Negative & Non-Significant" = "goldenrod1",
                  "Negative & Significant" = "firebrick4")

# subset the data
dat = coeff[coeff$var ==  "Species_Interaction",]

# Calculate the median value
median_value <- median(dat$mean)
 
# count how many values we have that will be graphed
pos = ddply(dat, .(neg_vs_pos), summarize,
            count = length(neg_vs_pos))
# make the label just below the largest count. 
y_position <- 0.75 * max(pos$count)

# test the plot 
p = 
ggplot(dat, aes(x = mean))+
  geom_histogram(aes(y = after_stat(count), fill = neg_vs_pos), bins = 50)+
  scale_fill_manual(values = color_values)+
  geom_vline(aes(xintercept = 0), linetype = "dashed")+
  geom_vline(aes(xintercept = median(mean)), color = "purple") + 
  annotate("text", x = median_value, y = y_position, label = paste("Median =", round(median(dat$mean), 2)),
           vjust = 1, hjust = -0.2, color = "black", size = 4) +
  theme_test()+
  labs(x = "Mean Species Interaction Value", y = "Frequency", fill = NULL, title = "All species pairs")

## save it!
# ggsave("figures/Histogram_species_interaction_value_ALL_species_pairs_20230626.png", p, 
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
  
  # count how many values we have that will be graphed
  pos = ddply(dat, .(neg_vs_pos), summarize,
              count = length(neg_vs_pos))
  # make the label just below the largest count. 
  y_position <- 0.75 * max(pos$count)
  
  # make the plot 
  p = 
    ggplot(dat, aes(x = mean))+
    geom_histogram(aes(y = after_stat(count), fill = neg_vs_pos), bins = 50)+
    scale_fill_manual(values = color_values)+
    geom_vline(aes(xintercept = 0), linetype = "dashed")+
    geom_vline(aes(xintercept = median(mean)), color = "purple") + 
    annotate("text", x = median_value, y = y_position, label = paste("Median =", round(median(dat$mean), 2)),
             vjust = 1, hjust = -0.2, color = "black", size = 4) +
    theme_test()+
    labs(x = "Mean Species Interaction Value", y = "Frequency", fill = NULL, title = gp)
  
  ## create a relevant file directory to save it
  # grab today's date
  day<-str_sub(Sys.Date(),-2)
  month<-str_sub(Sys.Date(),-5,-4)
  year<-str_sub(Sys.Date(),-10,-7)
  date = paste(year,month,day, sep = "")
  
  # construct a saving path 
  path = paste("figures/", gp, "/Histogram_species_interaction_values_", gp, "_", date, ".png", sep = "")
  
  ## save it 
  # ggsave(path, p, width = 12, height = 8, units = "in")
  
}
rm(p,i,gp,pos, y_position, dat, median_value, 
   day, month, year, date,path)


### Now make two more histos --> one w/ sub large carn and dom large carn
## Pt1 --> dominant
# subset the data
dat = coeff[coeff$var ==  "Species_Interaction" & grepl("DOM-Large_Carnivore", coeff$guild_pair),]

# Calculate the median value
median_value <- median(dat$mean)

# count how many values we have that will be graphed
pos = ddply(dat, .(neg_vs_pos), summarize,
            count = length(neg_vs_pos))
# make the label just below the largest count. 
y_position <- 0.5 * max(pos$count)

# test the plot 
p = 
  ggplot(dat, aes(x = mean))+
  geom_histogram(aes(y = after_stat(count), fill = neg_vs_pos), bins = 50)+
  scale_fill_manual(values = color_values)+
  geom_vline(aes(xintercept = 0), linetype = "dashed")+
  geom_vline(aes(xintercept = median(mean)), color = "purple") + 
  annotate("text", x = median_value, y = y_position, label = paste("Median =", round(median(dat$mean), 2)),
           vjust = 1, hjust = -0.2, color = "black", size = 4) +
  theme_test()+
  labs(x = "Mean Species Interaction Value", y = "Frequency", fill = NULL, title = "Large carnivores are dominant")
# # save it!
# ggsave("figures/Histogram_species_interaction_value_large_carnivores_dominant_20230626.png", p,
#        width = 12, height = 8, units = "in")

 
 ## Pt2 --> subordinate
# subset the data
dat = coeff[coeff$var ==  "Species_Interaction" & grepl("SUB-Large_Carnivore", coeff$guild_pair),]

# Calculate the median value
median_value <- median(dat$mean)

# count how many values we have that will be graphed
pos = ddply(dat, .(neg_vs_pos), summarize,
            count = length(neg_vs_pos))
# make the label just below the largest count. 
y_position <- 0.75 * max(pos$count)

# test the plot 
p = 
  ggplot(dat, aes(x = mean))+
  geom_histogram(aes(y = after_stat(count), fill = neg_vs_pos), bins = 50)+
  scale_fill_manual(values = color_values)+
  geom_vline(aes(xintercept = 0), linetype = "dashed")+
  geom_vline(aes(xintercept = median(mean)), color = "purple") + 
  annotate("text", x = median_value, y = y_position, label = paste("Median =", round(median(dat$mean), 2)),
           vjust = 1, hjust = -0.2, color = "black", size = 4) +
  theme_test()+
  labs(x = "Mean Species Interaction Value", y = "Frequency", fill = NULL, title = "Large carnivores are subordinate")
# # save it!
# ggsave("figures/Histogram_species_interaction_value_large_carnivores_subordinate_20230626.png", p,
#        width = 12, height = 8, units = "in")




#
###### Visualize interaction across community via matrix ######

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







###### Quick trouble shooting pt1, 20230614 ####

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

###### Quick trouble shooting pt2, 20230615 ####### 

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
