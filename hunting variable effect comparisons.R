######## Examine & Visualize hunting variable tests in co-abundance models 

####### Due to counter-intuitve positive results from the hunting variable 
###### Showing generally positive affects on wildlife abundance in co-abundance mods, 
##### We ran 2 tests with the hunting variable on 13 models 
### 1) Remove hunting variable completly
### 2) Use the log10 transformation of the hunting variable 
## All tests were run with HALF_LONG settings to ensure valid comparisons. 

## The goal is to compare the SIV between these different tests to examine
## if hunting drives species interations
### OR 
## If we see a positve result with the negative log10 version b/c values are negative 

### AND, a secondary goal is to see how wildlife respond to this hunting variable 
## both on their own and how the variable mediates species interactions. 


##### Zachary Amir, Z.Amir@uq.edu.au
#### Code made/finished on January 18th, 2024

## start fresh
rm(list = ls())

## load libs 
library(tidyverse)
library(plyr)
library(unmarked)
library(ggridges)

## unsure if this will go on github or not, just set WD to be safe 
setwd("/Users/zachary_amir/Dropbox/Zach PhD/Ch3 Trophic release project/SEA_trophic_cascades_co-abundance")


####### Import various types of models #####

### I have run ALL 66 preferred prey models with and without the regular hunting variable, 
## so start here

#### first import model preformance dataframe to know which mods are preferred or not 
preform = read.csv("results/summarized_results/model_performance_and_SIVs_20231215.csv")
pref_mods = preform$Species_Pair[preform$preference == "preferred"]

#### First import models with hunting variable 

## List all the relevant coefficent files 
files = list.files("results/HALF_LONG_testing_Oct_2023/coefficent_dataframes/")
files = files[!grepl("OLD", files)] # dont want old files

# store results here
coeff.res = list()

# loop thru each file to import into list
for(i in 1:length(files)){
  
  # import 
  d = read.csv(paste("results/HALF_LONG_testing_Oct_2023/coefficent_dataframes/", files[i], sep = ""))
  
  ## if there are pesky row.names, remove em!
  d$X = NULL
  
  # save it in the list 
  coeff.res[[i]] = d
  
  if(ncol(d) != 10){
    print(paste("The file", files[i], "with the index value of", i, 
                "has the wrong number of columns and wont combine well. INVESTIGTE!"))
  }
}
rm(d,i, files)

## bind together and inspect
hunt_coeff = do.call(rbind, coeff.res)

## thin to only the preferred prey mods
hunt_coeff = hunt_coeff[hunt_coeff$Species_Pair %in% pref_mods,]

## quick inspection
head(hunt_coeff)
str(hunt_coeff)
sort(unique(hunt_coeff$Species_Pair)) # 66 mods, good! 
rm(coeff.res)


#### Next, import models WITHOUT the hunting variable 

## List all the relevant coefficent files 
files = list.files("results/HALF_LONG_final_Jan_2024/coefficent_dataframes/")
files = files[!grepl("OLD", files)] # dont want old files

# store results here
coeff.res = list()

# loop thru each file to import into list
for(i in 1:length(files)){
  
  # import 
  d = read.csv(paste("results/HALF_LONG_final_Jan_2024/coefficent_dataframes/", files[i], sep = ""))
  
  ## if there are pesky row.names, remove em!
  d$X = NULL
  
  # save it in the list 
  coeff.res[[i]] = d
  
  if(ncol(d) != 10){
    print(paste("The file", files[i], "with the index value of", i, 
                "has the wrong number of columns and wont combine well. INVESTIGTE!"))
  }
}
rm(d,i, files)

## bind together and inspect
no_hunt_coeff = do.call(rbind, coeff.res)

## thin to only the preferred prey mods
no_hunt_coeff = no_hunt_coeff[no_hunt_coeff$Species_Pair %in% pref_mods,]

## quick inspection
head(no_hunt_coeff)
str(no_hunt_coeff)
sort(unique(no_hunt_coeff$Species_Pair)) # 66 mods, good! 
setdiff(no_hunt_coeff$Species_Pair, hunt_coeff$Species_Pair)
setdiff(hunt_coeff$Species_Pair, no_hunt_coeff$Species_Pair) # all present and accounted for! 
# keep it clean! 
rm(coeff.res)


### Finally, I also ran a test where I log transformed the variable
## but I messed up when saving data for this test and 
# only ran models on kind random species pairs 

### import the log hunting results 
files = list.files("results/Hunting_test_Dec_2023/coefficent_dataframes/")
files = files[grepl("LOG_HUNT", files)] # only want the logged results 

# store results here
coeff.res = list()

# loop thru each file to import into list
for(i in 1:length(files)){
  
  # import 
  d = read.csv(paste("results/Hunting_test_Dec_2023/coefficent_dataframes/", files[i], sep = ""))
  
  ## if there are pesky row.names, remove em!
  d$X = NULL
  
  # save it in the list 
  coeff.res[[i]] = d
  
  if(ncol(d) != 10){
    print(paste("The file", files[i], "with the index value of", i, 
                "has the wrong number of columns and wont combine well. INVESTIGTE!"))
  }
}
rm(d,i, files)

## bind together and inspect
log_hunt_coeff = do.call(rbind, coeff.res)

## quick inspection
head(log_hunt_coeff)
str(log_hunt_coeff)
sort(unique(log_hunt_coeff$Species_Pair)) # 12 pairs,
# check if these match any preferred mods?
log_hunt_coeff = log_hunt_coeff[log_hunt_coeff$Species_Pair %in% pref_mods,]
unique(log_hunt_coeff$Species_Pair) # only 2! Leave these for now, may or may not be useful
rm(coeff.res)


#### Combine all data into a single data frame for easier comparison
## but first make sure we can distinguish them! 
hunt_coeff$test = "regular_hunting"
no_hunt_coeff$test = "no_hunting"
log_hunt_coeff$test = "log_hunting"

# and combine 
coeff = rbind(hunt_coeff, no_hunt_coeff, log_hunt_coeff)
# keep it clean! 
rm(hunt_coeff, no_hunt_coeff, log_hunt_coeff)

## split apart species pair into sub vs dom species
coeff <- within(coeff, {
  sub_sp <- sapply(strsplit(Species_Pair, "~"), function(x) x[1])
  dom_sp <- sapply(strsplit(Species_Pair, "~"), function(x) x[2])
})
# clean it a bit 
coeff$dom_sp = gsub("DOM-", "", coeff$dom_sp)
coeff$sub_sp = gsub("SUB-", "", coeff$sub_sp)

## Add position to denote if the variable affects the dominant or subordinate species
coeff$position[grepl("SUB", coeff$species)] = "subordinate"
coeff$position[grepl("DOM", coeff$species)] = "dominant"


####### Compare SIV with and without hunting variable ######

## first pull out just the SIVs 
d = coeff[coeff$var == "Species_Interaction",]
table(d$test) # should be equal, except for the log mods 

## just remove log mods for now, not useful
d = d[d$test != "log_hunting",]

## check which models produced substantially different results 
# either in terms of significance, parameter convergence, or effect size < or > than 0.3
d$difference = as.character(NA)
for(i in 1:length(unique(d$Species_Pair))){
  
  # subset data for the relevant species pair 
  d1 = d[d$Species_Pair == unique(d$Species_Pair)[i],]
  
  # if the parameter significance is different, 
  if(d1$sig[d1$test == "regular_hunting"] !=
     d1$sig[d1$test == "no_hunting"]){
    
    # let us know! 
    d1$difference = "different"
    
    ## update dataset with new results 
    d$difference[d$Species_Pair == unique(d$Species_Pair)[i]] = unique(d1$difference)
    
    # Then skip the rest of the code b/c they are different and we dont want to overwrite it! 
    next
    
  }else{
    # but if not, also let us know
    d1$difference = "no_difference"
    ## update dataset with new results 
    d$difference[d$Species_Pair == unique(d$Species_Pair)[i]] = unique(d1$difference)
  } # end sig 
  
  # if either parameter did not converge,
  if(d1$Rhat[d1$test == "regular_hunting"] > 1.2 |
     d1$Rhat[d1$test == "no_hunting"] > 1.2){
    
    # let us know! 
    d1$difference = "different"
    
    ## update dataset with new results 
    d$difference[d$Species_Pair == unique(d$Species_Pair)[i]] = unique(d1$difference)
    
    # Then skip the rest of the code b/c they are different and we dont want to overwrite it! 
    next
    
  }else{
    # but if not, also let us know
    d1$difference = "no_difference"
    ## update dataset with new results 
    d$difference[d$Species_Pair == unique(d$Species_Pair)[i]] = unique(d1$difference)
    } # end convergence
  
  # if the effect sizes are very different,
  if(abs(d1$mean[d1$test == "regular_hunting"] - d1$mean[d1$test == "no_hunting"]) > 0.3){
    
    # let us know! 
    d1$difference = "different"
    
    ## update dataset with new results 
    d$difference[d$Species_Pair == unique(d$Species_Pair)[i]] = unique(d1$difference)
    
    # Then skip the rest of the code b/c they are different and we dont want to overwrite it! 
    next
    
  }else{
    # but if not, also let us know
    d1$difference = "no_difference"
    ## update dataset with new results 
    d$difference[d$Species_Pair == unique(d$Species_Pair)[i]] = unique(d1$difference)
    } # end effect size 
    
 
} # end per pair 
rm(i,d1)

### BUT, lets also make sure were not examining different but not relevant models 
## e.g. they are different, but both non-significant or poor convergence. 

## first save which pairs are differnet b/c it will change as the loop goes! 
diff_pairs = unique(d$Species_Pair[d$difference == "different"])

for(i in 1:length(diff_pairs)){
  
  # subset for 1 species
  d1 = d[d$Species_Pair == diff_pairs[i],]
  
  ## if both have large Rhat values, 
  if(d1$Rhat[1] > 1.2 & d1$Rhat[2] > 1.2){
    ## update dataset with new results 
    d$difference[d$Species_Pair == diff_pairs[i]] = "different_but_irrelevant"
  }
  
  ## if both have non-significant SIVs,
  if(d1$sig[1] == "Non-Significant" & d1$sig[2] == "Non-Significant"){
    ## update dataset with new results 
    d$difference[d$Species_Pair == diff_pairs[i]] = "different_but_irrelevant"
  }
  
}
rm(d1, i)

## inspect
anyNA(d$difference) #MUST BE F
table(d$difference) # good! minority are different, but still a third of results... 
## but only half of the differnet results are relevant, which should make the resutls easier to understand. 

## which species pairs were different?
sort(unique(d$Species_Pair[d$difference == "different"])) # quite different species! 11 pairs

## update diff_pairs now that we've excluded irrelevant ones
diff_pairs = unique(d$Species_Pair[d$difference == "different"])

## save which species showed different responses 
diff_sp = unique(c(d$dom_sp[d$difference == "different"],
                   d$sub_sp[d$difference == "different"]))
sort(diff_sp) # 12 species, will be useful later when looking at affect of hunting on species 

### visualize the differences

# Add a column for significance marker
d$marker <- ifelse(d$sig == "Significant", "*", "")

# set spacing 
dodge_width = 1

## set a logical order of species pairs based on most frequent large carns, then prey 
sort(table(d$dom_sp[d$difference == "different"]))
order = unique(c(d$Species_Pair[grepl("DOM-Neofelis_genus", d$Species_Pair) & d$difference == "different"],
                 d$Species_Pair[grepl("DOM-Panthera_tigris", d$Species_Pair) & d$difference == "different"],
                 d$Species_Pair[grepl("DOM-Panthera_pardus", d$Species_Pair) & d$difference == "different"],
                 d$Species_Pair[grepl("DOM-Cuon_alpinus", d$Species_Pair) & d$difference == "different"],
                 d$Species_Pair[grepl("DOM-Macaca_nemestrina", d$Species_Pair) & d$difference == "different"],
                 d$Species_Pair[grepl("DOM-Arctictis_binturong", d$Species_Pair) & d$difference == "different"]))

## create a new df w/ only different results and ordered factors
d_plot = d[d$difference == "different",]
d_plot$Species_Pair = factor(d_plot$Species_Pair, levels = order)

## make the plot
p =
  ggplot(d_plot, aes(x = test, y = mean, fill = test))+
  geom_bar(stat = "identity", position = position_dodge(width = dodge_width), 
           color = "black", width = dodge_width) +
  geom_errorbar(aes(ymin = mean-sd, ymax = mean+sd), 
                position = position_dodge(width = dodge_width), 
                width = 0, color = "black") +
  geom_text(aes(y = mean+.3*sign(mean), label = marker), position = position_dodge(width = dodge_width), size = 8) +
  geom_text(aes(y = .04, label = round(Rhat, 2)), position = position_dodge(width = dodge_width), size = 3, vjust = -.6, hjust = -.8) +
  theme_test()+
  facet_wrap(~Species_Pair)+
  labs(y = "Species Interaction Value", x = NULL)+
  theme(axis.text.y = element_text(size = 16),
        axis.title.x = element_text(size = 16),
        axis.text = element_text(color = "black"),
        text = element_text(family = "Helvetica"))
p
### Very interesting! Lots of differences to unpack here... 

# grab today's date  
day<-str_sub(Sys.Date(),-2)
month<-str_sub(Sys.Date(),-5,-4)
year<-str_sub(Sys.Date(),-10,-7)
date = paste(year,month,day, sep = "")
rm(year,month,day)
# And save it! 
# ggsave(paste("explore/SIV_differences_hunting_vs_no_hunting_pref_prey_11_mods_HALF_LONG_", date, ".png", sep = ""), p,
#        height = 8, width = 12, units = "in")

## clean up the enviro, were done here. 
rm(p, dodge_width, order, d_plot, d)


####### Examine hunting effect sizes #######

### There are two approaches here--
# 1) examine the effect size from co-abundance mods from same species, which is what we do in the MS
# 2) run single-species N-mixture models w/ hunting to get a single effect size. 
# 3) compare them! 

### for simplicity sake, we will only look at the species with differing results 
diff_sp
##BUT, after a bit of playing below, I know we dont have enough observations of all, so remove em and simplify! 
diff_sp = diff_sp[! diff_sp %in% c("Macaca_fascicularis",
                                   "Viverra_zibetha",
                                   "Macaca_nemestrina",
                                   "Arctictis_binturong")]

### first, extract hunting coefficent values for these species 
hunt_coeff = coeff[coeff$var == "Hunting" & coeff$dom_sp %in% diff_sp & coeff$sub_sp %in% diff_sp, ]

# inspect
head(hunt_coeff)
length(unique(hunt_coeff$Species_Pair)) # 40 mods here. 
length(unique(hunt_coeff$dom_sp));length(unique(hunt_coeff$sub_sp)) # 12 species both here, all represented. 
table(hunt_coeff$position) # equal 

## goal is to visualize like a ridgeline plot, where each species is the y axis, 
## and each species has a dominant and subordinate range
### AND, I will run a N-mixture model for each species and plot the single value from there. 
# Ideally the ranges per species overlap substantially! 

## quickly summarize the data to look for any species that may be more variabel than others
ddply(hunt_coeff, .(species), summarise,
      mean_ES = mean(mean),
      SD_ES = sd(mean), 
      n_obs = length(species))
## leopards are the crazy ones! But maybe too crazy, just start w/ capricornis 

## make a new column to combine dom and sub 
hunt_coeff$rel_species = hunt_coeff$species
hunt_coeff$rel_species = gsub("DOM-", "", hunt_coeff$rel_species)
hunt_coeff$rel_species = gsub("SUB-", "", hunt_coeff$rel_species)

## set factor order for nice plot
hunt_coeff$rel_species = factor(hunt_coeff$rel_species, 
                                levels = c("Muntiacus_genus", "Sus_scrofa","Capricornis_genus",
                                           "Rusa_unicolor",  "Neofelis_genus","Cuon_alpinus",
                                           "Panthera_pardus", "Panthera_tigris"))
                                           

#### Visualize the results 
p =
ggplot(hunt_coeff[hunt_coeff$rel_species != "Panthera_pardus",], 
       aes(x = mean, y = rel_species, fill = position)) +
  geom_density_ridges(scale=1, alpha = .8, jittered_points = T, point_alpha=.5)+
  geom_vline(aes(xintercept = 0), color = "black") +
  theme_test() +
  labs(x = "Hunting Effect Size")
## save em! 
# ggsave(paste("explore/hunting_non-transformed_effect_size_7_species_HALF_LONG_", date, ".png", sep = ""), p,
#        height = 8, width = 7, units = "in")


#### This is giving errors, and this is the right way to do this anyway
#### Skip for now, but leaving unhashed if worth revisiting. 
# #
# ##
# ### Before we make the plot, we will run a quick N-mixture model model on each rel species
# 
# ### Import data bundles 
# bdata = readRDS("data/send_to_HPC/OLD/Bundled_data_for_Bayes_co-abundance_mods_260_species_pairs_20231106.RDS")
# 
# ## extract relevant info for each species
# res = list() # store results here
# for(i in 1:length(diff_sp)){
#   
#   # thin to bundles w/ relevant species as dom
#   dat = bdata[grepl(paste("DOM-", diff_sp[i], sep = ""), names(bdata))]
#   # and just take the first for simplcity sake 
#   dat = dat[[1]]
#   
#   # now only want to retain sites where species is present (z.dom =1)
#   pos = which(dat$Z.dom == 1)
#   
#   # and thin matrix and obs covs to match 
#   y = dat$y.dom[pos,]
#   o_cov = list("num_cams" = dat$cams[pos,])
#   
#   ## make site covs based on dat
#   s_cov = data.frame("SU" = rownames(dat$y.dom[pos, ]),
#                      "hunting" = dat$hunt[pos],
#                      "flii" = dat$flii[pos],
#                      "hfp" = dat$hfp[pos],
#                      "elev" = dat$elev[pos])
#   
#   ## combine into a umf 
#   umf = unmarkedFramePCount(y = y, siteCovs = s_cov, obsCovs = o_cov)
#   
#   ## save it! 
#   res[[i]] = umf
#   names(res)[i] = diff_sp[i]
#   
# }
# rm(bdata, dat, o_cov, s_cov, umf, y, i, pos)
# 
# # inspect
# summary(res$Neofelis_genus) # seems like it worked! 
# ## vars arent well centered, but not suprising because we removed sites. 
# 
# ### Run quick N-mixture loop to compare all covaraites 
# for(i in 1:length(res)){
#   
#   ## select one umf 
#   umf = res[[i]]
#   
#   ## null mod
#   null = pcount(~num_cams ~ 1, umf)
#   
#   ## hunt mod 
#   m1 = pcount(~num_cams ~ hunting, umf)
#   ## flii mod
#   ## hfp mod
#   ## elev mod
#   
#   # rank em
#   
# }



###### Examine the other variables w/ and w/out hunting ######

## Matthew wants to know how the other variables included in the abundance formula differ 
## with and without the hunting variable.
# So reproduce the code from SIVs, but change to other vars 

## BUT, using a larger threshold for effect sizes than SIV b/c SIV is usually so small
# using 0.8 instead of 0.3 here, a 0.5 increase. 

## save relevant vars 
unique(coeff$var)
vars= c("FLII","HFP","Elevation")

## make it a loop! 
plots = list() # save plots here 
diffs = list() # save difference dataframes here 
for(v in 1:length(vars)){
  
  ## first pull out just the relevant var  
  d = coeff[coeff$var == vars[v],]
  
  ## just remove log mods for now, not useful
  d = d[d$test != "log_hunting",]
  
  ## check which models produced substantially different results 
  # either in terms of significance, parameter convergence, or effect size < or > than 0.3
  d$difference = as.character(NA)
  for(i in 1:length(unique(d$Species_Pair))){
    
    # subset data for the relevant species pair 
    d1 = d[d$Species_Pair == unique(d$Species_Pair)[i],]
    
    ### b/c both dom and sub species have each variable, we are going to assess each individually
    ## via a nested loop per positon 
    for(p in 1:length(unique(d1$position))){
      
      ## subset data for the relevant position
      dp = d1[d1$position == unique(d1$position)[p], ]
      
      # if the parameter significance is different, 
      if(dp$sig[dp$test == "regular_hunting"] !=
         dp$sig[dp$test == "no_hunting"]){
        
        # let us know! 
        dp$difference = "different"
        
        ## update dataset with new results 
        d$difference[d$Species_Pair == unique(d$Species_Pair)[i] & 
                       d$position == unique(d1$position)[p]] = unique(dp$difference)
        
        # Then skip the rest of the code b/c they are different and we dont want to overwrite it! 
        next
        
      }else{
        # but if not, also let us know
        dp$difference = "no_difference"
        ## update dataset with new results 
        d$difference[d$Species_Pair == unique(d$Species_Pair)[i] & 
                       d$position == unique(d1$position)[p]] = unique(dp$difference)
      } # end sig 
      
      # if either parameter did not converge,
      if(dp$Rhat[dp$test == "regular_hunting"] > 1.2 |
         dp$Rhat[dp$test == "no_hunting"] > 1.2){
        
        # let us know! 
        dp$difference = "different"
        
        ## update dataset with new results 
        d$difference[d$Species_Pair == unique(d$Species_Pair)[i] & 
                       d$position == unique(d1$position)[p]] = unique(dp$difference)
        
        # Then skip the rest of the code b/c they are different and we dont want to overwrite it! 
        next
        
      }else{
        # but if not, also let us know
        dp$difference = "no_difference"
        ## update dataset with new results 
        d$difference[d$Species_Pair == unique(d$Species_Pair)[i] & 
                       d$position == unique(d1$position)[p]] = unique(dp$difference)
      } # end convergence
      
      # if the effect sizes are very different,
      if(abs(dp$mean[dp$test == "regular_hunting"] - 
             dp$mean[dp$test == "no_hunting"]) > 0.8){
        
        # let us know! 
        dp$difference = "different"
        
        ## update dataset with new results 
        d$difference[d$Species_Pair == unique(d$Species_Pair)[i] & 
                       d$position == unique(d1$position)[p]] = unique(dp$difference)
        
        # Then skip the rest of the code b/c they are different and we dont want to overwrite it! 
        next
        
      }else{
        # but if not, also let us know
        dp$difference = "no_difference"
        ## update dataset with new results 
        d$difference[d$Species_Pair == unique(d$Species_Pair)[i] & 
                       d$position == unique(d1$position)[p]] = unique(dp$difference)
      } # end effect size 
      
    } # end per position 
    
  } # end per pair 
  rm(i,d1, dp)
  
  ### BUT, lets also make sure were not examining different but not relevant models 
  ## e.g. they are different, but both non-significant or poor convergence. 
  
  ## first save which pairs are different b/c it will change as the loop goes! 
  diff_pairs = unique(d$Species_Pair[d$difference == "different"])
  
  for(i in 1:length(diff_pairs)){
    
    # subset for 1 species
    d1 = d[d$Species_Pair == diff_pairs[i],]
    
    ### b/c both dom and sub species have each variable, we are going to assess each individually
    ## via a nested loop per positon 
    for(p in 1:length(unique(d1$position))){
      
      ## subset data for the relevant position
      dp = d1[d1$position == unique(d1$position)[p], ]
      
      ## if both have large Rhat values, 
      if(dp$Rhat[1] > 1.2 & dp$Rhat[2] > 1.2){
        ## update dataset with new results 
        d$difference[d$Species_Pair == diff_pairs[i] & 
                       d$position == unique(d1$position)[p]] = "different_but_irrelevant"
      }
      
      ## if both have non-significant SIVs,
      if(dp$sig[1] == "Non-Significant" & dp$sig[2] == "Non-Significant"){
        ## update dataset with new results 
        d$difference[d$Species_Pair == diff_pairs[i] & 
                       d$position == unique(d1$position)[p]] = "different_but_irrelevant"
      }
      
    } # end per position 
    
  } # end per diff_pairs
  rm(d1,dp, i)
  
  ## verify we have no NAs in differences-
  if(anyNA(d$difference)){
    print(paste("Error in code for var:", vars[v], "so loop is stopping now."))
    next
  }
  
  ## save info about differences
  dif = data.frame(table(d$difference))
  dif$var = vars[v]
  names(dif)[1] = "difference_level"
  
  ## save it in the list 
  diffs[[v]] = dif
  names(diffs)[v] = vars[v]
  
  # Add a column for significance marker
  d$marker <- ifelse(d$sig == "Significant", "*", "")
  
  # set spacing 
  dodge_width = 1
  
  ## create a new df w/ only species_pairs that have different results
  d_plot = d[d$Species_Pair %in% d$Species_Pair[d$difference == "different"] ,]

  ## make the plot
  p =
    ggplot(d_plot, aes(x = test, y = mean, fill = position))+
    geom_bar(stat = "identity", position = position_dodge(width = dodge_width), 
             color = "black", width = dodge_width) +
    geom_errorbar(aes(ymin = mean-sd, ymax = mean+sd), 
                  position = position_dodge(width = dodge_width), 
                  width = 0, color = "black") +
    scale_fill_manual(values = c("tan3","wheat3"))+
    geom_text(aes(y = mean+.3*sign(mean), label = marker), position = position_dodge(width = dodge_width), size = 8) +
    geom_text(aes(y = .04, label = round(Rhat, 1)), position = position_dodge(width = dodge_width), size = 3, vjust = -.6, hjust = -.5) +
    theme_test()+
    facet_wrap(~Species_Pair)+
    labs(y = paste(vars[v], "Effect Size"), x = NULL)+
    theme(axis.text.y = element_text(size = 16),
          axis.title.x = element_text(size = 16),
          axis.text = element_text(color = "black"),
          text = element_text(family = "Helvetica"))
  
  ## save it in the list
  plots[[v]] = p
  names(plots)[v] = vars[v]

  
} # end per var 
rm(p, dodge_width, d_plot, d, v, diff_pairs, dif)

## inspect 
diffs$FLII
48/(48+52+164) * 100 # 18% produced meaningfully different results 
plots$FLII

diffs$HFP
32/(32+62+170) * 100 # 12% produced different results 
plots$HFP

diffs$Elevation
28/(28+38+198) * 100 # 10% prodcued different results 
plots$Elevation

# grab today's date  
day<-str_sub(Sys.Date(),-2)
month<-str_sub(Sys.Date(),-5,-4)
year<-str_sub(Sys.Date(),-10,-7)
date = paste(year,month,day, sep = "")
rm(year,month,day)

# Save em!
# ggsave(paste("explore/FLII_differences_hunting_vs_no_hunting_pref_prey_HALF_LONG_", date, ".png", sep = ""), plots$FLII,
#        height = 10, width = 14, units = "in")
# ggsave(paste("explore/HFP_differences_hunting_vs_no_hunting_pref_prey_HALF_LONG_", date, ".png", sep = ""), plots$HFP,
#        height = 10, width = 14, units = "in")
# ggsave(paste("explore/ELEV_differences_hunting_vs_no_hunting_pref_prey_HALF_LONG_", date, ".png", sep = ""), plots$Elevation,
#        height = 8, width = 10, units = "in")

 ## clean up the enviro, were done here. 
rm(vars, date, plots, diffs)
