#### Compare different HPC runs of the co-abundance models
### Examine how SHORT, MEDIUM, HALF-LONG, and LONG MCMC settings influence:
## parameter preformance (effect size + SD)
## Bayes p-values + Chat values
## Abundance per SU + CI

## Zachary Amir, Z.Amir@uq.edu.au

## start fresh
rm(list = ls())

## load libs 
library(tidyverse)
library(plyr)

## set WD (should be automatic b/c RProj, but being safe anyway)
setwd("/Users/zachary_amir/Dropbox/Zach PhD/Trophic release project/SEA_trophic_cascades_co-abundance")

###### import data ######

## will create a nested list of all models per setting, 
## importing coefficents for now (SHORT is missing BPV values :( )

## first, import all short models (except missing PPCs)
files = list.files("results/SHORT_testing_Jul_2023/coefficent_dataframes/")
# store results here
coef = list()
for(i in 1:length(files)){
  path = paste("results/SHORT_testing_Jul_2023/coefficent_dataframes/", files[i], sep = "")
  d = read.csv(path)
  coef[[i]] = d
}

# save in a nested list 
short = list("coef" = coef)
rm(coef)

## next, import all medium models
file_coef = list.files("results/MIDDLE_testing_Jul_2023/coefficent_dataframes/")
file_ppc = list.files("results/MIDDLE_testing_Jul_2023/PPC_dataframes/")
file_ppc = file_ppc[! grepl("plotdata", file_ppc)] # dont want plotting data 
files = list("coef" = file_coef, "ppc" = file_ppc)
rm(file_coef, file_ppc)
# store results here
coef = list()
ppc = list()

for(i in 1:length(files)){
  
  # select coef or PPC files
  file = files[[i]]
  
  # and read them in seperately  
  if(names(files)[i] == "coef"){
    
    for(l in 1:length(file)){
      
      path = paste("results/MIDDLE_testing_Jul_2023/coefficent_dataframes/", file[l], sep = "")
      d = read.csv(path)
      coef[[l]] = d
      
    } # end per l 

  } # end coef condition
  
  if(names(files)[i] == "ppc"){
    
    for(l in 1:length(file)){
      
      path = paste("results/MIDDLE_testing_Jul_2023/PPC_dataframes/", file[l], sep = "")
      d = read.csv(path)
      ppc[[l]] = d
      
    } # end per l 
    
  } # end ppc condition  

}
rm(i,l,d,file,files)

# save in a nested list 
med = list("coef" = coef, "ppc" = ppc)


## next import all half long models 
file_coef = list.files("results/HALF_LONG_testing_Oct_2023/coefficent_dataframes/")
file_ppc = list.files("results/HALF_LONG_testing_Oct_2023/PPC_dataframes/")
file_ppc = file_ppc[! grepl("plotdata", file_ppc)] # dont want plotting data 
files = list("coef" = file_coef, "ppc" = file_ppc)
rm(file_coef, file_ppc)
# store results here
coef = list()
ppc = list()

for(i in 1:length(files)){
  
  # select coef or PPC files
  file = files[[i]]
  
  # and read them in seperately  
  if(names(files)[i] == "coef"){
    
    for(l in 1:length(file)){
      
      path = paste("results/HALF_LONG_testing_Oct_2023/coefficent_dataframes/", file[l], sep = "")
      d = read.csv(path)
      coef[[l]] = d
      
    } # end per l 
    
  } # end coef condition
  
  if(names(files)[i] == "ppc"){
    
    for(l in 1:length(file)){
      
      path = paste("results/HALF_LONG_testing_Oct_2023/PPC_dataframes/", file[l], sep = "")
      d = read.csv(path)
      ppc[[l]] = d
      
    } # end per l 
    
  } # end ppc condition  
  
}
rm(i,l,d,file,files)

# save in a nested list 
half = list("coef" = coef, "ppc" = ppc)


## last, import all full long models 
file_coef = list.files("results/LONG_testing_Sep_2023/coefficent_dataframes/")
file_ppc = list.files("results/LONG_testing_Sep_2023/PPC_dataframes/")
file_ppc = file_ppc[! grepl("plotdata", file_ppc)] # dont want plotting data 
files = list("coef" = file_coef, "ppc" = file_ppc)
rm(file_coef, file_ppc)
# store results here
coef = list()
ppc = list()

for(i in 1:length(files)){
  
  # select coef or PPC files
  file = files[[i]]
  
  # and read them in seperately  
  if(names(files)[i] == "coef"){
    
    for(l in 1:length(file)){
      
      path = paste("results/LONG_testing_Sep_2023/coefficent_dataframes/", file[l], sep = "")
      d = read.csv(path)
      coef[[l]] = d
      
    } # end per l 
    
  } # end coef condition
  
  if(names(files)[i] == "ppc"){
    
    for(l in 1:length(file)){
      
      path = paste("results/LONG_testing_Sep_2023/PPC_dataframes/", file[l], sep = "")
      d = read.csv(path)
      ppc[[l]] = d
      
    } # end per l 
    
  } # end ppc condition  
  
}
rm(i,l,d,file,files, path)

# save in a nested list 
full = list("coef" = coef, "ppc" = ppc)
rm(coef, ppc)

#
##
### combine all of the results 
res = list("SHORT" = short, "MEDIUM" = med, "HALF_LONG" = half, "FULL_LONG" = full)
rm(short, med, half, full)


###### Visualize SIV #####

## interested in tracking the effect size, standard error, and Rhat scores

## grab this information from matching species pairs! 
siv = list()
for(i in 1:length(res)){ # repeat for each MCMC setting
  
  # subset for one setting
  r = res[[i]]
  
  # make sure we are only working w/ coefficents, not PPC
  r = r[["coef"]]
  
  temp = list()
  for(l in 1:length(r)){ # repeat for each model 
    
    # subset for one dataframe/model
    dat = r[[l]]
    # only want SIV info 
    dat = dat[dat$var == "Species_Interaction",]
    # and add the MCMC setting
    dat$MCMC = names(res)[i]
    # save the results 
    temp[[l]] = dat
    names(temp)[l] = unique(dat$Species_Pair)
  }
  
  ## combine all results into one df
  df = do.call(rbind, temp)
  rownames(df) = NULL
  
  # save combined df in list
  siv[[i]] = df
  names(siv)[i] = names(res)[i]
  
}
rm(i,l,temp,dat,df,r)

## convert siv into one dataframe
siv = do.call(rbind, siv)

## remove any dog models
siv = siv[!grepl("Canis_lupus_familiaris", siv$Species_Pair),]

## set factor order of MCMC settings
siv$MCMC = factor(siv$MCMC, levels = c("SHORT","MEDIUM","HALF_LONG","FULL_LONG"))

## group based on species pairs
siv = siv %>% arrange(Species_Pair)

## add a column if the mod is top-down (preds dominant or bottom up)
siv$direction = "bottom-up"
siv$direction[grepl("DOM-Panthera", siv$Species_Pair)] = "top-down"
siv$direction[grepl("DOM-Cuon_alp", siv$Species_Pair)] = "top-down"
siv$direction[grepl("DOM-Neofelis", siv$Species_Pair)] = "top-down"


## only grab species_pairs that have full-long completed 
## (assuming they have lower levels completed too)
b = siv$Species_Pair[siv$MCMC == "FULL_LONG"]
hl = unique(b)
rm(b)

## chose a dodge width to use in the plots 
dw = 0.2

## make the plot 
p = 
ggplot(siv[siv$Species_Pair %in% hl,], aes(x = MCMC, y = mean, group = Species_Pair, color = direction))+
  geom_point(aes(size = Rhat), position=position_dodge(width=dw))+
  geom_line(position=position_dodge(width=dw))+
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd),position=position_dodge(width=dw), width = 0)+
  theme_minimal()+
  labs(title = "Species Interaction Value Preformance", y = "Mean Effect Size", x = "MCMC Settings")
# ggsave("explore/SIV_preformance_4_settings_20231013.png", p, width = 10, height = 7, units = "in")
### This is no longer a nice graph w/ 27 completed models that match here... but will save anyway 
# ggsave("explore/SIV_preformance_4_settings_20231025.png", p, width = 10, height = 7, units = "in")


## create a facet column for the large carnivore species as dom or sub --> 8 levels 
siv$group[grepl("DOM-Panthera_t", siv$Species_Pair)] = "Dom-tiger"
siv$group[grepl("DOM-Panthera_p", siv$Species_Pair)] = "Dom-leopard"
siv$group[grepl("DOM-Neofelis", siv$Species_Pair)] = "Dom-clouded_leopard"
siv$group[grepl("DOM-Cuon_al", siv$Species_Pair)] = "Dom-dhole"
siv$group[grepl("SUB-Panthera_t", siv$Species_Pair)] = "Sub-tiger"
siv$group[grepl("SUB-Panthera_p", siv$Species_Pair)] = "Sub-leopard"
siv$group[grepl("SUB-Neofelis", siv$Species_Pair)] = "Sub-clouded_leopard"
siv$group[grepl("SUB-Cuon_al", siv$Species_Pair)] = "Sub-dhole"


## try w/ a facet wrap instead
p = 
ggplot(siv[siv$Species_Pair %in% hl,], aes(x = MCMC, y = mean, group = Species_Pair, color = direction))+
  geom_point(aes(size = Rhat), position=position_dodge(width=dw))+
  geom_line(position=position_dodge(width=dw))+
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd),position=position_dodge(width=dw), width = 0)+
  facet_wrap(~group)+
  theme_test()+
  labs(title = "Species Interaction Value Preformance", y = "Mean Effect Size", x = "MCMC Settings")
## very useful to split this way!! 
## missing any dominant results for tigers and leopards! 
# ggsave("explore/SIV_preformance_4_settings_facet_by_Carnivore_20231025.png", p, width = 10, height = 7, units = "in")

## re-run on Nov 2nd --> missing dom dhole, probably last mods in named list sent to HPC 
# ggsave("explore/SIV_preformance_4_settings_facet_by_Carnivore_20231102.png", p, width = 10, height = 7, units = "in")



 ## remake for all species w/ short, med, and half-long settings
fl = siv$Species_Pair[siv$MCMC == "HALF_LONG"]

# subset the data
test = siv[siv$Species_Pair %in% fl,]
# clean up labels
test = test[test$MCMC != "FULL_LONG",]
test$MCMC = as.character(test$MCMC)
test$MCMC[test$MCMC == "HALF_LONG"] = "LONG"
test$MCMC = factor(test$MCMC, levels = c("SHORT","MEDIUM","LONG"))

## make the plot 
p=
ggplot(test, aes(x = MCMC, y = mean, group = Species_Pair, color = direction))+
  geom_jitter(aes(size = Rhat), position=position_dodge(width=dw))+
  geom_line(position=position_dodge(width=dw))+
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd),position=position_dodge(width=dw), width = 0)+
  facet_wrap(~group)+
  theme_test()+
  labs(title = "Species Interaction Value Preformance", y = "Mean Effect Size", x = "MCMC Settings")
# ggsave("explore/SIV_preformance_3_settings_and_facet_20231013.png", p, width = 14, height = 7, units = "in")
# ggsave("explore/SIV_preformance_3_settings_and_facet_per_carnivore_20231025.png", p, width = 14, height = 7, units = "in")
# ggsave("explore/SIV_preformance_3_settings_and_facet_per_carnivore_20231102.png", p, width = 14, height = 7, units = "in")

### try re-making the plot w/ facet per species_pair 
p = 
ggplot(test[test$Species_Pair %in% unique(test$Species_Pair)[1:71],], 
       aes(x = MCMC, y = mean, group = Species_Pair, color = direction))+
  geom_jitter(aes(size = Rhat), position=position_dodge(width=dw))+
  geom_line(position=position_dodge(width=dw))+
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd),position=position_dodge(width=dw), width = 0)+
  facet_wrap(~Species_Pair)+
  theme_test()+
  labs(title = "Species Interaction Value Preformance", y = "Mean Effect Size", x = "MCMC Settings")+
  theme(strip.text = element_text(size = 6.5)) # make the facet text smaller
## can only easily view via zooming off R 
ggsave("explore/SIV_preformance_3_settings_facet_per_pair_pt1_20231102.png", p, width = 20, height = 10, units = "in")

## do the same for pt2
p = 
  ggplot(test[test$Species_Pair %in% unique(test$Species_Pair)[72:142],], 
         aes(x = MCMC, y = mean, group = Species_Pair, color = direction))+
  geom_jitter(aes(size = Rhat), position=position_dodge(width=dw))+
  geom_line(position=position_dodge(width=dw))+
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd),position=position_dodge(width=dw), width = 0)+
  facet_wrap(~Species_Pair)+
  theme_test()+
  labs(title = "Species Interaction Value Preformance", y = "Mean Effect Size", x = "MCMC Settings")+
  theme(strip.text = element_text(size = 6.5)) # make the facet text smaller
## can only easily view via zooming off R 
ggsave("explore/SIV_preformance_3_settings_facet_per_pair_pt2_20231102.png", p, width = 20, height = 10, units = "in")


## keep the enviro clean! 
rm(test,hl,fl,dw,p)


#### Examine which species pairs are consistently producing poor SIVs 
check = siv[siv$mean > 3 | siv$mean < -3,]
## Visual comparison 

## save all of these problem pairs!
prob_pairs = unique(check$Species_Pair)

siv[siv$Species_Pair == "SUB-Lophura_ignita~DOM-Panthera_pardus",] # this is a crazy one, consistent across settings. 
siv[siv$Species_Pair == "SUB-Muntiacus_genus~DOM-Panthera_tigris",] # this is a good one, consistent across settings. 



###### Visualize PPC #####

## interested in BPV & Chat 

## grab this information from matching species pairs! 
ppc = list()
for(i in 1:length(res)){ # repeat for each MCMC setting
  
  # subset for one setting
  r = res[[i]]
  
  # make sure we are only working w/ coefficents, not PPC
  r = r[["ppc"]]
  
  ## skip short mods w/out PPC data 
  if(is.null(r)){
    next
  }
  
  temp = list()
  for(l in 1:length(r)){ # repeat for each model 
    
    # subset for one dataframe/model
    dat = r[[l]]
    # only want PPC info 
    dat = dat[, c("BPV.dom","BPV.sub","Chat.dom","Chat.sub","Species_Pair")]
    
    # Convert the wide data frame to long format
    long_df <- dat %>%
      pivot_longer(
        cols = -Species_Pair,
        names_to = "value.type",
        values_to = "value"
      )
    # Extract 'dom' and 'sub' from the 'value.type' column
    long_df <- long_df %>%
      separate(value.type, into = c("value.type", "sub.value"), sep = "\\.")
    # Rename the columns
    colnames(long_df) <- c("Species_Pair", "value_type", "position", "value")
    
    # and add the MCMC setting
    long_df$MCMC = names(res)[i]
    # save the results 
    temp[[l]] = long_df
    names(temp)[l] = unique(long_df$Species_Pair)
    
  }
  
  ## combine all results into one df
  df = do.call(rbind, temp)
  rownames(df) = NULL
  
  # save combined df in list
  ppc[[i]] = df
  names(ppc)[i] = names(res)[i]
  
}
rm(i,l,temp,dat,long_df,df,r)

## convert ppc into one dataframe
ppc = do.call(rbind, ppc)

## remove any dog models, from short & medium runs 
ppc = ppc[!grepl("Canis_lupus_familiaris", ppc$Species_Pair),]

## set factor order of MCMC settings
ppc$MCMC = factor(ppc$MCMC, levels = c("MEDIUM","HALF_LONG","FULL_LONG"))

## group based on species pairs
ppc = ppc %>% arrange(Species_Pair)

## add a column if the mod is top-down (preds dominant or bottom up)
ppc$direction = "bottom-up"
ppc$direction[grepl("DOM-Panthera", ppc$Species_Pair)] = "top-down"
ppc$direction[grepl("DOM-Cuon_alp", ppc$Species_Pair)] = "top-down"
ppc$direction[grepl("DOM-Neofelis", ppc$Species_Pair)] = "top-down"

## create a facet column for the large carnivore species as dom or sub --> 8 levels 
ppc$group[grepl("DOM-Panthera_t", ppc$Species_Pair)] = "Dom-tiger"
ppc$group[grepl("DOM-Panthera_p", ppc$Species_Pair)] = "Dom-leopard"
ppc$group[grepl("DOM-Neofelis", ppc$Species_Pair)] = "Dom-clouded_leopard"
ppc$group[grepl("DOM-Cuon_al", ppc$Species_Pair)] = "Dom-dhole"
ppc$group[grepl("SUB-Panthera_t", ppc$Species_Pair)] = "Sub-tiger"
ppc$group[grepl("SUB-Panthera_p", ppc$Species_Pair)] = "Sub-leopard"
ppc$group[grepl("SUB-Neofelis", ppc$Species_Pair)] = "Sub-clouded_leopard"
ppc$group[grepl("SUB-Cuon_al", ppc$Species_Pair)] = "Sub-dhole"


### Examine bpv for now 
bpv = ppc[ppc$value_type == "BPV", ]

## determine if values are good or bad fit
bpv$outcome = "bad_fit"
bpv$outcome[bpv$value >= 0.25 & bpv$value <= 0.75] = "good_fit"

## determine if each species pair is valid or not
for(i in 1:length(unique(bpv$Species_Pair))){
  
  # check one pair
  d = bpv[bpv$Species_Pair == unique(bpv$Species_Pair)[i], ]
  
  ## repeat for each MCMC setting
  for(l in 1:length(unique(d$MCMC))){
    
    dl = d[d$MCMC == unique(d$MCMC)[l], ]
    
    # if dom AND sub dont agree with each other, 
    if(length(unique(dl$outcome)) != 1){
      # mark as not valid, 
      bpv$valid[bpv$Species_Pair == unique(bpv$Species_Pair)[i] &
                  bpv$MCMC == unique(d$MCMC)[l]] = "NOT_valid_model"
      # and move on to the next
      next
    }
    
    # but if they both match and are good,
    if(unique(dl$outcome) == "good_fit"){
      # assign them as good
      bpv$valid[bpv$Species_Pair == unique(bpv$Species_Pair)[i] &
                  bpv$MCMC == unique(d$MCMC)[l]] = "valid_model"
    }
    
    # and the same for bad 
    # but if they both match and are good,
    if(unique(dl$outcome) == "bad_fit"){
      # assign them as good
      bpv$valid[bpv$Species_Pair == unique(bpv$Species_Pair)[i] &
                  bpv$MCMC == unique(d$MCMC)[l]] = "NOT_valid_model"
    }
    
  } # end per setting 
  
}
rm(i,d)

## now grab species pair, setting, and validity
dat = distinct(select(bpv, Species_Pair, MCMC, valid))

## for each setting, determine number of valid vs invalid mods 
bpv_percent = ddply(dat, .(MCMC), summarize,
                    num_good = length(valid[valid == "valid_model"]),
                    num_bad = length(valid[valid == "NOT_valid_model"]),
                    total = length(Species_Pair))
## so we can then calculate the percent
bpv_percent$perc_good = round(bpv_percent$num_good / bpv_percent$total * 100, 2)
bpv_percent$perc_bad = round(bpv_percent$num_bad / bpv_percent$total * 100, 2) 
bpv_percent # seems right! 

## restructure so its long 


### visualize BPV overall models 
p= 
ggplot(bpv, aes(y = value, x = MCMC, color = position))+
  geom_jitter(width = .1)+
  geom_text(data = bpv_percent, aes(label = paste(perc_good, "%", sep = ""), x = MCMC, y = 0.9) , position = position_dodge(width = 0.2), inherit.aes = F)+
  geom_hline(yintercept = 0.5)+
  geom_hline(yintercept = 0.75, linetype = "dashed", color = "red")+
  geom_hline(yintercept = 0.25, linetype = "dashed", color = "red")+
  theme_test()+
  labs(y = "Bayesian P-Value", x = NULL, fill = NULL)+
  coord_cartesian(ylim = c(0,1))
  # facet_wrap(~group, nrow = 4, ncol = 4)

## save it! 
# ggsave("explore/BPV_preformance_3_settings_overall_20231102.png", p, width = 8, height = 6, units = "in")

## now examine per grouping 
p= 
  ggplot(bpv, aes(y = value, x = MCMC, color = position))+
  geom_jitter(width = .1)+
  # geom_text(data = bpv_percent, aes(label = paste(perc_good, "%", sep = ""), x = MCMC, y = 0.9) , position = position_dodge(width = 0.2), inherit.aes = F)+
  geom_hline(yintercept = 0.5)+
  geom_hline(yintercept = 0.75, linetype = "dashed", color = "red")+
  geom_hline(yintercept = 0.25, linetype = "dashed", color = "red")+
  theme_test()+
  labs(y = "Bayesian P-Value", x = NULL, fill = NULL)+
  coord_cartesian(ylim = c(0,1))+
  facet_wrap(~group, nrow = 4, ncol = 4)

## save it! 
# ggsave("explore/BPV_preformance_3_settings_facet_per_group_20231102.png", p, width = 8, height = 6, units = "in")


## thin to species that made it thru half long
hl = unique(bpv$Species_Pair[bpv$MCMC == "HALF_LONG"])
bpv_2 = bpv[bpv$Species_Pair %in% hl, ]

## now examine per species pair, pt1
p=
  ggplot(bpv_2[bpv_2$Species_Pair %in% unique(bpv_2$Species_Pair)[1:71], ], 
         aes(y = value, x = MCMC, fill = position))+
  geom_col(color = "black", position = "dodge", color = 'black')+
  geom_hline(yintercept = 0.5)+
  geom_hline(yintercept = 0.75, linetype = "dashed", color = "red")+
  geom_hline(yintercept = 0.25, linetype = "dashed", color = "red")+
  theme_test()+
  labs(y = "Bayesian P-Value", x = NULL, fill = NULL)+
  coord_cartesian(ylim = c(0,1))+
  scale_fill_manual(values = c("tan3","wheat3"))+
  facet_wrap(~Species_Pair)+
  theme(strip.text = element_text(size = 6.5)) # make the facet text smaller
## can only easily view via zooming off R 
ggsave("explore/BPV_preformance_3_settings_facet_per_pair_pt1_20231102.png", p, width = 21, height = 10, units = "in")

## and pt2
p=
  ggplot(bpv_2[bpv_2$Species_Pair %in% unique(bpv_2$Species_Pair)[72:142], ], 
         aes(y = value, x = MCMC, fill = position))+
  geom_col(color = "black", position = "dodge", color = 'black')+
  geom_hline(yintercept = 0.5)+
  geom_hline(yintercept = 0.75, linetype = "dashed", color = "red")+
  geom_hline(yintercept = 0.25, linetype = "dashed", color = "red")+
  theme_test()+
  labs(y = "Bayesian P-Value", x = NULL, fill = NULL)+
  coord_cartesian(ylim = c(0,1))+
  scale_fill_manual(values = c("tan3","wheat3"))+
  facet_wrap(~Species_Pair)+
  theme(strip.text = element_text(size = 6.5)) # make the facet text smaller
## can only easily view via zooming off R 
ggsave("explore/BPV_preformance_3_settings_facet_per_pair_pt2_20231102.png", p, width = 21, height = 10, units = "in")


## keep it clean! 
rm(p, l, dl,hl, dat, bpv, bpv_2, bpv_percent)



##### Now do the same for Chat! 
chat = ppc[ppc$value_type == "Chat", ]

## determine if values are good or bad fit
chat$outcome = "bad_fit"
chat$outcome[chat$value >= 0.98 & chat$value <= 1.1] = "good_fit"

## determine if each species pair is valid or not
for(i in 1:length(unique(chat$Species_Pair))){
  
  # check one pair
  d = chat[chat$Species_Pair == unique(chat$Species_Pair)[i], ]
  
  ## repeat for each MCMC setting
  for(l in 1:length(unique(d$MCMC))){
    
    dl = d[d$MCMC == unique(d$MCMC)[l], ]
    
    # if dom AND sub dont agree with each other, 
    if(length(unique(dl$outcome)) != 1){
      # mark as not valid, 
      chat$valid[chat$Species_Pair == unique(chat$Species_Pair)[i] &
                  chat$MCMC == unique(d$MCMC)[l]] = "NOT_valid_model"
      # and move on to the next
      next
    }
    
    # but if they both match and are good,
    if(unique(dl$outcome) == "good_fit"){
      # assign them as good
      chat$valid[chat$Species_Pair == unique(chat$Species_Pair)[i] &
                  chat$MCMC == unique(d$MCMC)[l]] = "valid_model"
    }
    
    # and the same for bad 
    # but if they both match and are good,
    if(unique(dl$outcome) == "bad_fit"){
      # assign them as good
      chat$valid[chat$Species_Pair == unique(chat$Species_Pair)[i] &
                  chat$MCMC == unique(d$MCMC)[l]] = "NOT_valid_model"
    }
    
  } # end per setting 
  
}
rm(i,d)

## now grab species pair, setting, and validity
dat = distinct(select(chat, Species_Pair, MCMC, valid))

## for each setting, determine number of valid vs invalid mods 
chat_percent = ddply(dat, .(MCMC), summarize,
                    num_good = length(valid[valid == "valid_model"]),
                    num_bad = length(valid[valid == "NOT_valid_model"]),
                    total = length(Species_Pair))
## so we can then calculate the percent
chat_percent$perc_good = round(chat_percent$num_good / chat_percent$total * 100, 2)
chat_percent$perc_bad = round(chat_percent$num_bad / chat_percent$total * 100, 2) 
chat_percent # seems right! 

## restructure so its long 


### visualize chat overall models 
p= 
  ggplot(chat, aes(y = value, x = MCMC, color = position))+
  geom_jitter(width = .1)+
  geom_text(data = chat_percent, aes(label = paste(perc_good, "%", sep = ""), x = MCMC, y = 0.9) , position = position_dodge(width = 0.2), inherit.aes = F)+
  # geom_hline(yintercept = 1)+
  geom_hline(yintercept = 1.1, linetype = "dashed", color = "red")+
  theme_test()+
  labs(y = "Magnitude of Over-Dispersion (C-hat)", x = NULL, fill = NULL)+ 
  coord_cartesian(ylim = c(0.9,1.2))
  # facet_wrap(~group, nrow = 4, ncol = 4)

## save it! 
# ggsave("explore/Chat_preformance_3_settings_overall_20231102.png", p, width = 8, height = 6, units = "in")

## now examine per grouping 
p= 
  ggplot(chat, aes(y = value, x = MCMC, color = position))+
  geom_jitter(width = .1)+
  # geom_text(data = chat_percent, aes(label = paste(perc_good, "%", sep = ""), x = MCMC, y = 0.9) , position = position_dodge(width = 0.2), inherit.aes = F)+
  # geom_hline(yintercept = 1)+
  geom_hline(yintercept = 1.1, linetype = "dashed", color = "red")+
  theme_test()+
  labs(y = "Magnitude of Over-Dispersion (C-hat)", x = NULL, fill = NULL)+ 
  coord_cartesian(ylim = c(0.9,1.2))+
  facet_wrap(~group, nrow = 4, ncol = 4)

## save it! 
# ggsave("explore/Chat_preformance_3_settings_facet_per_group_20231102.png", p, width = 8, height = 6, units = "in")


## thin to species that made it thru half long
hl = unique(chat$Species_Pair[chat$MCMC == "HALF_LONG"])
chat_2 = chat[chat$Species_Pair %in% hl, ]

## now examine per species pair, pt1
p=
  ggplot(chat_2[chat_2$Species_Pair %in% unique(chat_2$Species_Pair)[1:71], ], 
         aes(y = value, x = MCMC, fill = position))+
  geom_col(color = "black", position = "dodge", color = 'black')+
  geom_hline(yintercept = 1.1, linetype = "dashed", color = "red")+
  theme_test()+
  labs(y = "Magnitude of Over-Dispersion (C-hat)", x = NULL, fill = NULL)+ 
  coord_cartesian(ylim = c(0.9,1.2))+
  scale_fill_manual(values = c("tan3","wheat3"))+
  facet_wrap(~Species_Pair)+
  theme(strip.text = element_text(size = 6.5)) # make the facet text smaller
## can only easily view via zooming off R 
# ggsave("explore/chat_preformance_3_settings_facet_per_pair_pt1_20231102.png", p, width = 21, height = 10, units = "in")

## and pt2
p=
  ggplot(chat_2[chat_2$Species_Pair %in% unique(chat_2$Species_Pair)[72:142], ], 
         aes(y = value, x = MCMC, fill = position))+
  geom_col(color = "black", position = "dodge", color = 'black')+
  geom_hline(yintercept = 1.1, linetype = "dashed", color = "red")+
  theme_test()+
  labs(y = "Magnitude of Over-Dispersion (C-hat)", x = NULL, fill = NULL)+ 
  coord_cartesian(ylim = c(0.9,1.2))+
  scale_fill_manual(values = c("tan3","wheat3"))+
  facet_wrap(~Species_Pair)+
  theme(strip.text = element_text(size = 6.5))  # make the facet text smaller
## can only easily view via zooming off R 
# ggsave("explore/chat_preformance_3_settings_facet_per_pair_pt2_20231102.png", p, width = 21, height = 10, units = "in")





######## Inspect problematic species pairs ######

### import the captures to check back w/ OG data
caps = read.csv("data/send_to_HPC/clean_captures_to_make_UMFs_20230821.csv")
## need to bind landscapes here too, import the meta
meta = read.csv("data/send_to_HPC/clean_metadata_to_make_UMFs_20230821.csv")
add = select(meta, cell_id_3km, Landscape)
caps = merge(caps, add, by = "cell_id_3km")
rm(add)

#### Examine which species pairs are consistently producing poor SIVs 
check = siv[siv$mean > 3 | siv$mean < -3,]
## save all of these problem pairs!
prob_pairs = unique(check$Species_Pair)

## inspect pair of good vs bad pairs for testing. 
siv[siv$Species_Pair == "SUB-Lophura_ignita~DOM-Panthera_pardus",] # this is a crazy one, consistent across settings. 
siv[siv$Species_Pair == "SUB-Muntiacus_genus~DOM-Panthera_tigris",] # this is a good one, consistent across settings. 


## I have a suspicion that most of these problematic pairs have a small amount of spaital overlap, 
## thus producing crazy estimates. 
### BUT, best to check with all species pairs and compare 

## store results here
land_inspection = data.frame(matrix(NA, ncol = 8, nrow = length(unique(siv$Species_Pair))))
names(land_inspection) = c("Species_Pair", "dom_species", "sub_species", 
                           "num_intersecting_landscapes", 
                           "dom_det_intersect", "sub_det_intersect",
                           "dom_indiv_intersect", "sub_indiv_intersect")
land_inspection$Species_Pair = unique(siv$Species_Pair)

## check every single pair! 
for(i in 1:length(unique(siv$Species_Pair))){
  
  ## split apart pairs
  pair = str_split(unique(siv$Species_Pair)[i], "~")[[1]]
  dom_sp = pair[grepl("DOM", pair)]
  dom_sp = gsub("DOM-", "", dom_sp) # for easier indexing
  sub_sp = pair[grepl("SUB", pair)]
  sub_sp = gsub("SUB-", "", sub_sp) # for easier indexing
  
  ## grab dominant landscapes
  dom_land = unique(caps$Landscape[caps$Species == dom_sp])
  ## grab subordinate landscapes
  sub_land = unique(caps$Landscape[caps$Species == sub_sp])
  
  ## grab only the matching landscapes
  lands = intersect(dom_land, sub_land)
  
  ## how many detections for each species across all shared landscapes?
  dom_det = nrow(caps[caps$Species == dom_sp &
                        caps$Landscape %in% lands,])
  sub_det = nrow(caps[caps$Species == sub_sp &
                        caps$Landscape %in% lands,])
  
  ## how many individuals for each species across all shared landscapes?
  dom_indiv = sum(caps$total_indiv_records[caps$Species == dom_sp &
                        caps$Landscape %in% lands])
  sub_indiv = sum(caps$total_indiv_records[caps$Species == sub_sp &
                        caps$Landscape %in% lands])

  ## save all this info 
  land_inspection$dom_species[land_inspection$Species_Pair == unique(siv$Species_Pair)[i]] = dom_sp
  land_inspection$sub_species[land_inspection$Species_Pair == unique(siv$Species_Pair)[i]] = sub_sp
  land_inspection$num_intersecting_landscapes[land_inspection$Species_Pair == unique(siv$Species_Pair)[i]] = length(lands)
  land_inspection$dom_det_intersect[land_inspection$Species_Pair == unique(siv$Species_Pair)[i]] = dom_det
  land_inspection$sub_det_intersect[land_inspection$Species_Pair == unique(siv$Species_Pair)[i]] = sub_det
  land_inspection$dom_indiv_intersect[land_inspection$Species_Pair == unique(siv$Species_Pair)[i]] = dom_indiv
  land_inspection$sub_indiv_intersect[land_inspection$Species_Pair == unique(siv$Species_Pair)[i]] = sub_indiv
  
}
rm(i, pair, dom_sp, dom_land, dom_det, dom_indiv, 
   sub_indiv, sub_det, sub_land, sub_sp, lands)

## account for crazy mods
land_inspection$crazy_mod = F
land_inspection$crazy_mod[land_inspection$Species_Pair %in% prob_pairs] = T

## remove rows w/ non relevant species
land_inspection = land_inspection[! land_inspection$sub_species %in% c("Bos_taurus", "Macaca_genus",
                                                                       "Hystrix_genus", "Presbytis_rubicunda",
                                                                       "Herpestes_genus", "Rheithrosciurus_macrotis",
                                                                       "Macaca_genus"),]
land_inspection = land_inspection[! land_inspection$dom_species %in% c("Bos_taurus", "Macaca_genus",
                                                                       "Hystrix_genus", "Presbytis_rubicunda",
                                                                       "Herpestes_genus", "Rheithrosciurus_macrotis",
                                                                       "Macaca_genus"),]
## total number of detections-
land_inspection$total_det = land_inspection$dom_det_intersect + land_inspection$sub_det_intersect
# and individuals! 
land_inspection$total_indiv = land_inspection$dom_indiv_intersect + land_inspection$sub_indiv_intersect

## remove data w/ no detections
land_inspection = land_inspection[land_inspection$total_det != 0, ] # for removing old data for non relevant species

## inspect detections
summary(land_inspection$total_det)
summary(land_inspection$total_det[land_inspection$crazy_mod == T]) # much lower values 

## insepct individuals
summary(land_inspection$total_indiv)
summary(land_inspection$total_indiv[land_inspection$crazy_mod == T]) # higher values than detections (obviously)


## Inspect mods that should have failed... how did these work?!
land_inspection[land_inspection$Species_Pair == "SUB-Cuon_alpinus~DOM-Sus_barbatus",]
siv[siv$Species_Pair == "SUB-Cuon_alpinus~DOM-Sus_barbatus",] # this should have failed! Not significant tho. Maybe ample sample size of dom saved the day. 
ppc[ppc$Species_Pair == "SUB-Cuon_alpinus~DOM-Sus_barbatus",] # Amazed that this passed the thresholds too. 
# check the inverse here
land_inspection[land_inspection$Species_Pair == "SUB-Sus_barbatus~DOM-Cuon_alpinus",]
siv[siv$Species_Pair == "SUB-Sus_barbatus~DOM-Cuon_alpinus",] # parameter did not converge --> not enough dominant individuals 
ppc[ppc$Species_Pair == "SUB-Sus_barbatus~DOM-Cuon_alpinus",] # good fit tho... go figure. 

land_inspection[land_inspection$Species_Pair == "SUB-Lophura_diardi~DOM-Panthera_pardus",] 
siv[siv$Species_Pair == "SUB-Lophura_diardi~DOM-Panthera_pardus",] # Huuuge SE --> not good. 
ppc[ppc$Species_Pair == "SUB-Lophura_diardi~DOM-Panthera_pardus",] # would have failed here. 

land_inspection[land_inspection$Species_Pair == "SUB-Viverra_tangalunga~DOM-Panthera_tigris",]
siv[siv$Species_Pair == "SUB-Viverra_tangalunga~DOM-Panthera_tigris",] # large SE
ppc[ppc$Species_Pair == "SUB-Viverra_tangalunga~DOM-Panthera_tigris",] # Its a decent fit! 


## inspect various conditions in the data 
select(land_inspection[land_inspection$total_det < 100,], Species_Pair, dom_det_intersect, sub_det_intersect, crazy_mod)
# at both 100 and 200, there are still true and false crazy mods, tho more true the lower you go
#  make 103 shared detection the relevant threshold --> essentially 100, but captures crazy mongoose-tiger mod. 

nrow(land_inspection[land_inspection$total_det < 100,]) # 20 mods
nrow(land_inspection[land_inspection$total_indiv < 100,]) # 18 mods
## probs better to be more inclusive and use total_det


## grab all pairs w/ less than 100 shared detections
low_det = land_inspection[land_inspection$total_det < 100,]
## but check which low_det are good
low_det_good = low_det[low_det$crazy_mod == F,]
low_det_good # no preferred prey. Check SIV and BPV next 

siv[siv$Species_Pair %in% low_det_good$Species_Pair,]
# no pairs are significant for either hypothesis. 
## SUB-Hemigalus_derbyanus~DOM-Cuon_alpinus is significant and a good fit, but positive relationship --> unsupportive of top-down hypothesis. 
# most have fairly high SEs as well. 

as.data.frame(ppc[ppc$Species_Pair %in% low_det_good$Species_Pair, ])
## most, but not all, do not have a good model fit. 

## now check the opposite, enough detections but still crazy mods
high_det_bad = land_inspection[land_inspection$total_det >= 100 & land_inspection$crazy_mod == T, ]
high_det_bad # tapir and dhole should be considered, but they are not a preferred pred-prey pair. 

siv[siv$Species_Pair %in% high_det_bad$Species_Pair, ] # the only mods w/out crazy estimates have failed parameters
ppc[ppc$Species_Pair %in% high_det_bad$Species_Pair, ] # surprisingly, sub-tapir, dom-dhole have a good fit! Doesnt matter tho, parameter didnt converge. 

## inspect SUB-Panthera_pardus~DOM-Paradoxurus_hermaphroditus b/c it became signficant and converged and (almost) a good fit. 
table(caps$Landscape[caps$Species == "Panthera_pardus"])
table(caps$Landscape[caps$Species == "Paradoxurus_hermaphroditus" & 
                       caps$Landscape %in% unique(caps$Landscape[caps$Species == "Panthera_pardus"])])
## quite complementary! 
## I say we leave it for now and see how we go. 

## andother inspection: SUB-Hemigalus_derbyanus~DOM-Cuon_alpinus --> should have failed but didnt and is a good fit. 
table(caps$Landscape[caps$Species == "Cuon_alpinus" & caps$Landscape %in% unique(caps$Landscape[caps$Species == "Hemigalus_derbyanus"])])
table(caps$Landscape[caps$Species == "Hemigalus_derbyanus" & 
                       caps$Landscape %in% unique(caps$Landscape[caps$Species == "Cuon_alpinus"])])
## very comparable for Kerinci seblat. 


###### MOVING FORWARD-->
## Only generate species pairs if there are at least 100 detections for both species across landscapes where both species are detected! 

summary(siv$mean[siv$Species_Pair %in% land_inspection$Species_Pair[land_inspection$total_det > 100] &
                   siv$MCMC == "HALF_LONG"]) ## this is a good range! but these mods are only half done! 

## end it on a clean note. 
rm(check, high_det_bad, land_inspection, low_det, low_det_good, prob_pairs)

