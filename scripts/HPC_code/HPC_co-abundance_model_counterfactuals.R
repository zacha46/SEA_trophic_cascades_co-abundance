###### Zachary Amir's co-abundance model, Z.Amir@uq.edu.au
##### Using 1) count matricies, 2) site covarites, & 3) observation covariates. 
#### Currently examining 260 pairwise interactions between 44 species based on trophic guilds

### Data has already been formatted to into the proper data bundle to run the co-abundance model on the HPC 
## But because models are too large to save ALL of them, 
# this script will extract only the relevant dataframes. 

### September 2023 Update-
## Will no longer generate prediction plots or monitor site-level abundance
# in an effort to save RAM and complete LONG settings. 

## Submitted to HPC on December 10th, 2024
# Zachary Amir, z.amir@uq.edu.au

####### Set up #####

library(jagsUI)
library(tidyverse)

## specify the number of cores to be uses, should be the same as the model 
options(mc.cores = 3)

#### Read Job Array index value into R
slurm = Sys.getenv("SLURM_ARRAY_TASK_ID")
slurm = as.numeric(slurm) #imports as character var, not numeric

#### Read external values from SLURM into R
setting = Sys.getenv("SETTING")  # MCMC setting
# pref = Sys.getenv("PREF")        # preferred prey models or not
# gb = Sys.getenv("GB")            # RAM

#### Also read in counterfactual settings! 
counter = Sys.getenv("COUNTER")


#### List all possible bundled data files 
files = list.files("data/co-abundance")[grepl("Bundled_data", list.files("data/co-abundance/"))]

# #for local testing
# files = list.files("data_GitHub_CoA_bundles/counterfactual_testing/")[grepl("Bundled_data", list.files("data_GitHub_CoA_bundles/counterfactual_testing/"))]

# ## Subset files for proper pref setting
# files = files[grepl(pref, files)]
# ## and for proper GB setting
# f = files[grepl(gb, files)]

## Subset files for proper counterfactual test 
f = files[grepl(counter, files)]

### Import the formatted data
dat = readRDS(paste("data/co-abundance/", f ,sep = ""))

# #for local testing
# dat = readRDS(paste("data_GitHub_CoA_bundles/counterfactual_testing/", f ,sep = ""))



## Thin to a single species pair
bdata = dat[[slurm]]

## and save the name of the species pair 
n = names(dat)[slurm]

### Verify that this species pair at this setting has not already been completed-
# First, create file name that mimics results to check if already present
res =  paste(setting, "_", counter, "_co-abundance_coefficents_", n, sep = "")

# Second, list all completed results
res_search = list.files("results/co-abundance/dec2024_counterfactual_coefficent_dataframes/")

# If the newly constructed file name matches ANY values already present in results,
if(any(grepl(res, res_search))){

  ## Print a message stating so
  print(paste("The species pair:", n, "run with MCMC settings:", setting, 
              "and counterfactual test:", counter, 
              "is already present in the results folder. This R script is terminating now."))

  ## and terminate the R script fully
  stop("The script was terminated.")
}

## matricies saved as char but need to be numeric 
bdata$y.dom = apply(bdata$y.dom, 2, as.numeric)
bdata$y.sub = apply(bdata$y.sub, 2, as.numeric)
bdata$cams = apply(bdata$cams, 2, as.numeric)

## new variable got saved as a 2D matrix, but needs to be a vector, so convert
bdata$comm_det = as.vector(bdata$comm_det)

## if results are not already present, start the models! 
print(paste("Begining to run co-abundance model for: ", n, " with MCMC settings: ", setting ,
            " and is starting at ", Sys.time(), sep = ""))


####### Make the model in BUGS language and run it ####

cat(file = "ZDA_Co_Abundance_Model_final_20241201.jags", 
    
    "model{

  # Priors for both species
  for (i in 1:2) {
      
    # Intercept
    a0[i] ~ dnorm(0, 0.01)
      
    # FLII
    a1[i] ~ dnorm(0, 0.01)
      
    # HFP
    a2[i] ~ dnorm(0, 0.01)
    
    # Elev
    a3[i] ~ dnorm(0, 0.01)
    
    # community detections
    a4[i] ~ dnorm(0, 0.01)
  
  
    # Det intercept
    b0[i] ~ dnorm(0, 0.01)
    
    # Cams
    b2[i] ~ dnorm(0, 0.01)
      
    ## OD params for observation model
    tau.p[i] <- pow(sd.p[i], -2) 
    sd.p[i] ~ dunif(0, 1)  #not sure how to define variance here... leaving it the same as K&R
      
    }
    
  ## Species interaction
  a5 ~ dnorm(0, 0.01)
    
  # Landscape RE hyper prior --> define it's variance
  sigma.a6 ~ dunif(0,5)
  var.a6 <- 1/(sigma.a6*sigma.a6) 
    
  sigma.a7 ~ dunif(0,5)
  var.a7 <- 1/(sigma.a7*sigma.a7)
    
  for (k in 1:narea) {
      
    a6[k] ~ dnorm(0,var.a6)
    a7[k] ~ dnorm(0,var.a7)
      
  }
  
  # year RE hyper prior --> define it's variance
  sigma.a8 ~ dunif(0,5)
  var.a8 <- 1/(sigma.a8*sigma.a8)

  sigma.a9 ~ dunif(0,5)
  var.a9 <- 1/(sigma.a9*sigma.a9)

  for (k in 1:nyear) {

    a8[k] ~ dnorm(0,var.a8)
    a9[k] ~ dnorm(0,var.a9)

  }
  
    # source RE hyper prior --> define it's variance
  sigma.b3 ~ dunif(0,5)
  var.b3 <- 1/(sigma.b3*sigma.b3)

  sigma.b4 ~ dunif(0,5)
  var.b4 <- 1/(sigma.b4*sigma.b4)

  for (k in 1:nsource) {

    b3[k] ~ dnorm(0,var.b3)
    b4[k] ~ dnorm(0,var.b4)

  }
    
    
  # Likelihood
  # Ecological model for true abundance per site
  for (j in 1:nsites) {
    
    # Abundance of Subordinate Species w/ iZIP
    N.sub[j] ~ dpois(lambda.sub[j] * Z.sub[j])

      log(lambda.sub[j]) <- a0[1] + a1[1]*flii[j] + a2[1]*hfp[j] + a3[1]*elev[j] + a4[1]*comm_det[j] + a5*N.dom[j] + a6[area[j]] + a8[year[j]]
    
    # Abundance of Dominant Species w/ iZIP 
    N.dom[j] ~ dpois(lambda.dom[j] * Z.dom[j])
    
      log(lambda.dom[j]) <- a0[2] + a1[2]*flii[j] + a2[2]*hfp[j] + a3[2]*elev[j] + a4[2]*comm_det[j] + a7[area[j]] + a9[year[j]]
                      
                      
    # Observation model for counts per replicated observation with OD params 
    for (k in 1:nreps) {
    
      ## Subordinate species
      y.sub[j,k] ~ dbin(p.sub[j,k], N.sub[j])
      
        p.sub[j,k] <- 1 / (1 + exp(-lp.lim.sub[j,k]))
        
          lp.lim.sub[j,k]<- min(250, max(-250, lp.sub[j,k])) #stabilize logit
          
            lp.sub[j,k] <- b0[1] +  b2[1]*cams[j,k] + b3[source[j]] + eps.p.sub[j,k]
        
              eps.p.sub[j,k] ~ dnorm(0, tau.p[1])
          
      
      ## Dominant species
      y.dom[j,k] ~ dbin(p.dom[j,k], N.dom[j])
      
        p.dom[j,k] <- 1 / (1 + exp(-lp.lim.dom[j,k]))
        
          lp.lim.dom[j,k]<- min(250, max(-250, lp.dom[j,k])) #stabilize logit
          
            lp.dom[j,k] <- b0[2] + b2[2]*cams[j,k] + b4[source[j]] + eps.p.dom[j,k]
        
              eps.p.dom[j,k] ~ dnorm(0, tau.p[2])
         
         
    ### PPC- Subordinate
      
    ## Expected count at site j, occasion k
    exp.sub[j,k]<- N.sub[j] * p.sub[j,k]
      
    ## Discrepency from Real Data
    # small denominator added to prevent dividing by zero
    E.sub[j,k] <- pow((y.sub[j,k] - exp.sub[j,k]), 2) / (exp.sub[j,k] + 0.5) 
      
    ## Simulate new sub counts from model 
    y.rep.sub[j,k] ~ dbin(p.sub[j,k], N.sub[j])
      
    ## Discrepency from Simulation
    E.rep.sub[j,k] <- pow((y.rep.sub[j,k] - exp.sub[j,k]), 2) / (exp.sub[j,k] + 0.5)
      
      
      
    ### PPC- Dominant
      
    ## Expected count at site j, occasion k
    exp.dom[j,k]<- N.dom[j] * p.dom[j,k]
      
    ## Discrepency from Real Data
    # small denominator added to prevent dividing by zero
    E.dom[j,k] <- pow((y.dom[j,k] - exp.dom[j,k]), 2) / (exp.dom[j,k] + 0.5) 
      
    ## Simulate new dom counts from model 
    y.rep.dom[j,k] ~ dbin(p.dom[j,k], N.dom[j])
      
    ## Discrepency from Simulation
    E.rep.dom[j,k] <- pow((y.rep.dom[j,k] - exp.dom[j,k]), 2) / (exp.dom[j,k] + 0.5)
      

    } #k
  } #j

  ## Derived Parameters
  
  #Chi-Square Test Statistic- Subordinate
  fit.sub = sum(E.sub[,])
  fit.rep.sub = sum(E.rep.sub[,])

  #Chi-Square Test Statistic- Dominant
  fit.dom = sum(E.dom[,])
  fit.rep.dom = sum(E.rep.dom[,])

} 
    ",fill = TRUE)


### Specify the parameters to be monitored.
params = c('a0', 'a1', 'a2', 'a3', 'a4', 'a5',    # Abundance parameters
           'b0',  'b2',                           # Detection parameters
           'var.a6','var.a7','a6', 'a7',          # Landscape random effects
           'var.a8','var.a9','a8', 'a9',          # Year random effects
           'var.b3','var.b4','b3', 'b4',          # source random effects
           "tau.p","eps.p.dom","eps.p.sub",       # OD params
           "fit.sub", "fit.rep.sub",              # Chi2 stat for sub, real then simulated
           "fit.dom", "fit.rep.dom"               # Chi2 stat for dom, real then simulated
           # "E.dom", "E.rep.dom",
           # "E.sub", "E.rep.sub",
           # "N.dom", "N.sub"                     # NOT monitoring abundance per site to save RAM on LONG mods. 
)                         


# Specify the initial values
inits = function() {
  list(a0=rnorm(2), a1=rnorm(2), a2=rnorm(2), a3=rnorm(2), a4=rnorm(2), a5=rnorm(1),
       b0=rnorm(2), b2=rnorm(2),
       sd.p = runif(2, .4, .8), # Experimented and medium values seem to produce better convergence.
       a6=rnorm(bdata$narea), a7=rnorm(bdata$narea),
       a8=rnorm(bdata$nyear), a9=rnorm(bdata$nyear),
       b3=rnorm(bdata$nsource), b4=rnorm(bdata$nsource),
       sigma.a6= runif(1, 3, 5), sigma.a7= runif(1, 3, 5), #crashes w/ inits >= 6, and has better convergence > 2
       sigma.a8= runif(1, 3, 5), sigma.a9= runif(1, 3, 5), #havent explored these values at all.
       sigma.b3= runif(1, 3, 5), sigma.b4= runif(1, 3, 5), #havent explored these values at all.
       N.sub = as.vector(apply(bdata$y.sub ,1,max, na.rm=T)),   
       N.dom = as.vector(apply(bdata$y.dom ,1,max, na.rm=T)),
       Z.dom = as.vector(rep(as.numeric(NA), length.out = length(bdata$Z.dom))),
       Z.sub = as.vector(rep(as.numeric(NA), length.out = length(bdata$Z.sub))) 
  )  # Royle recommends leaving Z as NA, and JAGS says it needs to be numeric          
}

# MCMC settings, based on assignment above
## Want burn-in to be ~20% of iterations and then thin = (ni - nb) / ideal n.eff (per chain), ideally 30000 in the long one. 
### Assess n.eff via (ni - nb)/nt * nc 
if(setting == "SHORT"){
  ni <- 3000;  nt <- 5; nb <- 60; nc <- 3; na = NULL      #quick test to make sure code works, 2.5 hr per mod
}
if(setting == "MIDDLE"){
  ni = 50000;  nt = 20; nb = 10000 ; nc <- 3; na = NULL   #examine parameter values --> use this for prelim testing. 36-49 hrs per mod
}
if(setting == "LONG"){
  ni = 250000;  nt = 20; nb = 50000 ; nc <- 3; na = NULL  #publication quality run --> ~160 hours per mod, all finish in < 2 weeks! 
}

# take the start time 
start = Sys.time()

### Run the model 
mod = jags(bdata, inits, params, "ZDA_Co_Abundance_Model_final_20241201.jags",
           ## MCMC settings
           n.chains = nc, n.adapt = na, n.thin = nt, 
           n.iter = ni, n.burnin = nb, parallel = T)

# take the end time 
end = Sys.time()


### DONT SAVE THE MODEL, they are too big 

print(paste("Finished running co-abundance model for: ", n, " at ", Sys.time(),
            ". This model took ", round(difftime(end, start, units = "hours"), 4), " hours to be completed. Beginning dataframe extractions now.", sep = ""))



####### Generate dataframe for coefficent bar plots ######

# Extract summary table from your model 
s = data.frame(mod$summary[,c(1:3,7:10)]) ## Only need relevant parameters

## move row names to new column var to begin re-labelling them. 
s$var = row.names(s)

## remove variables we dont need here to be concise 
s = s[!grepl("E", s$var),]
s = s[!grepl("eps", s$var),] 
s = s[!grepl("fit", s$var),] # pretty sure this is all PPC stuff anyway. 
s = s[!grepl("N.", s$var),] # dont need b/c its estimated later

# Add sub species name
s$species[grepl("\\[1]", s$var)] = str_split(n, "~")[[1]][1]  
s$species[grepl("a6", s$var)] = str_split(n, "~")[[1]][1] #landscape RE
s$species[grepl("a8", s$var)] = str_split(n, "~")[[1]][1] #year RE
s$species[grepl("b3", s$var)] = str_split(n, "~")[[1]][1] #source RE
s$species[grepl("sub", s$var)] = str_split(n, "~")[[1]][1]  
s$species[s$var == "a5"] = str_split(n, "~")[[1]][1] # interaction 

# Add dom species name
s$species[grepl("\\[2]", s$var)] = str_split(n, "~")[[1]][2]
s$species[grepl("a7", s$var)] = str_split(n, "~")[[1]][2] # landscape RE
s$species[grepl("a9", s$var)] = str_split(n, "~")[[1]][2] # year RE
s$species[grepl("b4", s$var)] = str_split(n, "~")[[1]][1] #source RE
s$species[grepl("dom", s$var)] = str_split(n, "~")[[1]][2]

# Clean up variables --> translate BUGS code to english! 
s$var[grepl("a0", s$var)]= "Abundance_intercept"
s$var[grepl("a1", s$var)]= "FLII"
s$var[grepl("a2", s$var)]= "HFP"
s$var[grepl("a3", s$var)]= "Elevation"
s$var[grepl("a4", s$var)]= "Community_detections"
s$var[grepl("a5", s$var)]= "Species_Interaction"

s$var[grepl("b0", s$var)]= "Detection_intercept"
s$var[grepl("b2", s$var)]= "Active_cams"

s$sig[s$overlap0 == 0] = "Significant"
s$sig[s$overlap0 == 1] = "Non-Significant"
#dont need overlap0 anymore
s$overlap0 = NULL

## REs
s$var[grepl("a6", s$var)] = gsub("a6", "Landscape", s$var[grepl("a6", s$var)])
s$var[grepl("a7", s$var)] = gsub("a7", "Landscape", s$var[grepl("a7", s$var)])
s$var[grepl("a8", s$var)] = gsub("a8", "Year", s$var[grepl("a8", s$var)])
s$var[grepl("a9", s$var)] = gsub("a9", "Year", s$var[grepl("a9", s$var)])
s$var[grepl("b3", s$var)] = gsub("b3", "Source", s$var[grepl("b3", s$var)])
s$var[grepl("b4", s$var)] = gsub("b4", "Source", s$var[grepl("b4", s$var)])

# Clean up col names
colnames(s)[c(3,4)] = c("lower", "upper")
# and remove rownames
rownames(s) = NULL

# Add the species pair
s$Species_Pair = n

## save it! 
day<-substr(Sys.Date(),9, 10)
month<-substr(Sys.Date(),6,7)
year<-substr(Sys.Date(),1,4)

path = paste(paste(paste(paste("results/co-abundance/dec2024_counterfactual_coefficent_dataframes/", slurm, "_", setting, "_", counter, "_co-abundance_coefficents_", n, "_", year,sep=""),month,sep=""),day,sep=""),".csv",sep="")
write.csv(s, path, row.names = F)


## Give us an update please! 
print(paste("Finished generating coefficent dataframe for: ", n, " at ", Sys.time(),
            ". Beginning PPC dataframes now.", sep = ""))


# ####### Generate dataframe for prediction plots ######
# 
# ## grab a fresh version of bdata that hasnt been transformed to match JAGS format
# bdata = dat[[slurm]]
# 
# ## We will need the metadata here to grab accurate landscape ID's 
# covars = read.csv(paste("data/ZDA_UMF/", 
#                       list.files("data/ZDA_UMF/")[grepl("clean_metadata", list.files("data/ZDA_UMF/"))],
#                       sep = ""))
# 
# ## Thin covars to match sampling units in species matrixs
# covars = covars[covars$cell_id_3km %in% rownames(bdata$y.dom),]
# 
# ## First we extract estimated abundance of both species per site 
# ## Then we predict prey abundance as a function of predator abundance
# 
# # Create empty df to fill in estimates abundance of both species
# est.dat = data.frame(matrix(NA, nrow = 0, ncol = 9))
# names(est.dat) = c("Sub_abundance", "Dom_abundance", 
#                    "lower_sub", "upper_sub", 
#                    "lower_dom", "upper_dom",
#                    "Sampling_Unit", "Landscape", 
#                    "Species_Pair")
# # Then fill it in!
# est.dat[1:length(colMeans(mod$sims.list$N.sub)), 1] = colMeans(mod$sims.list$N.sub)
# est.dat[1:length(colMeans(mod$sims.list$N.dom)), 2] = colMeans(mod$sims.list$N.dom)
# est.dat[1:length(colMeans(mod$sims.list$N.sub)), 3] = apply(mod$sims.list$N.sub, 2, quantile, probs = 0.025)
# est.dat[1:length(colMeans(mod$sims.list$N.sub)), 4] = apply(mod$sims.list$N.sub, 2, quantile, probs = 0.975)
# est.dat[1:length(colMeans(mod$sims.list$N.dom)), 5] = apply(mod$sims.list$N.dom, 2, quantile, probs = 0.025)
# est.dat[1:length(colMeans(mod$sims.list$N.dom)), 6] = apply(mod$sims.list$N.dom, 2, quantile, probs = 0.975)
# est.dat[1:length(colMeans(mod$sims.list$N.sub)), 7] = covars$cell_id_3km
# est.dat[1:length(colMeans(mod$sims.list$N.sub)), 8] = covars$Landscape
# est.dat[1:length(colMeans(mod$sims.list$N.sub)), 9] = n
# 
# 
# #### Prepare covariates included in BUGS model to make predictions
# # Extract covariates from data bundle 
# flii = bdata$flii
# hfp = bdata$hfp
# elev = bdata$elev
# hunt = bdata$hunt
# 
# # extract the dominant species abundance from the model and vary it across the observed range
# int = colMeans(mod$sims.list$N.dom)
# int = seq(from = min(int), to = max(int), length.out = length(colMeans(mod$sims.list$N.dom))) 
# # Important to keep the length.out the same at est.dat to merge together later! 
# 
# # Extract average overall landscape random effect
# # for the subordinate species  across all landscapes
# land = mod$sims.list$a6 # all RE values for each simulation of MCMC
# land = colMeans(land) # average RE value per landscape
# 
# # do the same for year
# yr = mod$sims.list$a8
# yr = colMeans(yr)
# 
# #Create empty df to fill predictions per row
# pred.dat = matrix(NA, nrow = length(int), ncol = 7)
# 
# for(l in 1:length(int)){ 
#   
#   ## Species interaction prediction
#   p = mod$sims.list$a0[,1] +             #Let abundance intercept vary
#     mod$sims.list$a1[,1] * mean(flii) +  #hold FLII constant 
#     mod$sims.list$a2[,1] * mean(hfp) +   #hold HFP constant
#     mod$sims.list$a3[,1] * mean(elev) +  #hold elev constant
#     mod$sims.list$a4[,1] * mean(hunt) +  #hold hunting constant
#     mod$sims.list$a5 * int[l] +          #vary predator abundance
#     land + yr                            #hold REs constant
#   
#   ## Back-transform from Poisson distribution
#   p = exp(p)
#   
#   ## Fill in the dataframe
#   pred.dat[l,1] = mean(p) #mean predicted abundance
#   pred.dat[l,2] = quantile(p, 0.025) # lower CI
#   pred.dat[l,3] = quantile(p, 0.975) # upper CI
#   pred.dat[l,4] = int[l] # var (Already back-transformed? --> yes!)
#   pred.dat[l,5] = "Sp_Interaction" #var name
#   pred.dat[l,6] = "Subordinate" #Dominant vs Subordinate Abundance Identifier
#   pred.dat[l,7] = n #species pair
#   
# }
# ## Will get warnings about multiplying different length objects, its ok. 
# 
# # Convert to dataframe
# pred.dat = as.data.frame(pred.dat)
# names(pred.dat) = c("Predicted.N", "lower", "upper", "var", 
#                     "var.name", "Position", "Species_Pair") #names need to match before rbinding
# # Ensure numbers are numeric
# pred.dat[,1:4] = lapply(pred.dat[,1:4], as.numeric)
# 
# ### Save these two data frames 
# day<-substr(Sys.Date(),9, 10)
# month<-substr(Sys.Date(),6,7)
# year<-substr(Sys.Date(),1,4)
# 
# path = paste(paste(paste(paste("results/co-abundance/prediction_dataframes/", slurm, "_", setting, "_predicted_abundance_plotdata_", n, "_", year,sep=""),month,sep=""),day,sep=""),".csv",sep="")
# write.csv(pred.dat, path)
# 
# path = paste(paste(paste(paste("results/co-abundance/prediction_dataframes/", slurm, "_", setting, "_estimated_abundance_per_SU_plotdata_", n, "_", year,sep=""),month,sep=""),day,sep=""),".csv",sep="")
# write.csv(est.dat, path)
# 
# print(paste("Finished generating prediction dataframes for: ", n, " at ", Sys.time(),
#             ". Beginning PPC dataframes now.", sep = ""))
# 
# 
####### Generate dataframe for PPC plots ########

# # Data Prep for Plot
# ppc.dat = data.frame("sub.sim" = mod$sims.list$fit.rep.sub,
#                      "sub.real" = mod$sims.list$fit.sub,
#                      "dom.sim" = mod$sims.list$fit.rep.dom,
#                      "dom.real" = mod$sims.list$fit.dom)
# 
# ## Add a label to color the dots
# ppc.dat$sub.label[mod$sims.list$fit.rep.sub > mod$sims.list$fit.sub]= "Above"
# ppc.dat$sub.label[mod$sims.list$fit.rep.sub < mod$sims.list$fit.sub]= "Below"
# 
# ppc.dat$dom.label[mod$sims.list$fit.rep.dom > mod$sims.list$fit.dom]= "Above"
# ppc.dat$dom.label[mod$sims.list$fit.rep.dom < mod$sims.list$fit.dom]= "Below"
# 
# ## dont forget to add the species pair! 
# ppc.dat$Species_Pair = n
# 
# #### Save it! 
# day<-substr(Sys.Date(),9, 10)
# month<-substr(Sys.Date(),6,7)
# year<-substr(Sys.Date(),1,4)
# 
# path = paste(paste(paste(paste("results/co-abundance/PPC_dataframes/", slurm, "_", setting, "_PPC_plotdata_", n, "_", year,sep=""),month,sep=""),day,sep=""),".csv",sep="")
# write.csv(ppc.dat, path)
# 
# print(paste("Finished generating PPC plotdata dataframe for: ", n, " at ", Sys.time(),
#             ". Beginning PPC values dataframe now.", sep = ""))


##### Generate dataframe w/ only BPV and C-hat values

## Create a dataframe to store results
da = data.frame(matrix(NA, nrow = 0, ncol = 11))
names(da) = c("Interaction_Estimate", "SD", "lower", "upper", "Rhat",
               "Significance", "BPV.dom", "BPV.sub", "Chat.dom", 
               "Chat.sub", "Species_Pair")

# Extract posterior mean values and Rhat
a = data.frame(mod$summary)
a = a[rownames(a) == "a5", c("mean", "sd","X2.5.", "X97.5.", "Rhat", "overlap0")]
# start filling in da
da[1, 1:6] = a

# clean up sig
da$Significance[da$Significance == 0] = "Significant"
da$Significance[da$Significance == 1] = "Non-Significant"

# BPV
da[1, 7] = mean(mod$sims.list$fit.rep.dom > mod$sims.list$fit.dom)
da[1, 8] = mean(mod$sims.list$fit.rep.sub > mod$sims.list$fit.sub)

#C-hat sub
a = mod$sims.list$fit.sub
b = mod$sims.list$fit.rep.sub
da[1, 9] = mean(a / b)

#C-hat dom
a = mod$sims.list$fit.dom
b = mod$sims.list$fit.rep.dom
da[1, 10] = mean(a / b)

#species pair
da[1, 11] = n

## save it!
day<-substr(Sys.Date(),9, 10)
month<-substr(Sys.Date(),6,7)
year<-substr(Sys.Date(),1,4)

path = paste(paste(paste(paste("results/co-abundance/dec2024_counterfactual_PPC_dataframes/", slurm, "_", setting, "_", counter, "_BPV_and_Chat_values_", n, "_", year,sep=""),month,sep=""),day,sep=""),".csv",sep="")
write.csv(da, path)

print(paste("Finished generating PPC plotdata dataframe for: ", n, " at ", Sys.time(),
            ". All dataframe extractions are now complete.", sep = ""))



