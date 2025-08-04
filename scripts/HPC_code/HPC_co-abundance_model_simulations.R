###### Zachary Amir's co-abundance model, Z.Amir@uq.edu.au
##### Using 1) count matricies, 2) site covarites, & 3) observation covariates. 
#### Currently running simulated datasets thru here! 

### Data has already been formatted to into the proper data bundle to run the co-abundance model on the HPC 
## The code to see how data was simulated can be found here: scripts/Step0_co-abundance_simulation.R

## Submitted to HPC on August 4th, 2024
# Zachary Amir, z.amir@uq.edu.au

####### Set up #####

## start fresh for local testing
# rm(list = ls())

# load libraries
library(jagsUI)
library(tidyverse)

## specify the number of cores to be uses, should be the same as the model 
options(mc.cores = 3)

#### Read Job Array index value into R
slurm = Sys.getenv("SLURM_ARRAY_TASK_ID")
slurm = as.numeric(slurm) #imports as character var, not numeric
# local testing 
# slurm = 1

## and MCMC setting! 
setting = Sys.getenv("SETTING")  # MCMC setting
# setting = "SHORT"

#### List all possible bundled data files 
files = list.files("data/bundled_data")[grepl("Bundled_data", list.files("data/bundled_data"))]

# #for local testing
# setwd("~/Dropbox/Zach PhD/Ch3 Trophic release project/SEA_TC_GitHub_data_storage/data/step2_output_CoA_bundles")
# files = list.files()[grepl("Bundled_data", list.files())]

## Make sure the file is for the simulated data
f = files[grepl("SIMULATION", files)]

### Import the formatted data
dat = readRDS(paste("data/bundled_data/", f ,sep = ""))

# #for local testing
# dat = readRDS(f)


## Thin to a single species pair
bdata = dat[[slurm]]

## and save the name of the simulation test
n = names(dat)[slurm]

### Verify that this species pair at this setting has not already been completed-
# First, create file name that mimics results to check if already present
res =  paste(setting, "_co-abundance_coefficents_", n, sep = "")

# Second, list all completed results
res_search = list.files("results/simulations/coefficent_dataframes/")

# If the newly constructed file name matches ANY values already present in results,
if(any(grepl(res, res_search))){
  
  ## Print a message stating so
  print(paste("The species pair:", n, "run with MCMC settings:", setting, "is already present in the results folder. This R script is terminating now."))
  
  ## and terminate the R script fully
  stop("The script was terminated.")
}

## remove biases in list
bdata$bias_applied = NULL

## remove simulated abundances from list and save in a separate DF 
abund = data.frame("N.dom_simulated" = bdata$N.dom,
                   "N.sub_simulated" = bdata$N.sub)
bdata$N.dom = NULL
bdata$N.sub = NULL

## if results are not already present, start the models! 
print(paste("Begining to run co-abundance simulation test: ", n, " with MCMC settings: ", setting ,
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
           "fit.dom", "fit.rep.dom",               # Chi2 stat for dom, real then simulated
           # "E.dom", "E.rep.dom",
           # "E.sub", "E.rep.sub",
           "N.dom", "N.sub"                     
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
  ni <- 1000;  nt <- 2; nb <- 500; nc <- 2; na = NULL      #quick test to make sure code works, 2.5 hr per mod
}
if(setting == "MIDDLE"){
  ni = 10000;  nt = 10; nb = 2000 ; nc <- 3; na = NULL   #examine parameter values --> use this for prelim testing. 36-49 hrs per mod
}
if(setting == "LONG"){
  ni = 25000;  nt = 10; nb = 5000 ; nc <- 3; na = NULL  #publication quality run --> ~160 hours per mod, all finish in < 2 weeks! 
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

print(paste("Finished running co-abundance simulation test: ", n, " at ", Sys.time(),
            ". This model took ", round(difftime(end, start, units = "hours"), 4), " hours to be completed. Beginning dataframe extractions now.", sep = ""))


####### Generate dataframe for coefficient bar plots ######

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
s$species[grepl("\\[1]", s$var)] = "subordinate_sp"
s$species[grepl("a6", s$var)] = "subordinate_sp" #landscape RE
s$species[grepl("a8", s$var)] = "subordinate_sp" #year RE
s$species[grepl("b3", s$var)] = "subordinate_sp" #source RE
s$species[grepl("sub", s$var)] = "subordinate_sp"  
s$species[s$var == "a5"] = "subordinate_sp" # interaction 

# Add dom species name
s$species[grepl("\\[2]", s$var)] = "dominant_sp"
s$species[grepl("a7", s$var)] = "dominant_sp" # landscape RE
s$species[grepl("a9", s$var)] = "dominant_sp" # year RE
s$species[grepl("b4", s$var)] = "dominant_sp" #source RE
s$species[grepl("dom", s$var)] = "dominant_sp"

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

# Add the simulation test
s$sim_test = n

## save it! 
day<-substr(Sys.Date(),9, 10)
month<-substr(Sys.Date(),6,7)
year<-substr(Sys.Date(),1,4)

path = paste(paste(paste(paste("results/simulations/coefficent_dataframes/", slurm, "_", setting, "_", "co-abundance_coefficents_", n, "_", year,sep=""),month,sep=""),day,sep=""),".csv",sep="")
write.csv(s, path, row.names = F)


## Give us an update please! 
print(paste("Finished generating coefficent dataframe for: ", n, " at ", Sys.time(),
            ". Beginning PPC dataframe now.", sep = ""))

####### Generate PPC dataframes ########

##### Generate dataframe w/ only BPV and C-hat values

## Create a dataframe to store results
da = data.frame(matrix(NA, nrow = 0, ncol = 11))
names(da) = c("Interaction_Estimate", "SD", "lower", "upper", "Rhat",
              "Significance", "BPV.dom", "BPV.sub", "Chat.dom", 
              "Chat.sub", "sim_test")

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

#sim_test
da[1, 11] = n

## save it!
day<-substr(Sys.Date(),9, 10)
month<-substr(Sys.Date(),6,7)
year<-substr(Sys.Date(),1,4)

path = paste(paste(paste(paste("results/simulations/PPC_dataframes/", slurm, "_", setting, "_BPV_and_Chat_values_", n, "_", year,sep=""),month,sep=""),day,sep=""),".csv",sep="")
write.csv(da, path)

print(paste("Finished generating PPC plotdata dataframe for: ", n, " at ", Sys.time(),
            ". Beginning prediction dataframe extractions now.", sep = ""))

####### Generate abundance estimates and compare w/ simulation ######

## Only extracting estimated abundance of both species per site
# and comparing w/ simulated abundance 

# Create empty df to fill in estimates abundance of both species
est.dat = data.frame(matrix(NA, nrow = 0, ncol = 7))
names(est.dat) = c("Sub_abundance", "Dom_abundance",
                   "lower_sub", "upper_sub",
                   "lower_dom", "upper_dom",
                   # "Sampling_Unit", "Landscape",
                   "sim_test")
# Then fill it in!
est.dat[1:length(colMeans(mod$sims.list$N.sub)), 1] = colMeans(mod$sims.list$N.sub)
est.dat[1:length(colMeans(mod$sims.list$N.dom)), 2] = colMeans(mod$sims.list$N.dom)
est.dat[1:length(colMeans(mod$sims.list$N.sub)), 3] = apply(mod$sims.list$N.sub, 2, quantile, probs = 0.025)
est.dat[1:length(colMeans(mod$sims.list$N.sub)), 4] = apply(mod$sims.list$N.sub, 2, quantile, probs = 0.975)
est.dat[1:length(colMeans(mod$sims.list$N.dom)), 5] = apply(mod$sims.list$N.dom, 2, quantile, probs = 0.025)
est.dat[1:length(colMeans(mod$sims.list$N.dom)), 6] = apply(mod$sims.list$N.dom, 2, quantile, probs = 0.975)
est.dat[1:length(colMeans(mod$sims.list$N.sub)), 7] = n

## Now add in the true simulated abundance values 
est.dat = cbind(abund, est.dat)

## options to calculate the differences
# option 1, coverage w/in credible interval
coverage_sub <- mean(est.dat$N.sub_simulated >= est.dat$lower_sub & est.dat$N.sub_simulated <= est.dat$upper_sub)
coverage_dom <- mean(est.dat$N.dom_simulated >= est.dat$lower_dom & est.dat$N.dom_simulated <= est.dat$upper_dom)
# ideally coverage is greater than .90%
if(any(coverage_sub > .9 & coverage_dom > .9)){
  print("The estimated abundance values provide acceptable coverage within 95% CI of the estimated values")
}else{
  print("The estimated abundance values DO NOT provide acceptable coverage within 95% CI of the estimated values")
}

# option 2, root mean square deviation -> how far are the posteriors from the truth?
# RMSE summarizes both the bias and the variance in abundance estimates --> shouldnt be greater than 5
rmse_sub <- sqrt(mean((est.dat$N.sub_simulated - est.dat$Sub_abundance)^2))
rmse_dom <- sqrt(mean((est.dat$N.dom_simulated - est.dat$Dom_abundance)^2))
# lower rmse is better!
if(any(rmse_sub > 5 & rmse_sub > 5)){
  print("The estimated abundance values are NOT close to the simulated truth.")
}else{
  print("The estimated abundance values are close to the simulated truth.")
}

### Save these two data frames
day<-substr(Sys.Date(),9, 10)
month<-substr(Sys.Date(),6,7)
year<-substr(Sys.Date(),1,4)

path = paste(paste(paste(paste("results/simulations/prediction_dataframes/", slurm, "_", setting, "_estimated_abundance_SIMULATION_comparison_results_", n, "_", year,sep=""),month,sep=""),day,sep=""),".csv",sep="")
write.csv(est.dat, path)

print(paste("Finished generating prediction dataframes for: ", n, " at ", Sys.time(),
            ". All dataframes have been extracted and saved from this model. Script is terminating now.", sep = ""))
