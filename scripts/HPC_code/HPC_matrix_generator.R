##### Format camera trap data into almost UMFs via HPC
#### it just takes too long on a regular comp,
### when I know it'll work on the HPC 

## Last edited February 28th, 2024
# Z.amir@uq.edu.au

library(tidyverse)
library(vegan)

#### import the data --> use relative file names so updates will work easily! 
species = readRDS(paste("data/ZDA_UMF/", 
                        list.files("data/ZDA_UMF/")[grepl("species_vector", list.files("data/ZDA_UMF/"))],
                        sep = ""))
group_sp = readRDS(paste("data/ZDA_UMF/", 
                        list.files("data/ZDA_UMF/")[grepl("group_living", list.files("data/ZDA_UMF/"))],
                        sep = ""))
caps = read.csv(paste("data/ZDA_UMF/", 
                      list.files("data/ZDA_UMF/")[grepl("clean_captures", list.files("data/ZDA_UMF/"))],
                      sep = ""))
meta = read.csv(paste("data/ZDA_UMF/", 
                      list.files("data/ZDA_UMF/")[grepl("clean_metadata", list.files("data/ZDA_UMF/"))],
                      sep = ""))

# ### Local testing
# species = readRDS(paste("data_GitHub/",
#                         list.files("data_GitHub/")[grepl("species_vector", list.files("data_GitHub/"))],
#                         sep = ""))
# group_sp = readRDS(paste("data_GitHub/",
#                          list.files("data_GitHub/")[grepl("group_living", list.files("data_GitHub/"))],
#                          sep = ""))
# caps = read.csv(paste("data_GitHub/",
#                      list.files("data_GitHub/")[grepl("clean_captures", list.files("data_GitHub/"))],
#                      sep = ""))
# meta = read.csv(paste("data_GitHub/",
#                      list.files("data_GitHub/")[grepl("clean_metadata", list.files("data_GitHub/"))],
#                      sep = ""))

### Assign key variables a value

# set maximum sampling duration to 100 days, b/c this is how I split surveys 
dur = 100 

# set sampling occasion window length to 5 b/c that's what I used in the co-abundance MS
w = 5 

## Read Job Array index value into R
slurm = Sys.getenv("SLURM_ARRAY_TASK_ID")
slurm = as.numeric(slurm) #imports as character var, not numeric

## select a single species 
sp = (species)[slurm] # select a single species

## let us know whats going on
print(paste("this SLURM_ARRAY_INDEX is: ", slurm, 
            "and it's running species: ", sp, sep = " "))

## Add a conditional statement if there are already relevant results in the results folder
if(any(grepl(sp, list.files("results/ZDA_UMF")))) {
  
  ## If true
  ## Print a message stating so
  print(paste("The species:", sp, "already has a matrix generated"))
  
  ## Terminate the R script fully
  stop("Terminating the script.")
}

#### Right away, remove any SU's w/ 7 or more active cameras because it created problems in the model
meta = meta[meta$cams_included_count < 7,]

## and then thin the caps also! 
caps = caps[caps$cell_id_3km %in% meta$cell_id_3km, ]

## rename caps and meta to match old code 
## And for easier local testing in case I make a mistake!
m = meta 
c = caps

## make year a character var so it isn't standardized
m$year = as.character(m$year)

## standardize site covaraites to ensure variables are comparable across models later
m.num<- m[,sapply(m, is.numeric)] 
m.std<- decostand(m.num, method = "standardize", na.rm = TRUE)
m.std = m.std[,colSums(is.na(m.std)) < nrow(m.std)] #remove any columns with no data
m.char<- m[,sapply(m, is.character)]
m<- data.frame(m.char, m.std)


#
##
###
#### Begin Creating Detection/Count History Matrix 
###
##
#


## Outline the structure (i.e. no data) of the matrix and add col and row names
mat = matrix(NA, 
             nrow = length(unique(c$cell_id_3km)), # number of rows are our sampling locations
             ncol = length(seq(from =1, to= max(c$seq))), # number of columns are our sampling occasions
             dimnames = list(as.character(unique(c$cell_id_3km)), # row names, then column names
                             seq(from =1, to= max(c$seq))))

## Determine when each sampling unit was active-
for(j in 1:length(unique(c$cell_id_3km))){
  
  a= c[c$cell_id_3km == unique(c$cell_id_3km)[j],] #subset for a specific sampling unit out of all possible options
  
  indx = seq(from = 1, to = max(a$seq)) #determine the sequence it was operational for 
  
  mat[j,indx]=0 # at row J, across all sampling occasions, put a zero there
}

## Fill in the matrix based on sampling unit and sampling occasion
for(j in 1:length(unique(c$cell_id_3km))){ #repeat for each sampling unit
  
  su = unique(c$cell_id_3km)[j] #specify the sampling unit
  
  a = c[c$cell_id_3km == su & c$Species == sp,] #subset captures data for specific sampling unit and for specific species
  
  # Fill in matrix w/ count data 
  if(nrow(a)>0 ){ #Bypass cameras without a detection and leave them as zero
    
    for(s in 1:length(a$Date)){ #repeat for each Date detected at the sampling unit
      
      d = a$Date[s] #specify the sampling date
      
      indx = a$seq[a$cell_id_3km == su & a$Date == d] #specify the matching date index
      
      ### Add a conditional for species that appear in large groups 
      if(sp %in% group_sp){
        
        # Convert all individuals to groups
        a$total_indiv_records[a$Species %in% group_sp
                              & a$total_indiv_records > 0]=  1
        
        mat[su,indx]= unique(a$total_indiv_records[a$Date == d])
        
      }else{ #for all other Species, use counts
        
        mat[su,indx]= max(unique(a$total_indiv_records[a$Date == d])) ## add max to bypass a tragulus typo. 
        
      } # end group conditional
    } # end per date loop 
  } # end nrow conditional statement
}# end matrix filling loop 


#
##
### Compress Detection/Count History Matrix into multi-day Sampling Occasions
##
#


## If all sampling units are active for less than 100 days, use the maximum sequence length for dur instead
if(max(c$seq) < dur){ dur = max(c$seq)}


## Create a new and empty compressed matrix to fit sampling occasions
dh.mat = matrix(nrow = nrow(mat), # number of rows are our sampling locations (same as above)
                ncol = round(dur/w),  # number of cols is our maximum sampling duration (dur) divided by shink window (w)
                dimnames = list(as.character(rownames(mat)), # row names are our sampling unit names
                                seq(from = 1, to = round(dur/w)))) # col names are the number of shrunk observation windows. 

### Sampling occasion loop-
for(u in 1:nrow(mat)){ #Repeat for each row in the matrix
  
  for(p in 1:round(dur/w)){ # Repeat for each sampling occasion window 
    
    # Outline the start dates of sampling occasion, using values provided in for-loop
    starts<-seq(from=1, to=dur, by=w)
    
    # Select a single start date
    l<-starts[p]
    
    ## by pass to limit matricies going out of bounds
    if(max(l:(l+w)) > dur){
      
      # create a new sequence that terminates at dur
      seq = l:(l+w)
      seq = seq[seq < dur]
      
    }else{
      
      #if its fine, save the name seq
      seq = l:(l+w)
      
    } # end out of bounds conditional 

    # Make if-else statement 
    ifelse(all(mat[u,seq]== "NA", na.rm= TRUE) == "TRUE", 
           #if all values in matrix @ row u across the sampling occasion window are all NA,
           dh.mat[u,p]<-NA, # then leave the sampling occasion as NA,
           dh.mat[u,p]<- sum(as.numeric(mat[u,seq]), na.rm = TRUE)) 
    # But if FALSE, take the the sum of the detections

  } # End loop per sampling occasion
  
} # End loop per row in matrix 


#
##
### Format Observation covariates to match compressed matrix 
##
#


## Generate observation covariates to match sampling occasions 
obs = distinct(dplyr::select(c, cell_id_3km, seq, 
                             num_cams_active_at_date)) ## Come here and change to include more obs.covs if we get them!

## create empty obs dataframe
obs.mat = (matrix(NA, 
                  nrow = length(unique(obs$cell_id_3km)), 
                  ncol = length(seq(from =1, to= max(obs$seq))),
                  dimnames = list(as.character(unique(obs$cell_id_3km)), # row names, then column names
                                  seq(from = 1, to= max(obs$seq)))))

for(u in 1:length(unique(obs$cell_id_3km))){ #repeat for each sampling unit
  
  # Select a single sampling unit (i.e. row)
  su = unique(obs$cell_id_3km)[u]
  
  # Select data from a single sampling unit
  o = obs[obs$cell_id_3km == su,]
  
  for(x in 1:max(o$seq)){ #repeat for each sequence 
    
    # Select the sequence (i.e. column)
    indx = seq(from = 1, to = max(o$seq), by = 1)[x]
    
    # select num active cams
    n = obs$num_cams_active_at_date[obs$cell_id_3km == su & obs$seq == indx] ## COME HERE IF YOU HAVE MORE OBS COVS TO ADD! 
    
    # Add conditional to force n to be 1 if no cams detected a species, but were active in the date sequence
    if(length(n) == 0 & indx <= max(o$seq)){ n = 1 } 
    
    ## Fill in the obs dataframe, matching per row and column
    obs.mat[su,indx] = n
    
  } # End per sequence
} # End per sampling unit


### Compress observation matrix to match sampling occasion window 

## Create a new and empty compressed matrix to fit sampling occasions
obs2 = matrix(nrow = nrow(obs.mat), # number of rows are our sampling locations (same as above)
              ncol = round(dur/w),  # number of cols is our maximum sampling duration (dur) divided by shink window (w)
              dimnames = list(as.character(rownames(obs.mat)), # row names are our sampling unit names
                              seq(from = 1, to = round(dur/w)))) # col names are the number of shrunk observation windows. 

for(u in 1:nrow(obs.mat)){ #Repeat for each row in the matrix
  
  for(p in 1:round(dur/w)){ # Repeat for each sampling occasion window 
    
    # Outline the start dates of sampling occasion, using values provided in for-loop
    starts<-seq(from=1, to=dur, by=w)
    
    # Select a single start date
    l<-starts[p]
    
    ## by pass to limit matricies going out of bounds
    if(max(l:(l+w)) > dur){
      
      # create a new sequence that terminates at dur
      seq = l:(l+w)
      seq = seq[seq < dur]
      
    }else{
      
      #if its fine, save the name seq
      seq = l:(l+w)
      
    } # end out of bounds conditional 
    
    
    # Make if-else statement 
    ifelse(all(obs.mat[u,seq]== "NA", na.rm= TRUE) == "TRUE", 
           #if all values in matrix @ row u across the sampling occasion window are all NA,
           obs2[u,p]<-NA, # then leave the sampling occasion as NA,
           obs2[u,p]<- sum(as.numeric(obs.mat[u,seq]), na.rm = TRUE)) 
    # But if FALSE, take the the sum of observation covariate
    
  } # End loop per sampling occasion
  
} # End loop per row in matrix 

#Standardize the observation covaraite
me = mean(as.vector(obs2),na.rm=T)
s = sd(as.vector(obs2),na.rm=T)
cams1<-(obs2-me)/s


## Replace all NA values in obs cov w/ a zero. 
for (i in 1:dim(cams1)[1]) { #Repeat for each sampling unit
  
  # Select sampling occasions where any camera was active
  a <- which(obs2[i,]>0)
  
  # and make all other occasions (i.e. NA) equal to zero
  cams1[i,-a]=0
  
}
scaled.obs = cams1


# Verify the meta matches the order as the Detection history matrix 
m = m[order(match(m$cell_id_3km, rownames(dh.mat))),]


## and bind sampling unit name to obs and dh matricies
dh.mat = as.data.frame(cbind(m$cell_id_3km, dh.mat))
names(dh.mat)[1] = "SU"

scaled.obs = as.data.frame(cbind(m$cell_id_3km, scaled.obs))
names(scaled.obs)[1] = "SU"


## save in a list
umf = list("y" = dh.mat, "siteCovs" = m, "ObsCovs" = scaled.obs)

## grab the date
day<-str_sub(Sys.Date(),-2)
month<-str_sub(Sys.Date(),-5,-4)
year<-str_sub(Sys.Date(),-10,-7)
date = paste(year,month,day, sep = "")

# Save it! 
path = paste(slurm, "_", sp, "_UMF_list_", date, ".RDS", sep = "")
saveRDS(umf, paste("results/ZDA_UMF/", path, sep = ""))

