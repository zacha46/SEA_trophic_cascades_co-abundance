#### R Script for running co-abundance models on the HPC 

#### Script for creating the simulation scenario framework which is an input for the HPC analysis

# Read in packages from HPC folders
library(jagsUI, lib.loc = "/home/s4633921/R/x86_64-pc-linux-gnu-library/4.1")
library('dplyr', lib.loc = '/home/s4633921/R/x86_64-pc-linux-gnu-library/4.1')

# Declare functions

# A little custom function to add dates to output files
date.wrap <- function(string){
  paste0(string, "_", Sys.Date(), "_JMS")
}

# Negating the 'in' operator for useful functionality
`%!in%` <- Negate(`%in%`)

# Create initial values for JAGS model
make_inits <- function(data_list){
  
  # max values observed in matricies 
  max_y_dom <- apply(data_list$y.dom, 1, max)
  max_y_sub <- apply(data_list$y.sub, 1, max)
  
  # ensure N starts >= observed maxima, and >=1 so log(N+eps) is finite
  Ndom_init <- pmax(as.integer(max_y_dom + 1L), 1L)
  Nsub_init <- pmax(as.integer(max_y_sub + 1L), 1L)
  
  # Force zero where Z==0
  Ndom_init <- data_list$Z.dom * Ndom_init
  Nsub_init <- data_list$Z.sub * Nsub_init
  
  # bundle in a list 
  return(list(
    # key params
    a0 = rnorm(2, 0, 0.2),
    a5 = rnorm(1, 0, 0.2),
    # site covs 
    a1 = rnorm(2, 0, 0.2),
    a2 = rnorm(2, 0, 0.2),
    a3 = rnorm(2, 0, 0.2),
    a4 = rnorm(1, 0, 0.2),
    # det covs 
    alpha.p.sub = plogis(0),
    alpha.p.dom = plogis(0),
    b0 = rnorm(2, 0, 0.2),
    b2 = rnorm(2, 0, 0.2),
    sd.p = runif(2, 0, 0.3),
    # abundance inital values 
    N.dom = Ndom_init,
    N.sub = Nsub_init 
  ))
}

# Function for extracting and plotting summary statistics from jags objects
extract_true_and_jags <- function(true_list, jags_out, param_names) {
  
  # --- Extract true values ---
  true_df <- do.call(rbind, lapply(param_names, function(nm) {
    if (!nm %in% names(true_list)) {
      stop("Parameter ", nm, " not found in the list")
    }
    
    vals <- true_list[[nm]]
    
    if (length(vals) == 2) {
      names_out <- paste0(nm, c("_dom", "_sub"))
    } else {
      names_out <- nm
    }
    
    data.frame(parameter = names_out, true = vals, stringsAsFactors = FALSE)
  }))
  
  rownames(true_df) <- NULL
  
  # --- Extract JAGS posterior summaries ---
  summary_table <- jags_out$summary
  jags_df <- do.call(rbind, lapply(true_df$parameter, function(pn) {
    if (grepl("_dom$", pn)) {
      jags_name <- sub("_dom$", "[1]", pn)
    } else if (grepl("_sub$", pn)) {
      jags_name <- sub("_sub$", "[2]", pn)
    } else {
      jags_name <- pn
    }
    
    if (!jags_name %in% rownames(summary_table)) {
      stop("Parameter ", jags_name, " not found in JAGS output")
    }
    
    vals <- summary_table[jags_name, c("2.5%", "50%", "97.5%")]
    data.frame(
      parameter = pn,
      lower = vals["2.5%"],
      median = vals["50%"],
      upper = vals["97.5%"],
      stringsAsFactors = FALSE
    )
  }))
  
  rownames(jags_df) <- NULL
  
  # --- Merge true and estimated ---
  merged_df <- merge(true_df, jags_df, by = "parameter")
  return(merged_df)
}

###------ Set up HPC parameters

#### Read external values from SLURM into R
setting = Sys.getenv("SETTING")  # MCMC setting

# MCMC settings, based on assignment above
## Want burn-in to be ~20% of iterations and then thin = (ni - nb) / ideal n.eff (per chain), ideally 3000 in the long one. 
### Assess n.eff via (ni - nb)/nt * nc 
if(setting == "SIMULATION"){
  ni <- 10000;  nt <- 10; nb <- 2000; nc <- 4; na = NULL
}

## specify the number of cores to be uses, should be the same as the model 
options(mc.cores = nc)

#### Read Job Array index value into R
args = commandArgs(trailingOnly = TRUE)
slurm = args[1]
slurm = as.numeric(slurm) # imports as character var, not numeric

# Read in the simulation grid
index <- readRDS("data/simulation_parameter_grid.rds")[slurm,'index']

# Read in correct data frame from master list
sim <- readRDS("data/sim_data_list.rds")[[index]]

# Construct initial values
inits_list <- list(make_inits(sim$data), 
                   make_inits(sim$data),
                   make_inits(sim$data),
                   make_inits(sim$data))

# read in model file                                                                                                                                                                 in model file

# Run jagsUI with the file path
params.mod <- c("a0","b0", "a1","a2","a3","a4","a5",
                "sigma.a6", "sigma.b3", 
                "lambda.dom", "lambda.sub",
                "y.dom", "y.sub", 
                "p.dom", "p.sub", 
                "N.dom", "N.sub",
                "p.val.dom", "p.val.sub",
                "c.hat.dom", "c.hat.sub") 

# Call the model 
mod <- jagsUI::jags(data   = sim$data,
                    inits  = inits_list,
                    parameters.to.save = params.mod,
                    model.file = "code/Co_abundance_simulation_model.txt",
                    n.chains = nc,
                    n.burnin = nb,
                    n.iter   = ni,  
                    n.thin   = nt,
                    parallel = TRUE,  
                    DIC      = FALSE)

# Params to compare
params <- c("a0", "a1", "a2", "a3", "a5", 
            "sigma.a6", "sigma.b3"####
            )

df_compare <- extract_true_and_jags(sim$true, mod, params)

# PPC Extraction
PPC.df <- data.frame("parameter" = c("p.val.dom", "p.val.sub", "c.hat.dom", "c.hat.sub"),
                     "true" = c(NA,NA, NA,NA),
                     "mean" = c(mod$mean$p.val.dom, mod$mean$p.val.sub, mod$mean$c.hat.dom, mod$mean$c.hat.sub))

# Bundle up all results
output <- list('parameters' = df_compare,
               'PPC' = PPC.df)

saveRDS(output, paste0("results/simulation_output_", slurm, ".rds"))





