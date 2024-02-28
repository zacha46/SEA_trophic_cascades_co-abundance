# Large carnivores and trophic cascades in tropical forests
## Code & data repository
### Last updated: February 28th, 2024
********

# Background

This GitHub repository includes code to implement Amir et al's co-abundance analysis to simultaneously quantify top-down and bottom-up forces shaping tropical forest food-webs using camera trap data. The citation will be provided here once it is accepted in a journal. 

This analysis uses scripts that work in R but also runs R scripts on High Performance Computers (HPC). The camera trap data that was used in this analysis was prepared using a multi-step data cleaning pipeline that is backed up on GitHub, but is currently kept private based on the data sharing agreements with our collaborators. To learn more about this camera trap history data standardization pipeline, please contact Zachary Amir (z.amir@uq.edu.au) or Matthew Luskin (m.luskin@uq.edu.au) to request access to the following GitHub repository: https://github.com/EcologicalCascadesLab/AsianCaptureHistories 


# Data

Except for the pre-cleaned camera trap data that originates from the multi-step data cleaning pipeline, all data to reproduce our analysis is provided on this Github repository. Data are divided into several different folders. 

* data_GitHub is the cleaned up version of our camera trap data that is used to format data into count histories and covariates to be run in co-abundance models. This data is used in the R script called HPC_matrix_generator.R, which runs on the HPC. 
* data_GitHub_UMFs are the bundled data matricies and covariates per species that were produced in the R script called HPC_matrix_generator.R. This information is fed into step 1 to create pairwise data bundles between two species, which are then fed to the co-abudnance models. 
* data_GitHub_CoA_bundles is the final output of step 1 that combines relevant pairwise count histories and packages them in a way for the co-abundance model to properly interpret. 
* results_final is the .csv outputs from running co-abundance models on the HPC. This information will be imported in step 2 to visualize and inspect co-abundance models. 


# Scripts

The scripts are broken down into three categories: 1) interactive R scripts, 2) HPC_code, and 3) SLURM_code. 
* Interactive R scripts are comprised of two steps and labelled as such. 
  * The first step imports our standardized camera trap data from DropBox and formats the data for our specific analysis, generates count history matrices and site and observation covariates for each species, and ultimately bundles the data into pairwise species that will be used in the co-abundance model. This script is formatted as a R Markdown script and can easily be broken into clear chunks. 
  * The second step imports the results from running the co-abundance model on the HPC to assess the level of support from each model and visualize the results in a number of ways. This script is more free-form than step 1 and is currently a regular R script, but will be converted to R Markdown once the analysis is finalized.

* HPC_code runs R scripts on the High Performance Computers. Each script runs in batch mode, instead of a more normal interactive mode, so the user does not interact with the script at all. Various printed messages are left in the HPC code to provide external updates regarding the progression of the code. This code is not meant to be used for learning how to use the HPC, that requires additional training. 
  * HPC_matrix_generator.R takes the formatted data in step 1 to produce count history matricies and format observation and site-level covariates for each species individually. This code runs per species on the HPC and is quite fast to run (< 1 hour).
  * HPC_co-abundance_model_final.R runs the co-abundance model on the pairwise data bundles generated at the end of step 1. This code runs per species pair and can take a short or long time based on the MCMC settings selected. Short versions take a few hours, while long versions take ~8 days. 

* SLURM_code is the way we communicate to the HPC to tell it how to run a specific R script. This is the only code in this project *not* coded in R, and uses Bash language instead. While there are many different versions of the co-abundance SLURM code, there is only one version of the matrix generation code. The multiple co-abundance SLURM codes only relate to different settings that are imported into the HPC script to specify computational or data requirements. This code is not meant to be used for learning how to use the HPC, that requires additional training. 
  * SLURM_generate_matricies.txt specifies a job array with 44 jobs (one per species), where each job requires 24 GB of RAM, 1 CPU, and is allowed to run for 1 hour maximum. The last line of the R script calls the specific R script on the HPC. 
  * SLURM_co-abundance_array is split between "SHORT", "MIDDLE", and "HALF_LONG" MCMC settings, which correspond to the same labels in the HPC_co-abundance_model_final.R code. The "HALF_LONG" code is further split into different RAM requirements, starting with 100, 250, and ending on 500 GB of RAM. The R script has been coded to bypass any completed models with the same settings. The remaining SLURM scripts are the output from this testing phase to facilitate the fastest possible analysis where pairwise data bundles have been split into community vs preferred prey with known RAM requirements. 
  
  
  