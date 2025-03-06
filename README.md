# Large carnivores and predator-prey associations in tropical forests
## Last updated: March 6th, 2025

[![GitHub Pages](https://img.shields.io/badge/GitHub-Pages-blue?logo=github)](https://zacha46.github.io/SEA_trophic_cascades_co-abundance/)


---

## **Background**

This GitHub repository includes code to implement Amir et al.'s co-abundance analysis, which creates separate models to quantify top-down (*i.e.*, predators suppress prey) and bottom-up (*i.e.*, prey bolster predators*) forces shaping tropical forest food webs using camera trap data. The citation will be provided here once the paper is accepted in a journal. 

This analysis is implemented using three types of coding scripts:

- **RMarkdown scripts** make up most of the analysis and are organized into five clearly labeled steps in the `scripts/` directory. Each RMarkdown file generates a `.html` summary of the information and results, which can viewed on the Repository's [GitHub Pages](https://zacha46.github.io/SEA_trophic_cascades_co-abundance/) or be downloaded in the `docs/` directory.  
- **R scripts** that run on [High-Performance Computers (HPC)](https://rcc.uq.edu.au/systems/high-performance-computing/bunya) are located in `scripts/HPC_code`.  
- **SLURM scripts** that communicate with the HPC are found in `scripts/SLURM_code`.  

The camera trap data used in this analysis was prepared using a multi-step data cleaning pipeline that is backed up on GitHub but is currently private due to data-sharing agreements with our collaborators. The dataset provided in this repository has been **de-identified**, with all latitude and longitude coordinates removed to protect sampling locations.  

To learn more about this camera trap data standardization pipeline, please contact [Zachary Amir](mailto:z.amir@uq.edu.au) or [Matthew Luskin](mailto:m.luskin@uq.edu.au) to request access to the [Asian Capture Histories GitHub Repository](https://github.com/EcologicalCascadesLab/AsianCaptureHistories)

---

## **Reports, data, and results**  

### **Analysis Reports**
The following reports summarize data transformations, statistics, and results at each step of analysis. Click on a report to view the full details:

1. [Step 1: Data Preparation](https://zacha46.github.io/SEA_trophic_cascades_co-abundance/Step_1_Prepare_data_to_make_count_history_matrix.html)
2. [Step 2: Generate Co-Abundance Bundles](https://zacha46.github.io/SEA_trophic_cascades_co-abundance/Step_2_Generate_co-abundance_pairwise_data_bundles.html)
3. [Step 3: Combine Results From HPC](https://zacha46.github.io/SEA_trophic_cascades_co-abundance/Step_3_combine_results_from_HPC.html)
4. [Step 4: Visualize Histograms](https://zacha46.github.io/SEA_trophic_cascades_co-abundance/Step_4_visualize_histogram_results.html)
5. [Step 5: Meta-Regression](https://zacha46.github.io/SEA_trophic_cascades_co-abundance/Step_5_meta-regression_and_visualization.html)

- Alternatively, you can download the reports as `.html` files found in the `docs/` directory.  

### **Data Structure**
Except for the cleaned camera trap data that originates from the [Asian Capture Histories GitHub Repository](https://github.com/EcologicalCascadesLab/AsianCaptureHistories) and is imported in `scripts/Step_1_Prepare_data_to_make_count_history_matrix.Rmd`, all data to reproduce our analysis is provided on this Github repository. Data are divided into several different folders. 

- `data_GitHub/` is a directory containing the data cleaned in `scripts/Step_1_Prepare_data_to_make_count_history_matrix.Rmd` and a few additional small data files. This data is formatted into count history matrices and similarly structured covariates using the script `scripts/HPC_code/HPC_matrix_generator_SEA_TC.R`.
- `data_GitHub_UMFs/` are the count history matrices and site- and observation-level covariates per species that were produced in the R script script `scripts/HPC_code/HPC_matrix_generator_SEA_TC.R`, saved as `.RDS files`. This information is fed into `scripts/Step_2_Generate_co-abundance_pairwise_data_bundles.Rmd` to create pairwise data bundles between two species, which are then used in the co-abundance models. 
- `data_GitHub_CoA_bundles/` is the output from `scripts/Step_2_Generate_co-abundance_pairwise_data_bundles.Rmd` that combines relevant pairwise count histories and bundles them as a `.RDS file` to be implemented in co-abundance models using the script `scripts/HPC_code/HPC_co-abundance_model_new_var_comm_det.R` or `scripts/HPC_code/HPC_co-abundance_model_counterfactuals.R` for counterfactual testing.

### **Results and Figures**
- `results_final/` is a directory containing the `.csv` outputs from running co-abundance models on the HPC. Results from the HPC are stored in directories based on the MCMC settings used to implement the model and the approximate date they were submitted. Results extracted from the HPC are further sub-divided into directories for `/coefficent_dataframes` and `/PPC_dataframes`. There are sub-directories within this folder, like `results_final/step3_output_combined_results` that contains summarized spreadsheets of key results included in our manuscript. 
- `figures/` is a directory containing `.png` outputs of figures and `.csv` data to inform figures from both `scripts/Step_4_visualize_histogram_results.Rmd` and `scripts/Step_5_meta-regression_and_visualization.Rmd`, each saved in their own labelled directory. 

## **Scripts**

Scripts are organized into three categories:
1. **RMarkdown scripts (`.Rmd`)** (five steps)
2. **HPC batch scripts (`scripts/HPC_code/`)**
3. **SLURM submission scripts (`scripts/SLURM_code/`)**

### **RMarkdown analysis workflow**
- `scripts/Step_1_Prepare_data_to_make_count_history_matrix.Rmd` imports our standardized camera trap data from DropBox, assesses species detections, and formats the data to generates count history matricesand site and observation covariates for each species using the HPC script `scripts/HPC_matrix_generator.R`. 
- `scripts/Step_2_Generate_co-abundance_pairwise_data_bundles.Rmd` bundles the formatted data into species-pairs that are used in co-abundance models run on the HPC script `scripts/HPC_code/HPC_co-abundance_model_new_var_comm_det.R`. In addition to the original co-abundance test, this script bundles data into species pairs for all counterfactual tests run on the HPC script `scripts/HPC_code/HPC_co-abundance_model_counterfactuals.R`.
- `scripts/Step_3_combine_results_from_HPC.Rmd` imports the results extracted from co-abundance model on the HPC to assess the level of support from each model. The results and `.html` file generated from this script form the basis of the numerically written results of our manuscript. 
- `scripts/Step_4_visualize_histogram_results.Rmd` visualizes the results using histograms and saves simplified results `.csv files` to include on the visualization outside of GitHub. 
- `scripts/Step_5_meta-regression_and_visualization.Rmd` implements a meta-regression on the co-abundance coefficents and visualizes the resulting output. The `.html file` of this script contains key information about the effect sizes of the meta-regression. 
  
### **HPC scripts**
HPC scripts run R code on the [High-Performance Computers (HPC)](https://rcc.uq.edu.au/systems/high-performance-computing/bunya). Each script runs in batch mode so the user does not interact with the script. Printed messages are left in the HPC code to provide external updates regarding the progression of the code. This code is not meant to be used for learning how to use the HPC, that is beyond the scope of this GitHub Repository, but [see here for more information](https://github.com/UQ-RCC/hpc-docs/blob/main/guides/Bunya-User-Guide.md). 
- `scripts/HPC_code/HPC_matrix_generator.R` takes the formatted data from `scripts/Step_1_Prepare_data_to_make_count_history_matrix.Rmd` and saved in `data_GitHub/` to produce count history matrices and format observation and site-level covariates for each species individually. This code runs per species on the HPC and is fast to run (< 1 hour). The resulting data from this HPC script is saved here: `data_UMFs/`.
- `scripts/HPC_code/HPC_co-abundance_model_new_var_comm_det.R` runs the co-abundance model on the pairwise data bundles generated from `scripts/Step_2_Generate_co-abundance_pairwise_data_bundles.Rmd` and stored in `data_CoA_bundles/`. This code runs per species pair and can take a short or long time based on the MCMC settings selected. Short versions take a few hours, while long versions take ~8 days per model. Completed results from the HPC are stored here: `results_final/LONG_5km_final_and_counterfactual_Dec2024_community_detections`.
- `scripts/HPC_code/HPC_co-abundance_model_counterfactuals.R` runs the counterfactual co-abundance models using the pairwise data bundles generated from `scripts/Step_2_Generate_co-abundance_pairwise_data_bundles.Rmd` and stored in `data_CoA_bundles/counterfactual_testing/`. The bundled data has been formatted for each counterfactual test, so the analytical code to test each model remains consistent. Duration of the script follows the same logic as the original test. 

### **SLURM scripts**
SLURM scripts submit jobs to the HPC by providing details about how to run a specific R script. While there are many different versions of the co-abundance SLURM code, there is only one version of the matrix generation code. The multiple co-abundance SLURM codes only relate to different settings that are imported into the HPC script to specify computational or data requirements. 
- `scripts/SLURM_code/SLURM_generate_matrices.txt` specifies a job array with 44 jobs (one per species), where each job requires 24 GB of RAM, 1 CPU, and is allowed to run for 1 hour maximum. The last line of the R script calls the specific R script on the HPC. 
- The remaining `scripts/SLURM_code/` is split between `scripts/SLURM_code/SHORT`, `scripts/SLURM_code/MIDDLE`, `scripts/SLURM_code/LONG`, which correlate to MCMC settings used to run co-abundance models. The "LONG" code is further split into different RAM requirements, starting with 100, 250, and ending on 500 GB of RAM. The R script has been coded to bypass any completed models with the same settings. The remaining SLURM scripts are the output from this testing phase to facilitate the fastest possible analysis where pairwise data bundles have been split into community vs preferred prey with known RAM requirements. 