#### r script for extracting processed HPC models for processing and visualisation

# packages
library(tidyverse)
library(patchwork)
library(cowplot)

# set path (change to dropbox)
path <- "/Users/sassen/Desktop/01_HPC_COA_Simulations/"

# import necessary environment files
simulation_parameter_grid <- readRDS(paste0(path, "data/simulation_parameter_grid.rds"))

# Step 1 - Species Interaction Value Plots
# Gather input params
collected <- simulation_parameter_grid |>
  group_by(scenario_name, SIV) |>
  summarise(indices = list(index), .groups = "drop")

# Plotting utility
plot_scenario <- function(scenario, collected) {
  # filter collected to only the scenario
  scenario_df <- collected %>% 
    filter(scenario_name == scenario)
  
  df_all <- map2_dfr(
    scenario_df$indices, scenario_df$SIV,
    function(indices, siv) {
      map_dfr(indices, function(slurm) {
        df <- readRDS(
          paste0(path, "results/simulation_output_", slurm, ".rds")
        )$parameters
      })
    }
  )
  
  filtered.siv <- df_all |>
    filter(parameter == 'a5')
  
  # plot all SIVs
  p1 <- ggplot(filtered.siv, aes(x = true, y = median)) +
    geom_linerange(aes(ymin = lower, ymax = upper), color = "skyblue", alpha = 0.1, size=1.5) +
    geom_point(color = "blue", size = 1, alpha = 0.2) +
    geom_point(aes(y = true), color = "red", shape = 18, size = 3) +
    geom_line(aes(y=0)) +
    annotate("text", 
             x = .5, 
             y = -1.8, 
             hjust = -0, vjust = -1.6,
             label = paste("Replicates: 100"),
             size = 4) +
    labs(
      x = "True SIV",
      y = "Posterior Estimate",
      title = paste0("Posterior SIV | ",unique(scenario_df$scenario_name) )
    ) +
    theme_classic(base_size = 14) + theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 1))
  
  
  df_all <- map2_dfr(
    scenario_df$indices, scenario_df$SIV,
    function(indices, siv) {
      map_dfr(indices, function(slurm) {
        readRDS(
          paste0(path, "results/simulation_output_", slurm, ".rds")
        )$abundance %>%
          mutate(SIV = siv, replicate = slurm)
      })
    }
  )
  
  coverage_df <- df_all %>%
    mutate(covered = true >= lower & true <= upper) %>%
    group_by(parameter, SIV, site) %>%
    summarise(coverage = mean(covered), .groups = "drop")
  
  coverage_site <- coverage_df %>%
    group_by(parameter, site) %>%
    summarise(
      coverage_prop = mean(coverage),  # averaged over SIV & replicates
      .groups = "drop"
    )
  
  p2<-ggplot(coverage_site, aes(x = factor(site), y = coverage_prop)) +
    geom_point(alpha = 0.6, color = "steelblue", size = 2) +
    geom_hline(yintercept = 0.75, linetype = "dashed", color = "grey40") +
    facet_wrap(~parameter, ncol = 1, scales = "free_y") + 
    labs(y = "Coverage proportion", x = "Site", title ='Predicted Abundance Coverage') +
    theme_classic(base_size = 13) +
    theme(panel.grid.minor = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
          axis.text.x = element_blank()) +
    ylim(min(min(coverage_site$coverage_prop)-0.05,.5), 1)
  
  p1/p2
  
}

#### 1 Base
pdf("/Users/sassen/Dropbox/Co-abundance Simulations Project/Figures/Base.pdf",
    height = 15, width = 12)
p_base <- plot_scenario("Base", collected)
print(p_base)
dev.off()

#### 2a Unmeasured Variation - Low
pdf("/Users/sassen/Dropbox/Co-abundance Simulations Project/Figures/Unmeasured Variation - Low.pdf",
    height = 15, width = 12)
p_2a <- plot_scenario("Unmeasured Variation - Low", collected)
print(p_2a)
dev.off()

#### 2b Unmeasured Variation - High
pdf("/Users/sassen/Dropbox/Co-abundance Simulations Project/Figures/Unmeasured Variation - High.pdf",
    height = 15, width = 12)
p_2b <- plot_scenario("Unmeasured Variation - High", collected)
print(p_2b)
dev.off()

#### 3a State Process Overdispersion - Low
pdf("/Users/sassen/Dropbox/Co-abundance Simulations Project/Figures/State Process Overdispersion - Low.pdf",
    height = 15, width = 12)
p_3a <- plot_scenario("State Process Overdispersion - Low", collected)
print(p_3a)
dev.off()

#### 3b State Process Overdispersion - High
pdf("/Users/sassen/Dropbox/Co-abundance Simulations Project/Figures/State Process Overdispersion - High.pdf",
    height = 15, width = 12)
p_3b <- plot_scenario("State Process Overdispersion - High", collected)
print(p_3b)
dev.off()

#### 4a  Unmeasured Spatial Covariate - Low
pdf("/Users/sassen/Dropbox/Co-abundance Simulations Project/Figures/Unmeasured Spatial Covariate - Low.pdf",
    height = 15, width = 12)
p_4a <- plot_scenario("Unmeasured Spatial Covariate - Low", collected)
print(p_4a)
dev.off()

#### 4b  Unmeasured Spatial Covariate - High
pdf("/Users/sassen/Dropbox/Co-abundance Simulations Project/Figures/Unmeasured Spatial Covariate - High.pdf",
    height = 15, width = 12)
p_4b <- plot_scenario("Unmeasured Spatial Covariate - High", collected)
print(p_4b)
dev.off()

#### 5a  Double Counting - Low
pdf("/Users/sassen/Dropbox/Co-abundance Simulations Project/Figures/Double Counting - Low.pdf",
    height = 15, width = 12)
p_5a <- plot_scenario("Double Counting - Low", collected)
print(p_5a)
dev.off()

#### 5b  Double Counting - High
pdf("/Users/sassen/Dropbox/Co-abundance Simulations Project/Figures/Double Counting - High.pdf",
    height = 15, width = 12)
p_5b <- plot_scenario("Double Counting - High", collected)
print(p_5b)
dev.off()

#### 6a  Spatial Spillover - Low
pdf("/Users/sassen/Dropbox/Co-abundance Simulations Project/Figures/Spatial Spillover - Low.pdf",
    height = 15, width = 12)
p_6a <- plot_scenario("Spatial Spillover - Low", collected)
print(p_6a)
dev.off()

#### 6b  Spatial Spillover - Low
pdf("/Users/sassen/Dropbox/Co-abundance Simulations Project/Figures/Spatial Spillover - High.pdf",
    height = 15, width = 12)
p_6b <- plot_scenario("Spatial Spillover - High", collected)
print(p_6b)
dev.off()

# Step 2. Grid plots
library(tidyverse)
library(patchwork)

plot_scenarios <- function(scenarios, collected, path) {
  
  # Filter only selected scenarios
  scenario_df <- collected %>% 
    filter(scenario_name %in% scenarios)
  
  # Read, tag, and combine all scenario data
  df_all <- map_dfr(1:nrow(scenario_df), function(i) {
    this_scenario <- scenario_df$scenario_name[i]
    these_indices <- scenario_df$indices[[i]]
    
    map_dfr(these_indices, function(slurm) {
      df <- readRDS(paste0(path, "results/simulation_output_", slurm, ".rds"))$parameters
      df$scenario_name <- this_scenario
      df
    })
  })
  
  # Filter for parameter of interest
  filtered.siv <- df_all %>%
    filter(parameter == "a5") %>%
    mutate(scenario_name = factor(scenario_name, levels = scenarios))
  
  # Faceted plot
  ggplot(filtered.siv, aes(x = true, y = median)) +
    geom_linerange(aes(ymin = lower, ymax = upper),
                   color = "skyblue", alpha = 0.1, size = 1.2) +
    geom_point(color = "blue", size = 1, alpha = 0.25) +
    geom_point(aes(y = true), color = "red", shape = 18, size = 2.5) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
    facet_wrap(~ scenario_name, scales = "fixed") +
    labs(
      x = "True SIV",
      y = "Posterior Estimate",
      title = "Posterior SIV Across Scenarios"
    ) +
    theme_classic(base_size = 13) +
    theme(
      strip.background = element_rect(fill = "grey90", color = "black"),
      strip.text = element_text(face = "bold"),
      panel.border = element_rect(color = "black", fill = NA),
      plot.title = element_text(hjust = 0.5),
      plot.margin = margin(10, 10, 10, 10)
    )
}

pdf("/Users/sassen/Dropbox/Co-abundance Simulations Project/Figures/Grid Overview SIV.pdf",
    height = 15, width = 20)
plot_scenarios(c("Base",
                 "Unmeasured Variation - Low","Unmeasured Variation - High",
                 "State Process Overdispersion - Low","State Process Overdispersion - High",
                 "Unmeasured Spatial Covariate - Low","Unmeasured Spatial Covariate - High",
                 "Double Counting - Low","Double Counting - High",
                 "Spatial Spillover - Low",
                 "Spatial Spillover - High"), 
               collected, 
               path = path)
dev.off()

plot_coverage_grid <- function(scenarios, collected, path) {
  
  # Filter only selected scenarios
  scenario_df <- collected %>% 
    filter(scenario_name %in% scenarios)
  
  # Read abundance results for all scenarios
  df_all <- map2_dfr(scenario_df$indices, scenario_df$SIV,
                     function(indices, siv) {
                       map_dfr(indices, function(slurm) {
                         readRDS(paste0(path, "results/simulation_output_", slurm, ".rds"))$abundance %>%
                           mutate(SIV = siv, replicate = slurm)
                       })
                     })
  
  # Merge scenario names
  scenario_names <- rep(scenario_df$scenario_name, lengths(scenario_df$indices))
  df_all <- df_all %>%
    mutate(scenario_name = rep(scenario_names, each = nrow(df_all) / length(scenario_names)))
  
  # Compute coverage
  coverage_df <- df_all %>%
    mutate(covered = true >= lower & true <= upper) %>%
    group_by(scenario_name, parameter, SIV, site) %>%
    summarise(coverage = mean(covered), .groups = "drop")
  
  coverage_site <- coverage_df %>%
    group_by(scenario_name, parameter, site) %>%
    summarise(coverage_prop = mean(coverage), .groups = "drop")
  
  coverage_site <- coverage_site %>%
    mutate(scenario_name = factor(scenario_name, levels = scenarios))
  
  # Plot coverage grid
  ggplot(coverage_site, aes(x = factor(site), y = coverage_prop)) +
    geom_point(alpha = 0.6, color = "steelblue", size = 2) +
    geom_hline(yintercept = 0.95, linetype = "dashed", color = "grey40") +
    facet_wrap(~ scenario_name + parameter, scales = "free_y", ncol = 4) +
    labs(y = "Coverage proportion", x = "Site",
         title = "Coverage Across Scenarios and Parameters") +
    theme_classic(base_size = 13) +
    theme(
      panel.grid.minor = element_blank(),
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
      strip.background = element_rect(fill = "grey90", color = "black"),
      strip.text = element_text(face = "bold"),
      plot.title = element_text(hjust = 0.5)
    ) +
    ylim(min(min(coverage_site$coverage_prop)-0.05,.5), 1)
}

pdf("/Users/sassen/Dropbox/Co-abundance Simulations Project/Figures/Grid Overview Abundance.pdf",
    height = 15, width = 20)
plot_coverage_grid(c("Base",
                 "Unmeasured Variation - Low","Unmeasured Variation - High",
                 "State Process Overdispersion - Low","State Process Overdispersion - High",
                 "Unmeasured Spatial Covariate - Low","Unmeasured Spatial Covariate - High",
                 "Double Counting - Low","Double Counting - High",
                 "Spatial Spillover - Low",
                 "Spatial Spillover - High"), 
               collected, 
               path = path)
dev.off()

### PPC plots
plot_scenario_PPC <- function(scenarios, collected, path, PPC_metric) {
  
  # Filter only selected scenarios
  scenario_df <- collected %>% 
    filter(scenario_name %in% scenarios)
  
  # Read, tag, and combine all scenario data
  df_all <- map_dfr(1:nrow(scenario_df), function(i) {
    this_scenario <- scenario_df$scenario_name[i]
    this_SIV <- scenario_df$SIV[i]
    these_indices <- scenario_df$indices[[i]]
    
    map_dfr(these_indices, function(slurm) {
      df <- readRDS(paste0(path, "results/simulation_output_", slurm, ".rds"))$PPC
      df$scenario_name <- this_scenario
      df$SIV <- this_SIV
      df
    })
  })
  
  # Filter for parameter of interest
  if(PPC_metric == 'p'){
    filtered.siv <- df_all %>%
      filter(parameter %in% c("p.val.dom",
                              "p.val.sub")) %>%
      mutate(scenario_name = factor(scenario_name, levels = scenarios))
    
    p<-ggplot(filtered.siv, aes(x = factor(SIV), y = mean, fill = parameter)) +
      #geom_point(position = position_dodge(width = 0.5), size = 2, alpha = 0.3) +
      geom_boxplot()+
      geom_hline(yintercept = 0.25, linetype = "dashed", color = "grey40") +
      geom_hline(yintercept = 0.75, linetype = "dashed", color = "grey40") +
      facet_wrap(~scenario_name, scales = "fixed") +
      labs(
        x = "SIV",
        y = "Bayesian P-Value",
        title = "Posterior Predictive P-values Across SIV"
      ) +
      scale_fill_manual(values = c("p.val.dom" = "brown", "p.val.sub" = "darkgreen"),
                         name = "Parameter") +
      theme_classic(base_size = 13) +
      theme(
        strip.background = element_rect(fill = "grey90", color = "black"),
        strip.text = element_text(face = "bold"),
        panel.border = element_rect(color = "black", fill = NA),
        plot.title = element_text(hjust = 0.5),
        plot.margin = margin(10, 10, 10, 10)
      )
    
  }
  
  if(PPC_metric == 'c'){
    filtered.siv <- df_all %>%
      filter(parameter %in% c("c.hat.dom",
                              "c.hat.sub")) %>%
      mutate(scenario_name = factor(scenario_name, levels = scenarios))
  
  
  p<-ggplot(filtered.siv, aes(x = factor(SIV), y = mean, fill = parameter)) +
    #geom_point(position = position_dodge(width = 0.5), size = 2, alpha = 0.3) +
    geom_boxplot()+
    geom_hline(yintercept = 0.8, linetype = "dashed", color = "grey40") +
    geom_hline(yintercept = 1.2, linetype = "dashed", color = "grey40") +
    facet_wrap(~scenario_name, scales = "fixed") +
    labs(
      x = "SIV",
      y = "Bayesian C-hat",
      title = "Posterior Predictive C-hat Across SIV"
    ) +
    scale_fill_manual(values = c("c.hat.dom" = "brown", "c.hat.sub" = "darkgreen"),
                       name = "Parameter") +
    theme_classic(base_size = 13) +
    theme(
      strip.background = element_rect(fill = "grey90", color = "black"),
      strip.text = element_text(face = "bold"),
      panel.border = element_rect(color = "black", fill = NA),
      plot.title = element_text(hjust = 0.5),
      plot.margin = margin(10, 10, 10, 10)
    )
  
  }
  p
}

plot_scenario_PPC(c("Base",
                    "Unmeasured Variation - Low","Unmeasured Variation - High",
                    "State Process Overdispersion - Low","State Process Overdispersion - High",
                    "Unmeasured Spatial Covariate - Low","Unmeasured Spatial Covariate - High",
                    "Double Counting - Low","Double Counting - High",
                    "Spatial Spillover - Low",
                    "Spatial Spillover - High"), 
                  collected,
                  path,
                  PPC_metric = 'p')


# Plotting P-values

#### P-values
scen.iterate <- c("Base",
  "Unmeasured Variation - Low","Unmeasured Variation - High",
  "State Process Overdispersion - Low","State Process Overdispersion - High",
  "Unmeasured Spatial Covariate - Low","Unmeasured Spatial Covariate - High",
  "Double Counting - Low","Double Counting - High",
  "Spatial Spillover - Low",
  "Spatial Spillover - High")

for(i in scen.iterate){
  pdf(paste0("/Users/sassen/Dropbox/Co-abundance Simulations Project/Figures/PPC/P_Values/",i,".pdf"),
      height = 15, width = 12)
  print(plot_scenario_PPC(c(i), 
                    collected,
                    path,
                    PPC_metric = 'p'))
  dev.off()
}

# C-hats

for(i in scen.iterate){
  pdf(paste0("/Users/sassen/Dropbox/Co-abundance Simulations Project/Figures/PPC/C_Hat/",i,".pdf"),
      height = 15, width = 12)
  print(plot_scenario_PPC(c(i), 
                          collected,
                          path,
                          PPC_metric = 'c'))
  dev.off()
}


### Final task - a table neatly summarising our scenarios

# grab it from the design grid, just first 11
# be mindful that SIVs and covariates swap

simulation_parameter_grid[1:11,]

output_scenario_table <- data.frame(
  'Scenario' = c("Unmeasured Variation",
                 "State Process Overdispersion", "Unmeasured Spatial Covariate",
                 "Double Counting", "Spatial Spillover", "Spatial Spillover"),
  'Parameters' = c('Unmeasured SD',
                            'Theta (Overdispersion)', 'Spatial Covariate Beta (Slope)',
                            'Double Counting Rate', 'Individual Spillover Rate','Population Fraction Exposed'),
  'Parameter_Value_Low' = c(0.2, 2, 0.5, 0.2, 0.2, 0.15),
  'Parameter_Value_High' = c(0.4, 1.5, 1.5, 0.5, 0.9, 0.15)
  
)

simulation_parameter_grid[1:11,] |>
  select(scenario_name, SIV) -> temp
  

result_summary_table <- rbind(temp, transform(temp, SIV = 1))


rm(temp)

### END

### The SIV Comparison plot, showing the passing or failing of counterfactuals
scenarios <- c("Base",
               "Unmeasured Variation - Low","Unmeasured Variation - High",
               "State Process Overdispersion - Low","State Process Overdispersion - High",
               "Unmeasured Spatial Covariate - Low","Unmeasured Spatial Covariate - High",
               "Double Counting - Low","Double Counting - High",
               "Spatial Spillover - Low",
               "Spatial Spillover - High")

plot.all.counter.factuals <- function(scenarios, collected){
  # Filter only selected scenarios
  scenario_df <- collected %>% 
    filter(scenario_name %in% scenarios)
  
  # Read, tag, and combine all scenario data
  df_all <- map_dfr(1:nrow(scenario_df), function(i) {
    this_scenario <- scenario_df$scenario_name[i]
    these_indices <- scenario_df$indices[[i]]
    
    map_dfr(these_indices, function(slurm) {
      df <- readRDS(paste0(path, "results/simulation_output_", slurm, ".rds"))$parameters[9,]
      df$scenario_name <- this_scenario
      df
    })
  })
  
  # Read, tag, and combine all scenario data
  df_all.ppc <- map_dfr(1:nrow(scenario_df), function(i) {
    this_scenario <- scenario_df$scenario_name[i]
    this_SIV <- scenario_df$SIV[i]
    these_indices <- scenario_df$indices[[i]]
    
    map_dfr(these_indices, function(slurm) {
      df <- readRDS(paste0(path, "results/simulation_output_", slurm, ".rds"))$PPC
      df$scenario_name <- this_scenario
      df$SIV <- this_SIV
      df
    })
  })
  
  ppc.df <- df_all.ppc |>
    mutate(replicate = as.integer(gl(nrow(df_all.ppc) / 4, 4))) |>
    select(replicate, SIV, parameter, mean)|>
    pivot_wider(
      names_from = parameter,
      values_from = mean
    )
  
  final.result.df <- cbind(df_all, ppc.df) |> 
    mutate(result = ifelse(p.val.dom >= .85 | p.val.dom <= .15 | p.val.sub >= .85 | p.val.sub <= .15 | (true < 0 & median > 0) | (true > 0 & median < 0) | (upper > 0 & lower < 0 ) , 
                           "Failed Support Criteria",
                           "Passed Support Criteria")) 
  
 final.result.df$scenario_name <- factor(final.result.df$scenario_name,
                                          levels = scenarios)
  
  
  # Faceted plot
  p<-ggplot(final.result.df, aes(x = true, y = median)) +
    # Add coloured quadrants first
    geom_rect(aes(xmin = -Inf, xmax = 0, ymin = 0, ymax = Inf),
              fill = "grey80", alpha = 0.1) +   # top-left
    geom_rect(aes(xmin = 0, xmax = Inf, ymin = -Inf, ymax = 0),
              fill = "grey80", alpha = 0.1) +   # bottom-right
    #geom_linerange(aes(ymin = lower, ymax = upper,
    #               color = result), alpha = 0.2, size = 1.2) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
    geom_jitter(aes(color = result), size = 1, alpha = 0.8, shape = 16) +
    geom_point(aes(y = true, color = "Truth"), shape = 18, size = 1.5) +
    facet_wrap(~ scenario_name,nrow = 3, ncol = 4, scales = "fixed") +
    labs(
      x = "True SIV",
      y = "Posterior Estimate",
      title = "Posterior SIV Across Scenarios",
      colour = 'Result'
    ) + ylim(-1.5, 1.5)+   # Text annotations
    annotate("text", x = -Inf, y = Inf, label = "SIV not valid",
             hjust = -0.1, vjust = 2.1, size = 2.3, fontface = "bold") +
    annotate("text", x = Inf, y = -Inf, label = "SIV not valid",
             hjust = 1.1, vjust = -1.1, size = 2.3, fontface = "bold") +
    theme_classic(base_size = 13) +
    theme(
      strip.background = element_rect(fill = "grey90", color = "black"),
      strip.text = element_text(face = "bold"),
      panel.border = element_rect(color = "black", fill = NA),
      plot.title = element_text(hjust = 0.5),
      plot.margin = margin(10, 10, 10, 10),
      legend.position = "bottom",
      legend.text = element_text(size = 12),
      legend.title = element_text(size = 13, face = "bold"),
    ) + scale_colour_manual(values =c("Failed Support Criteria" = alpha("#B0B0B0",0.8), 
                                      "Passed Support Criteria" = alpha("#6BAE6B",0.8),
                                      "Truth" = 'red'),
                            guide = guide_legend(override.aes = list(size = 4)))
  
  print(p)
}

# Save to Dropbox
pdf(paste0("/Users/sassen/Dropbox/Co-abundance Simulations Project/Figures/Final SIV Figures/SIV_Grid.pdf"),
    height = 10, width = 15)
plot.all.counter.factuals(scenarios, collected)
dev.off()

# Individuals
for (scen in scenarios){
  pdf(paste0("/Users/sassen/Dropbox/Co-abundance Simulations Project/Figures/Final SIV Figures/",scen,".pdf"),
      height = 7.5, width = 12)
  plot.all.counter.factuals(scen, collected)
  dev.off()
}

# Additional conceptual graph, neatly showing the progression of counterfactual application
plot.counter.factuals.conceptual <- function(scenarios, collected, type = '1'){
  lower.y = -1
  upper.y = 1
  if(scenarios %in% c("Unmeasured Spatial Covariate - Low",
                      c("Unmeasured Spatial Covariate - High"))){
    lower.y = -2.5
    upper.y = 2.5
  }
  
  # Filter only selected scenarios
  scenario_df <- collected %>% 
    filter(scenario_name %in% scenarios)
  
  # Read, tag, and combine all scenario data
  df_all <- map_dfr(1:nrow(scenario_df), function(i) {
    this_scenario <- scenario_df$scenario_name[i]
    these_indices <- scenario_df$indices[[i]]
    
    map_dfr(these_indices, function(slurm) {
      df <- readRDS(paste0(path, "results/simulation_output_", slurm, ".rds"))$parameters[9,]
      df$scenario_name <- this_scenario
      df
    })
  })
  
  # Read, tag, and combine all scenario data
  df_all.ppc <- map_dfr(1:nrow(scenario_df), function(i) {
    this_scenario <- scenario_df$scenario_name[i]
    this_SIV <- scenario_df$SIV[i]
    these_indices <- scenario_df$indices[[i]]
    
    map_dfr(these_indices, function(slurm) {
      df <- readRDS(paste0(path, "results/simulation_output_", slurm, ".rds"))$PPC
      df$scenario_name <- this_scenario
      df$SIV <- this_SIV
      df
    })
  })
  
  ppc.df <- df_all.ppc |>
    mutate(replicate = as.integer(gl(nrow(df_all.ppc) / 4, 4))) |>
    select(replicate, SIV, parameter, mean)|>
    pivot_wider(
      names_from = parameter,
      values_from = mean
    )
  
  if(type == 'All'){
    final.result.df <- cbind(df_all, ppc.df) |> 
      mutate(result = ifelse(p.val.dom >= .85 | p.val.dom <= .15 | p.val.sub >= .85 | p.val.sub <= .15 | (true < 0 & median > 0) | (true > 0 & median < 0) | (upper > 0 & lower < 0 ) , 
                             "Excluded",
                             "Passed Support Criteria"),
             )
    title.name = 'Supported'
  }
  
  if(type == '3'){
    final.result.df <- cbind(df_all, ppc.df) |> 
      mutate(result = ifelse(p.val.dom >= .85 | p.val.dom <= .15 | p.val.sub >= .85 | p.val.sub <= .15 , 
                             "Failed Support Criteria",
                             "Passed Support Criteria"))
    title.name = 'Unsupported_wrong_direction'
  }
  
  if(type == '2'){
    final.result.df <- cbind(df_all, ppc.df) |> 
      mutate(result = ifelse(p.val.dom >= .85 | p.val.dom <= .15 | p.val.sub >= .85 | p.val.sub <= .15, "Excluded", 
                             ifelse(
                             (upper > 0 & lower < 0 ) , 
                             "Failed Support Criteria",
                             "Passed Support Criteria")))
    title.name = 'Unsupported_unclear_SIV'
  }

  if(type == '1'){
    final.result.df <- cbind(df_all, ppc.df) |> 
      mutate(result = ifelse(p.val.dom >= .85 | p.val.dom <= .15 | p.val.sub >= .85 | p.val.sub <= .15 | upper > 0 & lower < 0, "Excluded", 
                             ifelse(
                               ((true < 0 & median > 0) | (true > 0 & median < 0) ) , 
                               "Failed Support Criteria",
                               "Passed Support Criteria")))
    title.name = 'Unsupported_wrong_direction'
  }
  
  final.result.df$scenario_name <- factor(final.result.df$scenario_name,
                                          levels = scenarios)
  
  
  # Faceted plot
  p1<- ggplot(final.result.df, aes(x = true, y = median)) +
    # Add coloured quadrants first
    geom_rect(aes(xmin = -Inf, xmax = 0, ymin = 0, ymax = Inf),
              fill ="grey80", alpha = ifelse(type == '1', 0.1, 0)) +   # top-left
    geom_rect(aes(xmin = 0, xmax = Inf, ymin = -Inf, ymax = 0),
              fill ="grey80", alpha = ifelse(type == '1', 0.1, 0)) +   # bottom-right
    #geom_linerange(aes(ymin = lower, ymax = upper,
    #               color = result), alpha = 0.2, size = 1.2) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
    geom_jitter(aes(color = result), size = 1, shape = 16) +
    geom_point(aes(y = true, color = "Truth"), shape = 18, size = 1.5) +
    #facet_wrap(~ scenario_name,nrow = 3, ncol = 4, scales = "fixed") +
    labs(
      x = "True SIV",
      y = "Posterior Estimate",
      title = title.name,
      colour = 'Result'
    ) + ylim(lower.y, upper.y)+   # Text annotations
    annotate("text", x = -Inf, y = Inf, label = ifelse(type == '1',"SIV not valid", ""),
             hjust = -0.1, vjust = 2.1, size = 2.3, fontface = "bold") +
    annotate("text", x = Inf, y = -Inf, label = ifelse(type == '1',"SIV not valid", ""),
             hjust = 1.1, vjust = -1.1, size = 2.3, fontface = "bold") +
    theme_classic(base_size = 13) +
    theme(
      strip.background = element_rect(fill = "grey90", color = "black"),
      strip.text = element_text(face = "bold"),
      panel.border = element_rect(color = "black", fill = NA),
      plot.title = element_text(hjust = 0.5),
      plot.margin = margin(10, 10, 10, 10),
      legend.position = "bottom",
      legend.text = element_text(size = 12),
      legend.title = element_text(size = 13, face = "bold"),
    ) + scale_colour_manual(values =c('Excluded' = rgb(1, 0, 0, alpha = 0),
                                      "Failed Support Criteria" = alpha("#B0B0B0",0.8), 
                                      "Passed Support Criteria" = alpha("#6BAE6B",0.8),
                                      "Truth" = 'red'),
                            guide = guide_legend(override.aes = list(size = 4)))
  return(p1)
}

for (scen in scenarios){
  # Plot all scenarios
  p_U3 <- plot.counter.factuals.conceptual(scenarios = scen, collected, type = '3') + theme(legend.position = "none") + theme(plot.margin = margin(5, 0.3, 5, 2))
  p_U2 <- plot.counter.factuals.conceptual(scenarios= scen, collected, type = '2') + theme(axis.title.y = element_blank())+ theme(legend.position = "none") + theme(plot.margin = margin(5,  0.3, 5,  0.3))
  p_U1 <- plot.counter.factuals.conceptual(scenarios= scen, collected, type = '1') + theme(axis.title.y = element_blank())+ theme(legend.position = "none") + theme(plot.margin = margin(5,  0.3, 5,  0.3))
  p_S <- plot.counter.factuals.conceptual(scenarios= scen, collected, type = 'All') + theme(axis.title.y = element_blank())+ theme(legend.position = "none") + theme(plot.margin = margin(5, 2, 5,  0.3))
  
  legend <- get_legend(
    p_U3 + theme(legend.position = "bottom",
                 legend.box.margin = margin(0,0,0,0))
  )
  combined <- plot_grid(
    p_U3, p_U2, p_U1, p_S,
    nrow = 1, align = "v", axis = "l"
  )
  final_plot <- plot_grid(
    combined, legend,
    ncol = 1, rel_heights = c(1, 0.1)  # 0.1 = space for legend
  )
  
  pdf(paste0("/Users/sassen/Dropbox/Co-abundance Simulations Project/Figures/SIV Counterfactuals/",scen,".pdf"),
      height = 5, width = 15)
  print(final_plot)
  dev.off()
}
# END


