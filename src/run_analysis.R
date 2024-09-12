# SCRIPT HEADER ---------------------------------------------------------------------------------------------
# This script runs a simple SIS model for malaria

# Lydia Braunack-Mayer
#l.brauanckmayer@gmail.com
#lydia.braunack-mayer@swisstph.ch

# 2024


# SCRIPT ----------------------------------------------------------------------------------------------------

## Set up ----------------------------------------------------------------------------------------------------

# Clear working directory
rm(list = ls())

# Set working directory
setwd("/scicore/home/penny/brauna0000/M3TPP/InterventionLayering/")

# Define ID for analysis
analysis.id <- "manuscript"

# Define results directory structure
if (!dir.exists(paste0("results/", analysis.id, "/"))) dir.create(paste0("results/", analysis.id, "/"))
if (!dir.exists(paste0("results/", analysis.id, "/plots/"))) dir.create(paste0("results/", analysis.id, "/plots/"))
if (!dir.exists(paste0("results/", analysis.id, "/data/"))) dir.create(paste0("results/", analysis.id, "/data/"))


# Load required packages
library(purrr)
library(dplyr)
library(tidyr)

# Source helper functions
source("src/helper_functions.R")

# Define maximum time frame modelled
tmax <- 365*10


## Define interventions --------------------------------------------------------------------------------------

# Define intervention coverage sample space
coverage.range <- seq(0, 1, by = 0.1)

# Create empty list to store interventions
interventions <- vector(mode = "list", length = 0)

# Create list of interventions
for (coverage in coverage.range) {
  
  # Set up structure to capture interventions
  interventions.drug <- vector(mode = "list", length = 0)
  interventions.vaccine <- vector(mode = "list", length = 0)
  interventions.CM <- vector(mode = "list", length = 0)
  interventions.test <- vector(mode = "list", length = 0)
  
  # Generate drug-based interventions
  interventions.drug[[1]] <- intervention.pe(t.start = 365 + 1,
                                             t.end = tmax,
                                             FUN.decay = weibull,
                                             id = paste0("drug - no effect, drug coverage - ", coverage),
                                             init.pe = 0,
                                             L = 0,
                                             k = 0)
  
  interventions.drug[[2]] <- intervention.pe(t.start = c(365 + 1),
                                             t.end = tmax,
                                             id = paste0("long acting drug - pulsed, drug coverage - ", coverage),
                                             FUN.decay = weibull,
                                             init.pe = 1*coverage,
                                             L = 31.1,
                                             k = 5.4)
  
  interventions.drug[[3]] <- intervention.pe(t.start = c(365 + 1, 365*2 + 1, 365*3 + 1),
                                             t.end = tmax,
                                             id = paste0("long acting drug - routine, drug coverage - ", coverage),
                                             FUN.decay = weibull,
                                             init.pe = 1*coverage,
                                             L = 31.1,
                                             k = 5.4)
  
  # Generate vaccine-based interventions
  interventions.vaccine[[1]] <- intervention.pe(t.start = 365 + 1,
                                                t.end = tmax,
                                                FUN.decay = weibull,
                                                id = paste0("vaccine - no effect, vaccine coverage - ", coverage),
                                                init.pe = 0,
                                                L = 0,
                                                k = 0)
  
  interventions.vaccine[[2]] <- intervention.pe(t.start = c(365 + 1),
                                                t.end = tmax,
                                                FUN.decay = weibull,
                                                id = paste0("vaccine - pulsed, vaccine coverage - ", coverage),
                                                init.pe = .91*coverage,
                                                L = 222,
                                                k = 0.69)
  
  interventions.vaccine[[3]] <- intervention.pe(t.start = c(365 + 1, 365*2 + 1, 365*3 + 1),
                                                t.end = tmax,
                                                FUN.decay = weibull,
                                                id = paste0("vaccine - routine, vaccine coverage - ", coverage),
                                                init.pe = .91*coverage,
                                                L = 222,
                                                k = 0.69)
  
  
  # Generated case-management interventions
  interventions.CM[[1]] <- intervention.scalar(t.int = 0:tmax,
                                               t.end = tmax,
                                               base.effect = 0,
                                               int.effect = 0,
                                               id = "case management - no effect")
  
  interventions.CM[[2]] <- intervention.scalar(t.int = (365 + 1):tmax,
                                               t.end = tmax,
                                               base.effect = 0,
                                               int.effect = 1 - (1 - 0.05)^(1/14),
                                               id = "case management - increased 5%")
  
  interventions.CM[[3]] <- intervention.scalar(t.int = (365 + 1):tmax,
                                               t.end = tmax,
                                               base.effect = 0,
                                               int.effect = 1 - (1 - 0.15)^(1/14),
                                               id = "case management - increased 15%")
  
  
  # Generate mass test and treatment interventions
  interventions.test[[1]] <- intervention.scalar(t.int = 365 + 1,
                                                 t.end = tmax,
                                                 base.effect = 1,
                                                 int.effect = 1,
                                                 id = paste0("test and treat - no effect, test and treat coverage - ", coverage))
  
  interventions.test[[2]] <- intervention.scalar(t.int = 365 + 1,
                                                 t.end = tmax,
                                                 base.effect = 1,
                                                 int.effect = 1 - coverage,
                                                 id = paste0("test and treat - pulsed, test and treat coverage - ", coverage))
  
  interventions.test[[3]] <- intervention.scalar(t.int = c(365 + 1, 365*2 + 1, 365*3 + 1),
                                                 t.end = tmax,
                                                 base.effect = 1,
                                                 int.effect = 1 - coverage,
                                                 id = paste0("test and treat - routine, test and treat coverage - ", coverage))
  
  
  # Create list with all combinations of all interventions
  temp <- list(interventions.drug, interventions.vaccine, interventions.CM, interventions.test) %>% 
    cross() %>%
    map(setNames, c("drug", "vaccine", "CM", "test"))
  
  interventions <- append(interventions, temp)
}

# See total number of intervention combinations
length(interventions)


## Run full model outputs ------------------------------------------------------------------------------------

# Define model parameters
param.tmax <- 365*10
param.delta <- 0
param.gamma <- (1 - 0.20)^(1/14)
param.beta <- seq(1 - param.gamma, 0.06, length = 9)[-1]
param.threshold <- 0

# Define initial values
N <- 10000

# Create list of parameters
params.list <- vector(mode = "list", length = length(param.beta))

for (i in 1:length(param.beta)) {
  params.list[[i]] <- list("tmax" = param.tmax,
                           "beta" = param.beta[i], 
                           "delta" = param.delta,
                           "gamma" = param.gamma,
                           "threshold" = param.threshold)
}

# Full run intervention model
out.int <- vector(mode = "list", length = length(interventions)*length(params.list))

for (j in 1:length(params.list)) {
  for (i in (length(interventions)*(j - 1) + 1):(length(interventions)*j)) {
    
    # Define init to initialise model at steady state
    init <- list("I" = N / params.list[[j]]$beta * (params.list[[j]]$beta - 1 + params.list[[j]]$gamma), "N" = N)
    
    # Run model
    out.int[[i]] <- sis.intervention(init = init, params = params.list[[j]], int = interventions[[i - length(interventions)*(j - 1)]])
    
    # Add index for model run
    out.int[[i]]$int$idnum <- i - length(interventions)*(j - 1)
  }
}


## Run post-processing ---------------------------------------------------------------------------------------

out.pp <- t(sapply(out.int, postprocessing))
out.pp <- as.data.frame(out.pp)
#View(out.pp)


## Prepare postprocessing data for analysis ------------------------------------------------------------------

head(out.pp)

df <- out.pp %>%
  mutate(across(.cols = c(ID:prevRed_10), unlist))

# Split intervention
df <- df %>%
  separate(intervention,
           into = c("drug", "drug.coverage", "vaccine", "vaccine.coverage", "CM", "test", "test.coverage"),
           sep = ",")

# Clean up labels
for(col in c("drug", "drug.coverage", "vaccine", "vaccine.coverage", "CM", "test", "test.coverage")) {
  df[, col] <- gsub(".*- ", "", df[, col])
}

# Add summary statistics
df <- df %>%
  group_by(ID, R0) %>%
  mutate("num.interventions" = sum(across(drug:test, ~ .x != " - no effect")),
         "drug.indicator" = drug != "no effect",
         "vaccine.indicator" = vaccine != "no effect",
         "CM.indicator" = CM != "no effect",
         "test.indicator" = test != "no effect",
         "pulsed.indicator" = any(grepl("pulsed", drug), grepl("pulsed", vaccine), grepl("pulsed", test)),
         "routine.indicator" = any(grepl("routine", drug), grepl("routine", vaccine), grepl("routine", test)))

# Convert back to matrix
df <- as.matrix(df)


## Write outputs to file -------------------------------------------------------------------------------------

saveRDS(out.int, paste0("results/", analysis.id, "/data/", analysis.id, "_simulations.rds"))
write.csv(df, paste0("results/", analysis.id, "/data/", analysis.id, "_postprocessed_data.csv"), row.names = TRUE)
