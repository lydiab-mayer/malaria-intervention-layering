# SCRIPT HEADER ---------------------------------------------------------------------------------------------
# This script plots results of a simple SIS model for malaria

# Lydia Braunack-Mayer
#l.brauanckmayer@gmail.com
#lydia.braunack-mayer@swisstph.ch

# 2024

# SCRIPT ----------------------------------------------------------------------------------------------------

## Set up ---------------------------------------------------------------------------------------------------

# Clear working directory
rm(list = ls())

# Set working directory
setwd("/scicore/home/penny/brauna0000/M3TPP/InterventionLayering/")

# Define ID for analysis
analysis.id <- "manuscript"

# Load required packages
library(dplyr)
library(ggplot2)

# Load postprocessing data
df <- read.csv(paste0("results/", analysis.id, "/data/", analysis.id, "_postprocessed_data.csv"))

# Define plot id
save.id <- "figB"


## Prepare data for plotting --------------------------------------------------------------------------------

# Define simulations to plot
df <- df %>%
  mutate(int.id = paste0("drug ", drug, ", vaccine ", vaccine, ", CM ", CM, ", test ", test))
unique(df$int.id)
index <- c("drug no effect, vaccine pulsed, CM no effect, test no effect",
           "drug no effect, vaccine pulsed, CM increased 5%, test no effect",
           "drug no effect, vaccine pulsed, CM increased 15%, test no effect",
           "drug pulsed, vaccine no effect, CM no effect, test pulsed",
           "drug pulsed, vaccine pulsed, CM no effect, test pulsed",
           "drug pulsed, vaccine no effect, CM increased 5%, test pulsed",
           "drug pulsed, vaccine pulsed, CM increased 5%, test pulsed",
           "drug pulsed, vaccine no effect, CM increased 15%, test pulsed",
           "drug pulsed, vaccine pulsed, CM increased 15%, test pulsed")

# Retain only selected simulations
df <- df %>%
  filter(int.id %in% index)

# Add in summary statistics for plotting
df <- df %>%
  mutate(int.coverage = ifelse(drug.coverage == vaccine.coverage & drug.coverage == test.coverage, drug.coverage, NA),
         layer.id = case_when(drug.indicator & !vaccine.indicator ~ "Drug alone",
                              !drug.indicator & vaccine.indicator ~ "Vaccine alone",
                              TRUE ~ "Combination drug + vaccine"))

# Filter out simulations with zero coverage
df <- df %>%
  filter(int.coverage != 0)


## Add outcome measures -------------------------------------------------------------------------------------

# Calculate gain in cumulative cases for drug alone vs. combo
temp1 <- df %>%
  filter(layer.id != "Vaccine alone") %>%
  group_by(R0, int.coverage, CM) %>%
  mutate(cases1YearDiff = (cases_2[vaccine == "no effect"] - 
                             cases_2[vaccine == "pulsed"]) / 
           cases_1[vaccine == "no effect"],
         cases2YearDiff = (cases_2[vaccine == "no effect"] + cases_3[vaccine == "no effect"] -
                             cases_2[vaccine == "pulsed"] - cases_3[vaccine == "pulsed"]) /
           (cases_1[vaccine == "no effect"]*2),
         cases3YearDiff = (cases_2[vaccine == "no effect"] + cases_3[vaccine == "no effect"] + cases_4[vaccine == "no effect"] - 
                             cases_2[vaccine == "pulsed"] - cases_3[vaccine == "pulsed"] - cases_4[vaccine == "pulsed"]) / 
           (cases_1[vaccine == "no effect"]*3)) %>%
  select(R0, int.coverage, CM, int.id, layer.id, cases1YearDiff, cases2YearDiff, cases3YearDiff, prev_1) %>%
  filter(layer.id != "Combination drug + vaccine") %>%
  distinct()

# Calculate gain in cumulative cases for vaccine alone vs. combo
temp2 <- df %>%
  filter(layer.id != "Drug alone") %>%
  group_by(R0, int.coverage, CM) %>%
  mutate(cases1YearDiff = (cases_2[drug == "no effect"] - 
                             cases_2[drug == "pulsed"]) / 
           cases_1[drug == "no effect"],
         cases2YearDiff = (cases_2[drug == "no effect"] + cases_3[drug == "no effect"] -
                             cases_2[drug == "pulsed"] - cases_3[drug == "pulsed"]) /
           (cases_1[drug == "no effect"]*2),
         cases3YearDiff = (cases_2[drug == "no effect"] + cases_3[drug == "no effect"] + cases_4[drug == "no effect"] - 
                             cases_2[drug == "pulsed"] - cases_3[drug == "pulsed"] - cases_4[drug == "pulsed"]) / 
           (cases_1[drug == "no effect"]*3)) %>%
  select(R0, int.coverage, CM, int.id, layer.id, cases1YearDiff, cases2YearDiff, cases3YearDiff, prev_1) %>%
  filter(layer.id != "Combination drug + vaccine") %>%
  distinct()

# Combine
df.plot <- rbind(temp1, temp2)

# Reformat
df.plot <- df.plot %>%
  mutate(cases1YearDiffFactor = cut(cases1YearDiff, breaks = seq(0, 1, 0.2), labels = c("0% - 20%",
                                                                                        "20% - 40%",
                                                                                        "40% - 60%",
                                                                                        "60% - 80%",
                                                                                        "80% - 100%")),
         cases2YearDiffFactor = cut(cases2YearDiff, breaks = seq(0, 1, 0.2), labels = c("0% - 20%",
                                                                                        "20% - 40%",
                                                                                        "40% - 60%",
                                                                                        "60% - 80%",
                                                                                        "80% - 100%")),
         cases3YearDiffFactor = cut(cases3YearDiff, breaks = seq(0, 1, 0.2), labels = c("0% - 20%",
                                                                                        "20% - 40%",
                                                                                        "40% - 60%",
                                                                                        "60% - 80%",
                                                                                        "80% - 100%")))

# Reformat labels
df.plot$CM <- ifelse(df.plot$CM == "no effect", "Baseline treatment rate",
                     ifelse(df.plot$CM == "increased 5%", "5% treatment rate increase", "15% treatment rate increase"))
df.plot$CM <- factor(df.plot$CM, levels = c("Baseline treatment rate", "5% treatment rate increase", "15% treatment rate increase"))
df.plot$layer.id <- paste0("Counterfactual: ", tolower(df.plot$layer.id))


## Plot data ------------------------------------------------------------------------------------------------

df.plot <- df.plot[df.plot$layer.id == "Counterfactual: drug alone", ]

# Year reduction
p <- ggplot(df.plot, aes(x = R0, y = int.coverage, fill = cases1YearDiffFactor)) +
  geom_tile(colour = "white") +
  facet_grid(. ~ CM)

p <- p + theme(panel.grid = element_blank(),
               panel.background = element_blank(),
               strip.background = element_blank(),
               strip.text = element_text(family = "serif", face = "bold", size = 10),
               axis.ticks = element_blank(),
               axis.title.x = element_text(family = "serif", face = "bold", size = 10),
               axis.title.y = element_text(family = "serif", face = "bold", size = 10, vjust = 2),
               legend.title = element_text(family = "serif", face = "bold", size = 10),
               title = element_text(family = "serif", face = "bold"),
               legend.key = element_blank(),
               legend.position = "bottom")

p <- p + scale_y_continuous(breaks = seq(min(df$int.coverage), max(df$int.coverage), 0.2),
                            labels = paste0(seq(min(df$int.coverage), max(df$int.coverage), 0.2)*100, "%")) +
  scale_fill_manual(values = c("#4393C3", "#87bad9", "#d3e6f1"))

p <- p + labs(title = "B", x = expression("R"["0"]), y = "Intervention coverage", fill = "Gain in cases averted due to combination drug + vaccine")# +

p


## Save plot ------------------------------------------------------------------------------------------------

ggsave(paste0("results/", analysis.id, "/plots/", save.id, ".jpeg"),
       p,
       width = 8,
       height = 3)
