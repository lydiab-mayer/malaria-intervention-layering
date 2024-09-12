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

# Load simulation data
out.int <- readRDS(paste0("results/", analysis.id, "/data/", analysis.id, "_simulations.rds"))

# Select simulations to plot
save.id <- "figA"
R0.id <- 3564
plot.id <- R0.id + c(677, 652, 680, 686, 661, 689, 695, 670, 698, 705, 655, 711, 714, 664, 720, 723, 673, 729)



## Prepare data for plotting --------------------------------------------------------------------------------

# Generate dataframe
df <- data.frame()

for (i in 1:length(plot.id)) {
  temp <- data.frame(out.int[[plot.id[i]]]$df, 
                     "id" = plot.id[i],
                     "int.id" = out.int[[plot.id[i]]]$int$id[1])
  df <- rbind(df, temp)
}

# Add intervention labels
df$deploy.id <- NA
df$deploy.id[df$id %in% plot.id[1:9]] <- "Pulsed deployment"
df$deploy.id[df$id %in% plot.id[10:18]] <- "Annual deployment"
df$deploy.id <- factor(df$deploy.id, levels = c("Pulsed deployment", "Annual deployment"))

df$plot.id <- NA
df$plot.id[df$id %in% plot.id[c(1:3, 1:3 + 9)]] <- "Baseline treatment rate"
df$plot.id[df$id %in% plot.id[c(4:6, 4:6 + 9)]] <- "5% treatment rate increase"
df$plot.id[df$id %in% plot.id[c(7:9, 7:9 + 9)]] <- "15% treatment rate increase"
df$plot.id <- factor(df$plot.id, levels = c("Baseline treatment rate",
                                            "5% treatment rate increase",
                                            "15% treatment rate increase"))

# Add pulsed labels
df$pulsed.id <- ifelse(df$id %in% plot.id[c(1, 4, 7, 10, 13, 16)], "Drug alone",
                       ifelse(df$id %in% plot.id[c(2, 5, 8, 11, 14, 17)], "Vaccine alone", "Combination drug + vaccine"))
df$pulsed.id <- factor(df$pulsed.id, levels = c("Drug alone", "Vaccine alone", "Combination drug + vaccine"))

# add time in years
df$year <- (df$t + 1)/ 365

# add % of total population infected
df$prop <- df$I / 10000


## Plot full data -------------------------------------------------------------------------------------------

p <- ggplot(df[df$deploy.id == "Pulsed deployment", ], aes(x = year, y = prop, linetype = pulsed.id)) +
  geom_ribbon(aes(ymin = 0, ymax = prop), alpha = 0.3, guide = "none", fill = "#4393C3") +
  geom_line(size = 0.6, colour = "#4393C3")

p <- p + facet_grid(. ~ plot.id)

p <- p + theme(panel.grid.major = element_line(linetype = "dotted", colour = "grey75"),
               panel.grid.minor = element_blank(),
               panel.background = element_blank(),
               strip.background = element_blank(),
               strip.text = element_text(family = "serif", face = "bold", size = 10),
               axis.ticks = element_blank(),
               axis.text = element_text(family = "serif", size = 10),
               axis.title.x = element_text(family = "serif", face = "bold", size = 10, vjust = -1),
               axis.title.y = element_text(family = "serif", face = "bold", size = 10, vjust = 3),
               title = element_text(family = "serif", face = "bold"),
               legend.title = element_blank(),
               legend.key = element_blank(),
               legend.text = element_text(family = "serif", size = 10),
               legend.position = "bottom") +
  scale_x_continuous(breaks = 0:10,
                     limits = c(0, 7)) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.2),
                     labels = paste0(seq(0, 100, by = 20), "%")) +
  scale_linetype_manual(values = c("dotted", "dashed", "solid"))

p <- p + labs(title = "A", x = "Year", y = "% of total population infected")

p


## Plot single plot -----------------------------------------------------------------------------------------

q <- ggplot(df[df$deploy.id == "Pulsed deployment" & df$plot.id == "Baseline treatment rate", ], aes(x = year, y = prop, linetype = pulsed.id)) +
  geom_ribbon(aes(ymin = 0, ymax = prop), alpha = 0.3, fill = "#4393C3") +
  geom_line(size = 0.6, colour = "#4393C3")

q <- q + theme(panel.grid.major = element_line(linetype = "dotted", colour = "grey75"),
               panel.grid.minor = element_blank(),
               panel.background = element_blank(),
               strip.background = element_blank(),
               strip.text = element_text(family = "sans", face = "bold", size = 8),
               axis.ticks = element_blank(),
               axis.text = element_text(family = "sans", size = 8),
               axis.title.x = element_text(family = "sans", face = "bold", size = 8, vjust = -1),
               axis.title.y = element_text(family = "sans", face = "bold", size = 8, vjust = 3),
               title = element_text(family = "sans", face = "bold"),
               legend.title = element_blank(),
               legend.key = element_blank(),
               legend.text = element_text(family = "sans", size = 6),
               legend.position = "bottom") +
  scale_x_continuous(breaks = 0:10,
                     limits = c(0, 7)) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.2),
                     labels = paste0(seq(0, 100, by = 20), "%")) +
  scale_linetype_manual(values = c("dotted", "dashed", "solid"))

q <- q + labs(x = "Year", y = "% of total population infected")

q


## Save plots -----------------------------------------------------------------------------------------------

ggsave(paste0("results/", analysis.id, "/plots/", save.id, ".jpeg"),
       plot = p,
       width = 8,
       height = 3)

ggsave(paste0("results/", analysis.id, "/plots/", save.id, "single.jpeg"),
       plot = q,
       width = 4,
       height = 3.2)
