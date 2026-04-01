# Animal response to runnels - Figure Generation
# Code by Amanda Wik and Nathalie Sommer
# 

# Libraries
library(vegan)
library(ape)
library(dplyr)
library(tidyverse)
library(textshape)
library(knitr)
library(ggplot2)
library(lme4)
library(ggeffects)  
library(lmerTest)
library(ggfortify)
library(MuMIn)
library(pbapply)
library(AICcmodavg)
library(rsq)
library(gtable)
library(gt)
library(TMB)
library(glmmTMB)
library(DHARMa)
library(ggrepel)
library(car)
library(betareg)
library(patchwork)

# Load data (run R Markdown first!)
combined_invert_data <- read.csv('Invertebrate_Data_2024.csv', 
                      header=TRUE, stringsAsFactors=FALSE)

# Convert to factors
combined_invert_data$Treatment <- factor(combined_invert_data$Treatment)
combined_invert_data$Condition <- factor(combined_invert_data$Water_Level)

# Add row names
combined_invert_data <- combined_invert_data %>%
  mutate(Site_ID = paste(Site, ID, sep = ""))

combined_invert_data <- combined_invert_data %>% column_to_rownames("Site_ID")

# Calculate abundance, richness and diversity

#Calculating Abundance (total number of invertebrates)
combined_invert_data <- combined_invert_data %>%
  mutate(Total_Epifauna = rowSums(across(Amphipods:Stratiomydae), na.rm = FALSE))
combined_invert_data <- combined_invert_data %>%
  mutate(Total_Infauna = rowSums(across(Capitella:Unknown_Infauna), na.rm = FALSE))

#Calculating Richness
combined_invert_data$Epifauna_Richness <- as.integer(specnumber(select(combined_invert_data, Amphipods:Stratiomydae)))
combined_invert_data$Infauna_Richness <- as.integer(specnumber(select(combined_invert_data, Capitella:Unknown_Infauna)))

#Calculating Shannon diversity
combined_invert_data$Epifauna_Shannon <- diversity(select(combined_invert_data, Amphipods:Stratiomydae), index = "shannon")
combined_invert_data$Infauna_Shannon <- diversity(select(combined_invert_data, Capitella:Unknown_Infauna), index = "shannon")


# Load fish data

fish_diet_data <- read.csv('Fish_Gut_Data_2025.csv', 
                           header=TRUE, stringsAsFactors=FALSE)

fish_survey_data <- read.csv('Fish_Survey_Data_2024.csv', 
                             header=TRUE, stringsAsFactors=FALSE)

fish_size_data <- read.csv('Fish_Size_Data_2024.csv', 
                           header=TRUE, stringsAsFactors=FALSE)

# Convert to factors
fish_survey_data$Tide <- factor(fish_survey_data$Tide)
fish_survey_data$Trap_Location <- factor(fish_survey_data$Trap_Location)
fish_survey_data$Treatment <- factor(fish_survey_data$Treatment)
fish_survey_data$Trap_Type <- factor(fish_survey_data$Trap_Type)

# Combine datasets
combined_fish_data <- merge(
  fish_survey_data,
  fish_diet_data,
  by = c("Site", "ID", "Tide"),
  all.y = TRUE
)

#Creating a dataframe of mean weight, mean length, and number of fish by trap (plot)
fish_stats_by_plot <- fish_size_data %>%
  group_by(Site, ID, Tide) %>%
  mutate(
    Mean_Weight_by_Trap = mean(Fish_Weight, na.rm = TRUE),
    Mean_Length_by_Trap = mean(Fish_Length, na.rm = TRUE),
    Total_Number_Fish_by_Trap = n()
  ) %>%
  ungroup() %>%
  select(-Fish_ID:-Fish_Weight) %>%
  distinct() 

# Add zero counts for empty traps
fish_survey_data <- fish_survey_data %>%
  left_join(fish_stats_by_plot %>% select(Site, ID, Tide, Total_Number_Fish_by_Trap),
            by = c("Site", "ID", "Tide")) %>%
  mutate(Total_Number_Fish_by_Trap = replace_na(Total_Number_Fish_by_Trap, 0))

#Adding site data to fish size dataframe
fish_field_size_data <- fish_size_data %>%
  left_join(fish_survey_data %>% select(Site:Treatment),
            by = c("Site", "ID", "Tide")) 

# Data preparation - excluding traps in the runnels
combined_fish_data_no_traps_in_runnels <- combined_fish_data %>% filter(Trap_Location != "Runnel")
fish_survey_data_no_traps_in_runnels <- fish_survey_data %>% filter(Trap_Location != "Runnel")
fish_field_size_data_no_traps_in_runnels <- fish_field_size_data %>% filter(Trap_Location != "Runnel")

# Prepare fish diet data
combined_fish_data <- combined_fish_data %>%
  mutate(Total_Crab = Crabs_Whole + Crab_Parts)

# Species richness
combined_fish_data$Diet_Richness <- specnumber(select(combined_fish_data, c(Algae:Beetle_Larvae, Diptera:Total_Crab)))

# Shannon diversity
combined_fish_data$Diet_Diversity_Shannon <- combined_fish_data %>%
  select(Algae:Beetle_Larvae, Diptera:Total_Crab) %>%
  mutate(across(everything(), as.numeric)) %>%
  as.matrix() %>%
  diversity(index = "shannon")

# Add row names
combined_fish_data <- combined_fish_data %>%
  mutate(Site_ID_Tide_Version = paste(Site, ID, Tide, Fish_ID, sep = ""))

combined_fish_data <- combined_fish_data %>% column_to_rownames("Site_ID_Tide_Version")

# Filter for runnel sites only
runnel_sites_data <- fish_survey_data %>% filter(Treatment == "Runnel")
runnel_fish_data_abundance <- combined_fish_data %>% filter(Treatment == "Runnel")
runnel_fish_data_size <- fish_field_size_data %>% filter(Treatment == "Runnel")

cat("Fish data processed successfully!\n")


# Generate figures

# Color palettes
pal_nature <- c("Control" = "#6C568CFF", "Runnel" = "#607345FF")
pal_tide <- c("H" = "#9386A6FF", "L" = "#583885")  
pal_location <- c("Platform" = "#6C568CFF", "Runnel" = "#607345FF")

pub_theme <- theme_bw(base_size = 10) +
  theme(
    plot.title = element_text(size = 11, face = "bold"),
    axis.title = element_text(size = 10),
    axis.text = element_text(size = 9),
    legend.title = element_text(size = 9),
    legend.text = element_text(size = 8),
    legend.position = "bottom",
    panel.grid.minor = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5)
  )

# Define group color palette for fish plots
pal_fish <- c(
  "Control_L" = "#583885",   
  "Control_H" = "#9386A6FF",   
  "Runnel_L"  = "#507520",   
  "Runnel_H"  = "#7F8C72FF"    
)

# Define the palette for within-runnel plots
pal_within <- c(
  "Platform_L" = "#4C7FB5",
  "Platform_H" = "#7D96AFFF",
  "Runnel_L"   = "#507520",
  "Runnel_H"   = "#7F8C72FF"
)


# Figure 3: Infauna abundance and composition

# Individual plots
infauna_richness_plot <- ggplot(combined_invert_data, aes(x = Treatment, y = Infauna_Richness, fill = Treatment)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7, width = 0.6) +
  geom_jitter(width = 0.2, alpha = 0.6, size = 1.5, shape = 21, fill = "white") +
  scale_fill_manual(values = pal_nature) +
  labs(x = "Treatment", y = "Infauna Richness") +
  pub_theme +
  theme(legend.position = "none") +
  annotate("text", x = 0.5, y = Inf, label = "a", hjust = 0.5, vjust = 1.5, size = 4, fontface = "bold") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15)))

infauna_abundance_plot <- ggplot(combined_invert_data, aes(x = Treatment, y = Total_Infauna, fill = Treatment)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7, width = 0.6) +
  geom_jitter(width = 0.2, alpha = 0.6, size = 1.5, shape = 21, fill = "white") +
  scale_fill_manual(values = pal_nature) +
  labs(x = "Treatment", y = "Infauna Abundance") +
  pub_theme +
  theme(legend.position = "none") +
  annotate("text", x = 0.5, y = Inf, label = "b", hjust = 0.5, vjust = 1.5, size = 4, fontface = "bold") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15)))

# Recreate infauna NMDS data for the figure
infauna_spp <- combined_invert_data %>% select(Capitella:Unknown_Infauna) %>% na.omit() 
rows_to_remove <- c("LBSL3", "LBNE4", "LBNG6", "LBSC1", "LBSC2", "LBSC3")
infauna_spp <- infauna_spp[!rownames(infauna_spp) %in% rows_to_remove, ]

infauna_treatment <- combined_invert_data %>% 
  drop_na(Total_Infauna) %>% 
  select("Site","ID","Treatment") %>% 
  mutate(Site_ID = paste(Site, ID, sep = "")) %>% 
  select(-Site, -ID) 
infauna_treatment <- infauna_treatment[!rownames(infauna_treatment) %in% rows_to_remove, ]

# Run NMDS
infauna_mds <- metaMDS(infauna_spp, distance = "bray", autotransform = FALSE)

# Create site scores
infauna_site_scrs <- as.data.frame(scores(infauna_mds, display = "sites"))
infauna_site_scrs <- cbind(infauna_site_scrs, Treatment = infauna_treatment$Treatment)
infauna_site_scrs_plot <- rownames_to_column(infauna_site_scrs, var = "Site_ID")

infauna_nmds_plot <- ggplot(infauna_site_scrs_plot, aes(x = NMDS1, y = NMDS2, color = Treatment, shape = Treatment)) +
  geom_point(size = 2.5, alpha = 0.8) +
  stat_ellipse(aes(color = Treatment), linewidth = 0.8, linetype = 2) +
  scale_color_manual(values = pal_nature) +
  scale_shape_manual(values = c(16, 17)) +
  labs(x = "NMDS1", y = "NMDS2") +
  pub_theme +
  theme(legend.position = "none") +
  annotate("text", x = -Inf, y = Inf, label = "c", hjust = -1, vjust = 1.5, size = 4, fontface = "bold") +
  scale_x_continuous(expand = expansion(mult = c(0.15, 0.05))) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15)))

# Combine plots
infauna_combined <- infauna_richness_plot / infauna_abundance_plot / infauna_nmds_plot

# Save as TIFF
 #ggsave("Figure_3_Infauna_Combined.tiff", infauna_combined, 
  #      width = 3.5, height = 9, dpi = 300, device = "tiff")


# Figure 4: Epifauna abundance and composition

# Individual plots
epifauna_richness_plot <- ggplot(combined_invert_data, aes(x = Treatment, y = Epifauna_Richness, fill = Treatment)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7, width = 0.6) +
  geom_jitter(width = 0.2, alpha = 0.6, size = 1.5, shape = 21, fill = "white") +
  scale_fill_manual(values = pal_nature) +
  labs(x = "Treatment", y = "Epifauna Richness") +
  pub_theme +
  theme(legend.position = "none") +
  annotate("text", x = 0.5, y = Inf, label = "a", hjust = 0.5, vjust = 1.5, size = 4, fontface = "bold") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15)))

epifauna_abundance_plot <- ggplot(combined_invert_data, aes(x = Treatment, y = Total_Epifauna, fill = Treatment)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7, width = 0.6) +
  geom_jitter(width = 0.2, alpha = 0.6, size = 1.5, shape = 21, fill = "white") +
  scale_fill_manual(values = pal_nature) +
  labs(x = "Treatment", y = "Epifauna Abundance") +
  pub_theme +
  theme(legend.position = "none") +
  annotate("text", x = 0.5, y = Inf, label = "b", hjust = 0.5, vjust = 1.5, size = 4, fontface = "bold") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15)))

# Recreate epifauna NMDS data for the figure
epifauna_spp <- combined_invert_data %>% select(Amphipods:Stratiomydae)
epifauna_spp <- epifauna_spp[!rownames(epifauna_spp) %in% "LBSC6", ]

epifauna_treatment <- combined_invert_data %>% 
  select("Site", "ID", "Treatment") %>% 
  mutate(Site_ID = paste(Site, ID, sep = "")) %>% 
  select(-Site, -ID)
epifauna_treatment <- epifauna_treatment[!rownames(epifauna_treatment) %in% "LBSC6", ]

# Run NMDS
epifauna_mds <- metaMDS(epifauna_spp, distance = "bray", autotransform = FALSE)

# Create site scores
epifauna_site_scrs <- as.data.frame(scores(epifauna_mds, display = "sites"))
epifauna_site_scrs <- cbind(epifauna_site_scrs, Treatment = epifauna_treatment$Treatment)
epifauna_site_scrs_plot <- rownames_to_column(epifauna_site_scrs, var = "Site_ID")

epifauna_nmds_plot <- ggplot(epifauna_site_scrs_plot, aes(x = NMDS1, y = NMDS2, color = Treatment, shape = Treatment)) +
  geom_point(size = 2.5, alpha = 0.8) +
  stat_ellipse(aes(color = Treatment), size = 0.8, linetype = 2) +
  scale_color_manual(values = pal_nature) +
  scale_shape_manual(values = c(16, 17)) +
  labs(x = "NMDS1", y = "NMDS2") +
  pub_theme +
  theme(legend.position = "none") +
  annotate("text", x = -Inf, y = Inf, label = "c", hjust = -1, vjust = 1.5, size = 4, fontface = "bold") +
  scale_x_continuous(expand = expansion(mult = c(0.15, 0.05))) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15)))

# Combine plots
epifauna_combined <- epifauna_richness_plot / epifauna_abundance_plot / epifauna_nmds_plot

# Save as TIFF
 #ggsave("Figure_4_Epifauna_Combined.tiff", epifauna_combined, 
        #width = 3.6, height = 9, dpi = 300, device = "tiff")


# Figure 5: Fish size and abundance between control and runnel sites

fish_survey_data_no_traps_in_runnels$Group <- with(fish_survey_data_no_traps_in_runnels, paste(Treatment, Tide, sep = "_"))

fish_abundance_plot <- ggplot(fish_survey_data_no_traps_in_runnels, aes(x = Treatment, y = Total_Number_Fish_by_Trap, fill = Group)) +
  geom_boxplot(position = position_dodge(width = 0.8), alpha = 0.7, outlier.shape = NA, width = 0.6) +
  geom_jitter(aes(color = Group), position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8),
              alpha = 0.6, size = 1.5, shape = 21, fill = "white") +
  scale_fill_manual(values = pal_fish, 
                    labels = c("Control_H" = "Control\nRising", "Control_L" = "Control\nReceding", 
                              "Runnel_H" = "Runnel\nRising", "Runnel_L" = "Runnel\nReceding"),
                    name = "") +
  scale_color_manual(values = pal_fish,
                     labels = c("Control_H" = "Control\nRising", "Control_L" = "Control\nReceding", 
                               "Runnel_H" = "Runnel\nRising", "Runnel_L" = "Runnel\nReceding"),
                     name = "") +
  labs(x = "Treatment", y = "Fish Count per Trap") +
  pub_theme +
  theme(legend.position = "top",
        legend.direction = "horizontal",
        legend.box = "horizontal",
        legend.margin = margin(0, 0, 0, 0),
        legend.key.size = unit(0.5, "cm"),
        legend.text = element_text(size = 8)) +
  annotate("text", x = 0.5, y = Inf, label = "a", hjust = 0.5, vjust = 1.5, size = 4, fontface = "bold") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15)))

fish_field_size_data_no_traps_in_runnels$Group <- with(fish_field_size_data_no_traps_in_runnels, paste(Treatment, Tide, sep = "_"))

fish_weight_plot <- ggplot(fish_field_size_data_no_traps_in_runnels, aes(x = Treatment, y = Fish_Weight, fill = Group)) +
  geom_boxplot(position = position_dodge(width = 0.8), alpha = 0.7, outlier.shape = NA, width = 0.6) +
  geom_jitter(aes(color = Group), position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8),
              alpha = 0.6, size = 1.5, shape = 21, fill = "white") +
  scale_fill_manual(values = pal_fish, guide = "none") +
  scale_color_manual(values = pal_fish, guide = "none") +
  labs(x = "Treatment", y = "Weight (g)") +
  pub_theme +
  theme(legend.position = "none") +
  annotate("text", x = 0.5, y = Inf, label = "b", hjust = 0.5, vjust = 1.5, size = 4, fontface = "bold") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15)))

fish_length_plot <- ggplot(fish_field_size_data_no_traps_in_runnels, aes(x = Treatment, y = Fish_Length, fill = Group)) +
  geom_boxplot(position = position_dodge(width = 0.8), alpha = 0.7, outlier.shape = NA, width = 0.6) +
  geom_jitter(aes(color = Group), position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8),
              alpha = 0.6, size = 1.5, shape = 21, fill = "white") +
  scale_fill_manual(values = pal_fish, guide = "none") +
  scale_color_manual(values = pal_fish, guide = "none") +
  labs(x = "Treatment", y = "Length (cm)") +
  pub_theme +
  theme(legend.position = "none") +
  annotate("text", x = 0.5, y = Inf, label = "c", hjust = 0.5, vjust = 1.5, size = 4, fontface = "bold") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15)))

# Combine plots
fish_combined <- fish_abundance_plot / fish_weight_plot / fish_length_plot

# Save as TIFF
 #ggsave("Figure_5_Fish_Combined.tiff", fish_combined, 
        #width = 3.6, height = 9, dpi = 300, device = "tiff")


# Figure 6: Fish size and abundance within runnel sites

runnel_sites_data$Group <- with(runnel_sites_data, paste(Trap_Location, Tide, sep = "_"))

within_runnel_abundance_plot <- ggplot(runnel_sites_data, aes(x = Trap_Location, y = Total_Number_Fish_by_Trap, fill = Group)) +
  geom_boxplot(position = position_dodge(width = 0.8), alpha = 0.7, outlier.shape = NA, width = 0.6) +
  geom_jitter(aes(color = Group), position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8),
              alpha = 0.6, size = 1.5, shape = 21, fill = "white") +
  scale_fill_manual(values = pal_within, 
                    labels = c("Platform_H" = "Platform\nRising", "Platform_L" = "Platform\nReceding", 
                              "Runnel_H" = "Runnel\nRising", "Runnel_L" = "Runnel\nReceding"),
                    name = "") +
  scale_color_manual(values = pal_within,
                     labels = c("Platform_H" = "Platform\nRising", "Platform_L" = "Platform\nReceding", 
                               "Runnel_H" = "Runnel\nRising", "Runnel_L" = "Runnel\nReceding"),
                     name = "") +
  labs(x = "Trap Location", y = "Fish Count per Trap") +
  pub_theme +
  theme(legend.position = "top",
        legend.direction = "horizontal",
        legend.box = "horizontal",
        legend.margin = margin(0, 0, 0, 0),
        legend.key.size = unit(0.5, "cm"),
        legend.text = element_text(size = 8)) +
  annotate("text", x = 0.5, y = Inf, label = "a", hjust = 0.5, vjust = 1.5, size = 4, fontface = "bold") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15)))

runnel_fish_data_size$Group <- with(runnel_fish_data_size, paste(Trap_Location, Tide, sep = "_"))

within_runnel_weight_plot <- ggplot(runnel_fish_data_size, aes(x = Trap_Location, y = Fish_Weight, fill = Group)) +
  geom_boxplot(position = position_dodge(width = 0.8), alpha = 0.7, outlier.shape = NA, width = 0.6) +
  geom_jitter(aes(color = Group), position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8),
              alpha = 0.6, size = 1.5, shape = 21, fill = "white") +
  scale_fill_manual(values = pal_within, guide = "none") +
  scale_color_manual(values = pal_within, guide = "none") +
  labs(x = "Trap Location", y = "Weight (g)") +
  pub_theme +
  theme(legend.position = "none") +
  annotate("text", x = 0.5, y = Inf, label = "b", hjust = 0.5, vjust = 1.5, size = 4, fontface = "bold") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15)))

runnel_fish_data_size$Group <- with(runnel_fish_data_size, paste(Trap_Location, Tide, sep = "_"))

within_runnel_length_plot <- ggplot(runnel_fish_data_size, aes(x = Trap_Location, y = Fish_Length, fill = Group)) +
  geom_boxplot(position = position_dodge(width = 0.8), alpha = 0.7, outlier.shape = NA, width = 0.6) +
  geom_jitter(aes(color = Group), position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8),
              alpha = 0.6, size = 1.5, shape = 21, fill = "white") +
  scale_fill_manual(values = pal_within, guide = "none") +
  scale_color_manual(values = pal_within, guide = "none") +
  labs(x = "Trap Location", y = "Length (cm)") +
  pub_theme +
  theme(legend.position = "none") +
  annotate("text", x = 0.5, y = Inf, label = "c", hjust = 0.5, vjust = 1.5, size = 4, fontface = "bold") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15)))

# Combine plots
within_runnel_combined <- within_runnel_abundance_plot / within_runnel_weight_plot / within_runnel_length_plot

# Save as TIFF
 #ggsave("Figure_6_Within_Runnel_Fish_Combined.tiff", within_runnel_combined, 
        #width = 3.5, height = 9, dpi = 300, device = "tiff")


# Figure 7: Fish diet NMDS
# Filter for receding tide data only
fish_diet_receding <- combined_fish_data %>% 
  filter(Tide == "H")

# Remove rows with 100% digested content (no diet data)
fish_diet_receding_no_100 <- fish_diet_receding %>%
  filter(Digested != "100")

# Remove any columns with all zeroes
fish_diet_receding_no_100 <- fish_diet_receding_no_100[, colSums(fish_diet_receding_no_100 != 0) > 0]

# Separate environmental data and species data
fish_diet_spp <- fish_diet_receding_no_100 %>% select(c(Algae:Beetles, Diptera:Total_Crab))

# Check for zero-variance columns and remove them
zero_var_cols <- sapply(fish_diet_spp, function(x) var(x, na.rm = TRUE) == 0)
if(any(zero_var_cols)) {
  cat("Removing zero-variance columns:\n")
  print(names(fish_diet_spp)[zero_var_cols])
  fish_diet_spp <- fish_diet_spp[, !zero_var_cols]
}

# Also check for very low variance columns
low_var_cols <- sapply(fish_diet_spp, function(x) {
  var_val <- var(x, na.rm = TRUE)
  return(var_val < 0.001)  # Very low variance threshold
})
if(any(low_var_cols)) {
  cat("Low-variance columns found:\n")
  print(names(fish_diet_spp)[low_var_cols])
}

# Environmental dataframe
fish_diet_env <- fish_diet_receding_no_100 %>% 
  select(c(Site:Tide, Fish_ID, Trap_Type:Treatment)) %>% 
  mutate(Site_ID_Tide_Version = paste(Site, ID, Tide, Fish_ID, sep = "")) %>% 
  select(-Site, -ID, -Tide, -Fish_ID) 

if(nrow(fish_diet_spp) < 3) {
  cat("WARNING: Too few samples for NMDS analysis\n")
} else if(ncol(fish_diet_spp) < 2) {
  cat("WARNING: Too few species for NMDS analysis\n")
} else {
  # Performing nmds with environmental (treatment) and species vectors
  fish_diet_mds <- metaMDS(fish_diet_spp, distance = "bray", autotransform = FALSE)
  
  # Check stress value
  cat("NMDS stress value:", round(fish_diet_mds$stress, 4), "\n")
}

# Only run envfit if NMDS was successful and we have species data
if(!inherits(fish_diet_mds, "try-error") && ncol(fish_diet_spp) > 0) {
  fish_diet_sppfit <- envfit(fish_diet_mds, fish_diet_spp, permutations = 999) # Fit species vectors
} else {
  cat("NMDS failed or no species data available. Check your data.\n")
}

# Saving initial NMDS results
fish_diet_site_scrs <- as.data.frame(scores(fish_diet_mds, display = "sites"))
fish_diet_site_scrs <- cbind(fish_diet_site_scrs, Treatment = fish_diet_env$Treatment)

# Removing outliers

# Outliers identified as samples with NMDS1 scores above the 99th percentile in the initial NMDS ordination of fish diet composition.
outlier_indices <- which(fish_diet_site_scrs$NMDS1 > quantile(fish_diet_site_scrs$NMDS1, 0.99))
outlier_samples <- rownames(fish_diet_site_scrs)[outlier_indices]
print(outlier_samples)  

# Remove rows with outliers
fish_diet_receding_no_100_no_outliers <- fish_diet_receding_no_100 %>%
  .[!(rownames(.) %in% outlier_samples), ]

fish_diet_receding_no_100_no_outliers <- fish_diet_receding_no_100_no_outliers[, colSums(fish_diet_receding_no_100 != 0) > 0]

# Re-running NMDS with outliers removed

# Separate environmental data and species data
fish_diet_spp <- fish_diet_receding_no_100_no_outliers %>% select(c(Algae:Beetles,Diptera:Total_Crab))

# Check for zero-variance columns and remove them (after outlier removal)
zero_var_cols_no_outliers <- sapply(fish_diet_spp, function(x) var(x, na.rm = TRUE) == 0)
if(any(zero_var_cols_no_outliers)) {
  cat("Removing zero-variance columns (after outlier removal):\n")
  print(names(fish_diet_spp)[zero_var_cols_no_outliers])
  fish_diet_spp <- fish_diet_spp[, !zero_var_cols_no_outliers]
}

# Also check for very low variance columns
low_var_cols_no_outliers <- sapply(fish_diet_spp, function(x) {
  var_val <- var(x, na.rm = TRUE)
  return(var_val < 0.001)  # Very low variance threshold
})
if(any(low_var_cols_no_outliers)) {
  cat("Low-variance columns found (after outlier removal):\n")
  print(names(fish_diet_spp)[low_var_cols_no_outliers])
}

# Environmental dataframe
fish_diet_env <- fish_diet_receding_no_100_no_outliers %>% 
  select(c(Site:Tide, Fish_ID, Trap_Type:Treatment)) %>% 
  mutate(Site_ID_Tide_Version = paste(Site, ID, Tide, Fish_ID, sep = "")) %>% 
  select(-Site, -ID, -Tide, -Fish_ID) 

# Check if we have enough data for meaningful NMDS
if(nrow(fish_diet_spp) < 3) {
  cat("WARNING: Too few samples for NMDS analysis\n")
} else if(ncol(fish_diet_spp) < 2) {
  cat("WARNING: Too few species for NMDS analysis\n")
} else {
  # Performing nmds with environmental (treatment) and species vectors
  fish_diet_mds <- metaMDS(fish_diet_spp, distance = "bray", autotransform = FALSE)
  
  # Check stress value
  cat("NMDS stress value:", round(fish_diet_mds$stress, 4), "\n")
}

# Only run envfit if NMDS was successful and we have species data
if(!inherits(fish_diet_mds, "try-error") && ncol(fish_diet_spp) > 0) {
  fish_diet_sppfit <- envfit(fish_diet_mds, fish_diet_spp, permutations = 999) # Fit species vectors
} else {
  cat("NMDS failed or no species data available. Check your data.\n")
}

# Saving new NMDS results
fish_diet_site_scrs <- as.data.frame(scores(fish_diet_mds, display = "sites"))
fish_diet_site_scrs <- cbind(fish_diet_site_scrs, Treatment = fish_diet_env$Treatment)

# Create the plot
fish_diet_nmds_plot <- ggplot(fish_diet_site_scrs, aes(x = NMDS1, y = NMDS2, color = Treatment, shape = Treatment)) +
  geom_point(size = 3, alpha = 0.8) +
  stat_ellipse(aes(color = Treatment), size = 1, linetype = 2) +
  scale_color_manual(values = pal_nature) +
  scale_shape_manual(values = c(16, 17)) +
  labs(x = "NMDS1", y = "NMDS2") +
  pub_theme +
  theme(legend.position = "bottom")

# Save as TIFF
 #ggsave("Figure_7_Fish_Diet_NMDS.tiff", fish_diet_nmds_plot, 
        #width = 7, height = 5, dpi = 300, device = "tiff")


# Supplementary figures

# S1: Trap type comparision (control sites only)

control_sites_data <- fish_survey_data %>% filter(Treatment == "Control")

trap_type_abundance_plot <- ggplot(control_sites_data, aes(x = Trap_Type, y = Total_Number_Fish_by_Trap, fill = Trap_Type)) +
  geom_boxplot(alpha = 0.7, width = 0.6) +
  geom_jitter(width = 0.2, alpha = 0.6, size = 1.5, shape = 21, fill = "white") +
  scale_fill_manual(values = c("Breder" = "#583885", "Bottle" = "#9386A6FF")) +
  labs(x = "Trap Type", y = "Fish Count per Trap") +
  pub_theme +
  theme(legend.position = "none")

ggsave("Supplementary_S1_Trap_Type_Comparison.png", trap_type_abundance_plot, 
width = 4, height = 3, dpi = 600, device = "png")

                                    
# S2: Percent digested content

digested_plot <- ggplot(combined_fish_data, aes(x = Tide, y = Digested, fill = Tide)) +
  geom_boxplot(alpha = 0.7, width = 0.6) +
  geom_jitter(width = 0.2, alpha = 0.6, size = 1.5, shape = 21, fill = "white") +
  scale_fill_manual(values = pal_tide) +
  scale_x_discrete(labels = c("L" = "Rising", "H" = "Receding")) +
  labs(x = "Tide", y = "Percent Digested Content") +
  pub_theme +
  theme(legend.position = "none")

#ggsave("Supplementary_S2_Digested_Content.tiff", digested_plot, 
#width = 4, height = 3, dpi = 300, device = "tiff")

# S2: Shannon diveristy for fish diet
diet_diversity_plot <- ggplot(combined_fish_data, aes(x = Treatment, y = Diet_Diversity_Shannon, fill = Treatment)) +
  geom_boxplot(alpha = 0.7, width = 0.6) +
  geom_jitter(width = 0.2, alpha = 0.6, size = 1.5, shape = 21, fill = "white") +
  scale_fill_manual(values = pal_nature) +
  labs(x = "Treatment", y = "Shannon Diversity") +
  pub_theme +
  theme(legend.position = "none")

#ggsave("Supplementary_S3_Diet_Diversity.tiff", diet_diversity_plot, 
#width = 4, height = 3, dpi = 300, device = "tiff")
