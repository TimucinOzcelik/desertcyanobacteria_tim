# Set the working directory
setwd("/Users/timozcelik/Desktop/Experimental")

# Load necessary libraries
library(readxl)
library(dplyr)
library(ggplot2)
library(tidyr)

# Read the data
data <- read_excel("temp_growth.xlsx")

# Get the initial ChlA measurement for each Location, Strain, and Temperature
initial_chla <- data %>%
  filter(Day == 0) %>%
  group_by(Location, Strain, Temperature) %>%
  summarise(InitialChlA = mean(ChlA), .groups = 'drop')

# Join the initial measurements with the main data frame and normalize as percentage
normalized_data <- data %>%
  left_join(initial_chla, by = c("Location", "Strain", "Temperature")) %>%
  mutate(PercentageChlA = (ChlA / InitialChlA) * 100) %>%
  select(-InitialChlA, -ChlA) %>%
  rename(ChlA = PercentageChlA)

# Calculate the mean and standard deviation for each group of triplicates
grouped_data <- normalized_data %>%
  group_by(Location, Strain, Temperature, Day) %>%
  summarise(
    MeanChlA = mean(ChlA),
    SdChlA = sd(ChlA),
    .groups = 'drop'
  )

# Function to plot the growth curves with error bars and customized shapes
plot_growth_curve <- function(data, strain, custom_colors) {
  ggplot(data, aes(x = Day, y = MeanChlA, group = interaction(Location, Temperature), color = Location)) +
    geom_line(size = 0.8) +
    geom_point(aes(shape = as.factor(Temperature)), size = 3) +
    geom_errorbar(aes(ymin = MeanChlA - SdChlA, ymax = MeanChlA + SdChlA), width = 0.2) +
    scale_shape_manual(values=c(16, 17)) +
    scale_color_manual(values=custom_colors) +
    scale_x_continuous(breaks=seq(0, max(data$Day), by = 7), name = "Time (days)") +
    labs(title = paste(strain, "Growth Curves"),
         y = "Biomass % (Chlorophyll A)",
         color = "Location",
         shape = "Temperature") +
    theme_minimal()
}

# Define custom colors for Phormidium and Chroococcidiopsis
phormidium_colors <- c("#93AA00", "#DB72FB")
chroococcidiopsis_colors <- c("#F8766D", "#00B9E3")

# Filter the data for Phormidium and Chroococcidiopsis strains
phormidium_data <- grouped_data %>% filter(Strain == 'P')
chroococcidiopsis_data <- grouped_data %>% filter(Strain == 'C')

# Plot the growth curves for Phormidium with custom colors
phormidium_plot <- plot_growth_curve(phormidium_data, "Phormidium", phormidium_colors)
print(phormidium_plot)

# Plot the growth curves for Chroococcidiopsis with custom colors
chroococcidiopsis_plot <- plot_growth_curve(chroococcidiopsis_data, "Chroococcidiopsis", chroococcidiopsis_colors)
print(chroococcidiopsis_plot)

#####
# Filter the data to only include the final biomass measurements (Day 21)
final_biomass_data <- normalized_data %>%
  filter(Day == 21)

library(broom)
library(writexl)

# Define a function to perform a t-test and return a data frame with the comparison and p-value
perform_test <- function(data, group_var, comparison_name) {
  test_result <- t.test(ChlA ~ get(group_var), data = data)
  data.frame(
    Comparison = comparison_name,
    P_Value = test_result$p.value
  )
}

# Create a list to store the comparison results
comparison_results <- list()

# Perform the comparisons
# Phormidium
comparison_results[["AT Phormidium 10 vs 25°C"]] <- perform_test(
  final_biomass_data %>% filter(Location == "AT", Strain == "P"),
  "Temperature",
  "AT Phormidium 10 vs 25°C"
)

comparison_results[["MV Phormidium 10 vs 25°C"]] <- perform_test(
  final_biomass_data %>% filter(Location == "MV", Strain == "P"),
  "Temperature",
  "MV Phormidium 10 vs 25°C"
)

comparison_results[["AT vs MV Phormidium at 10°C"]] <- perform_test(
  final_biomass_data %>% filter(Temperature == 10, Strain == "P"),
  "Location",
  "AT vs MV Phormidium at 10°C"
)

comparison_results[["AT vs MV Phormidium at 25°C"]] <- perform_test(
  final_biomass_data %>% filter(Temperature == 25, Strain == "P"),
  "Location",
  "AT vs MV Phormidium at 25°C"
)

# Chroococcidiopsis
comparison_results[["AT Chroococcidiopsis 10 vs 25°C"]] <- perform_test(
  final_biomass_data %>% filter(Location == "AT", Strain == "C"),
  "Temperature",
  "AT Chroococcidiopsis 10 vs 25°C"
)

comparison_results[["MV Chroococcidiopsis 10 vs 25°C"]] <- perform_test(
  final_biomass_data %>% filter(Location == "MV", Strain == "C"),
  "Temperature",
  "MV Chroococcidiopsis 10 vs 25°C"
)

comparison_results[["AT vs MV Chroococcidiopsis at 10°C"]] <- perform_test(
  final_biomass_data %>% filter(Temperature == 10, Strain == "C"),
  "Location",
  "AT vs MV Chroococcidiopsis at 10°C"
)

comparison_results[["AT vs MV Chroococcidiopsis at 25°C"]] <- perform_test(
  final_biomass_data %>% filter(Temperature == 25, Strain == "C"),
  "Location",
  "AT vs MV Chroococcidiopsis at 25°C"
)

# Combine all comparison results into one data frame
all_comparisons_df <- bind_rows(comparison_results)

# Write the results to an Excel file
write_xlsx(all_comparisons_df, "temp_growth_comparison_results.xlsx")
