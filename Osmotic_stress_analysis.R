library(tidyverse)
library(readxl)
library(ggplot2)
library(dplyr)

setwd("/Users/timozcelik/Desktop/Experimental")

data <- read_excel("AW_growth.xlsx")

# Assuming your data is loaded into a dataframe named 'data'
# Calculate the average day zero ChlA content for each strain and location
day_zero_avg <- data %>%
  filter(Day == 0) %>%
  group_by(Location, Strain) %>%
  summarise(Day_Zero_Avg_Chla = mean(ChlA))

# Merging the average day zero ChlA content back into the main data frame
data_merged <- merge(data, day_zero_avg, by = c("Location", "Strain"))

# Filtering out the day zero data
data_filtered <- data_merged %>%
  filter(Day == 21)

# Calculating the percent growth for each sample
data_filtered$Percent_Growth <- (data_filtered$ChlA / data_filtered$Day_Zero_Avg_Chla) * 100

# Calculating mean and standard deviation of percent growth
growth_stats <- data_filtered %>%
  group_by(Location, Strain, Osmotic_stress) %>%
  summarise(
    Mean = mean(Percent_Growth),
    SD = sd(Percent_Growth)
  )

# Adjusting the levels of Osmotic_stress to be in the order 0, -1, -2, -3
data_filtered$Osmotic_stress <- factor(data_filtered$Osmotic_stress, levels = c(0, -1, -2, -3))

# Calculating mean and standard deviation of percent growth
growth_stats <- data_filtered %>%
  group_by(Location, Strain, Osmotic_stress) %>%
  summarise(
    Mean = mean(Percent_Growth),
    SD = sd(Percent_Growth)
  )

# Creating a new column for the combined Location and Strain
growth_stats$Location_Strain <- with(growth_stats, paste(Location, Strain, sep = "_"))

# Setting colors for each Location_Strain combination
colors <- c("AT_P" = "#93AA00", "MV_P" = "#DB72FB", "AT_C" = "#F8766D", "MV_C" = "#00B9E3")

# Creating a combined plot with custom labels and bar labels
combined_plot <- ggplot(growth_stats, aes(x = Osmotic_stress, y = Mean, fill = Location_Strain)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_errorbar(aes(ymin = Mean - SD, ymax = Mean + SD), width = 0.2, position = position_dodge(0.9)) +
  geom_text(aes(label = round(Mean, 2)), vjust = -2, position = position_dodge(0.9), size = 2.5) +
  scale_fill_manual(values = colors) +
  facet_grid(Location ~ Strain, labeller = labeller(Strain = c(P = "Phormidium", C = "Chroococcidiopsis"))) +
  labs(x = "Osmotic Stress (MPa)", y = "% Growth") +
  theme_minimal()

# Print the combined plot
print(combined_plot)