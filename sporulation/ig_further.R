library(dplyr)
library(openxlsx)
library(tidyr)
library(ggplot2)
library(data.table)
library(tools)  # for file_path_sans_ext

table_path <- "outputs/rf_importance.xlsx"
# Check if the model file already exists
if (file.exists(table_path)) {
  # If the file exists, load the model
  combined_df2 <- read.xlsx(table_path)
  cat("Model loaded from file.\n")
}
names(combined_df2)[3:4] <- c("non_spore", "spore")
setDT(combined_df2)
# normalization
#combined_df2[, norm_spore_prop := spore_prop / sum(spore_prop)]
#combined_df2[, norm_non_spore_prop := non_spore_prop / sum(non_spore_prop)]
#combined_df2[, ratio := norm_spore_prop / norm_non_spore_prop]
combined_df2[, ratio := spore_prop / non_spore_prop]
# Add the 'class' column based on the ratio
combined_df2[, class := fifelse(ratio >= 0.8 & ratio <= 1.2, "equal",
                                fifelse(ratio > 1.1, "more in spore", "less in spore"))]
# Apply the new conditions to the 'class' column
combined_df2[, class := fifelse(prop > 0.75, "almost all in spore", class)]

# Apply the condition where spore_prop is approximately 0.0015 and non_spore_prop is 0
combined_df2[, class := fifelse(round(spore_prop, 4) > 0 & non_spore_prop == 0, "unique", class)]
table(combined_df2$class)
write.xlsx(combined_df2, "outputs/rf_importance.xlsx")

# Specify the directory containing the .xlsx files
dir_path <- "outputs/ranked_ig/"

# Get a list of all .xlsx files in the directory
file_list <- list.files(path = dir_path, pattern = "\\.xlsx$", full.names = TRUE)

# Read each file and store them in a list of data.tables, adding the 'binomial.name' column
data_list <- lapply(file_list, function(file) {
  dt <- as.data.table(read.xlsx(file))
  # Extract the base name without extension and add it as a new column
  dt[, binomial.name := file_path_sans_ext(basename(file))]
  return(dt)
})

# Optionally, combine all data.tables into one, if they have the same structure
ig.df <- rbindlist(data_list, use.names = TRUE, fill = TRUE)
setkey(combined_df2, gene)
setkey(ig.df, gene)

# Select only the necessary columns from combined_df2 for the join
combined_df2_subset <- combined_df2[, .(gene, spore_prop, non_spore_prop, prop, ratio, class)]

# Perform the join
ig.df <- merge(ig.df, combined_df2_subset, by = "gene", all.x = TRUE)

# Rank by 'binomial.name' and 'quantile_position' in descending order
setorder(ig.df, binomial.name, -quantile_position)
ig.df <- ig.df %>%
  dplyr::select(binomial.name, everything())
class_levels <- c("less in spore", "equal", "more in spore", "almost all in spore", "unique")

# Convert 'class' to a factor with the specified order
ig.df[, class := factor(class, levels = class_levels)]
# important sub
# Filter the data by quantile_position > 0.65
ig.sub <- ig.df[quantile_position > 0.65]

# Save ig.sub as a new data.table (optional, if you want to keep it as a separate object)
ig.sub <- copy(ig.sub)
ig.df[, quantile_group := ifelse(quantile_position > 0.65, "Above 65% Quantile", "Below 65% Quantile")]

# Create a 100% stacked bar plot using ggplot2
custom_colors <- c(
  "less in spore" = "#1f77b4",    # Blue (less in spore)
  "equal" = "#d3d3d3",            # Light grey (neutral/equal)
  "more in spore" = "#ff7f0e",    # Orange (more in spore)
  "almost all in spore" = "#d62728", # Red (almost all in spore)
  "unique" = "#800000"
)
ggplot(ig.df, aes(x = class, fill = quantile_group)) +
  geom_bar(position = "fill", stat = "count") +
  theme_bw() +
  labs(
    y = "Proportion",  # Adding y-axis title for clarity
    x = "Class",       # Label for the x-axis
    fill = "Quantile Group",  # Legend label for the color
    subtitle = "Proportion of genes by class and quantile group"
  ) +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values = c(
    "Above 65% Quantile" = "#ff7f0e",  # Orange for above 65%
    "Below 65% Quantile" = "#1f77b4"   # Blue for below 65%
  )) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.subtitle = element_text(hjust = 0.5)  # Center align the subtitle
  )
# Create the 100% stacked bar plot with the custom color palette
ggplot(ig.sub, aes(x = binomial.name, fill = class)) +
  geom_bar(position = "fill", stat = "count") +
  theme_bw() +
  labs(
    y = NULL,  # Remove y-axis title
    x = NULL,  # Remove x-axis title
    fill = "Class",
    subtitle = "Filtered for genes with average IG above 65%-quantile"  # Add subtitle
  ) +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values = custom_colors) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.subtitle = element_text(hjust = 0.5)  # Center align the subtitle
  )

ggplot(ig.sub, aes(x = class, y = quantile_position)) +
  geom_boxplot() +
  labs(
    subtitle = "Filtered for genes with average IG above the 65%-quantile",
    y = "Quantile Position",
    x = ""
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    plot.subtitle = element_text(size = 12),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10)
  )

model <- lm(quantile_position ~ binomial.name, data = ig.sub)
# Display the summary of the model
summary(model)
par(mfrow = c(2, 2))
plot(model)
