library(dplyr)
library(ggplot2)
library(openxlsx)
library(gridExtra)

table_path <- "outputs/rf_importance.xlsx"
# Check if the model file already exists
if (file.exists(table_path)) {
  # If the file exists, load the model
  combined_df2 <- read.xlsx(table_path)
  cat("Model loaded from file.\n")
}
regression <- lm(MeanDecreaseGini ~ spore_prop + non_spore_prop, data = combined_df2)
summary(regression)
par(mfrow = c(2, 2))  # Arrange plots in a 2x2 grid
plot(regression)

names(combined_df2)[3:4] <- c("non_spore", "spore")
p1 <- ggplot(combined_df2, aes(x = MeanDecreaseGini, y = spore_prop)) +
  geom_point(color = "blue", aes(size = spore), alpha = 0.2) +
  labs(y = "Prevalence in\nSporing Bacteria", x = "Mean Decrease Gini",
       size = "Spore\nPred Score") +
  scale_size_continuous(range = c(0.5, 7)) +
  theme_minimal()
p2 <- ggplot(combined_df2, aes(x = MeanDecreaseGini, y = non_spore_prop)) +
  geom_point(color = "red", aes(size = spore), alpha = 0.2) +
  labs(y = "Prevalence in\nNon-Sporing Bacteria", x = "Mean Decrease Gini",
       size = "Spore\nPred Score") +
  scale_size_continuous(range = c(0.5, 7)) +
  theme_minimal()
grid.arrange(
  p1 + theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank()),
  p2,
  ncol = 1,
  heights = c(1, 1)
)

# plot meandecrease gini as x and non-spore, spore on y
ggplot(combined_df2, aes(x = MeanDecreaseGini)) +
  geom_point(aes(y = spore, color = "Sporing"), size = 3, alpha = 0.7) +
  geom_point(aes(y = non_spore, color = "Non-Sporing"), size = 3, alpha = 0.7) +
  labs(y = "Contribution to prediction", x = "Mean Decrease Gini", color = "Prediction") +
  theme_minimal() +
  scale_color_manual(values = c("Sporing" = "blue", "Non-Sporing" = "red"))
