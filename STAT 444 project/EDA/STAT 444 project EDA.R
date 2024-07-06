# Load necessary libraries
library(ggplot2)
library(dplyr)
library(corrplot)

# Load the dataset
setwd('C:\\Users\\tyb_l\\Documents\\STAT 444 project\\EDA')

data <- read.csv("STAT 444 project data RAW.csv")

# Display the first few rows of the dataset
head(data)

# Summary statistics
summary(data)

# Data types and structure
str(data)

# Check for missing values
colSums(is.na(data))

# Convert date columns to Date type
data$date_1 <- as.Date(data$date_1, format = "%d-%m-%y")
data$date_2 <- as.Date(data$date_2, format = "%d-%m-%y")

# Plot histograms for numerical variables
numeric_vars <- sapply(data, is.numeric)
numeric_data <- data[, numeric_vars]

pdf("histograms.pdf")
par(mfrow = c(3, 3)) # Adjust the layout based on the number of numerical variables
for (col in names(numeric_data)) {
  hist(numeric_data[[col]], main = paste("Histogram of", col), xlab = col, col = "lightblue", border = "black")
}
dev.off()

# Correlation heatmap
cor_matrix <- cor(numeric_data, use = "complete.obs")
cor_matrix_filtered <- as.matrix(cor_matrix)
cor_matrix_filtered[abs(cor_matrix_filtered) < 0.9] <- 0
pdf("correlation_heatmap.pdf")
corrplot(cor_matrix, method = "color", type = "upper", tl.col = "black", tl.srt = 45)
corrplot(cor_matrix_filtered, method = "color", type = "upper", tl.col = "black", tl.srt = 45)
dev.off()

# based on new correltaion plot, we decided to remove:
# - temp2.c.
# - temp2_min.c.
# - wind_speed50_ave.m.s.
# - max_generation.mw.

# Scatter plots for some pairs of variables
pdf("scatter_plot_matrix.pdf")
pairs(numeric_data, main = "Scatter plot matrix")
dev.off()

# Boxplots for numerical variables grouped by month
pdf("boxplot_temp2_max_by_month.pdf")
ggplot(data, aes(x = as.factor(month), y = temp2_max.c.)) + 
  geom_boxplot() + 
  labs(title = "Boxplot of Max Temperature by Month", x = "Month", y = "Max Temperature (C)")
dev.off()

pdf("boxplot_wind_speed50_ave_by_month.pdf")
ggplot(data, aes(x = as.factor(month), y = wind_speed50_ave.m.s.)) + 
  geom_boxplot() + 
  labs(title = "Boxplot of Average Wind Speed by Month", x = "Month", y = "Average Wind Speed (m/s)")
dev.off()

# Time series plot for temperature and wind speed
pdf("time_series_temp2_ave.pdf")
ggplot(data, aes(x = date_1, y = temp2_ave.c.)) + 
  geom_line(color = "blue") + 
  labs(title = "Time Series of Average Temperature", x = "Date", y = "Average Temperature (C)")
dev.off()

pdf("time_series_wind_speed50_ave.pdf")
ggplot(data, aes(x = date_1, y = wind_speed50_ave.m.s.)) + 
  geom_line(color = "red") + 
  labs(title = "Time Series of Average Wind Speed", x = "Date", y = "Average Wind Speed (m/s)")
dev.off()

# Scatter plot of total demand vs. max generation
pdf("scatter_total_demand_vs_max_generation.pdf")

pdf("scatter_total_demand_vs_max_generation.pdf")

ggplot(data, aes(y = total_demand.mw., x = max_generation.mw.)) + 
  geom_point(color = "darkgreen") + 
  geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed") + 
  geom_smooth(method = "lm", color = "blue", se = FALSE) + 
  labs(title = "Scatter Plot of Total Demand vs. Max Generation", y = "Total Demand (MW)", x = "Max Generation (MW)")

dev.off()

