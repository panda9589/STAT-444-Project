# Load necessary library
library(readr)

setwd('C:\\Users\\tyb_l\\Documents\\STAT 444 project\\Initial fit')
# Load the data
data <- read_csv("STAT 444 project data RAW.csv")

normalize_column <- function(column) {
  mean_val <- mean(column, na.rm = TRUE)
  se_val <- sd(column, na.rm = TRUE) / sqrt(length(column[!is.na(column)]))
  normalized_column <- (column - mean_val) / se_val
  return(normalized_column)
}


# Apply the normalization to all numeric columns in the data frame
normalized_data <- data.frame(lapply(data, function(column) {
  if(is.numeric(column)) {
    normalize_column(column)
  } else {
    column
  }
}))

# View the normalized data
head(normalized_data)



# Remove rows with missing 'total_demand(mw)'
clean_data <- na.omit(normalized_data, cols = "total_demand(mw)")

# norm_data <- clean_data - 

# Fit the linear regression model using all other variables
# model <- lm(`total_demand(mw)` ~ . - date_1 - date_2, data = clean_data)
clean_data$temp2_max.c.
model <- lm(total_demand.mw. ~ . 
            - date_1 - date_2
            - temp2.c. - `temp2_min.c.`
            - `wind_speed50_ave.m.s.`
            - `max_generation.mw.`, data = clean_data)
# Summary of the model
summary(model)

# Plotting the diagnostics
par(mfrow = c(2, 2))
plot(model)

# Calculate residuals
residuals <- residuals(model)

# Calculate MSE
mse <- mean(residuals^2)

# Print MSE
print(mse)
print(residuals)
sum(data$`temp2_ave(c)`)
str(data$`temp2_ave(c)`)


# WLS ---------------------------------------------------------------------

#define weights to use
wt <- 1 / lm(abs(model$residuals) ~ model$fitted.values)$fitted.values^2
wls_model <- lm(total_demand.mw. ~ . 
                - date_1 - date_2
                - temp2.c. - `temp2_min.c.`
                - `wind_speed50_ave.m.s.`
                - `max_generation.mw.`, data = clean_data,
                weights = wt)
summary(wls_model)

# Plotting the diagnostics
par(mfrow = c(2, 2))
plot(wls_model)

# Calculate residuals
wls_residuals <- residuals(wls_model)

# Calculate MSE
wls_mse <- mean(wls_residuals^2)
wls_mse


# KNN ------------------------------------------------------------------------
library(FNN)
# Let's try a few values for k
#
# Drop unwanted columns
knn_data <- clean_data %>% select(-date_1, -date_2, -temp2.c., -temp2_min.c., -wind_speed50_ave.m.s., -max_generation.mw.)

# Prepare data for KNN regression
x <- as.matrix(knn_data %>% select(-total_demand.mw.))
y <- knn_data$total_demand.mw.

# Fit KNN models with different values of k
knn.fit5 <- knn.reg(train = x, y = y, k = 5)
knn.fit21 <- knn.reg(train = x, y = y, k = 21)
knn.fit51 <- knn.reg(train = x, y = y, k = 51)



