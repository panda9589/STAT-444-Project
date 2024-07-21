# List of required packages
packages <- c("grid", "gridExtra", "readr", "dplyr", "Matrix", "mgcv", 
              "splines", "gamair", "caret", "ggplot2", "gam", "vars", 
              "forecast", "vars", "rBayesianOptimization", "e1071")

# Function to check and install packages
install_if_missing <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
  }
}

# Check and install each package
invisible(lapply(packages, install_if_missing))

# Load the libraries
lapply(packages, library, character.only = TRUE)


linreg_kfold <- function(y,X,K) {
  
  # split the data into k folds
  n <- length(y) # n data points
  idx <- sample.int(n) # random permutation of 1...n
  kidx <- rep(1:K,each=n/K)
  foldidx <- split(idx,kidx)
  # check
  # all.equal(rep(n/K,K),Reduce(c,lapply(foldidx,length)))
  # all.equal(1:n,sort(Reduce(c,foldidx)))
  
  # fit the model and store prediction error
  errvec <- numeric(K)
  for (k in 1:K) {
    foldtmp <- foldidx[-k]
    trainidx <- Reduce(c,foldtmp)
    testidx <- foldidx[[k]]
    ytrain <- y[trainidx]
    Xtrain <- X[trainidx, ]
    ytest <- y[testidx]
    Xtest <- X[testidx, ]
    
    # Estimate beta from the training set
    betahat <- solve(crossprod(Xtrain), crossprod(Xtrain, ytrain))
    
    # Predict on the test set
    ytestpred <- Xtest %*% betahat
    
    # Compute error
    errvec[k] <- sum((ytestpred - ytest)^2)
  }
  # return CV(f)
  sum(errvec) /n
}
linreg_loo <- function(y, X) {
  # Leave-one-out CV using the one-step formula
  betahat <- solve(crossprod(X), crossprod(X, y))
  yhat <- X %*% betahat
  hii <- diag(X %*% solve(crossprod(X), t(X)))
  mean(((y - yhat) / (1 - hii))^2)
}

# put a copy of the dataset inside your "Documents" folder
# Load the data
data <- read_csv("STAT 444 project data RAW.csv")
colnames(data) <- c("date_1", "day", "month", "year", "temp2_c", "temp2_max_c",
                    "temp2_min_c", "temp2_ave_c", "surface_pressure_pa", 
                    "wind_speed50_max_m_s", "wind_speed50_min_m_s", 
                    "wind_speed50_ave_m_s", "prectotcorr", "total_demand_mw", 
                    "date_2", "max_generation_mw")
clean_data <- na.omit(data, cols = "total_demand_mw")
clean_data <- as.data.frame(clean_data)
data <- clean_data
data_scaled <- clean_data
data_scaled <- data_scaled %>% dplyr::select(-date_1, -date_2)
for (i in 1:(ncol(data_scaled))) data_scaled[ ,i] <- (data_scaled[ ,i] - mean(data_scaled[ ,i])) / sd(data_scaled[ ,i])
# Remove rows with missing 'total_demand(mw)'
clean_data <- as.data.frame(data_scaled)
clean_data <- clean_data %>% 
  dplyr::select(-temp2_c, -temp2_min_c, -wind_speed50_ave_m_s, -max_generation_mw)

# remove year and day too
#clean_data <- clean_data %>% 
#  dplyr::select(- year, -temp2_c, -temp2_min_c, -wind_speed50_ave_m_s, -max_generation_mw)

# Linear Regression --------------------------------------------------------------------------------------------------------------------------
# Fit the linear regression model using all other variables
# model <- lm(`total_demand(mw)` ~ . - date_1 - date_2, data = clean_data)
# model <- lm(total_demand.mw. ~ . 
#             - date_1 - date_2
#             - temp2.c. - `temp2_min.c.`
#             - `wind_speed50_ave.m.s.`
#             - `max_generation.mw.`, data = clean_data)
model <- lm(total_demand_mw ~ . , data = clean_data)

# Summary of the model
summary(model)

# 1. Residuals vs Fitted
plot1 <- ggplot(clean_data, aes_string(x = fitted(model), y = resid(model))) +
  geom_point() +
  geom_smooth(method = "loess", se = FALSE, color = "red") +
  labs(title = "Residuals vs Fitted", 
       x = "Fitted Values", 
       y = "Residuals")

# 2. Normal Q-Q
plot2 <- ggplot(clean_data, aes(sample = resid(model))) +
  stat_qq() +
  stat_qq_line() +
  labs(title = "Normal Q-Q", 
       x = "Theoretical Quantiles", 
       y = "Standardized Residuals")

# 3. Scale-Location
plot3 <- ggplot(clean_data, aes_string(x = fitted(model), y = sqrt(abs(resid(model))))) +
  geom_point() +
  geom_smooth(method = "loess", se = FALSE, color = "red") +
  labs(title = "Scale-Location", 
       x = "Fitted Values", 
       y = "Square Root of Absolute Residuals")

# 4. Residuals vs Leverage
plot4 <- ggplot(clean_data, aes_string(x = hatvalues(model), y = resid(model))) +
  geom_point() +
  labs(title = "Residuals vs Leverage", 
       x = "Leverage", 
       y = "Residuals") +
  geom_smooth(method = "loess", se = FALSE, color = "red")

library(gridExtra)
pdf("Linear_Regression_Residual_Analysis.pdf")
grid.arrange(plot1, plot2, plot3, plot4, ncol = 2, nrow = 2,
                           top = textGrob("Linear Regression Residual Analysis", 
                                          gp = gpar(fontsize = 15, fontface = "bold")))
dev.off()

# Calculate residuals
residuals <- residuals(model)

# Calculate MSE
mse <- mean(residuals^2)

# Print MSE
print(mse)
print(residuals)

# Prepare the data for cross-validation
# Ensure only relevant predictors are included as per the model specification
X <- model.matrix(total_demand_mw ~ ., data = clean_data)
y <- clean_data$total_demand_mw

# Specify the number of folds for k-fold cross-validation
# Get the number of observations (n)
n <- nrow(clean_data)

# Define the k values for CV
k_values <- c(5, 10, 20, n)  # Includes LOO as k=n

# Initialize a list or vector to store CV errors for each k
cv_errors <- numeric(length(k_values))

# Loop over the k values to perform CV
for (i in seq_along(k_values)) {
  K <- k_values[i]
  if (K == n) {
    # Leave-one-out cross-validation
    cv_errors[i] <- linreg_loo(y, X)  # Use the LOO function if different handling is needed
  } else {
    # k-fold cross-validation
    cv_errors[i] <- linreg_kfold(y, X, K)
  }
  cat("CV error for K =", K, ":", cv_errors[i], "\n")
}

# Optionally, you can print all CV errors together at the end
print(cv_errors)


# Saving the linear regression residual analysis plots


# WLS ---------------------------------------------------------------------

# Fit the initial linear regression model to get residuals
initial_model <- lm(total_demand_mw ~ ., data = clean_data)

# Define weights based on the residuals of the initial model
wt <- 1 / lm(abs(initial_model$residuals) ~ initial_model$fitted.values)$fitted.values^2

# Fit the WLS model using the calculated weights
wls_model <- lm(total_demand_mw ~ ., data = clean_data,
                weights = wt)

# Summary of the WLS model
summary(wls_model)

# Create diagnostic plots
# 1. Residuals vs Fitted
plot1 <- ggplot(clean_data, aes(x = fitted(wls_model), y = resid(wls_model))) +
  geom_point() +
  geom_smooth(method = "loess", se = FALSE, color = "red") +
  labs(title = "Residuals vs Fitted", 
       x = "Fitted Values", 
       y = "Residuals")

# 2. Normal Q-Q
plot2 <- ggplot(clean_data, aes(sample = resid(wls_model))) +
  stat_qq() +
  stat_qq_line() +
  labs(title = "Normal Q-Q", 
       x = "Theoretical Quantiles", 
       y = "Standardized Residuals")

# 3. Scale-Location
plot3 <- ggplot(clean_data, aes(x = fitted(wls_model), y = sqrt(abs(resid(wls_model))))) +
  geom_point() +
  geom_smooth(method = "loess", se = FALSE, color = "red") +
  labs(title = "Scale-Location", 
       x = "Fitted Values", 
       y = "Square Root of Absolute Residuals")

# 4. Residuals vs Leverage
plot4 <- ggplot(clean_data, aes(x = hatvalues(wls_model), y = resid(wls_model))) +
  geom_point() +
  labs(title = "Residuals vs Leverage", 
       x = "Leverage", 
       y = "Residuals") +
  geom_smooth(method = "loess", se = FALSE, color = "red")

# Arrange the plots in a 2x2 grid
library(gridExtra)
pdf("Weighted_Linear_Regression_Residual_Analysis.pdf")
grid.arrange(plot1, plot2, plot3, plot4, ncol = 2, nrow = 2,
                           top = textGrob("Weighted Linear Regression Residual Analysis", 
                                          gp = gpar(fontsize = 15, fontface = "bold")))

# Saving the weighted linear regression residual analysis plots
dev.off()
# Calculate residuals for the WLS model
wls_residuals <- residuals(wls_model)

# Calculate MSE for the WLS model
wls_mse <- mean(wls_residuals^2)
print(wls_mse)

# Fit the initial linear regression model to get residuals
initial_model <- lm(total_demand_mw ~ ., data = clean_data)

# Define weights based on the residuals of the initial model
wt <- 1 / lm(abs(initial_model$residuals) ~ initial_model$fitted.values)$fitted.values^2

# Create a new dataframe with weights
clean_data_wt <- clean_data
clean_data_wt$weights <- wt

# Define the trainControl object for k-fold cross-validation
train_control <- trainControl(method = "cv", number = 5)

# Define the formula for the model
formula <- total_demand_mw ~ .

# Train the WLS model using caret
wls_model_caret <- train(formula, data = clean_data_wt, method = "lm", 
                         weights = clean_data_wt$weights, trControl = train_control)
varImp(wls_model_caret)
# Print the results of cross-validation
print(wls_model_caret)
print(paste("Mean Squared Error: ", wls_model_caret$results$RMSE^2))
# KNN ------------------------------------------------------------------------
library(FNN)

# Prepare data for KNN regression
knn_data <- clean_data 

x <- as.matrix(knn_data %>% dplyr::select(-total_demand_mw))
y <- knn_data$total_demand_mw

# Define the trainControl object for k-fold cross-validation
train_control <- trainControl(method = "cv", number = 5)

# Create a data frame for the caret package
train_data <- data.frame(x, y)

# Perform k=5 k-fold cross-validation for KNN with different k values
set.seed(123) # for reproducibility

# Fit KNN model with k=5
knn_fit5 <- train(x, y, method = "knn",
                  tuneGrid = data.frame(k = 5),
                  trControl = train_control)

# Fit KNN model with k=21
knn_fit21 <- train(x, y, method = "knn",
                   tuneGrid = data.frame(k = 21),
                   trControl = train_control)

# Fit KNN model with k=51
knn_fit51 <- train(x, y, method = "knn",
                   tuneGrid = data.frame(k = 51),
                   trControl = train_control)

# Print the results of cross-validation
print(knn_fit5)
print(knn_fit21)
print(knn_fit51)

# Extract and print RMSE for each model
cat("MSE for k=5:", knn_fit5$results$RMSE^2, "\n")
cat("MSE for k=21:", knn_fit21$results$RMSE^2, "\n")
cat("MSE for k=51:", knn_fit51$results$RMSE^2, "\n")
# Combine the results into a single data frame for plotting
results <- rbind(knn_fit5$results, knn_fit21$results, knn_fit51$results)
results$k <- factor(c(rep(5, nrow(knn_fit5$results)),
                      rep(21, nrow(knn_fit21$results)),
                      rep(51, nrow(knn_fit51$results))),
                    levels = c(5, 21, 51))

# Plotting
ggplot(results, aes(x = k, y = RMSE, group = 1)) +
  geom_line() +
  geom_point() +
  labs(title = "Cross-Validation Performance Across Different k Values",
       x = "Number of Neighbors (k)", y = "Root Mean Squared Error (RMSE)")
# Saving the KNN cross-validation performance plot
pdf("KNN_CV_Performance.pdf")
ggplot(results, aes(x = k, y = RMSE, group = 1)) +
  geom_line() +
  geom_point() +
  labs(title = "Cross-Validation Performance Across Different k Values",
       x = "Number of Neighbors (k)", y = "Root Mean Squared Error (RMSE)")
dev.off()

# GAM -----------------------------------------------------------------------------------------------------------------------
knots <- 20

# Define the response variable
response_variable <- "total_demand_mw"

# Get all column names from the data except the response variable
covariates <- setdiff(names(clean_data), response_variable)

# Create the formula dynamically
smooth_terms <- paste0("s(", covariates, ", bs = 'bs', k = knots)", collapse = " + ")
formula <- as.formula(paste(response_variable, "~", smooth_terms))

# Fit the model using the dynamically created formula
gam_model <- mgcv::gam(formula, data = clean_data)
summary(gam_model)

cv_control <- trainControl(method = "cv", number = 5)
cv_gam <- train(total_demand_mw ~ temp2_max_c + surface_pressure_pa + wind_speed50_max_m_s + wind_speed50_min_m_s + prectotcorr,
                data = clean_data, method = "gamLoess",
                trControl = cv_control)

# Print the results of the cross-validation
print(cv_gam$results$RMSE^2)
cat("MSE for GAM:", cv_gam$results$RMSE^2, "\n")
par(mfrow=c(2, 3))  # Set layout to show multiple plots

# Residuals vs Fitted
plot(gam_model, residuals = TRUE, pch = 20, cex = 0.5, main = "Residuals vs. Fitted")

# Q-Q plot for residuals
qqnorm(resid(gam_model), main = "QQ Plot")
qqline(resid(gam_model), col = "red")

# Adding a common title across the top of all plots
mtext("GAM Residual vs Fitted + QQ Plot", side = 3, line = -2, outer = TRUE, cex = 1.3)
# Ridge regression -----------------------------------------------------------------------------------------
par(mfrow=c(1, 1))
X <- model.matrix(model)
glmnetridgecv <- cv.glmnet(X, clean_data$total_demand_mw, alpha = 0, nfolds = 5)
plot(glmnetridgecv)
minlambda <- glmnetridgecv$lambda.min
glmnetridge_nocv <- glmnet(X, clean_data$total_demand_mw, alpha = 0, nfolds = 5)
plot(glmnetridge_nocv, xvar = "lambda")
# Which variables do you think are those top curves?
round(t(glmnetridge_nocv$beta), 4)

glmnetridge_withcv <- glmnet(X, clean_data$total_demand_mw, alpha = 0, lambda = minlambda)
glmnetridge_withcv$beta # Coefficient estimates
cbind(glmnetridge_withcv$beta, coef(model))

ridge_model <- train(X, clean_data$total_demand_mw,
                     method = "glmnet",
                     trControl = train_control,
                     tuneGrid = expand.grid(.alpha = 0, .lambda = minlambda))
ridge_model$results$RMSE^2
# LASSO ---------------------------------------------------------------------------

glmnetlassocv <- cv.glmnet(X, clean_data$total_demand_mw, alpha = 1)
plot(glmnetlassocv)
minlambda <- glmnetlassocv$lambda.min
lambda1se <- glmnetlassocv$lambda.1se
glmnetlasso_nocv <- glmnet(X, clean_data$total_demand_mw, alpha = 1)
plot(glmnetlasso_nocv, xvar = "lambda")
cbind(
  coef(glmnetlasso_nocv, s = minlambda),
  coef(glmnetlasso_nocv, s = lambda1se)
)
lasso_model <- train(X, clean_data$total_demand_mw,
                     method = "glmnet",
                     trControl = train_control,
                     tuneGrid = expand.grid(.alpha = 1, .lambda = minlambda))
lasso_model$results$RMSE^2
# Elastic net ------------------------------------------------------------------------------------------------------------------------------------------------
# Define the function to optimize
elastic_net_cv <- function(alpha, lambda) {
  tryCatch({
    # Set up cross-validation control
    train_control <- trainControl(method = "cv", number = 5)
    
    # Train the model
    model <- train(
      X, y,
      method = "glmnet",
      trControl = train_control,
      tuneGrid = expand.grid(.alpha = alpha, .lambda = lambda)
    )
    
    # Return the negative RMSE as Bayesian Optimization maximizes the objective function
    list(Score = -min(model$results$RMSE), Pred = model$results$RMSE)
  }, error = function(e) {
    # Return a large negative score in case of an error to discourage this parameter set
    list(Score = -1e10, Pred = NA)
  })
}

# Define the bounds of the search space
bounds <- list(alpha = c(0, 1), lambda = c(0.001, 0.1))

# Run Bayesian Optimization
if (!exists(opt)){
  opt <- BayesianOptimization(
    FUN = elastic_net_cv,
    bounds = bounds,
    init_points = 10,
    n_iter = 50,
    acq = "ucb",
    kappa = 2.576,
    eps = 0.0,
    verbose = TRUE
  )
}

# Extract the best parameters
best_alpha <- opt$Best_Par["alpha"]
best_lambda <- opt$Best_Par["lambda"]

# Fit the final elastic net model using the best parameters
final_elastic_net_model <- glmnet(X, y, alpha = best_alpha, lambda = best_lambda)

# Coefficients of the final model
coefficients <- coef(final_elastic_net_model, s = best_lambda)
print(coefficients)


# Plotting the cross-validation results for elastic net
cv_results <- data.frame(lambda = elastic_net_model$results$lambda, RMSE = elastic_net_model$results$RMSE)
elastic_best_mse <- -opt$Best_Value
print(elastic_best_mse)

# PPR --------------------------------------------------------------------------------------------------------------------------------------------------------

ppr_model <- ppr(total_demand_mw~ ., data = clean_data,   nterms = 2, 
                 optlevel = 3, 
                 max.terms = 2,
                 sm.method = "gcvspline")

formula <- as.formula("total_demand_mw~ .")

# Function to perform cross-validation
cv_ppr <- function(data, formula, nterms, optlevel, max.terms, sm.method, k = 10) {
  n <- nrow(data)
  folds <- sample(rep(1:k, length.out = n))
  errors <- numeric(k)
  
  for (i in 1:k) {
    train_data <- data[folds != i, ]
    test_data <- data[folds == i, ]
    
    model <- ppr(formula, data = train_data, nterms = nterms, optlevel = optlevel, max.terms = max.terms, sm.method = sm.method)
    predictions <- predict(model, newdata = test_data)
    actuals <- test_data$total_demand_mw
    
    error <- mean((predictions - actuals)^2) # Mean Squared Error
    errors[i] <- error
  }
  
  return(mean(errors))
}

# Perform cross-validation
ppr_mse <- cv_ppr(clean_data, formula, nterms = 2, optlevel = 3, max.terms = 2, sm.method = "gcvspline", k = 10)
print(paste("Cross-Validation Error (MSE):", ppr_mse))
# multivariate time series ----------------------------------------------------------------------------------------------------------------------------------------
# Convert data to a time series object, here you should specify the frequency
# For example, if your data is daily, frequency would be 365, etc.
ts_data_orig <- clean_data
if("day" %in% colnames(ts_data_orig)) {
  # Remove the 'day' column from 'clean_data'
  ts_data_orig <- ts_data_orig[, !colnames(ts_data_orig) %in% "day"]
}
if("month" %in% colnames(ts_data_orig)) {
  # Remove the 'day' column from 'clean_data'
  ts_data_orig <- ts_data_orig[, !colnames(ts_data_orig) %in% "month"]
}
if("year" %in% colnames(ts_data_orig)) {
  # Remove the 'day' column from 'clean_data'
  ts_data_orig <- ts_data_orig[, !colnames(ts_data_orig) %in% "year"]
}

ts_data_orig$date_1 <- data$date_1

ts_data <- ts(ts_data_orig, frequency = 365, start=c(2018,1))

# Fit a VAR model, where 'p' is the order of the model
# You might need to determine the optimal 'p' using AIC or BIC, here I assume p=1 for simplicity
var_model <- VAR(ts_data, p = 1, type = "both")

# Calculate the initial length of the training set (50% of data)
initial_train_length <- floor(0.5 * nrow(ts_data))
test_length <- nrow(ts_data) - initial_train_length

# Function to convert index to (year, day of year)
index_to_year_day <- function(index) {
  year <- 2018 + (index - 1) %/% 365
  day <- (index - 1) %% 365 + 1
  return(c(year, day))
}

# Array to store forecast errors
errors <- numeric(test_length)

# Loop over each point in the test set, expanding the training set each time
for (i in 1:test_length) {
  # Expanding training set
  train_end_index <- initial_train_length + i - 1
  train_end <- index_to_year_day(train_end_index)
  train_set <- window(ts_data, end = train_end)
  
  # Next point to forecast
  test_index <- initial_train_length + i
  test_start_end <- index_to_year_day(test_index)
  test_set <- window(ts_data, start = test_start_end, end = test_start_end)
  
  # Fit VAR model to the training set
  var_model <- VAR(train_set, p = 1, type = "both")
  
  # Forecast the next point
  forecast_result <- predict(var_model, n.ahead = 1)
  
  # Assuming 'total_demand_mw' is the target variable; adjust index as needed based on your dataset
  forecasted_value <- forecast_result$fcst$total_demand_mw[1, "fcst"]
  actual_value <- test_set[1, "total_demand_mw"]
  
  # Calculate error for the forecast, here using Mean Squared Error (MSE)
  errors[i] <- (forecasted_value - actual_value)^2
}

# Calculate Mean Squared Error (MSE)
time_series_mse <- mean(errors)
print(paste("Mean Squared Error: ", time_series_mse))
# neural networks --------------------------------------------------------------------------------------------------------------------------------------------------
# check regression_NN.ipynb for MSE:
neural_network_mse = 0.98896531
# xgBoost ----------------------------------------------------------------------------------------------------------------------------------------------------------
# check regression_xbg.ipynb for MSE:
xgBoost_mse = 0.1603
# CV results comparison (so far):-----------------------------------------------------------------------------------------------------------------------------------
CV_results <- list(linear = cv_errors[1], # Linear model CV for K=5
                   wls = wls_model_caret$results$RMSE^2,
                   knn5 = knn_fit5$results$RMSE^2,
                   knn21 = knn_fit21$results$RMSE^2,
                   knn51 = knn_fit51$results$RMSE^2,
                   gam = cv_gam$results$RMSE^2, 
                   ridge = ridge_model$results$RMSE^2,
                   lasso = lasso_model$results$RMSE^2,
                   ppr = ppr_mse,
                   # time_series = time_series_mse,
                   neural_network = neural_network_mse,
                   elastic_xgBoost = xgBoost_mse,
                   elastic = elastic_best_mse)

# Convert the list to a data frame
models <- names(CV_results)
mse_values <- unlist(CV_results)
data_for_plot <- data.frame(Model = models, MSE = mse_values)

# Sort the data frame by MSE in ascending order
sorted_data_for_plot <- arrange(data_for_plot, MSE)

# Adjust the factor levels of 'Model' to follow the sorted order
sorted_data_for_plot$Model <- factor(sorted_data_for_plot$Model, levels = sorted_data_for_plot$Model)

# Create the plot
mse_plot <- ggplot(sorted_data_for_plot, aes(x = Model, y = MSE, fill = Model)) +
  geom_bar(stat = "identity", color = "black", show.legend = FALSE) +
  theme_minimal() +
  labs(title = "Comparison of MSE Across Models",
       x = "Model",
       y = "Mean Squared Error (MSE)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Improve label readability
mse_plot




# Saving the comparison of MSE across models
pdf("Model_MSE_Comparison3.pdf")
mse_plot
dev.off()

  

