# List of required packages
packages <- c("grid", "gridExtra", "readr", "dplyr", "Matrix", "mgcv", 
              "splines", "gamair", "caret", "ggplot2", "gam", "moments", "imputeTS")

# Function to check, install if missing, and load packages
install_and_load <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
  }
  # Load the package
  library(pkg, character.only = TRUE)
}

# Check, install, and load each package
invisible(lapply(packages, install_and_load))

data <- read_csv("STAT 444 project data RAW.csv")
colnames(data) <- c("date_1", "day", "month", "year", "temp2_c", "temp2_max_c",
                    "temp2_min_c", "temp2_ave_c", "surface_pressure_pa", 
                    "wind_speed50_max_m_s", "wind_speed50_min_m_s", 
                    "wind_speed50_ave_m_s", "prectotcorr", "total_demand_mw", 
                    "date_2", "max_generation_mw")
# for (i in 1:(ncol(data_scaled))) data_scaled[ ,i] <- (data_scaled[ ,i] - mean(data_scaled[ ,i])) / sd(data_scaled[ ,i])
clean_data <- data
clean_data <- na_interpolation(clean_data, option = "spline")
powerfun <- function(y, alpha) {
  y = y + 1e-8
  if (alpha == 0){
    y = log(y)
  }
  else if (alpha > 0) {
    y = y^alpha
  } else {
    y = -y^alpha
  }
}

alpha_vals = seq(-10, 10, 0.1)
best_alpha = c()
for (i in 1:(ncol(clean_data))){
  if (is.character(unname(unlist(clean_data[1,i])))){
    best_alpha = c(best_alpha, 1)
    next
  } 
  min_skew = Inf
  cur_best = 1
  for (alpha in alpha_vals){
    cat("alpha: ", alpha, " i: ", i, "\n")
    new_col = unlist(lapply(clean_data[ ,i], function(x) { return(powerfun(x, alpha)) }))
    # new_col = powerfun(clean_data[ ,i], alpha)
    new_col <- (new_col - mean(new_col)) / sd(new_col)
    cur_skew = (mean(new_col) - median(new_col))/sd(new_col)
    if (abs(cur_skew) < min_skew) {
      min_skew = abs(cur_skew)
      cur_best = alpha
    }
  }
  best_alpha = c(best_alpha, cur_best)
}

best_alpha <- setNames(as.list(best_alpha), colnames(clean_data))

best_alpha

















# Function to apply the best alpha transformation to each column
apply_power_transformation <- function(data, alphas) {
  transformed_data <- data
  for (i in seq_along(alphas)) {
    alpha <- alphas[[i]]
    column_name <- names(alphas)[i]
    if (alpha == 1) {
      transformed_data[[column_name]] <- data[[column_name]]
    }
    else transformed_data[[column_name]] <- powerfun(data[[column_name]], alpha)
  }
  return(transformed_data)
}

# Applying the transformation to the dataset
transformed_data <- apply_power_transformation(clean_data, best_alpha)

# Saving the transformed data to a CSV file
write_csv(transformed_data, "STAT_444_project_data_power_transformed.csv")

