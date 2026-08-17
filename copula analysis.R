# Copula-based joint modeling of crash count and risky driving behavior count

# 1. Load packages
required_packages <- c("glmmTMB", "readxl", "copula", "fitdistrplus", "writexl")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))]
if (length(missing_packages) > 0) stop(paste0("The following packages are required but not installed: ", paste(missing_packages, collapse = ", ")))

library(glmmTMB)
library(readxl)
library(copula)
library(fitdistrplus)
library(writexl)

# 2. Set file paths

input_file <- file.path("data", "example_data.xlsx")
output_dir <- "results"

if (!file.exists(input_file)) stop(paste0("The input file was not found: ", input_file))
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

output_excel <- file.path(output_dir, "Copula_joint_CDF_results.xlsx")
output_summary_FC <- file.path(output_dir, "model_summary_FC.txt")
output_summary_FB <- file.path(output_dir, "model_summary_FB.txt")
output_session_info <- file.path(output_dir, "sessionInfo.txt")

# 3. Read and validate data

data_raw <- read_excel(input_file)

required_columns <- c("TIME", "road", "FC", "FB", "V", "S", "LH", "S_std", "EX", "LE")
missing_columns <- setdiff(required_columns, names(data_raw))
if (length(missing_columns) > 0) stop(paste0("The following required columns are missing: ", paste(missing_columns, collapse = ", ")))

data <- data_raw[, required_columns]
data$TIME <- as.POSIXct(data$TIME, tz = "Asia/Shanghai")

numeric_columns <- c("FC", "FB", "V", "S", "LH", "S_std")
for (column in numeric_columns) data[[column]] <- as.numeric(data[[column]])

data$road <- as.factor(data$road)
data$EX <- as.factor(data$EX)
data$LE <- as.factor(data$LE)

number_before_removing_missing <- nrow(data)
data <- data[complete.cases(data[, required_columns]), ]
number_after_removing_missing <- nrow(data)

message("Number of observations before removing missing values: ", number_before_removing_missing)
message("Number of observations after removing missing values: ", number_after_removing_missing)

if (nrow(data) == 0) stop("No complete observations remain after removing missing values.")
if (any(data$FC < 0)) stop("FC contains negative values.")
if (any(data$FB < 0)) stop("FB contains negative values.")

# 4. Standardize continuous variables

variables_to_standardize <- c("V", "S", "LH", "S_std")
constant_variables <- variables_to_standardize[vapply(data[variables_to_standardize], function(x) is.na(sd(x)) || sd(x) == 0, FUN.VALUE = logical(1))]
if (length(constant_variables) > 0) stop(paste0("The following variables have zero variance: ", paste(constant_variables, collapse = ", ")))

for (column in variables_to_standardize) data[[column]] <- as.numeric(scale(data[[column]]))

# 5. Construct temporal variables

data$weekday <- factor(weekdays(data$TIME))
data$hour <- factor(format(data$TIME, "%H"))
data$weekday_hour <- interaction(data$weekday, data$hour, drop = TRUE)
data <- data[order(data$road, data$TIME), ]

time_levels <- sort(unique(data$TIME))
data$TIME_factor <- factor(data$TIME, levels = time_levels, ordered = TRUE)

road_time_check <- vapply(split(data$TIME, data$road), function(x) all(diff(as.numeric(x)) >= 0), FUN.VALUE = logical(1))
if (!all(road_time_check)) stop("TIME is not correctly ordered within one or more road segments.")

duplicated_road_time <- duplicated(data[, c("road", "TIME")])
if (any(duplicated_road_time)) stop(paste0("The input data contain ", sum(duplicated_road_time), " duplicated road-time observations."))

# 6. Fit the marginal model for crash count

message("Fitting the marginal model for crash count...")

fit_FC <- glmmTMB(FC ~ V + S + LH + EX + S_std + (1 | road) + (1 | road:weekday_hour) + (1 | road:hour) + ar1(TIME_factor + 0 | road), ziformula = ~ V + LE + EX + LH, data = data, family = nbinom1(link = "log"), na.action = na.fail)

print(summary(fit_FC))

# 7. Fit the marginal model for risky driving behavior count

message("Fitting the marginal model for risky driving behavior count...")

fit_FB <- glmmTMB(FB ~ V + S + LE + EX + LH + (1 | road) + ar1(TIME_factor + 0 | road) + (1 | road:weekday_hour) + (1 | road:hour), data = data, family = nbinom2(link = "log"), na.action = na.fail)

print(summary(fit_FB))

# 8. Save marginal-model summaries

capture.output(summary(fit_FC), file = output_summary_FC)
capture.output(summary(fit_FB), file = output_summary_FB)

# 9. Check model convergence

convergence_FC <- isTRUE(fit_FC$sdr$pdHess)
convergence_FB <- isTRUE(fit_FB$sdr$pdHess)

if (!convergence_FC) warning("The crash-count model may not have converged because the Hessian matrix is not positive definite.")
if (!convergence_FB) warning("The risky-driving-behavior model may not have converged because the Hessian matrix is not positive definite.")

# 10. Obtain predicted values

fitted_FC <- as.numeric(predict(fit_FC, type = "response"))
fitted_FB <- as.numeric(predict(fit_FB, type = "response"))

if (any(!is.finite(fitted_FC))) stop("The predicted crash-count values contain non-finite values.")
if (any(!is.finite(fitted_FB))) stop("The predicted risky-driving-behavior values contain non-finite values.")

# 11. Transform predicted values

log_offset <- 1e-4
log_fitted_FC <- log(fitted_FC + log_offset)
log_fitted_FB <- log(fitted_FB + log_offset)

fit_norm_FC <- fitdist(log_fitted_FC, "norm")
fit_norm_FB <- fitdist(log_fitted_FB, "norm")

mean_FC <- unname(fit_norm_FC$estimate["mean"])
sd_FC <- unname(fit_norm_FC$estimate["sd"])
mean_FB <- unname(fit_norm_FB$estimate["mean"])
sd_FB <- unname(fit_norm_FB$estimate["sd"])

if (!is.finite(sd_FC) || sd_FC <= 0) stop("The estimated standard deviation for crash predictions is invalid.")
if (!is.finite(sd_FB) || sd_FB <= 0) stop("The estimated standard deviation for risky-driving-behavior predictions is invalid.")

# 12. Calculate pseudo-observations

pseudo_FC <- pnorm(log_fitted_FC, mean = mean_FC, sd = sd_FC)
pseudo_FB <- pnorm(log_fitted_FB, mean = mean_FB, sd = sd_FB)

copula_epsilon <- 1e-6
pseudo_FC <- pmin(pmax(pseudo_FC, copula_epsilon), 1 - copula_epsilon)
pseudo_FB <- pmin(pmax(pseudo_FB, copula_epsilon), 1 - copula_epsilon)

u <- cbind(pseudo_FC, pseudo_FB)
colnames(u) <- c("pseudo_FC", "pseudo_FB")

if (any(!is.finite(u))) stop("The pseudo-observations contain non-finite values.")
if (any(u <= 0 | u >= 1)) stop("All pseudo-observations must lie strictly between 0 and 1.")

# 13. Define candidate Copula models

copula_models <- list(Gaussian = normalCopula(param = 0.1, dim = 2, dispstr = "un"), Clayton = claytonCopula(param = 1, dim = 2), Frank = frankCopula(param = 5, dim = 2), Gumbel = gumbelCopula(param = 1.5, dim = 2), Joe = joeCopula(param = 1.5, dim = 2))

start_values <- list(Gaussian = 0.1, Clayton = 1, Frank = 5, Gumbel = 1.5, Joe = 1.5)
optimization_methods <- c("L-BFGS-B", "Nelder-Mead", "BFGS", "CG")

# 14. Define the Copula fitting function

fit_one_copula <- function(copula_object, copula_name, observations, start_value, methods) {
  successful_fits <- list()
  
  for (optimization_method in methods) {
    message("Fitting ", copula_name, " Copula using ", optimization_method, "...")
    
    current_fit <- tryCatch(fitCopula(copula = copula_object, data = observations, method = "mpl", start = start_value, optim.method = optimization_method, hideWarnings = TRUE), error = function(e) {
      message("The ", copula_name, " Copula failed using ", optimization_method, ": ", conditionMessage(e))
      NULL
    })
    
    if (!is.null(current_fit)) {
      current_loglik <- tryCatch(as.numeric(logLik(current_fit)), error = function(e) NA_real_)
      if (is.finite(current_loglik)) successful_fits[[optimization_method]] <- current_fit
    }
  }
  
  if (length(successful_fits) == 0) {
    warning(paste0("All optimization methods failed for the ", copula_name, " Copula."))
    return(NULL)
  }
  
  loglik_values <- vapply(successful_fits, function(model) as.numeric(logLik(model)), FUN.VALUE = numeric(1))
  best_method <- names(which.max(loglik_values))
  list(fit = successful_fits[[best_method]], optimizer = best_method, logLik = unname(max(loglik_values)))
}

# 15. Fit candidate Copula models

copula_results <- lapply(names(copula_models), function(copula_name) fit_one_copula(copula_object = copula_models[[copula_name]], copula_name = copula_name, observations = u, start_value = start_values[[copula_name]], methods = optimization_methods))
names(copula_results) <- names(copula_models)

successful_models <- !vapply(copula_results, is.null, FUN.VALUE = logical(1))
copula_results <- copula_results[successful_models]

if (length(copula_results) == 0) stop("None of the candidate Copula models could be successfully fitted.")

# 16. Compare Copula models

number_of_observations <- nrow(u)

copula_comparison <- do.call(rbind, lapply(names(copula_results), function(copula_name) {
  current_result <- copula_results[[copula_name]]
  current_fit <- current_result$fit
  current_loglik <- current_result$logLik
  number_of_parameters <- length(coef(current_fit))
  current_AIC <- -2 * current_loglik + 2 * number_of_parameters
  current_BIC <- -2 * current_loglik + log(number_of_observations) * number_of_parameters
  estimated_parameters <- paste(round(coef(current_fit), 6), collapse = ", ")
  data.frame(Copula = copula_name, Parameter = estimated_parameters, Optimizer = current_result$optimizer, Log_likelihood = current_loglik, AIC = current_AIC, BIC = current_BIC, stringsAsFactors = FALSE)
}))

rownames(copula_comparison) <- NULL
copula_comparison <- copula_comparison[order(copula_comparison$AIC), ]

print(copula_comparison)

# 17. Select the best-fitting Copula

best_model_AIC <- copula_comparison$Copula[which.min(copula_comparison$AIC)]
best_model_BIC <- copula_comparison$Copula[which.min(copula_comparison$BIC)]

message("Best Copula based on AIC: ", best_model_AIC)
message("Best Copula based on BIC: ", best_model_BIC)

best_fit <- copula_results[[best_model_AIC]]$fit
best_copula <- best_fit@copula
best_parameter <- coef(best_fit)

message("Estimated parameter of the selected Copula: ", paste(round(best_parameter, 6), collapse = ", "))

# 18. Calculate dependence measures

kendall_tau <- tryCatch(as.numeric(tau(best_copula)), error = function(e) NA_real_)
tail_dependence <- tryCatch(lambda(best_copula), error = function(e) c(lower = NA_real_, upper = NA_real_))

lower_tail_dependence <- unname(tail_dependence[1])
upper_tail_dependence <- unname(tail_dependence[2])

selected_model_summary <- data.frame(Selection_criterion = c("AIC", "BIC"), Selected_Copula = c(best_model_AIC, best_model_BIC), stringsAsFactors = FALSE)

best_copula_parameters <- data.frame(Copula = best_model_AIC, Parameter = paste(round(best_parameter, 6), collapse = ", "), Kendalls_tau = kendall_tau, Lower_tail_dependence = lower_tail_dependence, Upper_tail_dependence = upper_tail_dependence, stringsAsFactors = FALSE)

# 19. Calculate the joint CDF

joint_CDF <- as.numeric(pCopula(u, copula = best_copula))

if (any(!is.finite(joint_CDF))) stop("The calculated joint CDF contains non-finite values.")
if (any(joint_CDF < 0 | joint_CDF > 1)) stop("The calculated joint CDF contains values outside [0, 1].")

# 20. Add results to the dataset

data$fitted_FC <- fitted_FC
data$fitted_FB <- fitted_FB
data$log_fitted_FC <- log_fitted_FC
data$log_fitted_FB <- log_fitted_FB
data$pseudo_FC <- pseudo_FC
data$pseudo_FB <- pseudo_FB
data$joint_CDF <- joint_CDF
data$best_copula <- best_model_AIC
data$weekday <- as.character(data$weekday)
data$hour <- as.character(data$hour)
data$weekday_hour <- as.character(data$weekday_hour)
data$TIME_factor <- as.character(data$TIME_factor)

# 21. Create model summaries

marginal_distribution_summary <- data.frame(
  Outcome = c("Crash count", "Risky driving behavior count"),
  Marginal_count_model = c("Zero-inflated negative binomial model, nbinom1", "Negative binomial model, nbinom2"),
  Transformation = rep(paste0("log(predicted value + ", log_offset, ")"), 2),
  Normal_mean = c(mean_FC, mean_FB),
  Normal_standard_deviation = c(sd_FC, sd_FB),
  Model_converged = c(convergence_FC, convergence_FB),
  stringsAsFactors = FALSE
)

print(summary(data$joint_CDF))
print(head(data[, c("TIME", "road", "FC", "FB", "fitted_FC", "fitted_FB", "pseudo_FC", "pseudo_FB", "joint_CDF", "best_copula")]))

# 22. Export results

write_xlsx(list(Joint_CDF_results = data, Copula_comparison = copula_comparison, Selected_model = selected_model_summary, Copula_parameters = best_copula_parameters, Marginal_distributions = marginal_distribution_summary), path = output_excel)

capture.output(sessionInfo(), file = output_session_info)

message("Analysis completed successfully.")
message("Selected Copula based on AIC: ", best_model_AIC)
message("Joint CDF results were saved to: ", output_excel)
message("Crash-count model summary was saved to: ", output_summary_FC)
message("Risky-driving-behavior model summary was saved to: ", output_summary_FB)
message("R session information was saved to: ", output_session_info)
