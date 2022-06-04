####### V1 Simulation study

######################################
# TO DO
# Change predictor variable to binary?
# Include different CI-levels
# Visualize results
######################################

######################################
### Steps
######################################
# Generate simulation data
# Split data
# Prediction model
# Calculate O:E ratio as reference
# Ampute data of validation set
# Multiple imputation of validaton set
# Predict values 
# Calculate O:E ratios for each imputed data set
# Pool
# Transform back
# Repeat x times and check for each round of simulation if ratio falls into true ratio

######################################
### Install and load required packages
######################################

## Specify required packages
packages = c("mice", 
             "lattice", 
             "metamisc",
             "dplyr",
             "ggplot2",
             "tidyverse",
             "caret",
             "MASS",
             "logit"
             )

## Install&load packages or load packages, if already installed
package.check <- lapply(
  packages,
  FUN = function(x) {
    if (!require(x, character.only = TRUE)) {
      install.packages(x, dependencies = TRUE)
      library(x, character.only = TRUE)
    }
  }
)

######################################
## Generate data
######################################

## Set seed to get reproducible results
set.seed(15)

## Set amount of study repetition
rep_amount = c(5)

## Initialize lists for results
result_per_sim = list()
result_table = list()
result_boolean_regular = vector()
result_boolean_log = vector()
result_boolean_sqrt = vector()

## Loop for repetition of study
for (j in 1:rep_amount){
  
  # Choose sample size
  n_sample = c(1000)

  # Continuous predictor variables
  continuous_var_1 = rnorm(n = n_sample, mean = 2, sd = 0.5)
  summary(continuous_var_1)
  #hist(continuous_var_1)
  
  continuous_var_2 = rnorm(n = n_sample, mean = 1.5, sd = 0.25)
  summary(continuous_var_2)
  #hist(continuous_var_2)
  
  # Create linear combination of predictor variables and intercept
  lin_combination_train = -4.5 + continuous_var_1 + continuous_var_2
  lin_combination_test = -4 + continuous_var_1 + continuous_var_2
  
  # Probability for response variable to be 1
  # Note: Due to application for logistic regression, use inverse logit function
  prob_outcome_train = 1/(1+exp(-lin_combination_train))
  prob_outcome_test = 1/(1+exp(-lin_combination_test))
  
  # Check that values are not approaching either 0 or 1 to avoid too deterministic approach
  # too close to zero?
  summary(prob_outcome_train)
  summary(prob_outcome_test)
  
  # Binary outcome variable as Bernoulli response variable
  # Note: Imbalanced Classification, desired probability for outcome_var = 1 between 20% and 30%
  outcome_var_train = rbinom(n = n_sample, size = 1, prob = prob_outcome_train)
  outcome_var_test = rbinom(n = n_sample, size = 1, prob = prob_outcome_test)
  #summary(outcome_var)
  #summary(outcome_var)
  #hist(outcome_var)
  
  # Combine it to a data frame
  train_data_complete = data.frame(
              continuous_var_1,
              continuous_var_2,
              outcome_var = outcome_var_train
              )
  
  test_data_complete = data.frame(
    continuous_var_1,
    continuous_var_2,
    outcome_var = outcome_var_test
  )

  ######################################
  ## Split data and develop prediction model
  ######################################
  
  # Split data into test and training data
  # split_prob = c(0.5)
  # training_samples_complete <- df_complete$outcome_var %>%
  #   createDataPartition(p = split_prob, list = FALSE)
  # 
  # train_data_complete = df_complete[training_samples_complete, ]
  # test_data_complete = df_complete[-training_samples_complete, ]

  # Create multiple logistic regression model based on complete data
  model_complete = glm(outcome_var ~ continuous_var_1 + continuous_var_2, 
                       data = train_data_complete,
                       family = binomial
  )
  summary(model_complete)$coef
  
  # # Fit model based on complete data
  # fit_model_complete <- model_complete %>% predict(test_data_complete,
  #                                                  type = "response"
  # )
  # # Predicted probability
  # summary(fit_model_complete)
  
  ######################################
  ## Calculate theoretical true O:E ratio
  ######################################
  
  # Define observed and expected events based on complete data
  sum(train_data_complete$outcome_var == 1)
  sum(test_data_complete$outcome_var == 1)
  n_observed_event = sum(test_data_complete$outcome_var == 1) #sum(test_data_complete$outcome_var == 1)
  n_expected_event = sum(train_data_complete$outcome_var == 1)#mean(fit_model_complete)*n_sample
  
  # O:E ratio with expected event calculated based on predicted probabilities
  reference_oe = oecalc(
    O = n_observed_event, # numeric vector of observed events
    E = n_expected_event, # numeric vector of expected events
    N = nrow(train_data_complete)
  )
  reference_oe
  
  ######################################
  ## Ampute test data 
  ######################################
  
  # Define parameters for amputation
  myprop = c(0.3) # Choose percentage of missing cells
  mypattern = rbind(c(0, 0, 1), c(1, 0 ,1), c(0, 1, 1)) # Choose missingness pattern (0 == missing, 1 == non-missing)
  myfreq = c(0.3, 0.35, 0.35) # Choose relative occurence of these patterns
  mymech = c("MAR") # Choose missingness mechanism (MAR based on literature)
  myweights = ampute.default.weights(mypattern, mymech) # Choose weights of weighted sum scores
  
  # Carry out amputation replacing it by NA
  # Note for interpretation: Proportion of missingness in terms of cells not in terms of cases
  amputation = ampute(data = test_data_complete,
                      patterns = mypattern,
                      prop = myprop, 
                      freq = myfreq, 
                      mech = mymech, 
                      weights = myweights, 
                      bycases = FALSE
                      )
  
  # save new data including amputations (NA) to data frame
  df_amputed = amputation$amp
  
  ######################################
  ## Evaluation of data amputation procedure
  ######################################
  
  # Show that amputation has worked
  head(amputation$amp)
  summary(amputation)
  
  # Visualization of missing data pattern
  md.pattern(amputation$amp)
  
  # Evaluation of amputation
  bwplot(amputation, which.pat = c(1, 3), descriptives = TRUE)
  xyplot(amputation, which.pat = 1)
  
  # Create multiple logistic regression model based on amputed data for comparison
  model_amputed = glm(outcome_var ~ continuous_var_1 + continuous_var_2, 
                      data = df_amputed, 
                      family = binomial
  )
  
  # Compare coefficients of amputed model vs. coefficients of complete model
  summary(model_amputed)$coef
  summary(model_complete)$coef
  
  ######################################
  ## Perform multiple imputation
  ######################################
  
  # Choose number of imputed data sets (reasonably between 5 and 10)
  imp_amount = c(10)
  
  # Choose number of iterations
  imp_iteration = c(20)
  
  # Choose imputation method
  # Interpretation: pmm default method for numeric data; logreg method for binary data
  # change to norm for smaller CI in results
  imp_method = c("pmm", "pmm", "logreg") 
  
  # Outcome variable as categorical variable for imp_method logreg
  df_amputed$outcome_var <- as.factor(df_amputed$outcome_var)
  
  # Impute data via mice
  imputation = mice(data = df_amputed, 
                m = imp_amount, 
                maxit = imp_iteration,
                method = imp_method, 
                print = TRUE,
                seed = 18
                )
  
  # Inspect convergence
  plot(imputation)
  
  # Save the imputed data sets
  # Note: alternative save as list of df
  imputed_all = complete(imputation, "long")
  
  ######################################
  ## Fit models and calculate (transformed) O:E ratios based on imputed data sets
  ######################################
  
  # Initialize lists for O:E ratios based on imputed data sets
  regular_oe_list = list()
  log_oe_list = list()
  sqrt_oe_list = list()
  
  # Loop for all imputed data sets
  for (i in 1:imp_amount){
    # Subset imputed data sets
    subset_imputed = subset(imputed_all, .imp == i)
    subset_imputed = dplyr::select(subset_imputed, 
                                   continuous_var_1, 
                                   continuous_var_2, 
                                   outcome_var)
    # Fit model on imputed data sets
    fit_imp = model_complete %>% predict(subset_imputed, 
                               type = "response")
    summary(fit_imp)
    # Calculate regular O:E and save it to list
    regular_oe_list[[i]] = oecalc(
      O = sum(subset_imputed$outcome_var == 1),
      E = mean(fit_imp)*n_sample, #*split_prob,
      N = nrow(subset_imputed)
    )
    # Calculate log O:E and save it to list
    log_oe_list[[i]] = oecalc(
      O = sum(subset_imputed$outcome_var == 1),
      E = mean(fit_imp)*n_sample, #*split_prob,
      N = nrow(subset_imputed),
      g = "log(OE)"
    )
    # Calculate square root O:E and save it to list
    sqrt_oe_list[[i]] = oecalc(
      O = sum(subset_imputed$outcome_var == 1),
      E = mean(fit_imp)*n_sample, #*split_prob,
      N = nrow(subset_imputed),
      g = "sqrt(OE)"
    )
  }
  
  ######################################
  ## Pooling procedure: Pool estimates
  ######################################
  
  # Initialize vectors for theta and theta.se
  regular_theta <- vector()
  log_theta <- vector()
  sqrt_theta <- vector()
  regular_theta_se <- vector()
  log_theta_se <- vector()
  sqrt_theta_se <- vector()
  
  # Save the o:e values to the vectors
  for (i in 1:imp_amount){
    # theta
    regular_theta[i] <- regular_oe_list[[i]]$theta
    log_theta[i] <- log_oe_list[[i]]$theta
    sqrt_theta[i] <- sqrt_oe_list[[i]]$theta
    # theta.se
    regular_theta_se[i] <- regular_oe_list[[i]]$theta.se
    log_theta_se[i] <- log_oe_list[[i]]$theta.se
    sqrt_theta_se[i] <- sqrt_oe_list[[i]]$theta.se
  }
  
  # Pool the o:e estimates (theta) using the mean
  pooled_theta_regular_oe = mean(regular_theta)
  pooled_theta_log_oe = mean(log_theta)
  pooled_theta_sqrt_oe = mean(sqrt_theta)
  
  ######################################
  ## Pooling procedure: Pool standard errors
  ######################################
  
  # Function to pool standard errors and receive respective t-values
  # Source: Epi and Big Data course
  # Note: Write new with less redundance as some variables already defined
  pooled_se <- function(est, se, n.imp){
    Qbar <- mean(est)
    # within-variance
    U <- sum(se**2)/n.imp 
    # between variance
    B <- sum((est - Qbar)**2)/(n.imp-1) 
    se_total <- sqrt(U + (1+1/n.imp)*B)
    r <- (1 + 1 / n.imp) * (B / U)
    v <- (n.imp - 1) * (1 + (1/r))^2
    t <- qt(0.975, v)
    res <- c(se_total, t)
    return(res)
  }
  
  # Pool the standard errors of the o:e estimates (theta.se)
  pooled_se_regular_oe = pooled_se(est = regular_theta, 
                                   se = regular_theta_se, 
                                   n.imp = imp_amount
                                     )
  
  pooled_se_log_oe = pooled_se(est = log_theta, 
                               se = log_theta_se,
                               n.imp = imp_amount
                               )
    
  pooled_se_sqrt_oe = pooled_se(est = sqrt_theta, 
                                se = sqrt_theta_se,
                                n.imp = imp_amount
                                )
  
  ######################################
  ## Pooling procedure: Pool Confidence intervals
  ######################################
  
  # Function to calculate confidence interval
  pooled_ci <- function(est, se, tvalue){
    theta.cilb <- est - tvalue * se
    theta.ciub <- est + tvalue * se
    ci <- c(theta.cilb, theta.ciub)
    return(ci)
  }
  
  # Pool the confidence intervals of the o:e estimates (theta.cilb and theta.ciub)
  pooled_ci_regular_oe = pooled_ci(est = pooled_theta_regular_oe,
                                     tvalue = pooled_se_regular_oe[2],
                                     se = pooled_se_regular_oe[1])
  
  pooled_ci_log_oe = pooled_ci(est = pooled_theta_log_oe,
                                   tvalue = pooled_se_log_oe[2],
                                   se = pooled_se_log_oe[1])
  
  pooled_ci_sqrt_oe = pooled_ci(est = pooled_theta_sqrt_oe,
                               tvalue = pooled_se_sqrt_oe[2],
                               se = pooled_se_sqrt_oe[1])
  
  ######################################
  ## Back-transformation of O:E ratios to original scale
  ######################################
  
  final_regular_oe = cbind(
    back_theta_regular_oe = pooled_theta_regular_oe,
    back_se_regular_oe = pooled_se_regular_oe[1],
    back_cilb_regular_oe = pooled_theta_regular_oe - pooled_se_regular_oe[2]*pooled_se_regular_oe[1],
    back_ciub_regular_oe = pooled_theta_regular_oe + pooled_se_regular_oe[2]*pooled_se_regular_oe[1]
    )
  colnames(final_regular_oe) <- c("theta", "theta.se", "theta.cilb", "theta.ciub")
  rownames(final_regular_oe) <- ("regular o:e ratio")
  
  final_log_oe = cbind(
    back_theta_log_oe = exp(pooled_theta_log_oe), 
    back_se_log_oe = pooled_se_log_oe[1], #exp(pooled_se_log_oe[1]),
    back_cilb_log_oe = exp(pooled_theta_log_oe - pooled_se_log_oe[2]*pooled_se_log_oe[1]),
    back_ciub_log_oe = exp(pooled_theta_log_oe + pooled_se_log_oe[2]*pooled_se_log_oe[1])
    )
  colnames(final_log_oe) <- c("theta", "theta.se", "theta.cilb", "theta.ciub")
  rownames(final_log_oe) <- ("log o:e ratio")

  final_sqrt_oe = cbind(
    back_theta_sqrt_oe =(pooled_theta_sqrt_oe)^2, 
    back_se_sqrt_oe = pooled_se_sqrt_oe[1], #(pooled_se_sqrt_oe[1])^2, 
    back_cilb_sqrt_oe = (pooled_theta_sqrt_oe - pooled_se_sqrt_oe[2]*pooled_se_sqrt_oe[1])^2,
    back_ciub_sqrt_oe = (pooled_theta_sqrt_oe + pooled_se_sqrt_oe[2]*pooled_se_sqrt_oe[1])^2
    )
  colnames(final_sqrt_oe) <- c("theta", "theta.se", "theta.cilb", "theta.ciub")
  rownames(final_sqrt_oe) <- ("square root o:e ratio")
  
  ######################################
  ## Compare performance of O:E measures with reference O:E measure
  ######################################
  
  # Results Table: Comparison of theta, theta.se and CI of different transformations
  true_oe_ratio <- as.list(reference_oe)[1:4]
  
  oe_comparison = rbind(true_oe_ratio, 
                        final_regular_oe, 
                        final_log_oe, 
                        final_sqrt_oe
                        )
  oe_comparison = as.data.frame(oe_comparison)
  
  # Results Boolean vector: Check if reference O:E lies in the 95% CI of the pooled O:E measure
  res_regular = isTRUE(final_regular_oe[3] < reference_oe$theta & reference_oe$theta < final_regular_oe[4])
  res_log = isTRUE(final_log_oe[3] < reference_oe$theta & reference_oe$theta < final_log_oe[4])
  res_sqrt = isTRUE(final_sqrt_oe[3] < reference_oe$theta & reference_oe$theta < final_sqrt_oe[4])
  
  # Results Visualization:
  # tbd
  
  ######################################
  ## Study repetition
  ######################################
  
  # Save results for study repetition
  result_table[[j]] = oe_comparison
  result_per_sim[[j]] = c(res_regular, res_log, res_sqrt)
  result_boolean_regular[[j]] = c(res_regular)
  result_boolean_log[[j]] = c(res_log)
  result_boolean_sqrt[[j]] = c(res_sqrt)
  
}

######################################
## Across study repetitions: Percentage of how often true O:E falls into 95% CI of respective O:E measure 
######################################

# Print tables of each round
print(result_table)
# print(result_per_sim)
# print(result_boolean_regular)
# print(result_boolean_log)
# print(result_boolean_sqrt)

# Interpretation: Percentage of how often true O:E falls into 95% CI of respective O:E measure 
percentage_coverage_regular = table(result_boolean_regular)["TRUE"]/rep_amount*100
percentage_coverage_log = table(result_boolean_log)["TRUE"]/rep_amount*100
percentage_coverage_sqrt = table(result_boolean_sqrt)["TRUE"]/rep_amount*100

# Print percentage of coverage
percentage_coverage_regular
percentage_coverage_log
percentage_coverage_sqrt

################################