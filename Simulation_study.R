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
rep_amount = c(50)

## Initialize lists for results
result_table_90 = list()
result_table_95 = list()
result_table_99 = list()

result_boolean_regular_90 = vector()
result_boolean_regular_95 = vector()
result_boolean_regular_99 = vector()

result_boolean_log_90 = vector()
result_boolean_log_95 = vector()
result_boolean_log_99 = vector()

result_boolean_sqrt_90 = vector()
result_boolean_sqrt_95 = vector()
result_boolean_sqrt_99 = vector()


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
  mypattern = rbind(c(0, 0, 1), c(1, 0 ,1), c(0, 1, 0)) # Choose missingness pattern (0 == missing, 1 == non-missing)
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
  ## Pooling procedure
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
    t_90 <- qt(0.950, v)
    t_95 <- qt(0.975, v)
    t_99 <- qt(0.995, v)
    res <- c(se_total, t_90, t_95, t_99)
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
  ## Back-transformation of O:E ratios to original scale
  ######################################
  
  # Initialize
  #[2] == 90%CI, [3] == 95%CI, [4] == 99%CI
  final_regular_oe = list()
  final_log_oe = list()
  final_sqrt_oe = list()
  
  for (k in 2:4){
    final_regular_oe[[k]] = cbind(
      theta = pooled_theta_regular_oe,
      theta.se = pooled_se_regular_oe[1],
      theta.cilb = pooled_theta_regular_oe - pooled_se_regular_oe[k]*pooled_se_regular_oe[1],
      theta.ciub = pooled_theta_regular_oe + pooled_se_regular_oe[k]*pooled_se_regular_oe[1]
      )
    
    final_log_oe[[k]] = cbind(
      theta = exp(pooled_theta_log_oe), 
      theta.se = pooled_se_log_oe[1],
      theta.cilb = exp(pooled_theta_log_oe - pooled_se_log_oe[k]*pooled_se_log_oe[1]),
      theta.ciub = exp(pooled_theta_log_oe + pooled_se_log_oe[k]*pooled_se_log_oe[1])
      )
  
    final_sqrt_oe[[k]] = cbind(
      theta =(pooled_theta_sqrt_oe)^2, 
      theta.se = pooled_se_sqrt_oe[1],
      theta.cilb = (pooled_theta_sqrt_oe - pooled_se_sqrt_oe[k]*pooled_se_sqrt_oe[1])^2,
      theta.ciub = (pooled_theta_sqrt_oe + pooled_se_sqrt_oe[k]*pooled_se_sqrt_oe[1])^2
      )
  }
  
  ######################################
  ## Compare performance of O:E measures with reference O:E measure
  ######################################
  
  # Subsample true oe ratio
  true_oe_ratio <- as.list(reference_oe)[1:4]
  
  # Results Table for 90% CI
  oe_comparison_90 = rbind(true_oe_ratio, 
                        final_regular_oe[[2]], 
                        final_log_oe[[2]], 
                        final_sqrt_oe[[2]]
                        )
  oe_comparison_90 = as.data.frame(oe_comparison_90)
  rownames(oe_comparison_90) <- c("true o:e", "regular o:e", "ln o:e", "sqrt o:e")
  
  # Results Table for 95% CI
  oe_comparison_95 = rbind(true_oe_ratio, 
                        final_regular_oe[[3]], 
                        final_log_oe[[3]], 
                        final_sqrt_oe[[3]]
  )
  oe_comparison_95 = as.data.frame(oe_comparison_95)
  rownames(oe_comparison_95) <- c("true o:e", "regular o:e", "ln o:e", "sqrt o:e")
  
  # Results Table for 99% CI
  oe_comparison_99 = rbind(true_oe_ratio, 
                           final_regular_oe[[4]], 
                           final_log_oe[[4]], 
                           final_sqrt_oe[[4]]
  )
  oe_comparison_99 = as.data.frame(oe_comparison_99)
  rownames(oe_comparison_99) <- c("true o:e", "regular o:e", "ln o:e", "sqrt o:e")
  
  # Boolean vector Results for 90% CI: Check if reference O:E lies in the 95% CI of the pooled O:E measure
  res_regular_90 = isTRUE(final_regular_oe[[2]][3] < reference_oe$theta & reference_oe$theta < final_regular_oe[[2]][4])
  res_log_90 = isTRUE(final_log_oe[[2]][3] < reference_oe$theta & reference_oe$theta < final_log_oe[[2]][4])
  res_sqrt_90 = isTRUE(final_sqrt_oe[[2]][3] < reference_oe$theta & reference_oe$theta < final_sqrt_oe[[2]][4])
  
  # Boolean vector Results for 95% CI: Check if reference O:E lies in the 95% CI of the pooled O:E measure
  res_regular_95 = isTRUE(final_regular_oe[[3]][3] < reference_oe$theta & reference_oe$theta < final_regular_oe[[3]][4])
  res_log_95 = isTRUE(final_log_oe[[3]][3] < reference_oe$theta & reference_oe$theta < final_log_oe[[3]][4])
  res_sqrt_95 = isTRUE(final_sqrt_oe[[3]][3] < reference_oe$theta & reference_oe$theta < final_sqrt_oe[[3]][4])
  
  # Boolean vector Results for 99% CI: Check if reference O:E lies in the 95% CI of the pooled O:E measure
  res_regular_99 = isTRUE(final_regular_oe[[4]][3] < reference_oe$theta & reference_oe$theta < final_regular_oe[[4]][4])
  res_log_99 = isTRUE(final_log_oe[[4]][3] < reference_oe$theta & reference_oe$theta < final_log_oe[[4]][4])
  res_sqrt_99 = isTRUE(final_sqrt_oe[[4]][3] < reference_oe$theta & reference_oe$theta < final_sqrt_oe[[4]][4])
  
  ######################################
  ## Study repetition
  ######################################
  
  # Save results for study repetition
  result_table_90[[j]] = oe_comparison_90
  result_table_95[[j]] = oe_comparison_95
  result_table_99[[j]] = oe_comparison_99

  result_boolean_regular_90[[j]] = c(res_regular)
  result_boolean_regular_95[[j]] = c(res_regular)
  result_boolean_regular_99[[j]] = c(res_regular)
  
  result_boolean_log_90[[j]] = c(res_log)
  result_boolean_log_95[[j]] = c(res_log)
  result_boolean_log_99[[j]] = c(res_log)
  
  result_boolean_sqrt_90[[j]] = c(res_sqrt)
  result_boolean_sqrt_95[[j]] = c(res_sqrt)
  result_boolean_sqrt_99[[j]] = c(res_sqrt)
}

######################################
## Across study repetitions: Percentage of how often true O:E falls into 95% CI of respective O:E measure 
######################################

# Print tables of each round
print(result_table_90)
print(result_table_95)
print(result_table_99)
# print(result_per_sim)
# print(result_boolean_regular)
# print(result_boolean_log)
# print(result_boolean_sqrt)

# Interpretation: Percentage of how often true O:E falls into 95% CI of respective O:E measure 
percentage_coverage_regular_90 = table(result_boolean_regular_90)["TRUE"]/rep_amount*100
percentage_coverage_log_90 = table(result_boolean_log_90)["TRUE"]/rep_amount*100
percentage_coverage_sqrt_90 = table(result_boolean_sqrt_90)["TRUE"]/rep_amount*100


percentage_coverage_regular_95 = table(result_boolean_regular_95)["TRUE"]/rep_amount*100
percentage_coverage_log_95 = table(result_boolean_log_95)["TRUE"]/rep_amount*100
percentage_coverage_sqrt_95 = table(result_boolean_sqrt_95)["TRUE"]/rep_amount*100

percentage_coverage_regular_99 = table(result_boolean_regular_99)["TRUE"]/rep_amount*100
percentage_coverage_log_99 = table(result_boolean_log_99)["TRUE"]/rep_amount*100
percentage_coverage_sqrt_99 = table(result_boolean_sqrt_99)["TRUE"]/rep_amount*100

# Print percentage of coverage
percentage_coverage_regular_90
percentage_coverage_log_90
percentage_coverage_sqrt_90

percentage_coverage_regular_95
percentage_coverage_log_95
percentage_coverage_sqrt_95

percentage_coverage_regular_99
percentage_coverage_log_99
percentage_coverage_sqrt_99

################################

# Results Visualization:
# tbd