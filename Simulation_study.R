####### V1 Simulation study 

######################################
### Steps
######################################
# Generate simulation data
# Split data
# Prediction model
# Calculate O/E ratio as reference
# Ampute data of validation set
# Multiple imputation of validaton set
# Predict values 
# Calculate O/E ratios for each imputed data set
# Pool
# Transform back
# Repeat x times and check for each round of simulation if ratio falls into true ratio
# Visualize results

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
             "ggplot2",
             "gridExtra"
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
set.seed(18)

## Set amount of study repetition
rep_amount = c(1000)

## Initialize lists for results
result_table = list()
all_true_oe = list()

all_regular_theta = list()
all_log_theta = list()
all_sqrt_theta = list()

distribution_prior_reg = list()
distribution_prior_log = list()
distribution_prior_sqrt = list()

## Loop for repetition of study
for (j in 1:rep_amount){
  
  # Choose sample size
  n_sample = c(1000)

  # Continuous predictor variables
  continuous_var_1 = rnorm(n = n_sample, mean = 2, sd = 0.5)
  #summary(continuous_var_1)
  #hist(continuous_var_1)
  
  continuous_var_2 = rnorm(n = n_sample, mean = 1.5, sd = 0.25)
  #summary(continuous_var_2)
  #hist(continuous_var_2)
  
  # Create linear combination of predictor variables and intercept
  lin_combination_train = -3.05 + continuous_var_1 + continuous_var_2
  lin_combination_test = -4.37 + continuous_var_1 + continuous_var_2
  
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
  
  # Additional noise
  # noise_param.1 = sample(seq_len(floor(length(outcome_var_train[outcome_var_train==1])*0.25)))
  # noise_param.0 = sample(seq_len(floor(length(outcome_var_train[outcome_var_train==0])*0.25)))
  # 
  # noise_param.1 = sample(seq_len(floor(length(outcome_var_test[outcome_var_test==1])*0.00)))
  # noise_param.0 = sample(seq_len(floor(length(outcome_var_test[outcome_var_test==0])*0.00)))
  
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
  
  # train_data_complete$outcome_var[noise_param.0] = 1
  # train_data_complete$outcome_var[noise_param.1] = 0
  # 
  # test_data_complete$outcome_var[noise_param.0] = 1
  # test_data_complete$outcome_var[noise_param.1] = 0

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
  ## Calculate true O/E ratio
  ######################################
  
  # Define observed and expected events based on complete data
  n_observed_event = sum(test_data_complete$outcome_var == 1) #sum(test_data_complete$outcome_var == 1)
  n_expected_event = sum(train_data_complete$outcome_var == 1)#mean(fit_model_complete)*n_sample
  
  # O/E ratio with expected event calculated based on predicted probabilities
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
  mypattern = rbind(c(1, 0, 1), c(0, 1, 1), c(0, 0, 1)) # Choose missingness pattern (0 == missing, 1 == non-missing)
  myfreq = c(1/3, 1/3, 1/3) # Choose relative occurence of these patterns
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
                print = TRUE
                )
  
  # Inspect convergence
  plot(imputation)
  
  # Save the imputed data sets
  # Note: alternative save as list of df
  imputed_all = complete(imputation, "long")
  
  ######################################
  ## Fit models and calculate (transformed) O/E ratios based on imputed data sets
  ######################################
  
  # Initialize lists for O/E ratios based on imputed data sets
  regular_oe_list = list()
  log_oe_list = list()
  sqrt_oe_list = list()
  
  regular_oe_theta = vector()
  log_oe_theta = vector()
  sqrt_oe_theta = vector()
  
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
    
    # Calculate regular O/E and save it to list
    regular_oe_list[[i]] = oecalc(
      O = sum(subset_imputed$outcome_var == 1),
      E = mean(fit_imp)*n_sample, #*split_prob,
      N = nrow(subset_imputed)
    )
    
    regular_oe_theta[[i]] = regular_oe_list[[i]][1]
    
    # Calculate log O/E and save it to list
    log_oe_list[[i]] = oecalc(
      O = sum(subset_imputed$outcome_var == 1),
      E = mean(fit_imp)*n_sample, #*split_prob,
      N = nrow(subset_imputed),
      g = "log(OE)"
    )
    
    log_oe_theta[[i]] = log_oe_list[[i]][1]
    
    # Calculate square root O/E and save it to list
    sqrt_oe_list[[i]] = oecalc(
      O = sum(subset_imputed$outcome_var == 1),
      E = mean(fit_imp)*n_sample, #*split_prob,
      N = nrow(subset_imputed),
      g = "sqrt(OE)"
    )
    
    sqrt_oe_theta[[i]] = sqrt_oe_list[[i]][1]
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
  
  # Save the O/E values to the vectors
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
  
  # Pool the O/E estimates (theta) using the mean
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
    res <- c(t_90, t_95, t_99, se_total)
    return(res)
  }
  
  # Pool the standard errors of the O/E estimates (theta.se)
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
  ## Back-transformation of O/E ratios to original scale
  ######################################
  
  # Loop over confidence intervals
  # Initialize [1] == 90%CI, [2] == 95%CI, [3] == 99%CI
  oe_comparison = list() 
  final_regular_oe = list()
  final_log_oe = list()
  final_sqrt_oe = list()
  
  for (k in 1:3){
    final_regular_oe[[k]] = cbind(
      theta = pooled_theta_regular_oe,
      theta.se = pooled_se_regular_oe[4],
      theta.cilb = pooled_theta_regular_oe - pooled_se_regular_oe[k]*pooled_se_regular_oe[4],
      theta.ciub = pooled_theta_regular_oe + pooled_se_regular_oe[k]*pooled_se_regular_oe[4]
      )
    
    final_log_oe[[k]] = cbind(
      theta = exp(pooled_theta_log_oe), 
      theta.se = pooled_se_log_oe[4],
      theta.cilb = exp(pooled_theta_log_oe - pooled_se_log_oe[k]*pooled_se_log_oe[4]),
      theta.ciub = exp(pooled_theta_log_oe + pooled_se_log_oe[k]*pooled_se_log_oe[4])
      )
  
    final_sqrt_oe[[k]] = cbind(
      theta =(pooled_theta_sqrt_oe)^2, 
      theta.se = pooled_se_sqrt_oe[4],
      theta.cilb = (pooled_theta_sqrt_oe - pooled_se_sqrt_oe[k]*pooled_se_sqrt_oe[4])^2,
      theta.ciub = (pooled_theta_sqrt_oe + pooled_se_sqrt_oe[k]*pooled_se_sqrt_oe[4])^2
      )
  
  ######################################
  ## Compare performance of O/E measures with reference O/E measure
  ######################################
  
  # Sub-sample of reference_oe to theta
  final_true_oe = as.list(reference_oe)[1:4]
      
  # Results Table with [[1]] == 90%, [[2]] == 95%, [[3]] == 99%, 
  oe_comparison[[k]] = rbind( 
                        final_true_oe,
                        final_regular_oe[[k]], 
                        final_log_oe[[k]], 
                        final_sqrt_oe[[k]]
                        )
  
  oe_comparison[[k]] = as.data.frame(oe_comparison[[k]])
  rownames(oe_comparison[[k]]) <- c("true o:e", "regular o:e", "ln o:e", "sqrt o:e")
  
  }  

  ######################################
  ## Study repetition
  ######################################

  # Save general results
  result_table[[j]] = oe_comparison
  all_true_oe[[j]] = final_true_oe$theta

  # Save theta to calculate bias
  all_regular_theta[[j]] = pooled_theta_regular_oe
  all_log_theta[[j]] = (exp(pooled_theta_log_oe))
  all_sqrt_theta[[j]] = ((pooled_theta_sqrt_oe)^2)
  
  # Save transformed O/E ratios based on imputed data set prior to pooling to check for improvement of normality shape
  distribution_prior_reg[[j]] = regular_oe_theta
  distribution_prior_log[[j]] = log_oe_theta
  distribution_prior_sqrt[[j]] = sqrt_oe_theta
}

######################################
## Calculations after the Simulation study
######################################

# Print tables of transformed O/E ratios of each iteration
print(result_table)

# Calculate true O/E ratio as mean of 1:rep_amount reference O/E ratio
mean_true_oe = (sum(unlist(all_true_oe)))/rep_amount

# Initialize vectors for loop
res_regular_90 = vector()
res_log_90 = vector()
res_sqrt_90 = vector()

res_regular_95 = vector()
res_log_95 = vector()
res_sqrt_95 = vector()

res_regular_99 = vector()
res_log_99 = vector()
res_sqrt_99 = vector()

res_regular_all = list()
res_log_all = list()
res_sqrt_all = list()

# Assign true value to vector, if respective estimate falls into respective CI
# Note: result_table[[rep_amount (j)]][[CI level (k)]][row, column] 
for (j in 1:rep_amount){
    # with respect to 90-CI
    res_regular_90[j] = isTRUE(result_table[[j]][[1]][2,3] < mean_true_oe & mean_true_oe < result_table[[j]][[1]][2,4])
    res_log_90[j] = isTRUE(result_table[[j]][[1]][3,3] < mean_true_oe & mean_true_oe < result_table[[j]][[1]][3,4])
    res_sqrt_90[j] = isTRUE(result_table[[j]][[1]][4,3] < mean_true_oe & mean_true_oe < result_table[[j]][[1]][4,4])
    # with respect to 95-CI
    res_regular_95[j] = isTRUE(result_table[[j]][[2]][2,3] < mean_true_oe & mean_true_oe < result_table[[j]][[2]][2,4])
    res_log_95[j] = isTRUE(result_table[[j]][[2]][3,3] < mean_true_oe & mean_true_oe < result_table[[j]][[2]][3,4])
    res_sqrt_95[j] = isTRUE(result_table[[j]][[2]][4,3] < mean_true_oe & mean_true_oe < result_table[[j]][[2]][4,4])
    # with respect to 99-CI
    res_regular_99[j] = isTRUE(result_table[[j]][[3]][2,3] < mean_true_oe & mean_true_oe < result_table[[j]][[3]][2,4])
    res_log_99[j] = isTRUE(result_table[[j]][[3]][3,3] < mean_true_oe & mean_true_oe < result_table[[j]][[3]][3,4])
    res_sqrt_99[j] = isTRUE(result_table[[j]][[3]][4,3] < mean_true_oe & mean_true_oe < result_table[[j]][[3]][4,4])
}

# Calculate percentage of how often true O/E falls into CI of respective O/E measure at level 0.9, 0.95, 0.99
percentage_coverage_regular_90 = table(res_regular_90)["TRUE"]/rep_amount
percentage_coverage_log_90 = table(res_log_90)["TRUE"]/rep_amount
percentage_coverage_sqrt_90 = table(res_sqrt_90)["TRUE"]/rep_amount

percentage_coverage_regular_95 = table(res_regular_95)["TRUE"]/rep_amount
percentage_coverage_log_95 = table(res_log_95)["TRUE"]/rep_amount
percentage_coverage_sqrt_95 = table(res_sqrt_95)["TRUE"]/rep_amount

percentage_coverage_regular_99 = table(res_regular_99)["TRUE"]/rep_amount
percentage_coverage_log_99 = table(res_log_99)["TRUE"]/rep_amount
percentage_coverage_sqrt_99 = table(res_sqrt_99)["TRUE"]/rep_amount

######################################
## Results: Table of coverage rate in %
######################################

coverage_regular = cbind(percentage_coverage_regular_90,
                         percentage_coverage_regular_95,
                         percentage_coverage_regular_99)

coverage_log = cbind(percentage_coverage_log_90,
                     percentage_coverage_log_95,
                     percentage_coverage_log_99)

coverage_sqrt = cbind(percentage_coverage_sqrt_90,
                      percentage_coverage_sqrt_95,
                      percentage_coverage_sqrt_99)

coverage_table = rbind(coverage_regular,
                       coverage_log,
                       coverage_sqrt)

coverage_table = as.data.frame(coverage_table)
rownames(coverage_table) <- c("Regular O:E Ratio", "Ln O:E Ratio", "Sqrt O:E Ratio")
names(coverage_table)[1] = "CI 0.90"
names(coverage_table)[2] = "CI 0.95"
names(coverage_table)[3] = "CI 0.99"

print(coverage_table)

######################################
## Results: Monte Carlo error of coverage rate
######################################
se_reg_90 = (sqrt((coverage_table[1,1]*(1-coverage_table[1,1])/rep_amount)))
se_reg_95 = (sqrt((coverage_table[1,2]*(1-coverage_table[1,2])/rep_amount)))
se_reg_99 = (sqrt((coverage_table[1,3]*(1-coverage_table[1,3])/rep_amount)))
  
se_log_90 = (sqrt((coverage_table[2,1]*(1-coverage_table[2,1])/rep_amount)))
se_log_95 = (sqrt((coverage_table[2,2]*(1-coverage_table[2,2])/rep_amount)))
se_log_99 = (sqrt((coverage_table[2,3]*(1-coverage_table[2,3])/rep_amount)))
  
se_sqrt_90 = (sqrt((coverage_table[3,1]*(1-coverage_table[3,1])/rep_amount)))
se_sqrt_95 = (sqrt((coverage_table[3,2]*(1-coverage_table[3,2])/rep_amount)))
se_sqrt_99 = (sqrt((coverage_table[3,3]*(1-coverage_table[3,3])/rep_amount)))

################################
# Results: Table of Bias
################################
diff_reg_true = vector()
diff_log_true = vector()
diff_sqrt_true = vector()

for (j in 1:rep_amount){
diff_reg_true[[j]] = all_regular_theta[[j]] - mean_true_oe
diff_log_true[[j]] = all_log_theta[[j]] - mean_true_oe
diff_sqrt_true[[j]] = all_sqrt_theta[[j]] - mean_true_oe
}

bias_regular = sum(diff_reg_true)/rep_amount
bias_log = sum(diff_log_true)/rep_amount
bias_sqrt = sum(diff_sqrt_true)/rep_amount

######################################
## Results: Monte Carlo error of bias
######################################
diff_reg_estimate = vector()
diff_log_estimate = vector()
diff_sqrt_estimate = vector()

for (j in 1:rep_amount){
  diff_reg_estimate[[j]] = (all_regular_theta[[j]] - mean(unlist(all_regular_theta)))^2
  diff_log_estimate[[j]] = (all_log_theta[[j]] - mean(unlist(all_regular_theta)))^2
  diff_sqrt_estimate[[j]] = (all_sqrt_theta[[j]] - mean(unlist(all_regular_theta)))^2
}

se_bias_reg = sqrt((1/(rep_amount*(rep_amount-1)))*sum(diff_reg_estimate))
se_bias_log = sqrt((1/(rep_amount*(rep_amount-1)))*sum(diff_log_estimate))
se_bias_sqrt = sqrt((1/(rep_amount*(rep_amount-1)))*sum(diff_sqrt_estimate))

################################
# Visualizations: Histograms of pooled estimates
################################
df_hist <- data.frame(theta_reg = unlist(all_regular_theta),
                      theta_ln = unlist(all_log_theta),
                      theta_sqrt = unlist(all_sqrt_theta))

plot_1 <- ggplot(df_hist, aes(x = theta_reg)) + 
  geom_histogram(aes(y = ..density..), bins=50, color="azure4", fill="white") +
  geom_vline(aes(xintercept=mean(theta_reg)), 
             color="darkred", linetype="solid", size=0.5)+
  geom_vline(aes(xintercept=mean_true_oe), 
             color="darkblue", linetype="solid", size=0.5)+
  theme_classic()+
  ylim(0, 15)+
  labs(x = "Regular O/E",
       y = "O/E ≈ 0.5")+
  theme(
    axis.title.x = element_text(color="black", size=12),
  )+
  geom_density(color="black")+
  theme(text=element_text(family="Times New Roman"))

plot_2 <- ggplot(df_hist, aes(x = theta_ln)) + 
  geom_histogram(aes(y = ..density..), bins=50, color="azure4", fill="white") +
  geom_vline(aes(xintercept=mean(theta_ln)),
             color="darkred", linetype="solid", size=0.5)+
  geom_vline(aes(xintercept=mean_true_oe),
             color="darkblue", linetype="solid", size=0.5)+
  theme_classic()+
  ylim(0, 15)+
  labs(x = "Ln O/E",
       y = "")+
  theme(
    axis.title.x = element_text(color="black", size=12),
  )+
  geom_density(color="black")+
  theme(text=element_text(family="Times New Roman"))

plot_3 <- ggplot(df_hist, aes(x = theta_sqrt)) + 
  geom_histogram(aes(y = ..density..), bins=50, color="azure4", fill="white") +
  geom_vline(aes(xintercept=mean(theta_sqrt)),
             color="darkred", linetype="solid", size=0.5)+
  geom_vline(aes(xintercept=mean_true_oe),
             color="darkblue", linetype="solid", size=0.5)+
  theme_classic()+
  ylim(0, 15)+
  labs(x = "Square root O/E",
       y = "")+
  theme(
    axis.title.x = element_text(color="black", size=12),
  )+
  geom_density(color = "black")+
  theme(text=element_text(family="Times New Roman"))

# size: 900*300
grid.arrange(plot_1, plot_2, plot_3,
             top = "",
             ncol = 3,
             nrow = 1)

################################
# Visualizations: Scatterplot on difference of theta to mean and to true oe
################################
df_diff = data.frame(diff_reg_true, 
                     diff_log_true, 
                     diff_sqrt_true)

df_diff$index = 1:nrow(df_diff)

sc_1 <- ggplot(df_diff, aes(x=diff_reg_true, y=index)) + 
  geom_point(color = "azure4", size=0.1, shape=20, alpha = 0.5)+
  labs(x = "Regular O/E",
       y = "O/E ≈ 0.5")+
  theme_classic()+
  geom_vline(aes(xintercept=mean(diff_reg_true)), # bias
             color="darkred", linetype="solid", size=0.5)+
  geom_vline(aes(xintercept=mean(diff_reg_true)-1.96*se_bias_reg), # monte carlo se of bias: bias-1.96*monte carlo se
             color="darkred", linetype="dashed", size=0.5)+
  geom_vline(aes(xintercept=mean(diff_reg_true)+1.96*se_bias_reg), # monte carlo se of bias: bias+1.96*monte carlo se
             color="darkred", linetype="dashed", size=0.5)+
  geom_vline(aes(xintercept=0), # true
             color="darkblue", linetype="solid", size=0.5)+
  theme(
    axis.title.x = element_text(color="black", size=12),
  )+
  theme(text=element_text(family="Times New Roman"))

sc_2 <- ggplot(df_diff, aes(x=diff_log_true, y=index)) + 
  geom_point(color = "azure4", size=0.1, shape=20, alpha = 0.5)+
  labs(x = "Ln O/E",
       y = "O/E ≈ 0.5")+
  theme_classic()+
  geom_vline(aes(xintercept=mean(diff_log_true)), # bias
             color="darkred", linetype="solid", size=0.5)+
  geom_vline(aes(xintercept=mean(diff_log_true)-1.96*se_bias_reg), # monte carlo se of bias: bias-1.96*monte carlo se
             color="darkred", linetype="dashed", size=0.5)+
  geom_vline(aes(xintercept=mean(diff_log_true)+1.96*se_bias_reg), # monte carlo se of bias: bias+1.96*monte carlo se
             color="darkred", linetype="dashed", size=0.5)+
  geom_vline(aes(xintercept=0), # true
             color="darkblue", linetype="solid", size=0.5)+
  theme(
    axis.title.x = element_text(color="black", size=12),
  )+
  theme(text=element_text(family="Times New Roman"))

sc_3 <- ggplot(df_diff, aes(x=diff_sqrt_true, y=index)) + 
  geom_point(color = "azure4", size=0.05, shape=20, alpha = 0.5)+
  labs(x = "Square root O/E",
       y = "O/E ≈ 0.5")+
  theme_classic()+
  geom_vline(aes(xintercept=mean(diff_sqrt_true)), # bias
             color="darkred", linetype="solid", size=0.5)+
  geom_vline(aes(xintercept=mean(diff_sqrt_true)-1.96*se_bias_reg), # monte carlo se of bias: bias-1.96*monte carlo se
             color="darkred", linetype="dashed", size=0.5)+
  geom_vline(aes(xintercept=mean(diff_sqrt_true)+1.96*se_bias_reg), # monte carlo se of bias: bias+1.96*monte carlo se
             color="darkred", linetype="dashed", size=0.5)+
  geom_vline(aes(xintercept=0), # true
             color="darkblue", linetype="solid", size=0.5)+
  theme(
    axis.title.x = element_text(color="black", size=12),
  )+
  theme(text=element_text(family="Times New Roman"))

# size: 600*350
grid.arrange(sc_1,
             sc_2,
             sc_3,
             top = "",
             ncol = 1,
             nrow = 3)

################################
# Visualizations: Lollipop plot of coverage rate and its monte carlo standard error
################################
# Note: Darkgreen == Regular, Darkorange == Log, Darkblue == Square root

# df Coverage at alpha = 0.1
df_cov_90 = data.frame(
  transformation = c("3", "2", "1"),
  cov_90 = c(unlist(coverage_table[,1])),
  se_90 = c(se_reg_90, se_log_90, se_sqrt_90))
                
# plot Coverage at alpha = 0.1
lolli_1 <- ggplot(df_cov_90) +
  geom_point(aes(x=cov_90, y=transformation, colour=transformation), size = 1) +
  geom_vline(xintercept = 0.9, color = "black")+
  xlim(0.7, 1)+
  scale_color_manual(name = "",
                     values = c("1" = "darkblue", "2" = "darkorange", "3" = "darkgreen"))+
  geom_segment( aes(x=cov_90, xend=0.9, y=transformation, yend=transformation, colour=transformation) , size=0.5) +
  geom_point( aes(x=cov_90-1.96*se_90, y=transformation, colour=transformation), shape=91, size = 5) +
  geom_point( aes(x=cov_90+1.96*se_90, y=transformation, colour=transformation), shape=93, size = 5) +
  theme_set(theme_bw())+
  labs(x = "α = 0.1", y = "O/E ≈ 0.5")+
  theme(legend.position="none", axis.text.y = element_blank())+
  theme(text=element_text(family="Times New Roman"))
  

# df Coverage at alpha = 0.05
df_cov_95 = data.frame(
  transformation = c("3", "2", "1"),
  cov_95 = c(unlist(coverage_table[,2])),
  se_95 = c(se_reg_95, se_log_95, se_sqrt_95))

# plot Coverage at alpha = 0.05
lolli_2 <- ggplot(df_cov_95) +
  geom_point(aes(x=cov_95, y=transformation, colour=transformation), size = 1) +
  geom_vline(xintercept = 0.95, color = "black")+
  xlim(0.7, 1)+
  scale_color_manual(name = "",
                     values = c("1" = "darkblue", "2" = "darkorange", "3" = "darkgreen"))+
  geom_segment( aes(x=cov_95, xend=0.95, y=transformation, yend=transformation, colour=transformation) , size=0.5) +
  geom_point( aes(x=cov_95-1.96*se_95, y=transformation, colour=transformation), shape=91, size = 5) +
  geom_point( aes(x=cov_95+1.96*se_95, y=transformation, colour=transformation), shape=93, size = 5) +
  theme_set(theme_bw())+
  theme(legend.position="none", axis.text.y = element_blank())+
  labs(x = "α = 0.05", y = "")+
  theme(text=element_text(family="Times New Roman"))

# df Coverage at alpha = 0.01
df_cov_99 = data.frame(
  transformation = c("3", "2", "1"),
  cov_99 = c(unlist(coverage_table[,3])),
  se_99 = c(se_reg_99, se_log_99, se_sqrt_99))

# plot Coverage at alpha = 0.05
lolli_3 <- ggplot(df_cov_99) +
  geom_point(aes(x=cov_99, y=transformation, colour=transformation), size = 1) +
  geom_vline(xintercept = 0.99, color = "black")+
  xlim(0.7, 1)+
  scale_color_manual(name = "",
                     values = c("1" = "darkblue", "2" = "darkorange", "3" = "darkgreen"))+
  geom_segment( aes(x=cov_99, xend=0.99, y=transformation, yend=transformation, colour=transformation) , size=0.5) +
  geom_point( aes(x=cov_99-1.96*se_99, y=transformation, colour=transformation), shape=91, size=5) +
  geom_point( aes(x=cov_99+1.96*se_99, y=transformation, colour=transformation), shape=93, size=5) +
  theme_set(theme_bw())+
  theme(legend.position="none", axis.text.y = element_blank())+
  labs(x = "α = 0.01", y = "")+
  theme(text=element_text(family="Times New Roman"))

# size: 800*200   
grid.arrange(lolli_1, lolli_2, lolli_3,
             top = "",
             ncol = 3,
             nrow = 1)
