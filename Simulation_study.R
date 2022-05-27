####### V1 Simulation study

######################################
# Open:
# Split ratio
# Reference O:E exactly 1?
# Convergence
# Reference O:E based on complete data or on split data?
######################################

######################################
# TO DO
# Write pool_se and pool_ci function
# Write function for CI in reference O:E
# Change predictor variable to binary
# Check for convergence. Update if full convergence. tbd. Pattern? No? Good.
# Loop over imputations
# Increase imputation to 10 as soon as loop implemented
# Check nrow of observed and expected events
# For loop simulation study repetition
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
# Repeat 1000 times
# Check for each round of simulation if ratio falls into true ratio

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
set.seed(18)

# Choose sample size 
n_sample = c(2000)

# Continuous predictor variables
continuous_var_1 = rnorm(n = n_sample, mean = 2.0, sd = 0.5)
summary(continuous_var_1)
hist(continuous_var_1)

continuous_var_2 = rnorm(n = n_sample, mean = 1.5, sd = 0.25)
summary(continuous_var_2)
hist(continuous_var_2)

# Create linear combination of predictor variables and intercept
lin_combination = -4.5 + continuous_var_1 + continuous_var_2

# Probability for response variable to be 1
# Note: Due to application for logistic regression, use inverse logit function
prob_outcome = 1/(1+exp(-lin_combination))

# Check that values are not approaching either 0 or 1 to avoid too deterministic approach
summary(prob_outcome)

# Binary outcome variable as Bernoulli response variable
# Note: Imbalanced Classification, desired probability for outcome_var = 1 between 20% and 30%
outcome_var = rbinom(n = n_sample, size = 1, prob = prob_outcome)
summary(outcome_var)
hist(outcome_var)

# Combine it to a data frame
df_complete = data.frame( 
           continuous_var_1,
           continuous_var_2,
           outcome_var
           )

######################################
## Split data and develop prediction model
######################################

## Split data into test and training data
split_prob = c(0.5)
training_samples_complete <- df_complete$outcome_var %>% 
  createDataPartition(p = split_prob, list = FALSE)

train_data_complete = df_complete[training_samples_complete, ]
test_data_complete = df_complete[-training_samples_complete, ]

# Create multiple logistic regression model based on complete data
model_complete = glm(outcome_var ~ continuous_var_1 + continuous_var_2, 
                     data = train_data_complete,
                     family = binomial
)
summary(model_complete)$coef

######################################
## Calculate reference O:E ratio based on logistic regression model with complete data
######################################

# Fit model based on complete data
fit_model_complete <- model_complete %>% predict(test_data_complete,
                                                 type = "response"
                                                 )
# Predicted probability
summary(fit_model_complete)

# Define observed and expected events based on complete data
n_observed_event = sum(df_complete$outcome_var == 1)
n_expected_event = mean(fit_model_complete)*n_sample

# O:E ratio with expected event calculated based on predicted probabilities
reference_oe = oecalc(
  O = n_observed_event, # numeric vector of observed events
  E = n_expected_event, # numeric vector of expected events
  N = 2000
)
reference_oe

######################################
## Ampute test data 
######################################

# Define parameters for amputation
myprop = c(0.15) # Choose percentage of missing cells (define reasonably between 0.1 and 0.2)
mypattern = rbind(c(0, 0, 1), c(1, 0 ,1), c(0, 1, 1)) # Choose missingness pattern (0 == missing, 1 == non-missing)
myfreq = c(0.1, 0.45, 0.45) # Choose relative occurence of these patterns
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

# Compare
summary(model_amputed)$coef
summary(model_complete)$coef

######################################
## Perform multiple imputation
######################################

# Choose number of imputed data sets (reasonably between 5 and 10)
imp_amount = c(3)

# Choose number of iterations
imp_iteration = c(1)

# Choose imputation method
# Interpretation: pmm default method for numeric data; logreg method for binary data
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

# Inspect and save the imputed data sets
imputed_data_1 = complete(imputation, 1)
imputed_data_2 = complete(imputation, 2)
imputed_data_3 = complete(imputation, 3)

######################################
## Fit the model with the imputed data sets
######################################

# Predict model based on imputed_data_1
fit_imp_1 <- model_complete %>% predict(imputed_data_1, 
                                             type = "response"
                                             )
summary(fit_imp_1)

# Predict model based on imputed_data_2
fit_imp_2 <- model_complete %>% predict(imputed_data_2, 
                                          type = "response"
)
summary(fit_imp_2)

# Predict model based on imputed_data_3
fit_imp_3 <- model_complete %>% predict(imputed_data_3, 
                                          type = "response"
)
summary(fit_imp_3)

######################################
## Calculate Regular O:E ratio
######################################

# Calculate regular O:E imp 1
n_observed_event_imp_1 = sum(imputed_data_1$outcome_var == 1)
n_expected_event_imp_1 = mean(fit_imp_1)*n_sample*split_prob

regular_oe_1 = oecalc(
  O = n_observed_event_imp_1,
  E = n_expected_event_imp_1,
  N = 1000
)
regular_oe_1

# Calculate regular O:E for imp 2
n_observed_event_imp_2 = sum(imputed_data_2$outcome_var == 1)
n_expected_event_imp_2 = mean(fit_imp_2)*n_sample*split_prob

regular_oe_2 = oecalc(
  O = n_observed_event_imp_2,
  E = n_expected_event_imp_2,
  N = 1000
)
regular_oe_2

# Calculate regular O:E for imp 3
n_observed_event_imp_3 = sum(imputed_data_3$outcome_var == 1)
n_expected_event_imp_3 = mean(fit_imp_3)*n_sample*split_prob

regular_oe_3 = oecalc(
  O = n_observed_event_imp_3,
  E = n_expected_event_imp_3,
  N = 1000
)
regular_oe_3

######################################
## Calculate Log(O:E) ratio
######################################

# Calculate log(O:E) for imp 1
log_oe_1 = oecalc(
  O = n_observed_event_imp_1,
  E = n_expected_event_imp_1,
  N = 1000,
  g = "log(OE)"
)
log_oe_1

# Calculate log(O:E) for imp 2
log_oe_2 = oecalc(
  O = n_observed_event_imp_2,
  E = n_expected_event_imp_2,
  N = 1000,
  g = "log(OE)"
)
log_oe_2

# Calculate log(O:E) for imp 3
log_oe_3 = oecalc(
  O = n_observed_event_imp_3,
  E = n_expected_event_imp_3,
  N = 1000,
  g = "log(OE)"
)
log_oe_3

######################################
## Calculate Square_root(O:E) ratio
######################################

# Calculate square root (O:E) for imp 1
sqrt_oe_1 = oecalc(
  O = n_observed_event_imp_1,
  E = n_expected_event_imp_1,
  N = 1000,
  g = "sqrt(OE)"
)
sqrt_oe_1

# Calculate square root (O:E) for imp 2
sqrt_oe_2 = oecalc(
  O = n_observed_event_imp_2,
  E = n_expected_event_imp_2,
  N = 1000,
  g = "sqrt(OE)"
)
sqrt_oe_2

# Calculate square root (O:E) for imp 3
sqrt_oe_3 = oecalc(
  O = n_observed_event_imp_3,
  E = n_expected_event_imp_3,
  N = 1000,
  g = "sqrt(OE)"
)
sqrt_oe_3

######################################
## Pooling procedure
######################################

# Bind values together
regular_oe = rbind(regular_oe_1, regular_oe_2, regular_oe_3)
log_oe = rbind(log_oe_1, log_oe_2, log_oe_3)
sqrt_oe = rbind(sqrt_oe_1, sqrt_oe_2, sqrt_oe_3)

# Pool the o:e estimates (theta)
pooled_theta_regular_oe = mean(regular_oe$theta)
pooled_theta_log_oe = mean(log_oe$theta)
pooled_theta_sqrt_oe = mean(sqrt_oe$theta)

# Pool the standard errors of the o:e estimates (theta.se)
# pooled_se_regular_oe = 
# pooled_se_log_oe = 
# pooled_se_sqrt_oe = 

# Pool the confidence intervals of the o:e estimates (theta.cilb and theta.ciub)
# pooled_cilb_regular_oe = 
# pooled_cilb_log_oe = 
# pooled_cilb_sqrt_oe = 
# 
# pooled_ciub_regular_oe =
# pooled_ciub_log_oe = 
# pooled_ciub_sqrt_oe = 

######################################
## Back-transformation of O:E ratios to original scale
######################################

back_theta_regular_oe = pooled_theta_regular_oe
# back_se_regular_oe = pooled_se_regular_oe
# back_cilb_regular_oe = pooled_cilb_regular_oe
# back_ciup_regular_oe = pooled_ciub_regular_oe
# 
# final_regular_oe = cbind(back_theta_regular_oe, 
#                          back_se_regular_oe, 
#                          back_cilb_regular_oe, 
#                          back_ciup_regular_oe
#                          )

back_theta_log_oe = exp(pooled_theta_log_oe)
#back_se_log_oe = exp(pooled_se_log_oe)
#back_cilb_log_oe = exp(pooled_cilb_log_oe)
#back_ciup_log_oe = exp(pooled_ciup_log_oe)
#
# final_log_oe = cbind(back_theta_log_oe, 
#                          back_se_log_oe, 
#                          back_cilb_log_oe, 
#                          back_ciup_log_oe
#                          )

back_theta_sqrt_oe =(pooled_theta_sqrt_oe)^2
#back_se_sqrt_oe = (pooled_se_sqrt_oe)^2
#back_cilb_sqrt_oe = (pooled_cilb_sqrt_oe)^2
#back_ciup_sqrt_oe = (pooled_ciup_sqrt_oe)^2
#
# final_sqrt_oe = cbind(back_theta_sqrt_oe, 
#                          back_se_sqrt_oe, 
#                          back_cilb_sqrt_oe, 
#                          back_ciup_sqrt_oe
#                          )

######################################
## Compare performance of O:E measures with reference O:E measure
######################################
# Check if reference O:E lies in the 95% CI of the O:E measure

######################################
## "Robustness check" with complete case analysis
######################################
# Testing: As data are missing MAR in the simulation study, they should be less biased compared to complete case analysis