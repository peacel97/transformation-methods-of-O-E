####### V1 Simulation study

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

# Install&load packages or load packages, if already installed
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
## Generate data for simulation: Define variables, distribution of data and sample size
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
summary(model_complete$fitted.values)
prob_fitted = as.data.frame(model_complete$fitted.values)

# Define observed and expected events based on complete data
n_observed_event = sum(df_complete$outcome_var == 1)
n_expected_event = round(mean(model_complete$fitted.values)*n_sample)

# O:E ratio with expected event calculated based on predicted probabilities
oe_complete_2 = oecalc(
  O = n_observed_event, # numeric vector of observed events
  E = n_expected_event, # numeric vector of expected events
  N = 1000
)
oe_complete_2

# # Predicted classes
# predicted_classes = ifelse(prob_fitted > 0.5, 1, 0)
# predicted_classes
# 
# # Rescale probability to [0;1] for predicting classes with threshold 0.5
# prob_fitted_scales <- data.frame(matrix(nrow = nrow(prob_fitted), ncol = ncol(prob_fitted)))
# for(i in seq_len(nrow(prob_fitted))){
#   column <- ((prob_fitted[i,1] - min(prob_fitted))/ (max(prob_fitted) - min(prob_fitted)))
#   prob_fitted_scales[i,1] <- column
# }
# names(prob_fitted_scales)[1] <- "classes"
# summary(prob_fitted_scales)
# 
# # Predicted classes
# predicted_classes = ifelse(prob_fitted_scales$classes > 0.5, 1, 0) 
# sum_predicted_classes = sum(predicted_classes == 1)
# sum_predicted_classes



# # O:E ratio with expected event calculated based on predicted classes
# n_expected_event_1 = round(sum_predicted_classes/(split_prob))
# oe_complete_1 = oecalc(
#              O = n_observed_event, # numeric vector of observed events
#              E = n_expected_event_1, # numeric vector of expected events
#              N = 1000
#              )
# oe_complete_1

######################################
## Ampute data from complete data set using MAR mechanism
######################################

# Define parameters for amputation
myprop = c(0.15) # Choose percentage of missing cells (define reasonably between 0.1 and 0.2)
mypattern = rbind(c(0, 0, 1), c(1, 0 ,1), c(0, 1, 1)) # Choose missingness pattern (0 == missing, 1 == non-missing)
myfreq = c(0.1, 0.45, 0.45) # Choose relative occurence of these patterns
mymech = c("MAR") # Choose missingness mechanism (MAR based on literature)
myweights = ampute.default.weights(mypattern, mymech) # Choose weights of weighted sum scores

# Carry out amputation replacing it by NA
# Note for interpretation: Proportion of missingness in terms of cells not in terms of cases
amputation = ampute(data = train_data_complete, # df_complete, 
                    patterns = mypattern,
                    prop = myprop, 
                    freq = myfreq, 
                    mech = mymech, 
                    weights = myweights, 
                    bycases = FALSE
                    )

# Show that amputation has worked
head(amputation$amp)
summary(amputation)

# Visualization of missing data pattern
md.pattern(amputation$amp)

# save new data including amputations (NA) to data frame
df_amputed = amputation$amp

# Evaluation of amputation
bwplot(amputation, which.pat = c(1, 3), descriptives = TRUE)
xyplot(amputation, which.pat = 1)

# check for convergence. Update if full convergence.
#tbd. Pattern? No? Good.

######################################
## Compare logistic regression model based on complete data and based on amputed data
######################################

# Create multiple logistic regression model based on complete data for comparison
model_complete = glm(outcome_var ~ continuous_var_1 + continuous_var_2, 
                    data = df_complete, 
                    family = binomial
)
summary(model_complete)$coef

# Create multiple logistic regression model based on amputed data for comparison
model_amputed = glm(outcome_var ~ continuous_var_1 + continuous_var_2, 
                    data = df_amputed, 
                    family = binomial
)
summary(model_amputed)$coef

######################################
## Split amputed data into training and test data
######################################

# Split amputed data into training and test data
split_prob_amputed = c(0.5)
training_samples_amputed <- df_amputed$outcome_var %>% 
  createDataPartition(p = split_prob_amputed, list = FALSE)

df_train_amputed = df_amputed[training_samples_amputed, ]
df_test_amputed = df_amputed[-training_samples_amputed, ]

######################################
## Perform multiple imputation
######################################

# Choose number of imputed data sets (reasonably between 5 and 10)
imp_amount = c(10)

# Choose number of iterations
imp_iteration = c(1)

# Choose imputation method
# Interpretation: pmm default method for numeric data; logreg method for binary data
imp_method = c("pmm", "pmm", "logreg") 

# Outcome variable as categorical variable for imp_method logreg
df_train_amputed$outcome_var <- as.factor(df_train_amputed$outcome_var)
df_test_amputed$outcome_var <- as.factor(df_test_amputed$outcome_var)

# Impute data via mice
imputation_train = mice(data = df_train_amputed, 
              m = imp_amount, 
              maxit = imp_iteration,
              method = imp_method, 
              print = TRUE,
              seed = 18
              )

imputation_test = mice(data = df_test_amputed, 
                        m = imp_amount, 
                        maxit = imp_iteration,
                        method = imp_method, 
                        print = TRUE,
                        seed = 18
                       )

# Inspect and save the imputed data sets
# Note: Due to imputation additional variables .imp and .id
imputed_data_train = complete(imputation_train, "long")
imputed_data_test = complete(imputation_test, "long")

######################################
## Fit model for each imputed data set
######################################

#### Here: Test with .imp == 1

# Subset data to relevant variables
imputed_data_train_1_sub = subset(imputed_data_train, .imp == 1)
imputed_data_train_1_sub = dplyr::select(imputed_data_train_1_sub, 
                              continuous_var_1, 
                              continuous_var_2,
                              outcome_var
                              )

# Build model based on imp_1
model_train_1 = glm(outcome_var ~ continuous_var_1 + continuous_var_2,
                    data = imputed_data_train_1_sub,
                    family = binomial 
                    )
summary(model_train_1)$coef

# Predict model based on imp_1
fit_model_imp_1 <- model_train_1 %>% predict(df_test_amputed, type = "response")
summary(model_train_1$fitted.values)

# Define observed and expected events for imp_1
n_observed_event_imp_1 = sum(imputed_data_train_1_sub$outcome_var == 1)
n_expected_event_imp_1 = round(mean(model_complete$fitted.values)*nrow(imputed_data_train_1_sub))

# O:E ratio with expected event calculated based on predicted probabilities
oe_complete_imp_1 = oecalc(
  O = n_observed_event_imp_1, # numeric vector of observed events
  E = n_expected_event_imp_1, # numeric vector of expected events
  N = nrow(imputed_data_train_1_sub)
)
oe_complete_imp_1

# log(O:E) ratio with expected event calculated based on predicted probabilities
log_oe_complete_imp_1 = oecalc(
  O = n_observed_event_imp_1, # numeric vector of observed events
  E = n_expected_event_imp_1, # numeric vector of expected events
  N = nrow(imputed_data_train_1_sub),
  g = "log(OE)"
)
log_oe_complete_imp_1

# square root O:E ratio with expected event calculated based on predicted probabilities
#tbd

######################################
## Pooling procedure
######################################

# Make object type compatible for pooling

# Pool the data sets to one
# pool.syn(imputation,  rule = "reiter2003")

######################################
## Back-transformation of O:E ratios
######################################

# Back-transformation O:E ratio to original scale
back_oe = exp()

# Back-transformation log(O:E) ratio to original scale

# Back-transformation square root (O:E) ratio to original scale

######################################
## Compare performance of different O:E measures by n = ??? repetitions of the simulation study
######################################

######################################
## Externally validate by new data set 
######################################

######################################
## "Robustness check" with complete case analysis
######################################
# Testing: As data are missing MAR in the simulation study, they should be less biased compared to complete case analysis