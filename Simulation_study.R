####### V0 Simulation study

######################################
### TO DO
######################################

# Threshold predicted class, normalize?
# Change one continuous variable to binary?
# External validation
# Shrinkage
# Pseudo code for chapter 3
# loop repetition

######################################
### Open questions
######################################

# Proportion for split
# Max iteration in imputation???
# Interpretation difference complete data and amputed data too little?
# Shrinkage?


######################################
### Steps
######################################
# Generate simulation data
# Calculate O:E ratio as reference
# Ampute data
# Split data in train and test data set
# Multiple imputation on train data set
# Calculate O:E ratio for each imputed data set
# Pool
# Check for each round of simulation if ratio falls into true ratio
# Externally validate

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
## Generate data for simulation
######################################

## Set seed to get reproducible results
set.seed(18)

## Define variables, distribution of data and sample size

# continuous variables could be e.g., blood pressure
continuous_var_1 = rnorm(n = 1000, mean = 2.0, sd = 0.5)
summary(continuous_var_1)
hist(continuous_var_1)

continuous_var_2 = rnorm(n = 1000, mean = 1.5, sd = 0.25)
summary(continuous_var_2)
hist(continuous_var_2)

# Create linear combination 
# with or without bias as intercept?
lin_combination = -4.5 + continuous_var_1 + continuous_var_2

# Probability for response variable to be 1
# Note: Due to application for logistic regression, use inverse logit function
prob_observed_outcome = 1/(1+exp(-lin_combination))

# Check that values are not approaching either 0 or 1 to avoid too deterministic approach
summary(prob_observed_outcome)

# Binary outcome variable as Bernoulli response variable
# e.g., diabetes positive (1) or negative (0)
# Desired probability for outcome_var = 1 between 20% and 30% to consider imbalance
outcome_var = rbinom(n = 1000, size = 1, prob = prob_observed_outcome)
summary(outcome_var)
hist(outcome_var)

# Combine it to a data frame
df_complete = data.frame( 
           continuous_var_1,
           continuous_var_2,
           outcome_var
           )

######################################
## Calculate reference O:E ratio based on logistic regression model with complete data
######################################

# Split data into test and training data
split_prob = c(0.7)
training_samples <- df_complete$outcome_var %>% 
  createDataPartition(p = split_prob, list = FALSE)

train_data = df_complete[training_samples, ]
test_data = df_complete[-training_samples, ]

# Create multiple logistic regression model based on train data of df_complete
model_complete = glm(outcome_var ~ continuous_var_1 + continuous_var_2, 
              data = train_data, 
              family = binomial
              )
summary(model_complete)$coef

# Predicted outcomes to take as values for expected outcomes in O:E ratio
prob_expected_outcome <- model_complete %>% predict(test_data, type = "response")
summary(model_complete$fitted.values)

# !!!! Predicted classes are not correct. Threshold? Normalize?
predicted_classes = ifelse(prob_expected_outcome > 0.5, 1, 0)
predicted_classes

# Value for O:E ratio by comparison of probabilities
oe_prob_based = mean(prob_observed_outcome) / mean(prob_expected_outcome)
oe_prob_based

# Calculate total O:E ratio & its standard error based on complete data
oe_complete <- oecalc(
             O = sum(df_complete$outcome_var == 1), # numeric vector of observed events
             E = sum(predicted_classes == 1)/(1-split_prob), # numeric vector of expected events
             N = 1000
             )

## Display result of total O:E ratio & its standard error based on complete data and logistic model
oe_complete

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
amputation = ampute(data = df_complete, 
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
split_prob_amputed = c(0.7)
training_samples_amputed <- df_amputed$outcome_var %>% 
  createDataPartition(p = split_prob_amputed, list = FALSE)

df_train_amputed = df_amputed[training_samples_amputed, ]
df_test_amputed = df_amputed[-training_samples_amputed, ]

######################################
## Perform multiple imputation on the training data
######################################

# Choose number of imputed data sets (reasonably between 5 and 10)
imp_amount = c(10)

# Choose number of iterations
imp_iteration = c(1)

# Outcome variable as categorical variable for imp_method logreg
df_train_amputed$outcome_var <- as.factor(df_train_amputed$outcome_var)

# Choose imputation method
# Interpretation: pmm default method for numeric data; logreg method for binary data
imp_method = c("pmm", "pmm", "logreg") 

# Impute data via mice
imputation = mice(data = df_train_amputed, 
              m = imp_amount, 
              maxit = imp_iteration,
              method = imp_method, 
              print = TRUE,
              seed = 18
              )

# Inspect and save the imputed data sets
#for (i in 1:imp_amount){
#  imp[i] = complete(imputation, i)
#  fit_imp[i] = glm(outcome_var ~ continuous_var_1 + continuous_var_2, 
#                   data = imp[i], 
#                   family = binomial)
#}

imp_1 = complete(imputation,1)
imp_2 = complete(imputation,2)
imp_3 = complete(imputation,3)
imp_4 = complete(imputation,4)
imp_5 = complete(imputation,5)
imp_6 = complete(imputation,6)
imp_7 = complete(imputation,7)
imp_8 = complete(imputation,8)
imp_9 = complete(imputation,9)
imp_10 = complete(imputation,10)

######################################
## Fit the prediction model for each imputed data set
######################################

# Fit the model for each imputed data set (here imp_1)
model_imputed_1 = glm(outcome_var ~ continuous_var_1 + continuous_var_2, 
                    data = imp_1, 
                    family = binomial
)

model_imputed_2 = glm(outcome_var ~ continuous_var_1 + continuous_var_2, 
                      data = imp_2, 
                      family = binomial
)

# Predict on test data (here imp_1)
fit_model_imp_1 <- model_imputed_1 %>% predict(df_test_amputed, type = "response")
summary(model_imputed_1$fitted.values)

fit_model_imp_2 <- model_imputed_2 %>% predict(df_test_amputed, type = "response")
summary(model_imputed_2$fitted.values)

prob_expected_outcome_imp_1 = mean(model_imputed_1$fitted.values)
prob_expected_outcome_imp_1

prob_expected_outcome_imp_2 = mean(model_imputed_2$fitted.values)
prob_expected_outcome_imp_2

######################################
## Calculate O:E, log(O:E) and square_root(O:E) for each imputed data set
######################################

# Probability for observed outcome for each imputed data set
prob_observed_outcome_imp_1 = sum(df_imp_1$outcome_var == 1) / nrow(df_imp_1)
prob_observed_outcome_imp_1

# Calculate regular O:E ratio for each imputed data set (here imp_1)
oe_prob_based_imp_1 = mean(prob_observed_outcome_imp_1) / prob_expected_outcome_imp_1
oe_prob_based_imp_1

# Probability for observed outcome for each imputed data set
prob_observed_outcome_imp_2 = sum(df_imp_2$outcome_var == 1) / nrow(df_imp_2)
prob_observed_outcome_imp_2

# Calculate regular O:E ratio for each imputed data set (here imp_1)
oe_prob_based_imp_2 = mean(prob_observed_outcome_imp_2) / prob_expected_outcome_imp_2
oe_prob_based_imp_2

# Calculate log(O:E) ratio for each imputed data set
log(oe_prob_based_imp_1)


######################################
## Pooling procedure
######################################

# Make object type compatible for pooling


# Pool the data sets to one
pool.syn(model_imputed_1,  rule = "reiter2003")


# Back-transformation total O:E ratio
back_oe = exp()















######################################
## Transformation of O:E and Pooling
######################################


#
#pooled <- pool.syn(fit, dfcom = NULL, rule = "reiter2003")
#
#summary(pooled)

# Calculate total O:E ratio & its standard error based on imputed model
oe_imputed <- oecalc(
  O = , # numeric vector of observed events
  E = , # numeric vector of expected events
  N = 1000
)

## Display result of O:E ratio
oe_imputed
plot(eo)

## Calculate log(O:E) ratio & its standard error based on imputed data
log_oe <- oecalc(
  O = , # numeric vector of observed events
  E = , # numeric vector of expected events
  N = 1000,
  g = "log(OE)"
)

## Display result of log(O:E) ratio
log_oe
plot(log_oe)

# Calculate square root transformation of O:E ratio based on imputed data

# Display results of square root transformation of O:E ratio

# Pooling procedure

# Back-transformation to the original scale

######################################
## Compare performance of different O:E measures
######################################



######################################
## Externally validate by new data set 
######################################

######################################
## Robustness check with complete case analysis
######################################

# Testing: As data are missing MAR in the simulation study, they should be less biased compared to complete case analysis

