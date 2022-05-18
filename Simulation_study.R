####### V0 Simulation study

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
             "MASS"
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
## Generate data
######################################

## Set seed to get reproducible results
set.seed(18)

## Define variables, distribution of data and sample size

# continuous variables could be e.g., blood pressure
continuous_var_1 = rnorm(n = 1000, mean = 80, sd = 5)
summary(continuous_var_1)
hist(continuous_var_1)

continuous_var_2 = rnorm(n = 1000, mean = 50, sd = 2)
summary(continuous_var_2)
hist(continuous_var_2)

# Create linear combination 
# with or without bias as intercept?
lin_combination = continuous_var_1 + continuous_var_2

# Probability for response variable to be 1
# Note: Due to application for logistic regression, use invlogit function
#prob_observed_outcome = 1/(exp(-lin_combination))
prob_observed_outcome = 1/(3+exp(-lin_combination))
prob_observed_outcome

# Add noise
# tbd

# Binary outcome variable as Bernoulli response variable
# e.g., diabetes positive (1) or negative (0)
# Desired probability for outcome_var = 1 between 20% and 30% to consider imbalance
#outcome_var = rbinom(n = 1000, size = 1, prob = prob_observed_outcome)
outcome_var = rbinom(n = 1000, size = 1, prob = 0.2)
summary(outcome_var)
hist(outcome_var)

# Combine it to a data frame
df_complete = data.frame( 
           continuous_var_1,
           continuous_var_2,
           outcome_var
           )

######################################
## Build logistic regression, predict outcomes & calculate original O:E ratio
######################################

# Split data into test and training data
training_samples <- df_complete$outcome_var %>% 
  createDataPartition(p = 0.7, list = FALSE)

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

# normalize
predicted_classes = ifelse(prob_expected_outcome > 0.5, 1, 0)
predicted_classes

# Assess model accuracy
mean(predicted_classes == test_data$outcome_var)

# Value for O:E ratio by comparison of probabilities
#prob_based_oe = mean(prob_observed_outcome) / mean(prob_expected_outcome)
prob_based_oe

# Calculate total O:E ratio & its standard error based on complete data
oe_complete <- oecalc(
             O = , # numeric vector of observed events
             E = , # numeric vector of expected events
             N = 1000
             )

## Display result of total O:E ratio & its standard error based on complete data and logistic model
oe_complete
plot(eo_complete)

######################################
## Ampute data from complete data set
######################################

# Define parameters for amputation
myprop = c(0.15) # Choose percentage of missing cells (define reasonably between 0.1 and 0.2)
mypattern = rbind(c(0, 0, 1), c(1, 0 ,1), c(0, 1, 1)) # Choose missingness pattern (0==missing, 1==non-missing)
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
#tbd

# Create multiple logistic regression model based on amputed data for comparison
model_amputed = glm(outcome_var ~ continuous_var_1 + continuous_var_2, 
                    data = df_amputed, 
                    family = binomial
)
summary(model_amputed)$coef

######################################
## Impute data
######################################

# Choose number of imputed data sets (reasonably between 5 and 10)
imp_amount = c(5)

# Outcome variable as categorical variable for imp_method logreg
df_amputed$outcome_var <- as.factor(df_amputed$outcome_var)

# Choose imputation method
# pmm default method for numeric data; logreg method for binary data
imp_method = c("pmm", "pmm", "logreg") 

# Impute data via mice
set.seed(18)
imputation = mice(data = df_amputed, 
              maxit = imp_amount, 
              method = imp_method, 
              print = TRUE
              )

# Inspect the imputed data sets
df_imp_1 = complete(imputation,1)
df_imp_2 = complete(imputation,2)
df_imp_3 = complete(imputation,3)
df_imp_4 = complete(imputation,4)
df_imp_5 = complete(imputation,5)

######################################
## Calculate O:E ratio and its transformations and visualize results based on imputed data
## Note: Implement correct O:E ratio, below just testing function
######################################

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

######################################
## Compare performance of different O:E measures
######################################

######################################
## Pooling procedure
######################################

######################################
## Evaluate performance of study
######################################

######################################
## Compare results to reference method
######################################



