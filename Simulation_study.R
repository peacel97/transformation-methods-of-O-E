####### V0 Simulation study

## Next steps:
#	Find suitable data generation mechanism
#	Check best imputation methods for mice
# Calculate o:e ratio in simulation study
# When to perform o:e ratio

######################################
### Install and load required packages
######################################

## Specify required packages
packages = c("mice", 
             "lattice", 
             "metamisc",
             "dplyr",
             "ggplot2"
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
# Note: vary sample between n = 500 and n = 1000 later on

# discrete variable could be e.g., age 
discrete_var = rpois(n = 500, lambda = 45)
summary(discrete_var)
hist(discrete_var)

# continuous variables could be e.g., blood pressure
continuous_var_1 = runif(n = 500, min = 60, max = 100) 
summary(continuous_var_1)
hist(continuous_var_1)

continuous_var_2 = runif(n = 500, min = 0, max = 1)
summary(continuous_var_2)
hist(continuous_var_2)

# binary outcome variable with imbalance (i.e., more 0 than 1 as defined in prob)
outcome_var = rbinom(n = 500, size = 1, prob = 0.2)
summary(outcome_var)
hist(outcome_var)

# Combine it to a data frame
df_complete = data.frame( 
           discrete_var,
           continuous_var_1,
           continuous_var_2,
           outcome_var
           )

# Create multiple logistic regression model based on complete data
model_complete = glm(outcome_var ~ discrete_var + continuous_var_1 + continuous_var_2, 
              data = df_complete, 
              family = binomial
              )
summary(model_complete)$coef

######################################
## Ampute data
######################################

# Choose percentage of missing cells (define reasonably between 0.1 and 0.2)
myprop = c(0.15)

# Choose missingness pattern
mypattern = c(1, 1, 1, 1)

# Choose relative occurence of these patterns
myfreq = c(0.25, 0.25, 0.25, 0.25)

# Choose missingness mechanism (MAR also default in ampute)
mymech = c("MAR")

# Choose weights of weighted sum scores
#myweights = 

# Carry out amputation replacing it by NA
amputation = ampute(data = df_complete, 
                       prop = myprop,
                       #freq = myfreq,
                       #mech = mymech
                       #weights = myweights
                       )

# Show that amputation has worked
head(amputation$amp)
summary(amputation)

# Visualization of missing data pattern
md.pattern(amputation$amp)

# Evaluation of amputation: 
bwplot(amputation, which.pat = c(1, 3), descriptives = TRUE)

# Evaluation of amputation:
xyplot(amputation, which.pat = 1)

# save new data including amputations (NA) to data frame
df_amputed = amputation$amp

# Create multiple logistic regression model based on amputed data
model_amputed = glm(outcome_var ~ discrete_var + continuous_var_1 + continuous_var_2, 
                    data = df_amputed, 
                    family = binomial
)
summary(model_amputed)$coef

######################################
## Impute data
######################################

# Choose number of imputed data sets (reasonably between 5 and 10)
imp_amount = c(5)

# Choose imputation method
# Note: Look for better imputation method (listed via methods(mice))
imp_method = c("pmm")

# Impute data via mice
imputation = mice(data = df_amputed, 
              maxit = imp_amount, 
              method = imp_method, 
              print = TRUE
              )

# Inspect the imputed data sets
# Note: Save to dataframe for each iteration
#for(i in 1:imp_amount) {
#  complete(imputation, i)
#}  

df_imp_1 = complete(imputation,1)
df_imp_2 = complete(imputation,2)
df_imp_3 = complete(imputation,3)
df_imp_4 = complete(imputation,4)
df_imp_5 = complete(imputation,5)

######################################
## View performance measures of imputed data sets
######################################

######################################
## Calculate O:E ratio and its transformations
######################################

## Calculate total O:E ratio & its standard error
oe <- oecalc(O = df_imp_1$outcome_var, E = df_complete$outcome_var, N = 1000, data = df_imp_1)

## Display result of O:E ratio
oe

## Calculate log(O:E) ratio & its standard error
log_oe <- oecalc(O = df_imp_1$outcome_var, E = df_complete$outcome_var, N = 1000, data = df_imp_1, g = "log(OE)")

## Display result of log(O:E) ratio
log_oe

######################################
## Pooling procedure
######################################

######################################
## Visualize results
######################################
