####### V0 Simulation study

######################################
### Install and load required packages
######################################

## Specify required packages
packages = c("mice", 
             "lattice", 
             "metamisc"
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

# Create a new dataframe which amputes values replacing it by NA
amputation_df = ampute(data = df_complete, 
                       prop = myprop,
                       bycases = FALSE
                       #freq = myfreq,
                       #mech = mymech
                       #weights = myweights
                       )

df_ammmmp = as.data.frame(amputation_df, row.names = NULL, optional = FALSE)

# Show that amputation has worked
head(amputation_df$amp)
summary(amputation_df)

# Visualization of missing data pattern
md.pattern(amputation_df$amp)

# Create multiple logistic regression model based on amputed data
model_amputed = glm(outcome_var ~ discrete_var + continuous_var_1 + continuous_var_2, 
                    data = amputed_df, 
                    family = binomial
                    )
summary(model_amputed)$coef

######################################
## Impute data
######################################

# Choose number of imputed data sets (reasonably between 5 and 10)
imp_amoung = c(5)

# Choose imputation method
# Note: Look for better imputation method (listed via methods(mice))
imp_method = c("pmm")

# Impute data via mice
imp_df = mice(amputation_df, 
              m = imp_amount, 
              method = imp_method
              )

######################################
## Perform transformation on O:E ratio
######################################

## Define parameters
#n = 
#n.events = 
#e.events = 
#df = 
#
## Calculate total O:E ratio & its standard error
#oe <- oecalc(O = n.events, E = e.events, N = n, data = df)
#
## Display result of O:E ratio
#oe
#
## Display total result in forest plot
#plot(oe)
#
## Calculate log O:E ratio & its standard error
#log_oe <- oecalc(O = n.events, E = e.events, N = n, data = df, g = "log(OE)")
#log_oe
#
## Display result in forest plot
#plot(log_oe)

######################################
## View performance measures
######################################

######################################
## Visualize results
######################################
