install.packages("readxl")
install.packages("tidyverse")
install.packages("plm")
install.packages("gplots")
install.packages("foreign")
install.packages("imputeTS")
install.packages("panelView")
install.packages("zoo")
install.packages("dplyr")
install.packages("nlme")
install.packages("lmtest")

library(readxl)
library(plm)
library(tidyverse)
library(gplots)
library(foreign)
library(panelView)
library(dplyr)
library(nlme)
library(lmtest)
library(tseries)


# data parametrage ________________________________________________________
# Load Data
data <- read_excel("data-science/lending_borrowing_data.xlsx")
head(data)

# Create log-transformed variables
# Compute log-transformed variables and log return of TVL
datanew <- data %>%
  arrange(Name, Date) %>%
  group_by(Name) %>%
  mutate(
    log_fees = log(Fees + 1),
    log_tokinc = log(token_incentives + 1), 
    log_coredev = log(core_dev + 1), 
    log_treasury = log(treasury + 1),
    log_tvl = log(TVL + 1),
    log_return_tvl = (log_tvl - lag(log_tvl)) / lag(log_tvl)  # Corrected relative return calculation
  ) %>%
  ungroup()

# Filter rows with NA or invalid values in log_return_tvl (this is optional in case of error)
datanew <- datanew %>% filter(!is.na(log_return_tvl) & log_return_tvl > 0)

# Convert data to panel format
data_pd <- pdata.frame(datanew, index = c("Name", "Date"), drop.index = FALSE)
data_pd <- make.pbalanced(data_pd, balance.type = "shared.times")

# Verify structure of the panel data
head(data_pd)
head(index(data_pd))

# Plot log-transformed dependent variable 
ggplot(datanew, aes(x = Date, y = log_return_tvl, group = Name, color = Name)) +
  geom_line(size = 1.2) +
  labs(title = "Log-Transformed Return TVL Over Time", x = "Time", y = "Log(Return TVL)") +
  theme_minimal() +
  theme(legend.position = "none")

#plot the data to get an overview _______________________________________________
plotmeans(log_return_tvl ~ Name, main = "Average Logarithmic Return of TVL by Name", data = datanew)
plotmeans(log_tvl ~ Date, main = "Average Logarithmic of TVL by Date", data = datanew)
plotmeans(log_return_tvl ~ Date, main = "Average Logarithmic of TVL by Date", data = datanew)
plotmeans(log_tvl ~ Name, main = "Average Logarithmic of TVL by Name", data = datanew)

# Here we are ploting avery independent variable against the dependent variable to check for patern and if high variation in the raw data are observable
datanew %>% filter (log_fees != "0") %>% ggplot(mapping= aes(x=log_fees, y= log_return_tvl, color = Name))+
  geom_point()+theme_bw()+theme(legend.position="bottom")+labs(title="return tvl against fees")

datanew %>% filter (log_treasury != "0") %>% ggplot(mapping= aes(x=log_treasury, y= log_return_tvl, color = Name))+
  geom_point()+theme_bw()+theme(legend.position="bottom")+labs(title="return tvl against treasury")

datanew %>%filter (log_tokinc != "0") %>% ggplot(mapping= aes(x=log_tokinc, y= log_return_tvl, color = Name))+
  geom_point()+theme_bw()+theme(legend.position="bottom")+labs(title="return tvl against token incentives")

datanew %>%filter (log_coredev != "0") %>% ggplot(mapping= aes(x=log_coredev, y= log_return_tvl, color = Name))+
  geom_point()+theme_bw()+theme(legend.position="bottom")+labs(title="return tvl against core dev")


# Perform Levin-Lin-Chu (LLC) stationarity test on log_return_tvl and log_tvl
llc_test_returntvl <- purtest(data_pd$log_return_tvl, test = "levinlin")
llc_test_tvl <- purtest(data_pd$log_tvl, test = "levinlin")

# Additional LLC tests for other variables
llc_test_fees <- purtest(data_pd$log_fees, test = "levinlin")
llc_test_treasury <- purtest(data_pd$log_treasury, test = "levinlin")
llc_test_core_dev <- purtest(data_pd$log_coredev, test = "levinlin")
llc_test_tokinc <- purtest(data_pd$log_tokinc, test = "levinlin")

# Display summary of the results
summary(llc_test_returntvl)
summary(llc_test_fees)
summary(llc_test_treasury)
summary(llc_test_core_dev)
summary(llc_test_tokinc)
summary(llc_test_tvl)

#here we are using adf test in order to double check for stationarity
sum(is.na(data_pd$log_return_tvl))
str(data_pd$log_return_tvl)
colnames(data_pd)
adf_test_tvl <- adf.test(data_pd$log_tvl, alternative = "stationary")
adf_test_return_tvl <- adf.test(data_pd$log_return_tvl, alternative = "stationary")
adf_test_incentives <- adf.test(data_pd$log_tokinc, alternative = "stationary")
adf_test_treasury <- adf.test(data_pd$log_treasury, alternative = "stationary")
adf_test_Fees <- adf.test(data_pd$log_fees, alternative = "stationary")
adf_test_core_dev <- adf.test(data_pd$log_coredev, alternative = "stationary")

# Display summary of the results
print(adf_test_return_tvl)
print(adf_test_tvl)
print(adf_test_incentives)
print(adf_test_treasury)
print(adf_test_Fees)
print(adf_test_core_dev)


#correct stationarity 
# Create differenced variables by grouping with dplyr
data_pd <- na.omit(data_pd) #optional in case of error by processing the code
data_pd <- data_pd %>%
  group_by(Name) %>%
  mutate(
    diff_log_return_tvl = log_return_tvl - lag(log_return_tvl),
    diff_log_tvl = log_tvl - lag(log_tvl),
    diff_log_fees = log_fees - lag(log_fees),
    diff_log_tokinc = log_tokinc - lag(log_tokinc),
    diff_log_coredev = log_coredev - lag(log_coredev),
    diff_log_treasury = log_treasury - lag(log_treasury)
  ) %>%
  ungroup()

# Remove any NAs from the differenced series
diff_log_treasury <- na.omit(diff_log_treasury)
diff_log_return_tvl <- na.omit(diff_log_return_tvl)
diff_log_tvl <- na.omit(diff_log_tvl)
diff_log_coredev <- na.omit(diff_log_coredev)
diff_log_tokinc <- na.omit(diff_log_tokinc)
diff_log_fees <- na.omit(diff_log_fees)


#Run the ADF test again to recheck if stationarity is corrected on all independent and dependent variables
adf_test_diff_treasury <- adf.test(diff_log_treasury, alternative = "stationary")
adf_test_diff_return_tvl <- adf.test(diff_log_return_tvl, alternative = "stationary")
adf_test_diff_tvl <- adf.test(diff_log_tvl, alternative = "stationary")
adf_test_diff_coredev <- adf.test(diff_log_coredev, alternative = "stationary")
adf_test_diff_fees <- adf.test(diff_log_fees, alternative = "stationary")
adf_test_diff_tokinc <- adf.test(diff_log_tokinc, alternative = "stationary")

#display results
print(adf_test_diff_return_tvl)
print(adf_test_diff_tvl)
print(adf_test_diff_treasury)
print(adf_test_diff_fees)
print(adf_test_diff_tokinc)
print(adf_test_diff_coredev)

# Here we are choosing the best model which fit the best our data
# Panel Data Models (OLS, FE, RE)


data_pd <- na.omit(data_pd) # again if there is arror 
# pooled (OLS) model:
ols <- plm(log_return_tvl ~ diff_log_fees + diff_log_tokinc + diff_log_coredev + diff_log_treasury, data = data_pd)
summary(ols)

# Fixed Effects Model (FE):
fe <- plm(log_return_tvl ~ diff_log_fees + diff_log_tokinc + diff_log_coredev + diff_log_treasury, 
          data = data_pd, model = "within", index = c("Name", "Date"))
summary(fe)

# test which model is the best between "OLS" and "FE":
pooltest(log_return_tvl ~ diff_log_fees + diff_log_tokinc + diff_log_coredev + diff_log_treasury, data = data_pd, index= c("Name","Date"), model="within")


# Random Effects Model (RE):
re <- plm(log_return_tvl ~ diff_log_fees + diff_log_tokinc + diff_log_coredev + diff_log_treasury, 
          data = data_pd, index = c("Name", "Date"), model = "random", random.method = "amemiya")
summary(re)

#checking for stationarity again after applying fixed effect (optional as well)
# Extract the residuals from the model
residuals_re <- residuals(re)
adf_test_residuals <- adf.test(residuals_re, alternative = "stationary")
print(adf_test_residuals)

# Hausman Test to compare "FE" and "RE":
phtest(fe, re)

#Check for bias in the model :
pdwtest(re) # Serial correlation
bptest(re, data = data_pd, studentize = FALSE) # Heteroscedasticity
pcdtest(re, test = "lm")  # Lagrange Multiplier test
pcdtest(re, test = "cd")  # Pesaran test

# Compute Driscol l-Kraay standard errors in order to reduce cross sectionall dependencies and heteroscedasticity 
cluster_se <- vcovHC(re, type = "HC1", cluster = "group")

# Display coefficients with cluster standard errors
coeftest(re, vcov = cluster_se)
r_squared <- summary(re)$r.squared
cat("R-squared:", r_squared, "\n")

#information criteria
# Extract necessary information
rss <- sum(resid(re)^2)  # Residual sum of squares
n <- nobs(re)            # Number of observations
k <- length(coef(re))    # Number of parameters in the model

# Compute Log-Likelihood
log_likelihood <- -n / 2 * (log(2 * pi) + log(rss / n) + 1)

# Compute AIC and BIC
aic_value <- -2 * log_likelihood + 2 * k
bic_value <- -2 * log_likelihood + k * log(n)

# Print results
cat("AIC:", aic_value, "\n")
cat("BIC:", bic_value, "\n")








#monte carlo _____________________________________________________________________

#get the input 
#get the mean and sd for each estimate based on the original sample 'data_pd'

mean_fees <- mean(data_pd$diff_log_fees, na.rm = TRUE)
sd_fees <- sd(data_pd$diff_log_fees, na.rm = TRUE)

mean_tokinc <- mean(data_pd$diff_log_tokinc, na.rm = TRUE)
sd_tokinc <- sd(data_pd$diff_log_tokinc, na.rm = TRUE)

mean_coredev <- mean(data_pd$diff_log_coredev, na.rm = TRUE)
sd_coredev <- sd(data_pd$diff_log_coredev, na.rm = TRUE)

mean_treasury <- mean(data_pd$diff_log_treasury, na.rm = TRUE)
sd_treasury <- sd(data_pd$diff_log_treasury, na.rm = TRUE)

# input the true values from the regression that we got before
true_intercept <- 0.000868
true_fees <- 0.00017
true_tokinc <- 0.000053
true_coredev <- 0.000332
true_treasury <- 0.000842

# Number of simulations
n_sim <- 10000

# Sample size for each simulation
n <- 500  

 #Storage for results
results <- data.frame(
  Intercept = numeric(n_sim),
  Fees = numeric(n_sim),
  Tokinc = numeric(n_sim),
  CoreDev = numeric(n_sim),
  Treasury = numeric(n_sim)
)


# Monte Carlo simulation loop
set.seed(123) 
for (i in 1:n_sim) {
  
  diff_log_fees <- rnorm(n, mean = mean_fees, sd = sd_fees)
  diff_log_tokinc <- rnorm(n, mean = mean_tokinc, sd = sd_tokinc)
  diff_log_coredev <- rnorm(n, mean = mean_coredev, sd = sd_coredev)
  diff_log_treasury <- rnorm(n, mean = mean_treasury, sd = sd_treasury)
  
  
  # Simulate dependent variable (log return)
  log_return <- true_intercept +
    true_fees * diff_log_fees +
    true_tokinc * diff_log_tokinc +
    true_coredev * diff_log_coredev +
    true_treasury * diff_log_treasury +
    rnorm(n, mean = 0, sd = 0.001)  
  
  # Fit the random effects model
  model <- lm(log_return ~ diff_log_fees + diff_log_tokinc + diff_log_coredev + diff_log_treasury)
  

  results[i, ] <- coef(model)
}

# Analyze and visualize the results
summary(results)

# Plot the distributions of the estimates
par(mfrow = c(2, 3))  
hist(results$Intercept, main = "Intercept", xlab = "Estimate", breaks = 30, col = "skyblue")
hist(results$Fees, main = "diff_log_fees", xlab = "Estimate", breaks = 30, col = "skyblue")
hist(results$Tokinc, main = "diff_log_tokinc", xlab = "Estimate", breaks = 30, col = "skyblue")
hist(results$CoreDev, main = "diff_log_coredev", xlab = "Estimate", breaks = 30, col = "skyblue")
hist(results$Treasury, main = "diff_log_treasury", xlab = "Estimate", breaks = 30, col = "skyblue")

#create a table which show how far from true estimate the simulated are 
# Calculate statistics for each coefficient
simulated_means <- colMeans(results)
simulated_medians <- apply(results, 2, median)

# True coefficients
true_coefficients <- c(
  Intercept = true_intercept,
  Fees = true_fees,
  Tokinc = true_tokinc,
  CoreDev = true_coredev,
  Treasury = true_treasury
)

# Calculate percentage differences
diff_mean_pct <- (simulated_means - true_coefficients) / true_coefficients * 100
diff_median_pct <- (simulated_medians - true_coefficients) / true_coefficients * 100

# Create a summary table
summary_table <- data.frame(
  Coefficient = names(true_coefficients),
  True_Estimate = true_coefficients,
  Simulated_Mean = simulated_means,
  Difference_From_Mean_Pct = diff_mean_pct,
  Simulated_Median = simulated_medians,
  Difference_From_Median_Pct = diff_median_pct
)

# Print the table
print(summary_table)
