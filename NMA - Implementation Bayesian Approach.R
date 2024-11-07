# Jihoon Lim
# Network Meta-Analysis Code
# Bayesian Approach

## Part A: Data Setup ####
# Install packages
install.packages("devtools")
devtools::install_github("MathiasHarrer/dmetar")
if (!require("remotes")) {
  install.packages("remotes")
}
remotes::install_github("MathiasHarrer/dmetar")
library(dmetar)
library(meta)
library(metafor)
install.packages("gemtc"); library(gemtc)
install.packages("rjags"); library(rjags)
library(grid)

# Read in Data
retention <- read.csv("C:/Users/limji/Desktop/McGill/Network Meta-Analysis/FINAL/retention.csv", header = TRUE, sep = ",")
trt_primary <- read.table(textConnection('id description
Buprenorphine "Buprenorphine"
Methadone "Methadone"
Naltrexone "Naltrexone"
Control "Control"
SROM "SROM"
                                            '), header = TRUE)

## Part B: Network Meta-Analysis ####
# | 1. Initialize analysis ####
## || a. Network Setup ####
network_primary <- mtc.network(data.ab = retention, treatments = trt_primary, description = "Bayesian NMA")
## || b. Network Evidence Map ####
plot(network_primary, use.description = TRUE, main = "Network Geometry for Retention")
## || c. Summary of Network Evidence ####
summary(network_primary)
print(network_primary)

# | 2. Set up Network Model ####
## || a. Fixed and Random effects model ####
model_primary_fe <- mtc.model(network_primary, likelihood = "binom", link = "log", linearModel = 'fixed', n.chain = 4)
model_primary_re <- mtc.model(network_primary, likelihood = "binom", link = "log", linearModel = 'random', n.chain = 4)

## || b. MCMC simulation - Burn-in 5000, iteration 10000, thin 100 ####
mcmc_primary_fe <- mtc.run(model_primary_re, n.adapt = 5000, n.iter = 100000, thin = 100); summary(mcmc_primary_fe)
mcmc_primary_re <- mtc.run(model_primary_re, n.adapt = 5000, n.iter = 100000, thin = 100); summary(mcmc_primary_re)
# Both models demonstrated similar DIC values. However, the FE model assumes that all studies share the same 
# common effect. However, it may not be reasonable to assume that there is one common effect size. On the other
# hand, the RE model assumes that the observed estimates of the treatment effect can vary across studies due
# to systematic differences in the treatment effect as well as random variations due to chance 
# (Riley, Higgins, & Deeks 2011, The BMJ).

## || c. Assess convergence of the MCMC algorithm ####
plot(mcmc_primary_re)
# Gelman-Rubin Plot
gelman.plot(mcmc_primary_re)
# Potential Scale Reduction Factor (PSRF) - Compares variation within each chain to variation between chain
gelman.diag(mcmc_primary_re)$mpsrf
## If PSRF < 1.05, we assume convergence of our MCMC simulations.

## || d. Assess consistency of the network model ####
# Node Splitting (This takes around 30 minutes!)
nodesplit_primary <- mtc.nodesplit(network_primary, 
                                   linearModel = "random",
                                   likelihood = "binom",
                                   link = "log",
                                   n.adapt = 5000, 
                                   n.iter = 100000, 
                                   thin = 100)
mtc.nodesplit.comparisons(network_primary)

# Summary of node splitting
summary(nodesplit_primary)
plot(summary(nodesplit_primary))
## If p-value < 0.05, there is inconsistency in our network.
# A. Check the direction of association
# B. Possible that studies comparing A and B included systematically different populations than other
# studies that also assessed A.
# C. Study design differences - Parallel vs. Crossover, etc.
# D. Direct comparison done on small sample with only a limited number of studies
# E. Warrants further studies comparing A and B

## Part C: Results ####
# | 1. Rankogram ####
rank.prob.primary <- rank.probability(mcmc_primary_re, preferredDirection = 1)
plot(rank.prob.primary, beside = TRUE, main = "Rankogram of the Preference of Treatment Options for Retention")

# | 2. Surface under the cumulative ranking (SUCRA) score ####
sucra.primary <- dmetar::sucra(rank.prob.primary, lower.is.better = FALSE); sucra.primary
plot(sucra.primary, main = "Surface Under the Cumulative Ranking Score of Treatment Options for Retention")

# | 3. Effect Size Estimates ####
results.primary <- round(exp(relative.effect.table(mcmc_primary_re)), 2); results.primary
write.csv(results.primary, "C:/Users/limji/Desktop/McGill/Network Meta-Analysis/results_primary.csv")

# | 4. Forest plot - Placebo as Reference ####
forest(relative.effect(mcmc_primary_re, t1 = "Control"), digits = 3)
grid.text("Risk Ratio for Retention for Each Medication against Control", .5, .7, gp = gpar(cex = 1))

## Part D: Network Meta-Regression ####
study.info <- read.csv("C:/Users/limji/Desktop/McGill/Network Meta-Analysis/meta-regression.csv", header = TRUE, sep = ",")
net.meta.reg <- list(retention = retention, study.info = study.info, treatments_bin = trt_primary)

# | 1. Meta-Regression Setup ####
network.mr <- mtc.network(data.ab = net.meta.reg$retention,
                          studies = net.meta.reg$study.info,
                          treatments = net.meta.reg$treatments_bin)
regressor1 <- list(coefficient = "shared", variable = "rob", control = "Control")
regressor2 <- list(coefficient = "shared", variable = "year_bin", control = "Control")

# | 2. Compile model ####
model.mr1 <- mtc.model(network.mr, likelihood = "binom", link = "log", type = "regression", regressor = regressor1)
model.mr2 <- mtc.model(network.mr, likelihood = "binom", link = "log", type = "regression", regressor = regressor2)

# | 3. MCMC simulation ####
mcmc.mr1 <- mtc.run(model.mr1, n.adapt = 5000, n.iter = 100000, thin = 100); summary(mcmc.mr1)
mcmc.mr2 <- mtc.run(model.mr2, n.adapt = 5000, n.iter = 100000, thin = 100); summary(mcmc.mr2)

# | 4. Effect Size Estimates ####
## || a. Risk of Bias ####
# Moderate to High
results.a.mr1 <- round(exp(relative.effect.table(mcmc.mr1, covariate = 1)), 2); results.a.mr1
# Low
results.a.mr0 <- round(exp(relative.effect.table(mcmc.mr1, covariate = 0)), 2); results.a.mr0

## || b. Year of Publication ####
# Before 2010
results.b.mr1 <- round(exp(relative.effect.table(mcmc.mr2, covariate = 0)), 2); results.b.mr1
# 2010 and After
results.b.mr0 <- round(exp(relative.effect.table(mcmc.mr2, covariate = 1)), 2); results.b.mr0

# | 5. Forest Plot ####
## || a. Risk of Bias ####
# Moderate to High
forest(relative.effect(mcmc.mr1, t1 = "Control", covariate = 1), digits = 3)
grid.text("Network Meta-Regression: Moderate to High Risk of Bias", .5, .7, gp = gpar(cex = 1))
# Low
forest(relative.effect(mcmc.mr1, t1 = "Control", covariate = 0), digits = 3)
grid.text("Network Meta-Regression: Low Risk of Bias", .5, .7, gp = gpar(cex = 1))

## || b. Year of Publication ####
# Before 2010
forest(relative.effect(mcmc.mr2, t1 = "Control", covariate = 0), digits = 3)
grid.text("Network Meta-Regression: Publication Year Before 2010", .5, .7, gp = gpar(cex = 1))
# 2010 and After
forest(relative.effect(mcmc.mr2, t1 = "Control", covariate = 1), digits = 3)
grid.text("Network Meta-Regression: Publication Year After 2010", .5, .7, gp = gpar(cex = 1))

# | 6. Compare Deviance Information Criterion (DIC) ####
## Network Meta-Regression vs. Original NMA Model
# Risk of Bias
summary(mcmc.mr1)$DIC; summary(mcmc_primary_re)$DIC
# Year of Publication
summary(mcmc.mr2)$DIC; summary(mcmc_primary_re)$DIC

## Part E: Secondary Outcome Analysis ####
# | 1. Read in Data ####
opioid <- read.csv("C:/Users/limji/Desktop/McGill/Network Meta-Analysis/opioid_use.csv", header = TRUE, sep = ",")
trt_secondary <- read.table(textConnection('id description
Buprenorphine "Buprenorphine"
Methadone "Methadone"
Naltrexone "Naltrexone"
Control "Control"
                                            '), header = TRUE)

# | 2. Initialize analysis ####
## || a. Network Setup ####
network_secondary <- mtc.network(data.ab = opioid, treatments = trt_secondary, description = "Bayesian NMA")
## || b. Network Evidence Map ####
plot(network_secondary, use.description = TRUE, main = "Network Geometry for Opioid Use")
## || c. Summary of Network Evidence ####
summary(network_secondary)
print(network_secondary)

# | 3. Set up Network Model ####
## || a. Fixed and random effect model ####
model_secondary_fe <- mtc.model(network_secondary, likelihood = "binom", link = "log", 
                                linearModel = 'fixed', n.chain = 4)
model_secondary_re <- mtc.model(network_secondary, likelihood = "binom", link = "log", 
                                linearModel = 'random', n.chain = 4)

## || b. MCMC simulation - Burn-in 5000, iteration 10000, thin 100 ####
mcmc_secondary_fe <- mtc.run(model_secondary_re, n.adapt = 5000, n.iter = 100000, thin = 100)
summary(mcmc_secondary_fe)
mcmc_secondary_re <- mtc.run(model_secondary_re, n.adapt = 5000, n.iter = 100000, thin = 100)
summary(mcmc_secondary_re)

## || c. Assess convergence of the MCMC algorithm ####
plot(mcmc_secondary_re)
# Gelman-Rubin Plot
gelman.plot(mcmc_secondary_re)
# Potential Scale Reduction Factor (PSRF) - Compares variation within each chain to variation between chain
gelman.diag(mcmc_secondary_re)$mpsrf
## If PSRF < 1.05, we assume convergence of our MCMC simulations.

## || d. Assess consistency of the network model ####
# Node Splitting
nodesplit_secondary <- mtc.nodesplit(network_secondary, 
                                     linearModel = "random",
                                     likelihood = "binom",
                                     link = "log",
                                     n.adapt = 5000, 
                                     n.iter = 100000, 
                                     thin = 100)
mtc.nodesplit.comparisons(network_secondary)
# Summary of node splitting
summary(nodesplit_secondary)
plot(summary(nodesplit_secondary), digits = 3)

# | 4. Results ####
## || a. Rankogram ####
rank.prob.secondary <- rank.probability(mcmc_secondary_re, preferredDirection = -1)
plot(rank.prob.secondary, beside = TRUE, main = "Rankogram of the Preference of Treatment Options for Opioid Use")

## || b. Surface under the cumulative ranking (SUCRA) score ####
sucra.secondary <- dmetar::sucra(rank.prob.secondary, lower.is.better = FALSE); sucra.secondary
plot(sucra.secondary)

## || c. Effect Size Estimates ####
results.secondary <- round(exp(relative.effect.table(mcmc_secondary_re)), 2); results.secondary
write.csv(results.secondary, "C:/Users/limji/Desktop/McGill/Network Meta-Analysis/results_secondary.csv")

## || d. Forest plot - Placebo as Reference ####
forest(relative.effect(mcmc_secondary_re, t1 = "Control"), digits = 2)
grid.text("RR for Illicit Opioid Use for Each Medication against Control", .5, .7, gp = gpar(cex = 1))
