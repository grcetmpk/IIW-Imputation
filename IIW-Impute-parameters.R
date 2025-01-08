source("IIW-Impute-functions.R")
require(parallel)
ncoresavailable <- parallel::detectCores()
nclusters <- makeCluster(ncoresavailable - 5)


# install.packages("beepr")
library(beepr)



##############################################################################################################
#~~~~~~~~~~~~~~~~~~~~~~ SCHEME B: ALL INTENSITY COVARIATES SUBJECT TO MISSINGNESS  ~~~~~~~~~~~~~~~~~~~~~~~~~~#
##############################################################################################################


######## ~~~~~~~~~~~~~Continuous OUTCOME, with splines~~~~~~~~~~~~~~~~~~~~ #######
set.seed(324)
n = 100
beta1 = 0.5
beta2 = 1
beta3 = 0.4
beta4 = 0.3
beta5 = 0.4
beta6 = 0 #0.3
gamma1 = 0.5
gamma2 = 0.4
gamma3 = 0.5
gamma4 = 0.3
gamma5 = -0.3
gamma6 = 0 #0.3
tau = 2
N = 1000
nimputations = 5 #number of imputations for the MI method
whichmissing <- c("V1t", "V2t", "V3t", "V4t") #which covariates we allow to be missing


schemes = c("Naive", "NoMissingness", "ObsTimesOnly", "MCAR", "MAR", "MNAR")
proportions = c(0.25, 0.5, 0.90)

outcome = "continuous"
usesplines = T

t1 <- Sys.time()
# results_continuous_n100 <- simulateResultsIIW(N, n, beta1, beta2, beta3, beta4, beta5, beta6, gamma1, gamma2, gamma3,
#                                           gamma4, gamma5, gamma6, tau, schemes, proportions, whichmissing, nimputations,
#                                           outcome, usesplines,
#                                           outputfulldatalist = TRUE, nclusters = nclusters)
# saveRDS(results_continuous_n100, "results_continuous_n100.rds")
t2 <- Sys.time()
t2 - t1 #7 hours

results_continuous_n100 <- readRDS("results_continuous_n100.rds")
results_continuous_n100$biasmat
beep()

results_continuous_n100$numevents
results_continuous_n100$newn



cleanUpResults(results_continuous_n100$biasmat, proportions, schemes, methodnames = c("CCA", "LOCF", "SI", "MI", "missForest"))
cleanUpResults(results_continuous_n100$msemat, proportions, schemes, methodnames = c("CCA", "LOCF", "SI", "MI", "missForest"))
cleanUpResults(results_continuous_n100$coverage_beta1, proportions, schemes, methodnames = c("CCA", "LOCF", "SI", "MI", "missForest"))





######## ~~~~~~~~~~~~~ Binary OUTCOME, with splines~~~~~~~~~~~~~~~~~~~~ #######
set.seed(2348)
n = 100
beta1 = 0.5
beta2 = 1
beta3 = 0.4
beta4 = 0.3
beta5 = 0.4
beta6 = 0 #0.3
gamma1 = 0.5
gamma2 = 0.4
gamma3 = 0.5
gamma4 = 0.3
gamma5 = -0.3
gamma6 = 0 #0.3
tau = 2
N = 1000
nimputations = 5 #number of imputations for the MI method

schemes = c("Naive", "NoMissingness", "ObsTimesOnly", "MCAR", "MAR", "MNAR")
proportions = c(0.25, 0.5, 0.90)

outcome = "binary"

t1 <- Sys.time()
results_binary_n100 <- simulateResultsIIW(N, n, beta1, beta2, beta3, beta4, beta5, beta6, gamma1, gamma2, gamma3,
                                          gamma4, gamma5, gamma6, tau, schemes, proportions, nimputations,
                                          outcome, usesplines,
                                          outputfulldatalist = TRUE, nclusters = nclusters)
saveRDS(results_binary_n100, "results_binary_n100.rds")
t2 <- Sys.time()
t2 - t1 #7 hours

results_binary_n100 <- readRDS("results_binary_n100.rds")
results_binary_n100$biasmat
beep()

results_binary_n100$numevents
results_binary_n100$newn



cleanUpResults(results_binary_n100$biasmat, proportions, schemes, methodnames = c("CCA", "LOCF", "SI", "MI", "missForest"))
cleanUpResults(results_binary_n100$msemat, proportions, schemes, methodnames = c("CCA", "LOCF", "SI", "MI", "missForest"))
cleanUpResults(results_binary_n100$coverage_beta1, proportions, schemes, methodnames = c("CCA", "LOCF", "SI", "MI", "missForest"))

cleanUpResults(results_binary_n100$extremitysums, proportions, schemes, methodnames = c("CCA", "LOCF", "SI", "MI", "missForest"))
cleanUpResults(results_binary_n100$bias_beta1_noextremity, proportions, schemes, methodnames = c("CCA", "LOCF", "SI", "MI", "missForest"))
cleanUpResults(results_binary_n100$mse_beta1_noextremity, proportions, schemes, methodnames = c("CCA", "LOCF", "SI", "MI", "missForest"))
cleanUpResults(results_binary_n100$, proportions, schemes, methodnames = c("CCA", "LOCF", "SI", "MI", "missForest"))





##############################################################################################################
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ SCHEME A: ONE COVARIATE MISSING ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
##############################################################################################################


######## ~~~~~~~~~~~~~Continuous OUTCOME, with splines~~~~~~~~~~~~~~~~~~~~ #######
set.seed(324)
n = 100
beta1 = 0.5
beta2 = 1
beta3 = 0.4
beta4 = 0.3
beta5 = 0.4
beta6 = 0 #0.3
gamma1 = 0.5
gamma2 = 0.4
gamma3 = 0.5
gamma4 = 0.3
gamma5 = -0.3
gamma6 = 0 #0.3
tau = 2
N = 1000
nimputations = 5 #number of imputations for the MI method
whichmissing <- c("V1t") #which covariates we allow to be missing


schemes = c("Naive", "NoMissingness", "ObsTimesOnly", "MCAR", "MAR", "MNAR")
proportions = c(0.25, 0.5, 0.90)

outcome = "continuous"
usesplines = T

t1 <- Sys.time()
results_continuous_n100_V1t <- simulateResultsIIW(N, n, beta1, beta2, beta3, beta4, beta5, beta6, gamma1, gamma2, gamma3,
                                          gamma4, gamma5, gamma6, tau, schemes, proportions, whichmissing, nimputations,
                                          outcome, usesplines,
                                          outputfulldatalist = TRUE, nclusters = nclusters)
saveRDS(results_continuous_n100_V1t, "results_continuous_n100_V1t.rds")
t2 <- Sys.time()
t2 - t1 #7 hours

results_continuous_n100_V1t <- readRDS("results_continuous_n100_V1t.rds")
results_continuous_n100_V1t$biasmat
beep()

results_continuous_n100_V1t$numevents
results_continuous_n100_V1t$newn



cleanUpResults(results_continuous_n100_V1t$biasmat, proportions, schemes, methodnames = c("CCA", "LOCF", "SI", "MI", "missForest"))
cleanUpResults(results_continuous_n100_V1t$msemat, proportions, schemes, methodnames = c("CCA", "LOCF", "SI", "MI", "missForest"))
cleanUpResults(results_continuous_n100_V1t$coverage_beta1, proportions, schemes, methodnames = c("CCA", "LOCF", "SI", "MI", "missForest"))









