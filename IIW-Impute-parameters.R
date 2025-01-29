source("IIW-Impute-functions.R")
require(parallel)
ncoresavailable <- parallel::detectCores()
nclusters <- makeCluster(ncoresavailable - 5)


# install.packages("beepr")
library(beepr)



##############################################################################################################
#~~~~~~~~~~~~~~~~~~~~~~ SCHEME A: 4 INTENSITY COVARIATES SUBJECT TO MISSINGNESS  ~~~~~~~~~~~~~~~~~~~~~~~~~~#
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


missingnesstypes = c("Naive", "NoMissingness", "ObsTimesOnly", "MCAR", "MAR", "MNAR")
proportions = c(0.25, 0.5, 0.90)

outcome = "continuous"
usesplines = T

t1 <- Sys.time()
# results_continuous_n100 <- simulateResultsIIW(N, n, beta1, beta2, beta3, beta4, beta5, beta6, gamma1, gamma2, gamma3,
#                                           gamma4, gamma5, gamma6, tau, missingnesstypes, proportions, whichmissing, nimputations,
#                                           outcome, usesplines, scheme = "A",
#                                           outputfulldatalist = TRUE, nclusters = nclusters, scheme)
# saveRDS(results_continuous_n100, "results_continuous_n100.rds")
t2 <- Sys.time()
t2 - t1 #7 hours

results_continuous_n100 <- readRDS("results_continuous_n100.rds")
results_continuous_n100$biasmat
beep()

results_continuous_n100$numevents
results_continuous_n100$newn



cleanUpResults(results_continuous_n100$biasmat, proportions, missingnesstypes, methodnames = c("CCA", "LOCF", "SI", "MI", "missForest"))
cleanUpResults(results_continuous_n100$msemat, proportions, missingnesstypes, methodnames = c("CCA", "LOCF", "SI", "MI", "missForest"))
cleanUpResults(results_continuous_n100$coverage_beta1, proportions, missingnesstypes, methodnames = c("CCA", "LOCF", "SI", "MI", "missForest"))





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
gamma4 = 0.4
gamma5 = -0.3
gamma6 = 0 #0.3
tau = 2
N = 1000
nimputations = 5 #number of imputations for the MI method
whichmissing <- c("V1t", "V2t", "V3t", "V4t") #which covariates we allow to be missing

missingnesstypes = c("Naive", "NoMissingness", "ObsTimesOnly", "MCAR", "MAR", "MNAR")
proportions = c(0.25, 0.5, 0.90)

outcome = "binary"

t1 <- Sys.time()
# results_binary_n100 <- simulateResultsIIW(N, n, beta1, beta2, beta3, beta4, beta5, beta6, gamma1, gamma2, gamma3,
#               gamma4, gamma5, gamma6, tau, missingnesstypes, proportions, whichmissing, nimputations,
#               outcome, usesplines = T, scheme = "A",
#               outputfulldatalist = TRUE, nclusters = nclusters)
# saveRDS(results_binary_n100, "results_binary_n100.rds")
results_binary_n100 <- readRDS("results_binary_n100.rds")

t2 <- Sys.time()
t2 - t1 #7 hours

results_binary_n100 <- readRDS("results_binary_n100.rds")
results_binary_n100$biasmat
beep()

results_binary_n100$numevents
results_binary_n100$newn



cleanUpResults(results_binary_n100$biasmat, proportions, missingnesstypes, methodnames = c("CCA", "LOCF", "SI", "MI", "missForest"))
cleanUpResults(results_binary_n100$msemat, proportions, missingnesstypes, methodnames = c("CCA", "LOCF", "SI", "MI", "missForest"))
cleanUpResults(results_binary_n100$coverage_beta1, proportions, missingnesstypes, methodnames = c("CCA", "LOCF", "SI", "MI", "missForest"))

cleanUpResults(results_binary_n100$extremitysums, proportions, missingnesstypes, methodnames = c("CCA", "LOCF", "SI", "MI", "missForest"))
cleanUpResults(results_binary_n100$bias_beta1_noextremity, proportions, missingnesstypes, methodnames = c("CCA", "LOCF", "SI", "MI", "missForest"))
cleanUpResults(results_binary_n100$mse_beta1_noextremity, proportions, missingnesstypes, methodnames = c("CCA", "LOCF", "SI", "MI", "missForest"))
cleanUpResults(results_binary_n100$coverage_beta1_noextremity, proportions, missingnesstypes, methodnames = c("CCA", "LOCF", "SI", "MI", "missForest"))










##############################################################################################################
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ SCHEME B: less variable ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
##############################################################################################################


######## ~~~~~~~~~~~~~Binary Outcome, scheme B~~~~~~~~~~~~~~~~~~~~ #######
set.seed(9756)
n = 100
beta1 = 0.5
beta2 = 1
beta3 = 0.4
beta4 = 0.5
beta5 = 0.6
beta6 = 0 #0.3
gamma1 = 0.6
gamma2 = 0.6
gamma3 = 0.5
gamma4 = 0.6
gamma5 = -0.2
gamma6 = 0 #0.3
tau = 2
N = 1000
scheme = "B"
nimputations = 5 #number of imputations for the MI method
whichmissing <- c("V1t", "V2t", "V3t", "V4t") #which covariates we allow to be missing


missingnesstypes = c("Naive", "NoMissingness", "ObsTimesOnly", "MCAR", "MAR", "MNAR")
proportions = c(0.25, 0.5, 0.9)

outcome = "binary"
usesplines = T

t1 <- Sys.time()
# results_binary_n100_B <- simulateResultsIIW(N, n, beta1, beta2, beta3, beta4, beta5, beta6, gamma1, gamma2, gamma3,
#                                           gamma4, gamma5, gamma6, tau, missingnesstypes, proportions, whichmissing, nimputations,
#                                           outcome, usesplines = F, scheme = "B",
#                                           outputfulldatalist = TRUE, nclusters = nclusters)
# saveRDS(results_binary_n100_B, "results_binary_n100_B.rds")
t2 <- Sys.time()
t2 - t1 #7 hours

results_binary_n100_B <- readRDS("results_binary_n100_B.rds")
results_binary_n100_B$biasmat
beep()

results_binary_n100_B$numevents
results_binary_n100_B$newn




cleanUpResults(results_binary_n100_B$biasmat, proportions, missingnesstypes, methodnames = c("CCA", "LOCF", "SI", "MI", "missForest"))
cleanUpResults(results_binary_n100_B$msemat, proportions, missingnesstypes, methodnames = c("CCA", "LOCF", "SI", "MI", "missForest"))
cleanUpResults(results_binary_n100_B$coverage_beta1, proportions, missingnesstypes, methodnames = c("CCA", "LOCF", "SI", "MI", "missForest"))

cleanUpResults(results_binary_n100_B$extremitysums, proportions, missingnesstypes, methodnames = c("CCA", "LOCF", "SI", "MI", "missForest"))
cleanUpResults(results_binary_n100_B$bias_beta1_noextremity, proportions, missingnesstypes, methodnames = c("CCA", "LOCF", "SI", "MI", "missForest"))
cleanUpResults(results_binary_n100_B$mse_beta1_noextremity, proportions, missingnesstypes, methodnames = c("CCA", "LOCF", "SI", "MI", "missForest"))
cleanUpResults(results_binary_n100_B$coverage_beta1_noextremity, proportions, missingnesstypes, methodnames = c("CCA", "LOCF", "SI", "MI", "missForest"))






#### CoNTINUOUS OUTCOME

set.seed(76)
n = 100
beta1 = 0.5
beta2 = 1
beta3 = 0.4
beta4 = 0.5
beta5 = 0.6
beta6 = 0 #0.3
gamma1 = 0.6
gamma2 = 0.6
gamma3 = 0.5
gamma4 = 0.7
gamma5 = 0.6
gamma6 = 0 #0.3
tau = 2
N = 1000
scheme = "B"
nimputations = 5 #number of imputations for the MI method
whichmissing <- c("V1t", "V2t", "V3t", "V4t") #which covariates we allow to be missing


missingnesstypes = c("Naive", "NoMissingness", "ObsTimesOnly", "MCAR", "MAR", "MNAR")
proportions = c(0.25, 0.5, 0.9)

outcome = "continuous"
usesplines = T

t1 <- Sys.time()
# results_continuous_n100_B <- simulateResultsIIW(N, n, beta1, beta2, beta3, beta4, beta5, beta6, gamma1, gamma2, gamma3,
#                                           gamma4, gamma5, gamma6, tau, missingnesstypes, proportions, whichmissing, nimputations,
#                                           outcome, usesplines = F, scheme = "B",
#                                           outputfulldatalist = TRUE, nclusters = nclusters)
# saveRDS(results_continuous_n100_B, "results_continuous_n100_B.rds")
t2 <- Sys.time()
t2 - t1 #4 hours

results_continuous_n100_B <- readRDS("results_continuous_n100_B.rds")
results_continuous_n100_B$biasmat
beep()

results_continuous_n100_B$numevents
results_continuous_n100_B$newn



cleanUpResults(results_continuous_n100_B$biasmat, proportions, missingnesstypes, methodnames = c("CCA", "LOCF", "SI", "MI", "missForest"))
cleanUpResults(results_continuous_n100_B$msemat, proportions, missingnesstypes, methodnames = c("CCA", "LOCF", "SI", "MI", "missForest"))
cleanUpResults(results_continuous_n100_B$coverage_beta1, proportions, missingnesstypes, methodnames = c("CCA", "LOCF", "SI", "MI", "missForest"))



