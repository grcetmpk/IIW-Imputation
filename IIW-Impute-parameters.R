source("IIW-Impute-functions.R")

set.seed(5976)
n = 100
beta1 = 0.5
beta2 = 1
beta3 = 0.9
beta4 = 0.5
beta5 = 0.5
beta6 = 0
gamma1 = 0.5
gamma2 = 0.5
gamma3 = 0.6
gamma4 = 0.6
gamma5 = 0.5
gamma6 = 0.3
tau = 2
N = 100
nimputations = 5 #number of imputations for the MI method

schemes = c("Naive", "NoMissingness", "ObsTimesOnly", "MCAR", "MAR", "MNAR")
proportions = c(0.1, 0.5, 0.75)


# simulateOneIIW(n, beta1, beta2, beta3, beta4, beta5, gamma1, gamma2, gamma3,
#                           gamma4, gamma5, gamma6, tau, schemes, proportions, nimputations)

time1 <- Sys.time()
testresults <- simulateResultsIIW(N, n, beta1, beta2, beta3, beta4, beta5, gamma1, gamma2, gamma3, 
                               gamma4, gamma5, gamma6, tau, schemes, proportions, nimputations,
                               outputfulldatalist = FALSE, inParallel = F, nclusters = NULL)
time2 <- Sys.time()
tottime_1sim <- time2 - time1
tottime_1sim # few outliers

testresults


# parallelization stuff
ncores <- detectCores() - 1
nclusters <- makeCluster(ncores)
inParallel = F


