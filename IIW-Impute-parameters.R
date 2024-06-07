source("IIW-Impute-functions.R")

set.seed(500)
n = 100
beta1 = 0.5
beta2 = 0.4
beta3 = 0.4
beta4 = 0.5
beta5 = 0.3
beta6 = 0
gamma1 = 0.2
gamma2 = 0.2
gamma3 = 0.3
gamma4 = 0.2
gamma5 = 0.3
gamma6 = 0.1
tau = 3
N = 1000

schemes = c("Naive", "NoMissingness", "ObsTimesOnly", "MCAR", "MAR", "MNAR")
proportions = c(0.05, 0.2, 0.5, 0.7, 0.9)


time.start <- Sys.time()
simulateOneIIW(n, beta1, beta2, beta3, beta4, beta5, gamma1, gamma2, gamma3, 
                          gamma4, gamma5, gamma6, tau, schemes, proportions)
time.end <- Sys.time()

simtime <- time.end - time.start

# parallelization stuff
ncores <- detectCores() - 1
nclusters <- makeCluster(ncores)
inParallel = F

#total estimated time for 1000 sim runs - 300 hours lol (2 weeks)


