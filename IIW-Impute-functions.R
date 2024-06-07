require(survival)
require(geepack)
require(nleqslv)
require(knitr)
require(kableExtra)
require(dplyr)
require(doParallel)
library(Rcpp)
require(mice)
library(missForest)

expit <- function(x){ return(exp(x)/(1+exp(x)))}

ncoresavailable <- parallel::detectCores()

##### This script contains the functions necessary for data generation and analysis/imputation
##### for the third project in Grace's thesis.


### Simulate irreg longitudinal data using bernoulli draws, discretizing the observation times
### from 0 to \tau, by 0.01. First, starting with a non-time-varying treatment


gendata_IIW <- function(n, beta1, beta2, beta3, beta4, beta5, gamma1, gamma2, gamma3, gamma4, gamma5, gamma6, tau){
  
  # Simulates observation times using bernoulli draws with probabilities proportional 
  # to the intensity. assume the treatment is time-invariant and randomized
  #
  # n: number of subjects
  # tau: maximum follow up
  
  
  mu_V1_D0 <- 1
  var_V1_D0 <- 0.25 
  mu_V1_D1 <- 2
  var_V1_D1 <- 0.5
  
  var_phi <- 0.25
  var_epsilon <- 1
  
  
  # generate subjects one by one and then combine
  id <- 1
  simdatafull <- data.frame(matrix(NA, ncol = 16))
  colnames(simdatafull) <- c("id", "times", "D", "V1t", "V2t", "V3t", "V4t", "V5t", "censortime", "eta", "cexp_V1_D", "cexp_V2_D", 
                             "cexp_V3_D", "cexp_V4_D", "cexp_V5_D", "y")
  
  
  disctimes <- seq(0, tau, by = 0.01)
  
  while(id <= n){
    
    ## Generate Covariates
    
    ## Generate baseline covariate
    
    # generate treatment assignment (time-invariant) at each time point
    D <- rep(rbinom(1, 1, 0.5), length(disctimes))
    
    
    # generate observation times confounder (time varying) V1t ~ N( mu_V_D0, var_V_D0) if D(t) = 0 and 
    # N( mu_V_D1, var_V_D1) if D = 1
    
    V1t <- ifelse(D == 0, rnorm(1, mu_V1_D0, sqrt(var_V1_D0)), rnorm(1, mu_V1_D1, sqrt(var_V1_D1)))
    
    # generate V2 (which will be used in V2(t) = V2log(t)
    
    V2 <- ifelse(D[1] == 0, runif(1, 0, 1), runif(1, 1, 2))
    
    # generate V3(t) which is categorical
    if(D[1] == 0){
      pV3t <- c(0.1, 0.4, 0.5)
    }else{
      pV3t <- c(0.7, 0.2, 0.1)
    }
    
    V3t <- sample(c(0,1,2), length(disctimes), pV3t, replace = TRUE)
    
    
    # generate V4(t)
    
    V4t <- rbinom(length(disctimes), 1, 0.3)
    
    # generate V5, which will be used in V5(t)
    
    V5 <- runif(1,0,1)
    
    
    # generate random effect
    phi <- rnorm(1, 0, var_phi)
    
    
    # simulate censoring time (uniform)
    censortime <- runif(1, tau/2, tau)
    
    # calculate eta
    eta <- rgamma(1, shape = 100, scale = 0.01) #eta is for random effect with mean 1 and sd 0.1 
    #(obs times within subject correlated if sd !0)
    
    simdata <- data.frame(rep(id, length(disctimes)), disctimes, D, V1t, V2*log(disctimes + 1), V3t, V4t, V5*sqrt(disctimes)/2, 
                          rep(censortime, length(disctimes)), rep(eta, length(disctimes)))
    colnames(simdata) <- c("id", "times", "D", "V1t", "V2t", "V3t", "V4t", "V5t", "censortime", "eta")
    
    # need conditional expectation and variance of Z | X
    
    simdata$cexp_V1_D <- ifelse(simdata$D == 0, mu_V1_D0, mu_V1_D1)
    simdata$cexp_V2_D <- ifelse(simdata$D == 0, 0.5, 1.5)*log(simdata$times + 1)
    simdata$cexp_V3_D <- ifelse(simdata$D == 0, 0.1*0 + 0.4*1 + 0.5*2, 0.7*0 + 0.2*1 + 0.1*2)
    simdata$cexp_V4_D <- 0.3 # not dependent on D
    simdata$cexp_V5_D <- 0.5 # not dependent on D (and not involved in the generation of Y) 
    
    # generate time-varying outcome at each possible time
    simdata$y <- (2 - simdata$times) + beta1*simdata$D + beta2*(simdata$V1 - simdata$cexp_V1_D) +
      beta3*(simdata$V2 - simdata$cexp_V2_D) + beta4*(simdata$V3 - simdata$cexp_V3_D) + beta5*(simdata$V4 - simdata$cexp_V4_D)+
      rep(phi, dim(simdata)[1]) + rnorm(dim(simdata)[1], 0, sqrt(var_epsilon))
    
    
    
    simdatafull <- rbind(simdatafull, simdata)
    
    
    
    
    id = id + 1
    
  }
  
  simdata <- simdatafull[-1,] #remove empty first row
  
  ### calculate intensities for each counterfactual observation time:
  
  
  simdata$intensities <- simdata$eta*sqrt(simdata$times)/2*exp(gamma1*simdata$D + gamma2*simdata$V1t + gamma3*simdata$V2t +
                                                         gamma4*simdata$V3t + gamma5*simdata$V4t + gamma6*simdata$V5t)
  simdata$prObs <- 0.01*simdata$intensities#simdata$intensities/max(simdata$intensities)
  simdata$Yobserved <- rbinom(nrow(simdata), 1, simdata$prObs)
  
  
  simdata <- simdata %>%
    filter(times < censortime)
  
  
  ### filter censored observations, unobserved observations (do this before or after calculating maximum intensity???)
  simdata_observed <- simdata %>%
    filter(Yobserved == 1)
  

  #check the generated data
  numevents <- summary(tapply(simdata_observed$Yobserved, simdata_observed$id, sum)) 
  
  
  # #also get data for individuals at baseline (all same here)
  newn <- length(unique(simdata_observed$id)) #number of people after censoring etc
  
  out <- list(simdata, simdata_observed, numevents, newn)
  names(out) <- c("simdata_full", "simdata_Yobserved", "numevents", "newn")
  
  return(out) 
  
}


makeMissing <- function(intensitydat, scheme, p){

  if(scheme == "NoMissingness"){
    
    intensitydat$Zobserved <- 1
    return(intensitydat) # return full dataset with no missingness
    
    }else if(scheme %in% c("ObsTimesOnly", "Naive")){
      
      intensitydat$Zobserved <-  intensitydat$Yobserved # Z only measured if it Y is observed
      
      # change any data with missing Z to NAs
      intensitydat$V1t[intensitydat$Zobserved == 0] <- NA
      intensitydat$V2t[intensitydat$Zobserved == 0] <- NA
      intensitydat$V3t[intensitydat$Zobserved == 0] <- NA
      intensitydat$V4t[intensitydat$Zobserved == 0] <- NA
      intensitydat$V5t[intensitydat$Zobserved == 0] <- NA
      
      
      return(intensitydat)
      
    }else if(scheme == "MCAR"){
      # randomly select p% of observations to be missing

      intensitydat$Zobserved <- rbinom(dim(intensitydat)[1], 1, 1 - p) 
      
      # change any data with missing = 1 to NAs
      intensitydat$V1t[intensitydat$Zobserved == 0] <- NA
      intensitydat$V2t[intensitydat$Zobserved == 0] <- NA
      intensitydat$V3t[intensitydat$Zobserved == 0] <- NA
      intensitydat$V4t[intensitydat$Zobserved == 0] <- NA
      intensitydat$V5t[intensitydat$Zobserved == 0] <- NA
      
      return(intensitydat)
      
    }else if(scheme == "MAR"){
      
      #determine the level of eta0 given eta1 to get the specified missingness level
      
      x <- as.matrix(intensitydat[, 3]) #matrix form of dataset for covariate in MAR model (D)
      eta1 <- 3 #set as anything for now
      
      f <- function(t) {            # Define a path through parameter space
        sapply(t, function(y) mean(1 / (1 + exp(-y -x %*% eta1))))
        # (sapply makes this function vectorizable)
      }
      
      # Find parameters (alpha, beta) yielding any specified proportions `p`.
      
                           
      results <- sapply(p, function(p) {
        alpha <- uniroot(function(t) f(t) - p, c(-1e6, 1e6), tol = .Machine$double.eps^0.5)$root
        c(alpha, f(alpha))})
      dimnames(results) <- list(c("alpha", "f(alpha)"), p=p)
      
      eta0 <- results[1,1]
      
      pmiss <- expit(eta0 + eta1*intensitydat$D)
      intensitydat$Zobserved <- rbinom(length(pmiss), 1, 1 - pmiss) #bernoulli draws to see which observed
      #sum(intensitydat$missing)/length(intensitydat$missing)
      
      # change any data with missing = 1 to NAs
      intensitydat$V1t[intensitydat$Zobserved == 0] <- NA
      intensitydat$V2t[intensitydat$Zobserved == 0] <- NA
      intensitydat$V3t[intensitydat$Zobserved == 0] <- NA
      intensitydat$V4t[intensitydat$Zobserved == 0] <- NA
      intensitydat$V5t[intensitydat$Zobserved == 0] <- NA
      
      
      return(intensitydat)
    
    }else if(scheme == "MNAR"){
      
      x <- as.matrix(intensitydat[, c(3, 8, 16)]) #matrix form of dataset for covariate in MAR model (D)
      eta <- c(3, 2, 2) #set as anything for now
      
      f <- function(t) {            # Define a path through parameter space
        sapply(t, function(y) mean(1 / (1 + exp(-y -x %*% eta))))
        # (sapply makes this function vectorizable)
      }
      
      # Find parameters (alpha, beta) yielding any specified proportions `p`.
      
      
      results <- sapply(p, function(p) {
        alpha <- uniroot(function(t) f(t) - p, c(-1e6, 1e6), tol = .Machine$double.eps^0.5)$root
        c(alpha, f(alpha))})
      dimnames(results) <- list(c("alpha", "f(alpha)"), p=p)
      
      eta0 <- results[1,1]
      
      pmiss <- expit(eta0 + eta[1]*intensitydat$D + eta[2]*intensitydat$V5t+eta[3]*intensitydat$y)
      intensitydat$Zobserved <- rbinom(length(pmiss), 1, 1 - pmiss) #bernoulli draws to see which observed
      #sum(intensitydat$Zobserved)/length(intensitydat$Zobserved)
      
      # change any data with missing = 1 to NAs
      intensitydat$V1t[intensitydat$Zobserved == 0] <- NA
      intensitydat$V2t[intensitydat$Zobserved == 0] <- NA
      intensitydat$V3t[intensitydat$Zobserved == 0] <- NA
      intensitydat$V4t[intensitydat$Zobserved == 0] <- NA
      intensitydat$V5t[intensitydat$Zobserved == 0] <- NA
      
      return(intensitydat)
      
    }else{
      
      return(NA)
    }
  }


IIW <- function(intensitydat, obsdat){
  
  ##### perform IIW on the imputed data set
  
  #include a variable counting observation number, to be used for lagging time for Surv function
  intensitydat$obsnumber <- with(intensitydat, ave(id, id, FUN = seq_along))
  
  # #create lagged time variable
  intensitydat$times.lag <- intensitydat$times[c(nrow(intensitydat ),1:(nrow(intensitydat )-1))]
  intensitydat$times.lag[intensitydat$obsnumber == 1] <- -0.01
  
  gamma.hat <- coxph(Surv(times.lag, times, Yobserved) ~ D  + V1t + V2t + V3t + V4t + V5t - 1, data = intensitydat)$coef
  delta.hat <- coxph(Surv(times.lag, times, Yobserved) ~ D - 1, data = intensitydat)$coef
  
  ### perform IIW on the "observed" X and Y data
  
  iiw <- exp(cbind(obsdat$D)%*%delta.hat)/
    exp(cbind(obsdat$D, obsdat$V1t, obsdat$V2t, 
              obsdat$V3t, obsdat$V4t, obsdat$V5t)%*%gamma.hat)
  
  beta1 <- summary(glm(y ~ D + offset(2-times) - 1, data=obsdat, weights = iiw))$coef[1,1]
  
  return(beta1)
  
  
}


CCA <- function(intensitydat, obsdat){
  
  #intensitydat -> data used for IIW which may or may not have missingness (indicated by Zmissing)
  #obsdat -> data on longitudinal outcome, only observed at "observation times" with no other missingness, used to fit the outcome model
  
  
  #complete case analysis
      
  intensitydat <- na.omit(intensitydat)
  out <- IIW(intensitydat, obsdat)
  return(out)
  
  # #include a variable counting observation number, to be used for lagging time for Surv function
  # intensitydat$obsnumber <- with(intensitydat, ave(id, id, FUN = seq_along))
  #     
  # # #create lagged time variable
  # intensitydat$times.lag <- intensitydat$times[c(nrow(intensitydat ),1:(nrow(intensitydat )-1))]
  # intensitydat$times.lag[intensitydat$obsnumber == 1] <- -0.01
  #     
  # gamma.hat <- coxph(Surv(times.lag, times, Yobserved) ~ D  + V1t + V2t + V3t + V4t + V5t - 1, data = intensitydat)$coef
  # delta.hat <- coxph(Surv(times.lag, times, Yobserved) ~ D - 1, data = intensitydat)$coef
  #     
  # 
  # #use gamma.hat and delta.hat to estimate IIW, one estimate per observation in the obsdat
  # 
  # iiw <- exp(cbind(obsdat$D)%*%delta.hat)/
  # exp(cbind(obsdat$D, obsdat$V1t, obsdat$V2t, obsdat$V3t, obsdat$V4t, 
  #           obsdat$V5t)%*%gamma.hat)
  #     
  # beta1 <- summary(glm(y ~ D + offset(2-times) - 1, data=obsdat, weights = iiw))$coef[1,1]
  # 
  # return(beta1)
}


LOCF <- function(intensitydat, obsdat){
  

    #do LOCF
    imputedintensitydata <- intensitydat %>%
      group_by(id) %>%
      tidyr::fill(names(intensitydat), .direction = "downup") %>%
      ungroup()
    
    out <- IIW(imputedintensitydata, obsdat)
    return(out)
    
    # ##### perform IIW on the imputed data set
    # 
    # #include a variable counting observation number, to be used for lagging time for Surv function
    # imputedintensitydata$obsnumber <- with(imputedintensitydata, ave(id, id, FUN = seq_along))
    # 
    # # #create lagged time variable
    # imputedintensitydata$times.lag <- imputedintensitydata$times[c(nrow(imputedintensitydata ),1:(nrow(imputedintensitydata )-1))]
    # imputedintensitydata$times.lag[imputedintensitydata$obsnumber == 1] <- -0.01
    # 
    # gamma.hat <- coxph(Surv(times.lag, times, Yobserved) ~ D  + V1t + V2t + V3t + V4t + V5t - 1, data = imputedintensitydata)$coef
    # delta.hat <- coxph(Surv(times.lag, times, Yobserved) ~ D - 1, data = imputedintensitydata)$coef
    # 
    # ### perform IIW on the "observed" X and Y data
    # 
    # iiw <- exp(cbind(obsdat$D)%*%delta.hat)/
    #   exp(cbind(obsdat$D, obsdat$V1t, obsdat$V2t, 
    #             obsdat$V3t, obsdat$V4t, obsdat$V5t)%*%gamma.hat)
    # 
    # beta1 <- summary(glm(y ~ D + offset(2-times) - 1, data=obsdat, weights = iiw))$coef[1,1]
    # 
    # return(beta1)
  }
  
  
SI <- function(intensitydat, obsdat){
 #regression switching imputation with PPM
  
  # only use observed Z covariates (D, V1t - v5t) for imputation
  
  intensitydat2 <- intensitydat[, 1:8]
  intensitydat2$D <- as.factor(intensitydat2$D)
  intensitydat2$V4t <- as.factor(intensitydat2$V4t)
  
  imp <- mice::mice(intensitydat2, m = 1, maxit = 5, meth = "pmm", seed = 100, printFlag = F)
  imputedintensitydat <- complete(imp)
  
  intensitydat_imputed_full <- intensitydat #copy original df with missingness
  
  #replace columns with imputed data
  intensitydat_imputed_full$V1t <- imputedintensitydat$V1t
  intensitydat_imputed_full$V2t <- imputedintensitydat$V2t
  intensitydat_imputed_full$V3t <- imputedintensitydat$V3t
  intensitydat_imputed_full$V4t <- imputedintensitydat$V4t
  intensitydat_imputed_full$V5t <- imputedintensitydat$V5t
  
  out <- IIW(intensitydat_imputed_full, obsdat)
  
  return(out)
}


MI <- function(intensitydat, obsdat){
  #regression switching imputation with PPM
  

  #include a variable counting observation number, to be used for lagging time for Surv function
  intensitydat$obsnumber <- with(intensitydat, ave(id, id, FUN = seq_along))

  # #create lagged time variable
  intensitydat$times.lag <- intensitydat$times[c(nrow(intensitydat ),1:(nrow(intensitydat )-1))]
  intensitydat$times.lag[intensitydat$obsnumber == 1] <- -0.01
  
  
  
  imp <- mice::mice(intensitydat, m = 5, maxit = 5, meth = "pmm", seed = 100, printFlag = F)
  imputedintensitydat <- complete(imp)
  
  intensitydat_imputed_full <- intensitydat #copy original df with missingness
  

  ##### perform IIW on the imputed data set
  
  gamma.hat <- summary(pool(with(imp, coxph(Surv(times.lag, times, Yobserved) ~ D  + V1t + V2t + V3t + V4t + V5t - 1, data = intensitydat))))[,2]
  delta.hat <- summary(pool(with(imp, coxph(Surv(times.lag, times, Yobserved) ~ D - 1, data = intensitydat))))[,2]
  
  ### perform IIW on the "observed" X and Y data
  
  iiw <- exp(cbind(obsdat$D)%*%delta.hat)/
    exp(cbind(obsdat$D, obsdat$V1t, obsdat$V2t, 
              obsdat$V3t, obsdat$V4t, obsdat$V5t)%*%gamma.hat)
  
  out <- summary(glm(y ~ D + offset(2-times) - 1, data=obsdat, weights = iiw))$coef[1,1]

  return(out)
}
  
  
NIRF <- function(intensitydat, obsdat){

  intensitydat2 <- intensitydat
  intensitydat2$D <- as.factor(intensitydat2$D)
  intensitydat2$V4t <- as.factor(intensitydat2$V4t)
  intensitydat2$Yobserved <- as.factor(intensitydat2$Yobserved)
  intensitydat2$V3t <- as.factor(intensitydat2$V3t)
  
  
  registerDoParallel(cores= ncoresavailable-2)
  start.time <-Sys.time()
  imp <- missForest(intensitydat2, maxiter = 5, ntree = 50, parallelize = "forests") #change to larger numbers later
  end.time <- Sys.time()
  time.taken <- end.time - start.time
  # time.taken
  imputedintensitydata <- imp$ximp
  
  #convert back to numeric data
  imputedintensitydata$D <- as.numeric(imputedintensitydata$D)
  imputedintensitydata$V4t <- as.numeric(imputedintensitydata$V4t)
  imputedintensitydata$Yobserved <- as.numeric(imputedintensitydata$Yobserved)
  imputedintensitydata$V3t <- as.numeric(imputedintensitydata$V3t)
  
  
  out <- IIW(imputedintensitydata, obsdat)
  return(out)
  
}  



simulateOneIIW<- function(n, beta1, beta2, beta3, beta4, beta5, gamma1, gamma2, gamma3, 
                          gamma4, gamma5, gamma6, tau, schemes, proportions){
  # Simulates one instance of the simulation, obtaining estimates for beta1 under various weighting
  # IIW uses stabilized weights
  
  
  # get required data from the data gen function 
  singlerun <- gendata_IIW(n, beta1, beta2, beta3, beta4, beta5, gamma1, gamma2, gamma3, gamma4, gamma5, gamma6, tau)
  simdata_full <- singlerun$simdata_full #this is the censored data, but shows data at each 0.01 increment (times at which X and Y are actualy
  #observed are indicated with observed = 1 in the dataset)
  simdata_observed <- singlerun$simdata_Yobserved
  numevents <- singlerun$numevents
  newn <- singlerun$newn
  
  
  # need a function that imputes the missingness mechanism (MCAR< MAR,  etc, and a proportion of missingness. It also 
  # applies all 8 imputation methods)
  
  betamat <- data.frame(matrix(NA, nrow = 1, ncol = 5*((length(schemes)-3)*length(proportions)+1) + 2))
  schemeno = 1
  schemenames <- c()
  
  for(scheme in schemes){
    
    if(scheme == "Naive"){
      
      #No IIW, no adjustments
      
      schemenames <- c(schemenames, scheme)
 
      beta1_naive <- summary(glm(y ~ D + offset(2-times) - 1, data=simdata_observed))$coef[1,1]
      
      betamat[1, schemeno] <- beta1_naive
      schemeno = schemeno + 1
      
    }else if(scheme == "NoMissingness"){
      
      #perform IIW, but not the adjustments for missingness since there is no missing data
      
      schemenames <- c(schemenames, scheme)
      
      intensitydata <- makeMissing(simdata_full, scheme = scheme, p = p)
      intensitydata <- intensitydata[ , c(1:8, 19)] #intensity data only includes id, times, Z, Yobserved indicator
      
      #include a variable counting observation number, to be used for lagging time for Surv function
      intensitydata$obsnumber <- with(intensitydata, ave(id, id, FUN = seq_along))
      
      
      # #create lagged time variable
      intensitydata$times.lag <- intensitydata$times[c(nrow(intensitydata ),1:(nrow(intensitydata )-1))]
      intensitydata$times.lag[intensitydata$obsnumber == 1] <- -0.01
      
      
      gamma.hat <- coxph(Surv(times.lag, times, Yobserved) ~ D  + V1t + V2t + V3t + V4t + V5t - 1, data = intensitydata)$coef
      delta.hat <- coxph(Surv(times.lag, times, Yobserved) ~ D - 1, data = intensitydata)$coef
      
      iiw <- exp(cbind(simdata_observed$D)%*%delta.hat)/
        exp(cbind(simdata_observed$D, simdata_observed$V1t, simdata_observed$V2t, simdata_observed$V3t, simdata_observed$V4t, 
                  simdata_observed$V5t)%*%gamma.hat)
      
      beta1_nm <- summary(glm(y ~ D + offset(2-times) - 1, data=simdata_observed, weights = iiw))$coef[1,1]
      
      betamat[1, schemeno] <- beta1_nm
      
      schemeno = schemeno + 1
      
    }else{
      
      
    for(p in proportions){
      
        betavec <- c() #stores betas for each of the 8 methods for a single scheme/proportion
        

        #invoke missingness
        intensitydata <- makeMissing(simdata_full, scheme = scheme, p = p)
        intensitydata <- intensitydata[, c(1:8, 19)] #intensity data only includes id, times, Z,  y (for MNAR), Yobserved indicator
        
        
        #### impute data using the different methods
        betavec <- c(betavec, CCA(intensitydata, simdata_observed), 
                     LOCF(intensitydata, simdata_observed), SI(intensitydata, simdata_observed),
                     MI(intensitydata, simdata_observed), NIRF(intensitydata, simdata_observed))
        
        
        if(scheme == "ObsTimesOnly"){
          schemenames <- c(schemenames, c(paste(scheme, "CCA", sep = "_"), 
                          paste(scheme, "LOCF", sep = "_"), paste(scheme, "SI", sep = "_"),
                          paste(scheme, "MI", sep = "_"), paste(scheme, "NIRF", sep = "_")))
        }else{
          schemenames <- c(schemenames, c(paste(scheme, p, "CCA", sep = "_"), 
                          paste(scheme, p, "LOCF", sep = "_"), paste(scheme, p, "SI", sep = "_"),
                          paste(scheme, p, "MI", sep = "_"), paste(scheme, p, "NIRF", sep = "_")))
        }
        
        
        
        betamat[1, schemeno:(schemeno+4)] <- betavec
        
        schemeno = schemeno + 5
        
        if(scheme == "ObsTimesOnly"){
          break
        }
      }
    }
  }
  

  
  colnames(betamat) <- schemenames
  
  out <- list(betamat, schemenames)
  names(out) <- c("betamat", "schemenames")
  
  
  return(out)
}



simulateResultsIIW<-  function(N, n, beta1, beta2, beta3, beta4, beta5, gamma1, gamma2, gamma3, 
                               gamma4, gamma5, gamma6, tau, schemes, proportions,
                               outputfulldatalist = FALSE, inParallel, nclusters = NULL){

  # Simulates N instances of the each scheme
  # computes in parallel if inParallel = T
  #

  nschemes <- 4*((length(schemes)-3)*length(proportions)+1) + 2


  results_beta1 <- matrix(data = NA, nrow = N, ncol = nschemes)

    for(i in 1:N){
      if(i%%1 == 0){print(i)}
      simrun <- simulateOneIIW(n, beta1, beta2, beta3, beta4, beta5, gamma1, gamma2, gamma3, gamma4,
                               gamma5, gamma6, tau, schemes, proportions)
      results_beta1[i, ] <- as.matrix(simrun$betamat)
    }
  



  colnames(results_beta1) <- simrun$schemenames

  bias_beta1 <- round(apply(results_beta1, FUN =  mean, MARGIN = 2)  - beta1 ,3)
  names(bias_beta1) <- simrun$schemenames

  var_beta1 <- round(apply(results_beta1, FUN = var, MARGIN = 2), 3)
  names(var_beta1) <- simrun$schemenames

  mse_beta1 <- round(apply((results_beta1 - beta1)^2, FUN = mean, MARGIN = 2), 3)
  names(mse_beta1) <- simrun$schemenames


# 
#   ## ## ##
# 
#   naive_beta1 <- c(bias_beta1[1], var_beta1[1], mse_beta1[1])
#   names(naive_beta1) <- c("Bias", "Var", "MSE")
# 
#   IIW_beta1_full <- c(bias_beta1[2], var_beta1[2], mse_beta1[2])
#   names(IIW_beta1_full) <- c("Bias", "Var", "MSE")
# 
#   IIW_beta1_obstimes <- c(bias_beta1[3], var_beta1[3], mse_beta1[3])
#   names(IIW_beta1_obstimes) <- c("Bias", "Var", "MSE")
# 
#   if(outputfulldatalist == TRUE){
#     out <- list(bias_beta1, var_beta1, mse_beta1,
#                 naive_beta1, IIW_beta1_full, IIW_beta1_obstimes, fullresults_beta1)
# 
#     names(out) <- c('bias_beta1', 'var_beta1', 'mse_beta1',
#                     'naive_beta1',  "IIW_beta1_full", "IIW_beta1_obstimes", "fullresults_beta1")
#   }else{
#     out <- list(bias_beta1, var_beta1, mse_beta1,
#                 naive_beta1,  IIW_beta1_full, IIW_beta1_obstimes)
# 
#     names(out) <- c('bias_beta1', 'var_beta1', 'mse_beta1',
#                     'naive_beta1',  "IIW_beta1_full", "IIW_beta1_obstimes")
#   }
# 
#   return(out)

}
# 
# 
# 
# simulateALLFIPTIW_bd <- function(N, n, beta1, beta2, beta3, gamma1, gamma2vec, gamma3vec, 
#                                  alpha0, alpha1vec, tau, outputfulldatalist = FALSE, censinform = F, 
#                                  eta1 = NULL, eta2 = NULL, eta3 = NULL, 
#                                  inParallel = T, nclusters = NULL){
#   #N: number of simulation runs
#   #n: vector of sample sizes
#   #beta1: coefficient for Xi(t) in logistic outcome model
#   #beta2: vector of coefficients to consider for outcome generation model
#   #gamma1, gamma2 parameters for intensity for Xi(t) and Zi, respectively
#   #tau: maximum follow-up time
#   #inParallel: runs in parallel with nclusters if true
#   
#   #This function aggregates simulation results for varying n and beta2
#   
#   
#   
#   resultsmat <- matrix(NA, nrow = 1, ncol = 15)
#   fulldatalist <- list()
#   fulldatalistnames <- c()
#   
#   i = 1
#   
#   
#   for(gamma2 in gamma2vec){
#     for(gamma3 in gamma3vec){
#       for(alpha1 in alpha1vec){
#         
#         
#         print(paste("Now on gamma2 = ", gamma2, ", gamma3 =", gamma3, ", alpha1 =  ", alpha1, sep = ""))
#         result <- simulateResultsFIPTIW_bd(N, n, beta1, beta2, beta3, gamma1, gamma2, gamma3, 
#                                            alpha0, alpha1, tau,timevarD, outputfulldatalist, inParallel, nclusters)
#         
#         resultsmat <- rbind(resultsmat, c(gamma2,  gamma3, alpha1,
#                                           result$naive_beta1, result$IIW_beta1, 
#                                           result$IPW_beta1, result$FIPTIW_beta1))
#         
#         # resultsmat_stabilized <- rbind(resultsmat_stabilized, c(gamma2,  gamma3, alpha1,
#         #                                                         result$naive_beta1, result$IIWstab_beta1, 
#         #                                                         result$IPWstab_beta1, result$FIPTIWstab_beta1))
#         # 
#         
#         listname <- paste("fulldata_gamma2_",gamma2, "_gamma3_", gamma3, "_alpha1_", alpha1, sep = "")
#         
#         if(outputfulldatalist == TRUE){
#           fulldatalist[[i]] <- result$fullresults_beta1
#           fulldatalistnames <- c(fulldatalistnames, listname)
#         }
#         
#         
#         i <- i + 1
#         
#       }
#     }
#   }
#   
#   if(outputfulldatalist == TRUE){
#     names(fulldatalist) <- fulldatalistnames
#   }
#   
#   
#   
#   resultsmat<- resultsmat[-1,]
#   # resultsmat_stabilized <- resultsmat_stabilized[-1,]
#   
#   colnames(resultsmat) <- c( "gamma2",  "gamma3", "alpha1", "Bias", "Var" , "MSE" , "Bias", "Var" , "MSE" , "Bias", "Var" , "MSE", "Bias", "Var" , "MSE" )
#   # colnames(resultsmat_stabilized) <- c( "gamma2", "gamma3", "alpha1", "Bias", "Var" , "MSE" , "Bias", "Var" , "MSE" , "Bias", "Var" , "MSE", "Bias", "Var" , "MSE" )
#   
#   if(outputfulldatalist == TRUE){
#     out <- list(resultsmat, fulldatalist)
#     names(out) <- c("resultsmat", "fulldatalist")
#   }else{
#     out <- list(resultsmat)
#     names(out) <- c("resultsmat")
#   }
#   
#   return(out)
# }

