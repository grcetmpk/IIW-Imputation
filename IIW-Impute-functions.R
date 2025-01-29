require(survival)
require(geepack)
require(nleqslv)
require(knitr)
require(kableExtra)
require(dplyr)
require(doParallel)
library(Rcpp)
require(mice)
require(missForest)
require(geepack)
require(splines)
require(NNMIS)

expit <- function(x){ return(exp(x)/(1+exp(x)))}



##### This script contains the functions necessary for data generation and analysis/imputation
##### for the third project in Grace's thesis.


### Simulate irreg longitudinal data using bernoulli draws, discretizing the observation times
### from 0 to \tau, by 0.01. First, starting with a non-time-varying treatment



gendata_IIW <- function(n, beta1, beta2, beta3, beta4, beta5, beta6, gamma1, 
                        gamma2, gamma3, gamma4, gamma5, gamma6, tau, outcome){
  
  # Simulates observation times using bernoulli draws with probabilities proportional 
  # to the intensity. assume the treatment is time-invariant and randomized
  
  # SCHEME A
  #
  # Simulates continuous outcome
  #
  # n: number of subjects
  # tau: maximum follow up
  
  
  mu_V1_D0 <- 1
  var_V1_D0 <- 1 #var_V1_D0 <- 0.25 
  mu_V1_D1 <- 3
  var_V1_D1 <- 2
  
  var_phi <- 0.25
  var_epsilon <- 1
  
  
  # generate subjects one by one and then combine
  id <- 1
  simdatafull <- data.frame(matrix(NA, ncol = 21))
  colnames(simdatafull) <- c("id", "times", "D", "V1t", "V2t", "V3t", "V4t", "V5t", "censortime", "eta", "cexp_V1_D", "cexp_V2_D", 
                             "cexp_V3_D", "cexp_V4_D", "cexp_V5_D", "cvar_V1_D", "cvar_V2_D", 
                             "cvar_V3_D", "cvar_V4_D", "cvar_V5_D", "y")
  
  
  disctimes <- seq(0, tau, by = 0.01)
  
  while(id <= n){
    
    ## Generate Covariates
    
    ## Generate baseline covariate
    
    # generate treatment assignment (time-invariant) at each time point
    D <- rep(rbinom(1, 1, 0.5), length(disctimes))
    
    

    # generate observation times confounder (time varying) V1t ~ N( mu_V_D0, var_V_D0) if D(t) = 0 and 
    # N( mu_V_D1, var_V_D1) if D = 1
      
    V1t <- round(ifelse(D == 0, rnorm(length(disctimes), mu_V1_D0, sqrt(var_V1_D0)), 
                          rnorm(length(disctimes), mu_V1_D1, sqrt(var_V1_D1))), 3)
      

    
    
    # generate V2 (which will be used in V2(t) = V2log(t)
    
    V2t <- round(ifelse(D == 0, runif(length(disctimes), 0, 1), runif(length(disctimes), 1, 2)), 3)

    
    # generate V3(t) which is binary
    # V3 <- ifelse(D[1] == 0, rbinom(1, 1, 0.2), rbinom(1, 1, 0.7))
    # V3t <- rep(V3, length(disctimes))
    
    V3t <- ifelse(D == 0, rbinom(length(disctimes), 1, 0.2), rbinom(length(disctimes), 1, 0.7)) #time varying version
    
    
    # generate V4(t)
    
    V4t <- round(ifelse(D == 0, runif(length(disctimes), 0.5, 2)*(disctimes-1)^2, 
                  runif(length(disctimes), 1.5, 4)*(disctimes-1)^2), 3)
    
    
    # generate V5, which will be used in V5(t)
    
    V5 <- ifelse(D == 0, rnorm(length(disctimes), 0, 1), rnorm(length(disctimes), 2, sqrt(0.5)))
    V5t <- round(V5*(-(disctimes)/2),3)
    
    
    # generate random effect
    phi <- rnorm(1, 0, var_phi)
    
    
    # simulate censoring time (uniform)
    censortime <- runif(1, tau/2, tau)
    
    # calculate eta
    eta <- rgamma(1, shape = 100, scale = 0.01) #eta is for random effect with mean 1 and sd 0.1 
    #(obs times within subject correlated if sd !0)
    
    simdata <- data.frame(rep(id, length(disctimes)), disctimes, D, V1t, V2t, V3t, V4t, V5t, 
                          rep(censortime, length(disctimes)), rep(eta, length(disctimes)))
    colnames(simdata) <- c("id", "times", "D", "V1t", "V2t", "V3t", "V4t", "V5t", "censortime", "eta")
    
    # need conditional expectation and variance of Z | X
    

    

    
    simdata$cexp_V1_D <- ifelse(simdata$D == 0, mu_V1_D0, mu_V1_D1)
    simdata$cexp_V2_D <- ifelse(simdata$D == 0, 0.5, 1.5)
    simdata$cexp_V3_D <- ifelse(simdata$D == 0, 0.2, 0.7)
    simdata$cexp_V4_D <- ifelse(simdata$D == 0, 2.5/2*(simdata$times - 1)^2, 5.5/2*(simdata$times - 1)^2)
    simdata$cexp_V5_D <- ifelse(simdata$D == 0, 0, 2)*(-simdata$times/2)

    simdata$cvar_V1_D <- ifelse(simdata$D == 0, var_V1_D0, var_V1_D1)
    simdata$cvar_V2_D <- 1/12
    simdata$cvar_V3_D <- ifelse(simdata$D == 0, 0.2*0.8, 0.2*0.7)
    simdata$cvar_V4_D <- ifelse(simdata$D == 0, 1/12*(2 - 0.5)^2*((simdata$times - 1)^2)^2, 1/12*(4 - 1.5)^2*((simdata$times -1)^2)^2)
    simdata$cvar_V5_D <- ifelse(simdata$D == 0, (-simdata$times/2)^2, 0.5^2*(-simdata$times/2)^2)

    if(outcome == "binary"){
      M <- sqrt(beta2^2*simdata$cvar_V1_D + beta3^2*simdata$cvar_V2_D + beta4^2*simdata$cvar_V3_D + 
                  beta5^2*simdata$cvar_V4_D + beta6^2*simdata$cvar_V5_D + var_epsilon + var_phi)/1.7
      
      fstar <- (2-simdata$times)*M - beta2*simdata$cexp_V1_D - beta3*simdata$cexp_V2_D - 
        beta4*simdata$cexp_V3_D - beta5*simdata$cexp_V4_D - beta6*simdata$cexp_V5_D
      
      
      
      # generate time-varying outcome at each possible time
      simdata$y <- ifelse(fstar + beta1*M*simdata$D + beta2*simdata$V1t + beta3*simdata$V2t + 
                            beta4*simdata$V3t + beta5*simdata$V4t + beta6*simdata$V5t + rep(phi, dim(simdata)[1]) + 
                            rnorm(dim(simdata)[1], 0, sqrt(var_epsilon)) > 0, 1, 0)
      
    }else if (outcome == "continuous"){
      
      simdata$y <- round((2 - simdata$times) + beta1*simdata$D + beta2*(simdata$V1 - simdata$cexp_V1_D) +
        beta3*(simdata$V2 - simdata$cexp_V2_D) + beta4*(simdata$V3 - simdata$cexp_V3_D) + beta5*(simdata$V4 - simdata$cexp_V4_D)+
        + beta6*(simdata$V5 - simdata$cexp_V5_D)+ rep(phi, dim(simdata)[1]) + 
        rnorm(dim(simdata)[1], 0, sqrt(var_epsilon)),3)
      
      
    }else{
      return(print(paste("incorrect type of outcome providedd")))
    }
    
    
    simdatafull <- rbind(simdatafull, simdata)
    
    
    
    
    id = id + 1
    
  }
  
  
  
  
  simdata <- simdatafull[-1,] #remove empty first row
  
  ### calculate intensities for each counterfactual observation time:
  
  
  simdata$intensities <- simdata$eta*sqrt(simdata$times)/2*exp(gamma1*simdata$D + gamma2*simdata$V1t + gamma3*simdata$V2t +
                                                                 gamma4*simdata$V3t + gamma5*simdata$V4t + gamma6*simdata$V5t)
  simdata$prObs <- ifelse(0.05*simdata$intensities > 1, 1, 0.01*simdata$intensities) #simdata$intensities/max(simdata$intensities)
  simdata$Yobserved <- rbinom(nrow(simdata), 1, simdata$prObs)
  
  
  # simdata <- simdata %>%
  #   filter(times < censortime)

  
  ### filter censored observations, unobserved observations (do this before or after calculating maximum intensity???)
  simdata_observed <- simdata %>%
    filter(Yobserved == 1) %>%
    select(id, times, V1t, V2t, V3t, V4t, V5t, D, y, Yobserved)
  
  
  #check the generated data
  numevents <- summary(tapply(simdata_observed$Yobserved, simdata_observed$id, sum)) 
  
  
  # #also get data for individuals at baseline (all same here)
  newn <- length(unique(simdata_observed$id)) #number of people after censoring etc
  
  out <- list(simdata, simdata_observed, numevents, newn)
  names(out) <- c("simdata_full", "simdata_Yobserved", "numevents", "newn")
  
  return(out) 
  
}




gendata_IIWB <- function(n, beta1, beta2, beta3, beta4, beta5, beta6, gamma1, 
                        gamma2, gamma3, gamma4, gamma5, gamma6, tau, outcome){
  
  # Simulates observation times using bernoulli draws with probabilities proportional 
  # to the intensity. assume the treatment is time-invariant and randomized
  #
  # Simulates continuous outcome
  #
  # n: number of subjects
  # tau: maximum follow up
  
  
  var_phi <- 0.25
  var_epsilon <- 1
  
  
  # generate subjects one by one and then combine
  id <- 1
  simdatafull <- data.frame(matrix(NA, ncol = 21))
  colnames(simdatafull) <- c("id", "times", "D", "V1t", "V2t", "V3t", "V4t", "V5t", "censortime", "eta", "cexp_V1_D", "cexp_V2_D", 
                             "cexp_V3_D", "cexp_V4_D", "cexp_V5_D", "cvar_V1_D", "cvar_V2_D", 
                             "cvar_V3_D", "cvar_V4_D", "cvar_V5_D", "y")
  
  
  disctimes <- seq(0, tau, by = 0.01)
  
  while(id <= n){
    
    ## Generate Covariates
    
    ## Generate baseline covariate
    
    # generate treatment assignment (time-invariant) at each time point
    D <- rep(rbinom(1, 1, 0.5), length(disctimes))
    
    

    V1t <- round(ifelse(D == 0, rnorm(1, 0, sqrt(0.5)), 
                          rnorm(1, 2, sqrt(0.5))), 3)*disctimes/3
 
    
    
    # generate V2(t) (time varying)
    
    V2t <- round(ifelse(D == 0, runif(length(disctimes), 0, 1), runif(length(disctimes), 1, 2)), 3)*log(disctimes + 1)/2
    
    
    # generate V3(t) , timeinvar
    
    V3t <- ifelse(D == 0, rbinom(length(disctimes), 1, 0.2), rbinom(length(disctimes), 1, 0.7))

    
    # generate V4(t), time-invariant
    
    V4t <- round(ifelse(D == 0, runif(1, 0.5, 2), 
                        runif(1, 1.5, 3)), 3)*(disctimes-1)/3
    
    
    # generate V5, time-var
    
    V5t <- ifelse(D == 0, rbinom(length(disctimes), 1, 0.8), rbinom(length(disctimes), 1, 0.4))
    
    # generate random effect
    phi <- rnorm(1, 0, var_phi)
    
    censortime <- runif(1, tau/2, tau)
  
    
    # calculate eta
    eta <- rgamma(1, shape = 100, scale = 0.01) #eta is for random effect with mean 1 and sd 0.1 
    #(obs times within subject correlated if sd !0)
    
    simdata <- data.frame(rep(id, length(disctimes)), disctimes, D, V1t, V2t, V3t, V4t, V5t, 
                          rep(censortime, length(disctimes)), rep(eta, length(disctimes)))
    colnames(simdata) <- c("id", "times", "D", "V1t", "V2t", "V3t", "V4t", "V5t", "censortime", "eta")
    
    # need conditional expectation and variance of Z | X

    
    simdata$cexp_V1_D <- ifelse(simdata$D == 0, 0, 2)*simdata$times/3
    simdata$cexp_V2_D <- ifelse(simdata$D == 0, 0.5, 1.5)*log(simdata$times + 1)/2
    simdata$cexp_V3_D <- ifelse(simdata$D == 0, 0.2, 0.7)
    simdata$cexp_V4_D <- ifelse(simdata$D == 0, 2.5/2, 5.5/2)*(simdata$times-1)/3
    simdata$cexp_V5_D <- ifelse(simdata$D == 0, 0.8, 0.4)
    
    simdata$cvar_V1_D <- 0.5*(simdata$times/3)^2
    simdata$cvar_V2_D <- 1/12*(log(simdata$times + 1)/2)^2
    simdata$cvar_V3_D <- ifelse(simdata$D == 0, 0.2*0.8, 0.2*0.7)
    simdata$cvar_V4_D <- ifelse(simdata$D == 0, 1/12*(2 - 0.5)^2*((simdata$times-1)/3)^2, 1/12*(4 - 1.5)^2)*((simdata$times-1)/3)^2
    simdata$cvar_V5_D <- ifelse(simdata$D == 0, 0.8*0.2, 0.4*0.6)
    
    if(outcome == "binary"){
      M <- sqrt(beta2^2*simdata$cvar_V1_D + beta3^2*simdata$cvar_V2_D + beta4^2*simdata$cvar_V3_D + 
                  beta5^2*simdata$cvar_V4_D + beta6^2*simdata$cvar_V5_D + var_epsilon + var_phi)/1.7
      
      fstar <- (2-simdata$times)*M - beta2*simdata$cexp_V1_D - beta3*simdata$cexp_V2_D - 
        beta4*simdata$cexp_V3_D - beta5*simdata$cexp_V4_D - beta6*simdata$cexp_V5_D
      
      
      
      # generate time-varying outcome at each possible time
      simdata$y <- ifelse(fstar + beta1*M*simdata$D + beta2*simdata$V1t + beta3*simdata$V2t + 
                            beta4*simdata$V3t + beta5*simdata$V4t + beta6*simdata$V5t + rep(phi, dim(simdata)[1]) + 
                            rnorm(dim(simdata)[1], 0, sqrt(var_epsilon)) > 0, 1, 0)
      
    }else if (outcome == "continuous"){
      
      simdata$y <- round((2 - simdata$times) + beta1*simdata$D + beta2*(simdata$V1 - simdata$cexp_V1_D) +
                           beta3*(simdata$V2 - simdata$cexp_V2_D) + beta4*(simdata$V3 - simdata$cexp_V3_D) + beta5*(simdata$V4 - simdata$cexp_V4_D)+
                           + beta6*(simdata$V5 - simdata$cexp_V5_D)+ rep(phi, dim(simdata)[1]) + 
                           rnorm(dim(simdata)[1], 0, sqrt(var_epsilon)),3)
      
      
    }else{
      return(print(paste("incorrect type of outcome providedd")))
    }
    
    
    simdatafull <- rbind(simdatafull, simdata)
    
    
    
    
    id = id + 1
    
  }
  
  
  
  
  simdata <- simdatafull[-1,] #remove empty first row
  
  ### calculate intensities for each counterfactual observation time:
  
  
  simdata$intensities <- simdata$eta*sqrt(simdata$times)/2*exp(gamma1*simdata$D + gamma2*simdata$V1t + gamma3*simdata$V2t +
                                                                 gamma4*simdata$V3t + gamma5*simdata$V4t + gamma6*simdata$V5t)
  simdata$prObs <- ifelse(0.01*simdata$intensities > 1, 1, 0.01*simdata$intensities) #simdata$intensities/max(simdata$intensities)
  simdata$Yobserved <- rbinom(nrow(simdata), 1, simdata$prObs)
  

  # simdata <- simdata %>%
  #   filter(times < censortime)
  
  
  ### filter censored observations, unobserved observations (do this before or after calculating maximum intensity???)
  simdata_observed <- simdata %>%
    filter(Yobserved == 1) %>%
    select(id, times, V1t, V2t, V3t, V4t, V5t, D, y, Yobserved)
  
  
  #check the generated data
  numevents <- summary(tapply(simdata_observed$Yobserved, simdata_observed$id, sum)) 
  
  
  # #also get data for individuals at baseline (all same here)
  newn <- length(unique(simdata_observed$id)) #number of people after censoring etc
  
  out <- list(simdata, simdata_observed, numevents, newn)
  names(out) <- c("simdata_full", "simdata_Yobserved", "numevents", "newn")
  
  return(out) 
  
}









makeMissing <- function(intensitydat, missingnesstype, p, whichmissing){
  
  
  # We are going to assume one additional covariate (V5t) is always known to impose MAR
  

  if(missingnesstype == "NoMissingness"){
    
    intensitydat$Zobserved <- 1
    intensitydat <- intensitydat %>% select(id, times, D, V1t, V2t, V3t, V4t, V5t, y, Zobserved, Yobserved)
    
    return(intensitydat) # return full dataset with no missingness
    
    }else if(missingnesstype %in% c("ObsTimesOnly", "Naive")){
      
      #special case of MNAR, any V that can be misisng is only known when Y is 
      
      intensitydat$Zobserved <-  intensitydat$Yobserved # Z only measured if it Y is observed
      
      if("V1t" %in% whichmissing){
        intensitydat$V1t[intensitydat$Zobserved == 0] <- NA
      }
      
      if("V2t" %in% whichmissing){
        intensitydat$V2t[intensitydat$Zobserved == 0] <- NA
      }
      
      if("V3t" %in% whichmissing){
        intensitydat$V3t[intensitydat$Zobserved == 0] <- NA
      }
      
      if("V4t" %in% whichmissing){
        intensitydat$V4t[intensitydat$Zobserved == 0] <- NA
      }
      
      
      intensitydat <- intensitydat %>% select(id, times, D, V1t, V2t, V3t, V4t, V5t, Zobserved, Yobserved)
      
      
      return(intensitydat)
      
    }else if(missingnesstype == "MCAR"){
      # randomly select p% of observations to be missing

      intensitydat$Zobserved <- rbinom(dim(intensitydat)[1], 1, 1 - p) 
      
      # change any data with missing Z to NAs
      if("V1t" %in% whichmissing){
        intensitydat$V1t[intensitydat$Zobserved == 0] <- NA
      }
      
      if("V2t" %in% whichmissing){
        intensitydat$V2t[intensitydat$Zobserved == 0] <- NA
      }
      
      if("V3t" %in% whichmissing){
        intensitydat$V3t[intensitydat$Zobserved == 0] <- NA
      }
      
      if("V4t" %in% whichmissing){
        intensitydat$V4t[intensitydat$Zobserved == 0] <- NA
      }
      
      intensitydat <- intensitydat %>% select(id, times, D, V1t, V2t, V3t, V4t, V5t, Zobserved, Yobserved)
      
      
      return(intensitydat)
      
    }else if(missingnesstype == "MAR"){
      
      #missingness dependent on fully observed D and V5t
      
      x <- as.matrix(intensitydat[, c(3, 8)]) #matrix form of dataset for covariate in MAR model (D)
      eta <- c(1, 2) #set as anything for now
      
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
      
      pmiss <- expit(eta0 + eta[1]*intensitydat$D + eta[2]*intensitydat$V5t)
      intensitydat$Zobserved <- rbinom(length(pmiss), 1, 1 - pmiss) #bernoulli draws to see which observed
      #sum(intensitydat$Zobserved)/length(intensitydat$Zobserved)
      
      # change any data with missing = 1 to NAs
      if("V1t" %in% whichmissing){
        intensitydat$V1t[intensitydat$Zobserved == 0] <- NA
      }
      
      if("V2t" %in% whichmissing){
        intensitydat$V2t[intensitydat$Zobserved == 0] <- NA
      }
      
      if("V3t" %in% whichmissing){
        intensitydat$V3t[intensitydat$Zobserved == 0] <- NA
      }
      
      if("V4t" %in% whichmissing){
        intensitydat$V4t[intensitydat$Zobserved == 0] <- NA
      }
      
      intensitydat <- intensitydat %>% select(id, times, D, V1t, V2t, V3t, V4t, V5t, Zobserved, Yobserved)

      
      return(intensitydat)
    
    }else if(missingnesstype == "MNAR"){
      
      #missingness dependent on V4t, which is not fully observed
      
      x <- as.matrix(intensitydat[, c(6, 7)])#, 21)])#c(3, 7, 8, 21)]) #matrix form of dataset for covariate in MAR model (D)
      eta <- c(1, 2)#, 3) #set as anything for now
      
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
      
      #pmiss <- expit(eta0 + eta[1]*intensitydat$D + eta[2]*intensitydat$V4t + eta[3]*intensitydat$V5t+eta[4]*intensitydat$y)
      pmiss <- expit(eta0 + eta[1]*intensitydat$V3t + eta[2]*intensitydat$V4t) # + eta[2]*intensitydat$y)
      intensitydat$Zobserved <- rbinom(length(pmiss), 1, 1 - pmiss) #bernoulli draws to see which observed
      #sum(intensitydat$Zobserved)/length(intensitydat$Zobserved)
      
      # change any data with missing = 1 to NAs
      if("V1t" %in% whichmissing){
        intensitydat$V1t[intensitydat$Zobserved == 0] <- NA
      }
      
      if("V2t" %in% whichmissing){
        intensitydat$V2t[intensitydat$Zobserved == 0] <- NA
      }
      
      if("V3t" %in% whichmissing){
        intensitydat$V3t[intensitydat$Zobserved == 0] <- NA
      }
      
      if("V4t" %in% whichmissing){
        intensitydat$V4t[intensitydat$Zobserved == 0] <- NA
      }
      
      intensitydat <- intensitydat %>% select(id, times, D, V1t, V2t, V3t, V4t, V5t, Zobserved, Yobserved)

      
      return(intensitydat)
      
    }else{
      
      return(NA)
    }
  }


IIW <- function(intensitydat, obsdat, outcome, usesplines){
  
  
  
  ##### perform IIW on the imputed data set
    
  # terti<-quantile(0:tau , c(0.3333, 0.66666), type = 1) 
  terti<-quantile(0:tau , c(0.5), type = 1) 
  
  #include a variable counting observation number, to be used for lagging time for Surv function
  intensitydat$obsnumber <- with(intensitydat, ave(id, id, FUN = seq_along))
  
  # #create lagged time variable
  intensitydat$times.lag <- intensitydat$times[c(nrow(intensitydat ),1:(nrow(intensitydat )-1))]
  intensitydat$times.lag[intensitydat$obsnumber == 1] <- -0.01
  
  gamma.hat <- coxph(Surv(times.lag, times, Yobserved) ~ D  + V1t + V2t + V3t + V4t - 1, data = intensitydat)$coef
  delta.hat <- coxph(Surv(times.lag, times, Yobserved) ~ D - 1, data = intensitydat)$coef
  
  ### perform IIW on the "observed" X and Y data
  

  iiw <- exp(cbind(obsdat$D)%*%delta.hat)/
    exp(cbind(obsdat$D, obsdat$V1t, obsdat$V2t, 
              obsdat$V3t, obsdat$V4t)%*%gamma.hat)
  
  
  beta1_results <- fitOutcomeModel(obsdat, outcome, usesplines, terti, iiwweights = iiw)
  
 
  
  return(beta1_results)
  
  
}


CCA <- function(intensitydat, obsdat, outcome, usesplines){
  
  #intensitydat -> data used for IIW which may or may not have missingness (indicated by Zmissing)
  #obsdat -> data on longitudinal outcome, only observed at "observation times" with no other missingness, used to fit the outcome model
  
  #complete case analysis
      
  intensitydat <- na.omit(intensitydat)
  out <- IIW(intensitydat, obsdat, outcome, usesplines)
  return(out)
}


LOCF <- function(intensitydat, obsdat, outcome, usesplines){

    #do LOCF
    imputedintensitydata <- intensitydat %>%
      group_by(id) %>%
      tidyr::fill(names(intensitydat), .direction = "downup") %>%
      ungroup()
    
    out <- IIW(imputedintensitydata, obsdat, outcome, usesplines)
    return(out)
    
  }
  
  
SI <- function(intensitydat, obsdat, outcome, usesplines, missingnesstype){
 #regression switching imputation with PPM
  
  # only use observed Z covariates (D, V1t - v5t) for imputation
  
  
  intensitydat2 <- intensitydat
  intensitydat2$D <- as.factor(intensitydat2$D)
  intensitydat2$V3t <- as.factor(intensitydat2$V3t)
  
  imp <- mice::mice(intensitydat2, m = 1, maxit = 5, meth = "pmm", seed = 100, printFlag = F)
  imputedintensitydat <- complete(imp)
  
  intensitydat_imputed_full <- intensitydat #copy original df with missingness
  
  #replace columns with imputed data
  intensitydat_imputed_full$V1t <- imputedintensitydat$V1t
  intensitydat_imputed_full$V2t <- imputedintensitydat$V2t
  intensitydat_imputed_full$V3t <- imputedintensitydat$V3t
  intensitydat_imputed_full$V4t <- imputedintensitydat$V4t
  intensitydat_imputed_full$V5t <- imputedintensitydat$V5t
  
  out <- IIW(intensitydat_imputed_full, obsdat, outcome, usesplines)
  
  return(out)
}


MI <- function(intensitydat, obsdat, nimputations, outcome, usesplines){
  #regression switching imputation with PPM
  

  
  
  if(outcome == "binary"){
    
    familymod = binomial(link = "logit")
    
  }else if(outcome == "continuous"){
    
    familymod = gaussian(link = "identity")
  }else{
    return(print("insufficient outcome type provided"))
  }
  
  intensitydat$D <- as.factor(intensitydat$D)
  intensitydat$V3t <- as.factor(intensitydat$V3t)
  
  # terti<-quantile(0:tau , c(0.3333, 0.66666), type = 1) 
  terti<-quantile(0:tau , c(0.5), type = 1) 

  #include a variable counting observation number, to be used for lagging time for Surv function
  intensitydat$obsnumber <- with(intensitydat, ave(id, id, FUN = seq_along))

  # #create lagged time variable
  intensitydat$times.lag <- intensitydat$times[c(nrow(intensitydat ),1:(nrow(intensitydat )-1))]
  intensitydat$times.lag[intensitydat$obsnumber == 1] <- -0.01
  
  
  
  imp <- mice::mice(intensitydat, m = nimputations, maxit = 5, meth = "pmm", seed = 100, printFlag = F)


  ##### perform IIW on the imputed data set
  
  gamma.hat <- summary(pool(with(imp, coxph(Surv(times.lag, times, Yobserved) ~ D  + V1t + V2t + V3t + V4t - 1))))[,2]
  delta.hat <- summary(pool(with(imp, coxph(Surv(times.lag, times, Yobserved) ~ D - 1))))[,2]
  
  ### perform IIW on the "observed" X and Y data
  
  iiw <- exp(cbind(obsdat$D)%*%delta.hat)/
    exp(cbind(obsdat$D, obsdat$V1t, obsdat$V2t, 
              obsdat$V3t, obsdat$V4t)%*%gamma.hat)
  
  beta1_results <- fitOutcomeModel(obsdat, outcome, usesplines, terti, iiwweights = iiw)

  return(beta1_results)
}
  
  
NIRF <- function(intensitydat, obsdat, outcome, usesplines){

  
  intensitydat2 <- intensitydat
  intensitydat2$D <- as.factor(intensitydat2$D)
  intensitydat2$V3t <- as.factor(intensitydat2$V3t)
  intensitydat2$Yobserved <- as.factor(intensitydat2$Yobserved)
  # 

  # registerDoParallel(cores= ncoresavailable-2)
  imp <- missForest(intensitydat2, maxiter = 5, ntree = 15, parallelize = "no") 

  imputedintensitydata <- imp$ximp
  
  imputedintensitydata$D <- as.numeric(as.character(imputedintensitydata$D))
  imputedintensitydata$V3t <- as.numeric(as.character(imputedintensitydata$V3t))
  imputedintensitydata$Yobserved <- as.numeric(as.character(intensitydat$Yobserved))

  out <- IIW(imputedintensitydata, obsdat, outcome, usesplines)
  return(out)
  
}  



fitOutcomeModel <- function(obsdat, outcome, usesplines, terti, iiwweights = NULL){
  
  #determine which family to use for geeglm based on outcome type
  
  if(outcome == "binary"){
    
    familymod = binomial(link = "logit")
    
  }else if(outcome == "continuous"){
    
    familymod = gaussian(link = "identity")
  }else{
    return(print("insufficient outcome type provided"))
  }
  
  
  # fit model based on spline usage
  

  if(usesplines == T){
    beta1mod <- summary(geeglm(y ~ D + bs(obsdat$times, degree=3,knots=c(terti)), 
                               family = familymod, data=obsdat, id = id,
                               weights = iiwweights))$coef
    
    beta1est <- beta1mod[2,1]
    se_beta1 <- beta1mod[2,2]
    
    ll_CI_beta1_95 <- beta1est - 1.96*se_beta1
    ul_CI_beta1_95 <- beta1est + 1.96*se_beta1
    
    beta1_covered <- ifelse(ll_CI_beta1_95 <= beta1 && beta1 <= ul_CI_beta1_95, 1, 0)
    
    
  }else{
    beta1mod <- summary(geeglm(y ~ D + offset(2 - times) - 1, 
                                  family = familymod, data=obsdat, id = id, weights = iiwweights))$coef
    
    beta1est <- beta1mod[1]
    se_beta1 <- beta1mod[2]
    
    ll_CI_beta1_95 <- beta1est - 1.96*se_beta1
    ul_CI_beta1_95 <- beta1est + 1.96*se_beta1
    
    beta1_covered <- ifelse(ll_CI_beta1_95 <= beta1 && beta1 <= ul_CI_beta1_95, 1, 0)
    
  }
  

  
  out <- list(beta1est, se_beta1, beta1_covered)
  names(out) <- c("beta1est", "se_beta1", "beta1_covered")
  
  return(out)
  
}




simulateOneIIW<- function(n, beta1, beta2, beta3, beta4, beta5, beta6, gamma1, gamma2, gamma3, 
                          gamma4, gamma5, gamma6, tau, missingnesstypes, proportions, whichmissing, nimputations, 
                          outcome, usesplines, scheme){
  # Simulates one instance of the simulation, obtaining estimates for beta1 under various weighting
  # IIW uses stabilized weights
  
  if(scheme == "A"){
    # get required data from the data gen function 
    singlerun <- gendata_IIW(n, beta1, beta2, beta3, beta4, beta5, beta6,
                             gamma1, gamma2, gamma3, gamma4, gamma5, gamma6, tau, outcome)
  }else{
    singlerun <- gendata_IIWB(n, beta1, beta2, beta3, beta4, beta5, beta6,
                             gamma1, gamma2, gamma3, gamma4, gamma5, gamma6, tau, outcome)
  }

  
  simdata_full <- singlerun$simdata_full #counterfactual but censored data. Has Y even if not observed. (used for intensity model)
  
  simdata_observed <- singlerun$simdata_Yobserved # only data that is observed at observation times (used for outcome model)
  numevents <- singlerun$numevents
  newn <- singlerun$newn
  
  # terti<-quantile(0:tau , c(0.3333, 0.66666), type = 1) 
  terti<-quantile(0:tau , c(0.5), type = 1) 

  
  # need a function that imputes the missingness mechanism (MCAR< MAR,  etc, and a proportion of missingness. It also 
  # applies all 5 missing data methods)

  betamat <- data.frame(matrix(NA, nrow = 1, ncol = 5*((length(missingnesstypes)-3)*length(proportions)+1) + 2))
  semat <- data.frame(matrix(NA, nrow = 1, ncol = 5*((length(missingnesstypes)-3)*length(proportions)+1) + 2))
  coveragemat <- data.frame(matrix(NA, nrow = 1, ncol = 5*((length(missingnesstypes)-3)*length(proportions)+1) + 2))

  
  missingnesstypeno = 1
  missingnesstypenames <- c()
  
  for(missingnesstype in missingnesstypes){
    
    if(missingnesstype == "Naive"){
      
      #No IIW, no adjustments
      
      missingnesstypenames <- c(missingnesstypenames, missingnesstype)
 

      beta1_naive_results <- fitOutcomeModel(simdata_observed, outcome, usesplines, terti, iiwweights = NULL)
      
      
      betamat[1, missingnesstypeno] <- beta1_naive_results$beta1est
      semat[1, missingnesstypeno] <- beta1_naive_results$se_beta1
      coveragemat[1, missingnesstypeno] <- beta1_naive_results$beta1_covered
      
      missingnesstypeno = missingnesstypeno + 1
      


      
    }else if(missingnesstype == "NoMissingness"){
      
      #perform IIW, but not the adjustments for missingness since there is no missing data
      
      missingnesstypenames <- c(missingnesstypenames, missingnesstype)
      
      intensitydata <- makeMissing(simdata_full, missingnesstype = missingnesstype, p, whichmissing)

      
      #include a variable counting observation number, to be used for lagging time for Surv function
      intensitydata$obsnumber <- with(intensitydata, ave(id, id, FUN = seq_along))
      
      
      # #create lagged time variable
      intensitydata$times.lag <- intensitydata$times[c(nrow(intensitydata ),1:(nrow(intensitydata )-1))]
      intensitydata$times.lag[intensitydata$obsnumber == 1] <- -0.01
      
      
      gamma.hat <- coxph(Surv(times.lag, times, Yobserved) ~ D  + V1t + V2t + V3t + V4t - 1, data = intensitydata)$coef
      delta.hat <- coxph(Surv(times.lag, times, Yobserved) ~ D - 1, data = intensitydata)$coef
      
      iiw <- exp(cbind(simdata_observed$D)%*%delta.hat)/
        exp(cbind(simdata_observed$D, simdata_observed$V1t, simdata_observed$V2t, simdata_observed$V3t, simdata_observed$V4t)%*%gamma.hat)
      
      
      beta1_nm_results <- fitOutcomeModel(simdata_observed, outcome, usesplines, terti, iiwweights = iiw)
      
      
    
      betamat[1, missingnesstypeno] <- beta1_nm_results$beta1est
      semat[1, missingnesstypeno] <- beta1_nm_results$se_beta1
      coveragemat[1, missingnesstypeno] <- beta1_nm_results$beta1_covered
      
      missingnesstypeno = missingnesstypeno + 1
      
      
    }else{
      
      
    for(p in proportions){
      
        betavec <- c() #stores betas for each of the methods for a single missingnesstype/proportion
        sevec <- c()
        coveragevec <- c()

        #invoke missingness
        intensitydata <- makeMissing(simdata_full, missingnesstype = missingnesstype, p = p, whichmissing)
        
        
        #### impute data using the different methods
        CCAresults <- CCA(intensitydata, simdata_observed, outcome, usesplines)
        LOCFresults <- LOCF(intensitydata, simdata_observed, outcome, usesplines)
        SIresults <- SI(intensitydata, simdata_observed, outcome, usesplines)
        MIresults <- MI(intensitydata, simdata_observed, nimputations, outcome, usesplines)
        NIRFresults <- NIRF(intensitydata, simdata_observed, outcome, usesplines)
        
        betavec <- c(betavec, CCAresults$beta1est,
                     LOCFresults$beta1est, 
                     SIresults$beta1est,
                     MIresults$beta1est,
                     NIRFresults$beta1est)
        
        sevec <- c(sevec, CCAresults$se_beta1, 
                     LOCFresults$se_beta1, 
                     SIresults$se_beta1,
                     MIresults$se_beta1,
                     NIRFresults$se_beta1)
        
        coveragevec <- c(coveragevec, CCAresults$beta1_covered,
                         LOCFresults$beta1_covered, 
                     SIresults$beta1_covered,
                     MIresults$beta1_covered,
                     NIRFresults$beta1_covered)
        
        
        if(missingnesstype == "ObsTimesOnly"){
          missingnesstypenames <- c(missingnesstypenames, c(paste(missingnesstype, "CCA", sep = "_"), paste(missingnesstype, "LOCF", sep = "_"), 
                                          paste(missingnesstype, "SI", sep = "_"),
                          paste(missingnesstype, "MI", sep = "_"), paste(missingnesstype, "NIRF", sep = "_")))
        }else{
          missingnesstypenames <- c(missingnesstypenames, c(paste(missingnesstype, p, "CCA", sep = "_"), 
                                          paste(missingnesstype, p, "LOCF", sep = "_"), 
                                          paste(missingnesstype, p, "SI", sep = "_"),
                          paste(missingnesstype, p, "MI", sep = "_"), paste(missingnesstype, p, "NIRF", sep = "_")))
        }
        
        
        
        betamat[1, missingnesstypeno:(missingnesstypeno+4)] <- betavec
        semat[1, missingnesstypeno:(missingnesstypeno+4)] <- sevec
        coveragemat[1, missingnesstypeno:(missingnesstypeno+4)] <- coveragevec

        missingnesstypeno = missingnesstypeno + 5
        
        

        
        if(missingnesstype == "ObsTimesOnly"){
          break
        }
      }
    }
  }
  

  
  colnames(betamat) <- missingnesstypenames
  
  out <- list(betamat, semat, coveragemat, missingnesstypenames, simdata_full, numevents, newn)
  names(out) <- c("betamat", "semat", "coveragemat", "missingnesstypenames", "simdata_full", "numevents", "newn")
  
  
  return(out)
}



comb <- function(...){
  mapply("rbind", ..., SIMPLIFY = F)
}

simulateResultsIIW <- function(N, n, beta1, beta2, beta3, beta4, beta5, beta6, gamma1, gamma2, gamma3, 
                               gamma4, gamma5, gamma6, tau, missingnesstypes, proportions, whichmissing, nimputations, 
                               outcome, usesplines, scheme,
                               outputfulldatalist = T, nclusters = 2){
  # Simulates N instances of the each missingnesstype
  # computes in parallel if inParallel = T
  #

  nmissingnesstypes <- 5*((length(missingnesstypes)-3)*length(proportions)+1) + 2


  registerDoParallel(nclusters)
  
  results_beta1<- foreach(i = 1:N, .combine = comb, .export = c("expit", "simulateOneIIW",
                                                     "geeglm", "coxph", "Surv", "makeMissing", "CCA", "LOCF",
                                                      "SI", "MI", "NIRF", "IIW",
                                                     "gendata_IIW", "gendata_IIWB", "bs", "fitOutcomeModel"
                                                             ),
               .packages = c("dplyr", "geepack", "mice", "missForest", "splines", "parallel")) %dopar% {
               simrun <- simulateOneIIW(n, beta1, beta2, beta3, beta4, beta5, beta6, gamma1, gamma2, gamma3, gamma4,
                                        gamma5, gamma6, tau, missingnesstypes, proportions, whichmissing, nimputations,
                                        outcome, usesplines, scheme)

               betamat <- as.matrix(simrun$betamat)
               colnames(betamat) <- simrun$missingnesstypenames
               
               semat <- as.matrix(simrun$semat)
               colnames(semat) <- simrun$missingnesstypenames
               
               coveragemat <- as.matrix(simrun$coveragemat)
               colnames(coveragemat) <- simrun$missingnesstypenames
               
               nobsmat <- as.matrix(t(as.numeric(simrun$numevents)))
               colnames(nobsmat) <- c("Min", "Q1", "Med", "Mean", "Q3", "Max")
               
               newn <- simrun$newn
              
               
               list(betamat, semat, coveragemat, nobsmat, newn)}
  
    
    



  # Results including extreme values
  
  bias_beta1 <- round(apply(results_beta1[[1]], FUN =  mean, MARGIN = 2)  - beta1 ,3)

  var_beta1 <- round(apply(results_beta1[[1]], FUN = var, MARGIN = 2), 3)

  mse_beta1 <- round(apply((results_beta1[[1]] - beta1)^2, FUN = mean, MARGIN = 2), 3)
  
  avgse_beta1 <- round(apply(results_beta1[[2]], FUN =  mean, MARGIN = 2) ,3)
  
  coverage_beta1 <- round(apply(results_beta1[[3]], FUN =  mean, MARGIN = 2) ,3)
  
  numevents <-  round(apply(results_beta1[[4]], FUN =  mean, MARGIN = 2) ,3)
  
  newn <- round(apply(results_beta1[[5]], FUN =  mean, MARGIN = 2) ,3)
  
  
  # Results EXcluding extreme values
  toremove <- which(rowSums(abs(results_beta1[[1]]) > 100) > 0) #extreme if beta > 100
  
  if(length(toremove) > 0){
    
    extremitysums <- colSums(results_beta1[[1]] > 100) + colSums(results_beta1[[1]] < -100) #count how many extreme
    
    
    results_beta1_noextremity <- results_beta1[[1]] #copy
    se_noextremity <- results_beta1[[2]] #copy
    covered_noextremity <- results_beta1[[3]] #copy
    numevents_noextremity <- results_beta1[[4]] #copy
    
    #remove flagged observations
    
    results_beta1_noextremity <- results_beta1_noextremity[-c(toremove), ] #remove extreme results (beta>100)
    se_noextremity <- se_noextremity[-c(toremove), ] #remove extreme results (beta>100)
    covered_noextremity <- covered_noextremity[-c(toremove), ] #remove extreme results (beta>100)
    numevents_noextremity <- numevents_noextremity[-c(toremove), ] #remove extreme results (beta>100)
    
    
    
    nremoved <- nrow(results_beta1[[1]]) - nrow(results_beta1_noextremity)
    
    bias_beta1_noextremity <- round(apply(results_beta1_noextremity, FUN =  mean, MARGIN = 2)  - beta1 ,3)
    
    var_beta1_noextremity <- round(apply(results_beta1_noextremity, FUN = var, MARGIN = 2), 3)
    
    mse_beta1_noextremity <- round(apply((results_beta1_noextremity - beta1)^2, FUN = mean, MARGIN = 2), 3)
    
    avgse_beta1_noextremity <- round(apply(se_noextremity, FUN =  mean, MARGIN = 2) ,3)
    
    coverage_beta1_noextremity <- round(apply(covered_noextremity, FUN =  mean, MARGIN = 2) ,3)
    
    avgnumevents_noextremity <-  round(apply(numevents_noextremity, FUN =  mean, MARGIN = 2) ,3)
    
    
    out <- list(bias_beta1, var_beta1, mse_beta1, 
                avgse_beta1, coverage_beta1, numevents,
                bias_beta1_noextremity, var_beta1_noextremity, mse_beta1_noextremity,
                avgse_beta1_noextremity, coverage_beta1_noextremity, avgnumevents_noextremity,
                nremoved, extremitysums, newn,
                results_beta1)
    
    names(out) <- c("biasmat", "varmat", "msemat", 
                    "avgse_beta1", "coverage_beta1", "numevents",
                    'bias_beta1_noextremity', 'var_beta1_noextremity', 'mse_beta1_noextremity',
                    'avgse_beta1_noextremity', 'coverage_beta1_noextremity', 'avgnumevents_noextremity',
                    'nremoved', 'extremitysums', 'newn',
                    'results_beta1')
    
    
  } else{
    
    out <- list(bias_beta1, var_beta1, mse_beta1, 
                avgse_beta1, coverage_beta1, numevents, newn,
                results_beta1)
    names(out) <- c("biasmat", "varmat", "msemat", 
                    "avgse_beta1", "coverage_beta1", "numevents", 'newn',
                    'results_beta1')
    
  }
  

  return(out)
}




cleanUpResults <- function(resultsmat, proportions, missingnesstypes, methodnames){
  
  #create a clean table
  
  nproportions <- length(proportions)
  nmissingnesstypes <- length(missingnesstypes)
  nmethods <- length(methodnames)
  

   ######## ANY MAT ################################

   outmat <- matrix(NA, nrow = (nmissingnesstypes-3)*nproportions + 3, ncol = 2 + nmethods)
   colnames(outmat) <- c("Missingness Mechanism", "pmiss", methodnames)
   outmat[1,] <- c("No IIW/Naive", "-", resultsmat[1], "-", "-", "-", "-")
   outmat[2,] <- c("No Missingness", 0, resultsmat[2], "-", "-", "-", "-")

   #Obstimesonly
   outmat[3,] <- c(missingnesstypes[3], "-", resultsmat[3:7])

   #MCAR
   outmat[4,] <- c(missingnesstypes[4], proportions[1], resultsmat[8:12])
   outmat[5,] <- c("", proportions[2], resultsmat[13:17])
   outmat[6,] <- c("", proportions[3], resultsmat[18:22])

   #MAR
   outmat[7,] <- c(missingnesstypes[5], proportions[1], resultsmat[23:27])
   outmat[8,] <- c("", proportions[2], resultsmat[28:32])
   outmat[9,] <- c("", proportions[3], resultsmat[33:37])

   #MNAR
   outmat[10,] <- c(missingnesstypes[6], proportions[1], resultsmat[38:42])
   outmat[11,] <- c("", proportions[2], resultsmat[43:47])
   outmat[12,] <- c("", proportions[3], resultsmat[48:52])

   outtab <- kable(outmat, digits = 3, format = "latex", booktabs = T)
  
   return(outtab)

  
}



