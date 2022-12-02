#Defining Libraries
library(nimble)
library(MCMCvis)

#####################################################
# NIMBLE Model Code:-

code_m1 <- nimbleCode({
  #--------------------------------------
  
  # Priors
  phi[1:2] ~ ddirch(alpha[1:2])
  p[1:2] ~ ddirch(alpha[1:2])
  
  
  # Likelihood 1 (for first year capture)
  #Likelihood 2 (from 2nd year of capture)
  for (i in 1:N){
    for (t in first[i]:first[i]){
      a[i,t] <-1
    }
    for (t in (first[i]+1):nyears){
      
      mu1[i,t]<- (a[i,t-1]*phi[1])
      a[i,t] ~ dbern(mu1[i,t])
      
      mu2[i,t] <- a[i,t]*p[1]
      h[i,t] ~ dbern(mu2[i,t])
    }
  }
})

#########################################################
# Preparing Data:-

cap_data <- read.csv(file.choose(),header = F)
alive_data <- read.csv(file.choose(), header= F)

N = nrow(cap_data)
nyears = ncol(cap_data)

for (i in 1:7){
  names(cap_data)[i]<- i
}

m<-colnames(cap_data)[max.col(cap_data, ties.method = "first")]

first <- rep(NA, N)
for (i in 1:N){
  first[i]<- as.numeric(m[i])
}

##########################################################
# MCMC:-

model_m1 <- nimbleModel(code_m1, 
                        constants = list(N = N,
                                         nyears = nyears,
                                         alpha = c(1,1),
                                         first = first),
                        data = list(a= alive_data,
                                    h = cap_data), 
                        inits = list(phi = rep(1/2,2), p = rep(1/2,2)))

McMC_m1 <- buildMCMC(model_m1, enableWAIC = TRUE)
compile_model_m1<- compileNimble(model_m1)
compile_MCMC_m1 <- compileNimble(McMC_m1, project = model_m1)
samples_m1 <- runMCMC(compile_MCMC_m1, niter = 10000, nchains = 2,
                      nburnin = 500, thin = 10,summary = TRUE,
                      WAIC = TRUE)

waic <- samples_m1$WAIC
waic

samples<- samples_m1$samples

############################################################

#Model Summary:-

MCMCsummary(samples,params = "p", round = 5)

MCMCsummary(samples,params = "p", round = 5)$mean

MCMCtrace(samples,pdf = F, ind = TRUE,
          Rhat = TRUE, n.eff = TRUE)