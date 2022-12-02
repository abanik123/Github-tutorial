## number of individuals
N <- 800
## number of observation periods
T <- 8
## time period of first capture
first <- rep(1, N)
## length of observation history from first capture
len <- T - first + 1
## survival probability
phi <- 0.7
## detection probability
pVec <- rep(c(0.2, 0.8), each=N/2)

## simulate z (alive / dead status),
## and y (encounter histories)
set.seed(0)
z <- matrix(NA, nrow=N, ncol=T)
y <- matrix(NA, nrow=N, ncol=T)
for(i in 1:N) {
  z[i, first[i]] <- y[i, first[i]] <- 1
  for(t in (first[i]+1):T) {
    z[i,t] <- rbinom(1, 1, phi*z[i,t-1])
    y[i,t] <- rbinom(1, 1, pVec[i]*z[i,t])
  }
}

##################################################
##################################################

## Non-Parametric Capture-Recapture Model
code <- nimbleCode({
  phi ~ dunif(0, 1)
  alpha ~ dgamma(1, 1)
  xi[1:N] ~ dCRP(conc=alpha, size=N)
  for(i in 1:M)   p[i] ~ dunif(0, 1)
  for(i in 1:N) {
    y[i,first[i]:T] ~ dCJS_ss(phi, p[xi[i]], len=len[i])
  }
})
M <- 100   ## maximum number of subgroups
constants <- list(N=N, T=T, first=first, len=len, M=M)
data <- list(y=y)
inits <- list(phi=0.5, alpha=1, xi=rep(1,N), p=rep(0.5,M))
Rmodel <- nimbleModel(code, constants, data, inits)

McMC_m1 <- buildMCMC(Rmodel, enableWAIC = TRUE)
compile_model_m1<- compileNimble(Rmodel)
compile_MCMC_m1 <- compileNimble(McMC_m1, project = Rmodel)
samples_m1 <- runMCMC(compile_MCMC_m1, niter = 100, nchains = 2,
                      nburnin = 5, thin = 10,summary = TRUE,
                      WAIC = TRUE)