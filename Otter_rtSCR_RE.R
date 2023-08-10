#==============================================================================#
#                                                                              #
#                RANDOM THINNING SPATIAL CAPTURE-RECAPTURE (rtSCR)             #
#                            ~~~ Lutra lutra ~~~                               #
#                   Random effects on baseline detection rate (traps)          #
#                            09/08/2023 9:41:42                                #
#                                                                              #
#==============================================================================#

setwd('C:/.../BData')

library(nimble)
library(basicMCMCplots)
library(coda)
library(lattice)
library(raster)
library(secr)
library(jagsUI)
library(makeJAGSmask)
library(rgdal)
library(scrbook)
library(mcmcOutput)
library(MCMCvis)
library(tidyverse)
library(ggplot2)
source("Functions_SCR.R")

# Available habitat
tablas <- raster("PNTD_ras.tif")

otterID.ch <- read.capthist("OtterData.txt", "OtterTraps.txt", 
  detector='count', cov="sex", noccasions=1)
summary(otterID.ch)
otterID<-aperm(otterID.ch,c(1,3,2))
yobs<-apply(otterID,c(1,2),sum)
(nind<-dim(yobs)[1])   # Detected and ID

# Centroids of hexagons
traplocs<-traps(otterID.ch)
X<-data.matrix(traplocs)
rownames(X)<-1:nrow(traplocs)
colnames(X)<-c("X","Y")

# Plot of ID
my_window <- extent(c(min(traplocs[,1])-1000, max(traplocs[,1])+1000, 
                    c(min(traplocs[,2])-1000, max(traplocs[,2])+1000)))
plot(my_window, type="n", asp=TRUE)
tot<-apply(otterID, 2,sum)
symbols(traplocs, circles=tot*100, inches=F,bg="#228B2219", fg=NULL, add=T)
points(traplocs, pch="+", col="blue")
# Spiderplot of ID
spiderplotJJ5(otterID, traplocs, buffer=2000, lwd=1)

# Data augmentation
M<-50
# y detected with data augmentation
y.obs<-array(0,c(M,nrow(X)))
y.obs[1:nind,]<-yobs

# Non-ID detections
otterNonID.ch <- read.capthist("otter_NonID.txt", "OtterTraps.txt", 
  detector='count', noccasions=1)
summary(otterNonID.ch)
otterNonID<-aperm(otterNonID.ch,c(1,3,2))
nnidd<-apply(otterNonID,c(2,3),sum)

ras<-stack(tablas)
names(ras)<-c('tablas')
mymask <- convertRaster(ras, as.data.frame(X))
str(mymask)

# Sampling effort by cell
Eff<-read.table("Effort.txt", header=TRUE)[,2]
Eff<-(Eff-mean(Eff))/sd(Eff)

## rt-SCR model
code <- nimbleCode({

  sigma ~ dunif(0, 500)     # set up the priors for the parameters
  beta ~ dnorm(0,0.01)
  psi ~ dunif(0, 1)
  # Random effects hyperparameters
  sigma.p ~ dunif(0, 10)
  mu0 ~ dnorm(0, 0.01)
  id.prob ~ dunif(0,1)
  
  for(j in 1:nTraps){
    eps[j] ~ dnorm(mu0, sd=sigma.p)
    log(p0[j]) <- beta*Eff[j] + eps[j]
  }

  for (i in 1:M){             # loop through the augmented population
    z[i] ~ dbern(psi)         # state of individual 'i' (real or imaginary)
    S[i, 1] ~ dunif(1, upperLimit[1]) # uniform priors for the activity 
    S[i, 2] ~ dunif(1, upperLimit[2]) #    centres for each individual
    pOK[i] <- habMat[trunc(S[i,1]), trunc(S[i,2])] # habitat check
    OK[i] ~ dbern(pOK[i])     # OK[i] = 1, the ones trick
    Dsq[i,1:nTraps] <- (S[i,1]-trapMat[1:nTraps,1])^2 + 
                       (S[i,2]-trapMat[1:nTraps,2])^2   
    lam[i,1:nTraps] <- p0[1:nTraps]*exp(-Dsq[i,1:nTraps]/(2*sigma^2))*nOcc*z[i]
    # Vectorized Poisson
    y.full[i,1:nTraps] ~ dPoissonVector(lam[i,1:nTraps]) # vectorized Poisson
    # y.full.sim[i,1:nTraps] ~ dPoissonVector(lam[i,1:nTraps]) # GoF data simulation
    for(j in 1:nTraps) {
      y.obs[i,j] ~ dbin(id.prob, y.full[i,j])
    }
  }

  N <- sum(z[1:M])   # derive number (check that N << M)
  sigmaR <- sigma * pixelWidth   # sigma at real scale
})

# Vectorized Poisson
dPoissonVector <- nimbleFunction(
  run = function(x = double(1), lambda = double(1),
  log = integer(0, default = 0)) {
    J <- length(x)
    ans <- 0.0
    for(j in 1:J)
      ans <- ans + dpois(x[j], lambda[j], 1)
    returnType(double())
    if(log) return(ans)
    else return(exp(ans))
  })

rPoissonVector  <- nimbleFunction(
  run = function(n = integer(), lambda = double(1)) {
    J <- length(lambda)
    ans<- numeric(J)
    for(j in 1:J)
      ans[j] <- rpois(1, lambda[j])
    returnType(double(1))
    return(ans)
  })

registerDistributions(list(
  dPoissonVector = list(
    BUGSdist = "dPoissonVector(lambda)",
    Rdist = "dPoissonVector(lambda)",
    discrete = TRUE,
    range = c(0, Inf),
    types = c('value = double(1)', 'lambda = double(1)'))
))

# Using Metropolis-Hastings that obey: yNonID <- yTRUE - yID
source("ID_sampler.R")

nTraps <- nrow(traplocs)
# DATA
str(data   <-   list(habMat=mymask$habMat,
                     trapMat=mymask$trapMat,
                     upperLimit=mymask$upperLimit,
                     OK = rep(1, M),
                     Eff=Eff,  
                     y.obs=y.obs))

# CONSTANTS                     
str(constants <-list(M = M, 
                     nTraps = nTraps,
                     pixelWidth=mymask$pixelWidth,
                     nnidd=nnidd,
                     nOcc = 1))


## INITS
# Inits for si
set.seed(1)
trapMat<-mymask$trapMat
ys<-apply(y.obs,c(1,2),sum)
s.start <- trapMat[sample(1:nTraps, size=M, replace=TRUE),]
for(i in 1:nind){
if(sum(apply(ys,1,sum)>0)) next
  s.start [i,1]<- mean( trapMat[ys[i,]>0,1] )
  s.start [i,2]<- mean( trapMat[ys[i,]>0,2] )
}
d <- e2dist(s.start[1:M,], trapMat)
sigma.sst  <- 150  # sigma.sst: choose an appropriate value
lam0.sst <- 0.1    # lam0.sst: choose an appropriate value
lam <- lam0.sst * exp( -(d^2)/(2 * sigma.sst^2))
K<-1

# Init for nonID indiviuals
yi <- array(0, c(M, nTraps, K)) # resighting array
for (j in 1:nTraps) {
  for (k in 1:K) {
    if (nnidd[j, k] > 0) {
      probs <- lam[ ,j]
      probs <- probs / sum(probs)
      latent.id <- sample(1:M, nnidd[j,k], prob = probs, replace = FALSE)
      yi[latent.id , j, k] <- 1
    }
  } # end of k
}   # end of j

# Init for true (yT) population: yT = yNonID + yID
yT <- apply(yi,c(1,2),sum) + y.obs
zst<-apply(yT, 1, sum); zst[zst>0]<-1
(id.prob.sst<-sum(y.obs)/(sum(y.obs)+sum(nnidd)))


str(inits    <- list(z = rep(1, M),
                     sigma=sigma.sst,
                     beta=runif(1,-5,5),
                     psi=runif(1,0,1),
                     mu0=runif(1,-5,5),
                     sigma.p=runif(1,0,5),
                     y.full=yT,
                     #y.full.sim=yT,
                     id.prob=id.prob.sst,
                     S = s.start))

# Estimates to monitor:
params <- c('N','psi','beta','sigmaR','mu0','sigma.p','id.prob')

Rmodel <- nimbleModel(code=code, constants=constants, data=data, inits=inits)
# Function to initialize 'complex' nodes
InitNod<-function(simNodes){
  simNodeScalar <- Rmodel$expandNodeNames(simNodes)
  allNodes <- Rmodel$getNodeNames()
  nodesSorted <- allNodes[allNodes %in% simNodeScalar]
  set.seed(1) # to fix simulations
  for(n in nodesSorted) {
    Rmodel$simulate(n)
    depNodes <- Rmodel$getDependencies(n)
    Rmodel$calculate(depNodes)
  }
}
InitNod(simNodes = 'eps')
Rmodel$initializeInfo()
Rmodel$calculate()
Cmodel <- compileNimble(Rmodel)

conf<-configureMCMC(Rmodel,  monitors=params, thin = 5)

#### --- load custom samplers at end of script  --- ####
# replace with new sampler for y (sample without replacement with sum 
# to n[j,k] constraint)
conf$removeSampler("y.full")
for(j in 1:nTraps){
  conf$addSampler(target = paste(paste("y.full[1:",M,", ",j,"]"), sep=""),
                  type = 'IDSampler',control = list(nnidd = nnidd[j], j=j, M=M),
                  silent = TRUE)
}

# Rebuild and compile with new sampler
conf$removeSamplers("S")
ACnodes <- paste0("S[", 1:constants$M, ", 1:2]")
for(node in ACnodes) {
  conf$addSampler(target = node,
                  type = "AF_slice",
                  control = list(adaptScaleOnly = TRUE),
                  silent = TRUE)
}

MCMC <- buildMCMC(conf)

CompMCMC <- compileNimble(MCMC, project = Rmodel)

## Execute MCMC algorithm and extract samples
nb <- 10000       # Burnin
ni <- 50000 + nb  # Iters
nc <- 3           # Chains


start.time2<-Sys.time()
outNim <- runMCMC(CompMCMC, niter = ni , nburnin = nb , nchains = nc, 
                  inits=inits, setSeed = TRUE, progressBar = TRUE, 
                  samplesAsCodaMCMC = TRUE)
end.time<-Sys.time()
end.time-start.time2 # post-compilation run time

summary(outNim)
