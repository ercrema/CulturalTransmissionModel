### Read Data ###
obsData=read.csv("data/observedFrequencies.csv",row.names=1)

### Source Functions ###
source("src/utility.R")
source("src/equilibrium.R")

### Define Prior Ranges ###
prior_mu=c(0.001,0.01)
prior_b=c(-0.5,0.5)
prior_rho=c(1,5)
prior_wy=c(0.2,5)
prior_r=c(1/0.75,1/0.25)

### Define Constants and Settings ###
duration=20
sampleRange=apply(obsData,1,sum)
H=c(17,10,12,14,15,16,12,8)
warmUp=5000
nsim=5000 #number of simulations


### Create Parameter Space ###
parameters=data.frame(mu=runif(nsim,prior_mu[1],prior_mu[2]),
		      rho=runif(nsim,prior_rho[1],prior_rho[2]),
		      b=runif(nsim,prior_b[1],prior_b[2]),
		      wy=runif(nsim,prior_wy[1],prior_wy[2]),
		      w=NA,
		      r=runif(nsim,prior_r[1],prior_r[2]),
		      epsilon=NA)

### Execute Simulations ###
for (x in 1:nsim)
	{
		eta=ceiling(c(sampleRange*parameters$r[x]/(round(mean(H))*parameters$rho[x])))
		parameters$w[x]=ceiling(sum(eta)/(duration*length(sampleRange))*parameters$wy[x])
		tmp=equilibrium(eta=eta,
				w=parameters$w[x],
				rho=parameters$rho[x],
				mu=parameters$mu[x],
				b=parameters$b[x],
				warmUp=warmUp,
				sampleRange=sampleRange,
				Hhat=round(mean(H)),
				obsData=as.matrix(obsData),
				T=nrow(obsData),
				diagnostic=FALSE)
		parameters$epsilon[x]=tmp$raw.epsilon
	}

### Compute Posterior ### 
cutoff=0.1 #in this case a top-skim of the best 10%
threshold=cutoff*nsim
index=order(parameters$epsilon)[1:threshold]
posterior=parameters[index,]

### Plot Posterior ###
library(LaplacesDemon)
p.interval(posterior$b,plot=TRUE,xlim=c(prior_b))


















