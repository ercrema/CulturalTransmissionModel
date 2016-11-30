### Read Data ###
obsData=as.matrix(read.csv("data/observedFrequencies.csv",row.names=1))

### Source Functions ###
source("src/utility.R")
source("src/varpoptrans.R")

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
nsim=5000 #number of simulations
transitions=7

### Create Parameter Space ###
parameters=data.frame(mu=runif(nsim,prior_mu[1],prior_mu[2]),
		      rho=runif(nsim,prior_rho[1],prior_rho[2]),
		      b=runif(nsim,prior_b[1],prior_b[2]),
		      wy=runif(nsim,prior_wy[1],prior_wy[2]),
		      r=runif(nsim,prior_r[1],prior_r[2]))

### Loop index ###
tiv=data.frame(h1=1:7,
	       h2=2:8,
	       target=2:8,
	       ha=c(1,1:6),
	       hb=1:7,
	       hc=2:8,
	       startF=1:7,
	       endF=2:8)

epsilonMatrix=matrix(NA,nrow=nsim,ncol=7)

### Execute Simulations ###

for (x in 1:nsim)
	{
		target=sampleRange*parameters$r[x]
	for (y in 1:transitions)
		{
			eta=round(((target[tiv$target[y]]*2)/parameters$rho[x] -
				   H[tiv$h1[y]] + H[tiv$h2[y]])/
				(H[tiv$h1[y]] + H[tiv$h2[y]]))
			w=ceiling(parameters$wy[x]/duration*eta)
			Hseq=H[c(tiv$ha[y], tiv$hb[y],tiv$hc[y])]
			tmp=varpoptrans(startSeq=obsData[tiv$startF[y],],
					   endSeq=obsData[tiv$endF[y],],
					   H=Hseq,
					   rho=parameters$rho[x],w=w,
					   mu=parameters$mu[x],
					   b=parameters$b[x],
					   eta=eta,
					   alpha=1,diagnostic=FALSE)
			epsilonMatrix[x,y]=as.numeric(tmp)
		}
	}
parameters=as.data.frame(cbind(parameters,epsilonMatrix))
colnames(parameters)=c("mu","rho","b","wy","r","phase7to8","phase8to9","phase9to10","phase10to11","phase11to12","phase12to13","phase13to14")



### Compute Posterior ### 
cutoff=0.1 #in this case a top-skim of the best 10%
threshold=cutoff*nsim
posterior_b95=matrix(NA,nrow=2,ncol=7)
posterior_b50=matrix(NA,nrow=2,ncol=7)
phases=c("phase7to8","phase8to9","phase9to10","phase10to11","phase11to12","phase12to13","phase13to14")
library(LaplacesDemon)

for (x in 1:7)
	{
		index=order(parameters[phases[x]])[1:threshold]
		posterior=parameters[index,]
		posterior_b95[,x]=p.interval(posterior$b,prob=0.95)[1:2]
		posterior_b50[,x]=p.interval(posterior$b,prob=0.50)[1:2]
	}

### Plot Posterior ###
plot(1,xlim=c(0,8),ylim=prior_b,type="n",axes=F,xlab="Phases",ylab="b")
axis(side=1,at=1:7,label=as.roman(8:14))
axis(side=2)
for (x in 1:7)
	{
		rect(xleft=x-0.3,xright=x+0.3,
		     ybottom=posterior_b95[1,x],
		     ytop=posterior_b95[2,x])
		     
		rect(xleft=x-0.3,xright=x+0.3,
		     ybottom=posterior_b50[1,x],
		     ytop=posterior_b50[2,x],
		     col=1)
	}
abline(h=0,lty=2,col=2)






					




