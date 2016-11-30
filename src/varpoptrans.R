##################################################
#    Variable Population/Transmission Version    #
##################################################

### PARAMETERS:
## rho ... production rate per household
## mu  ... innovation rate
## b   ... frequency bias
## wy   ... memory (in calendar years)
## r   ... 1/sampling fraction
## H ... number of households for phase i-1, i, and i+1
## startSeq ... observed frequencies of cultural variants at phase i-1
## endSeq ... observed frequencies of cultural variants at phase i
## eta ... number of transmission events
## duration ... duration of each phase (in calendar years)
## alpha ... parameter for the Dirichlet sampling 
## diagnostic ... if set to TRUE returns predicted frequencies of cultural variants for each phase.

varpoptrans<-function(startSeq,endSeq,H,rho,w,eta,mu,b,alpha,diagnostic=FALSE)
{
    ##define number of possible traits that can be discovered
    k_initial=length(startSeq)
    newTraitCounter=k_initial+1
    index.check=which(startSeq>0) #traits that exist during the first phase

    ##Define number of individuals in the last m steps of the previous phase#
    lastHouses<-ceiling(seq(H[1],H[2],length.out=eta)[c(eta-w+1):eta]*rho)	
    ##Define change in the number of individuals in the current phase
    Nseq<-ceiling(seq(H[2],H[3],length.out=eta)*rho)
    ##Define the initial pool of artefacts#
    names(startSeq)=1:length(startSeq)
    ## Conduct Dirichlet sampling ##
    totalSamplePool=sample(x=1:k_initial,size=sum(lastHouses),prob=rdirichlet(1,startSeq+alpha),replace=TRUE)
    totalSamplePool[1:length(index.check)]=index.check #ensure all initial variants exist in the sample pool
    totalSamplePool=sample(totalSamplePool) #shuffle
    traitList<-vector("list",length=w)
    for (x in 1:length(traitList)) #redistribute
        {
            traitList[[x]]=totalSamplePool[1:lastHouses[x]]
            totalSamplePool=totalSamplePool[-c(1:lastHouses[x])]
        }
    initialTraitList=traitList



    ## Core Simulation Starts Here ##


    ##Settings:
    if(Nseq[1]<Nseq[length(Nseq)]){option="increase"}
    if(Nseq[1]>Nseq[length(Nseq)]){option="decrease"}
    timeCounter=1
    population<-numeric()
    for (t in 1:eta)
        {
            ##Change in the number of potters#	
            if(t==1){currentTraits=numeric(length=Nseq[1])}	
            if (t>1)
                {
                    ##population decrease:        
                    if (option=="decrease") {currentTraits<-sample(currentTraits,size=Nseq[t])}
                    ##population increase:        
                    if (option=="increase") {currentTraits<-c(currentTraits,sample(currentTraits,size=Nseq[t]-length(currentTraits),replace=TRUE))}
                }
            
            ##transmission
            if (length(unique(unlist(traitList)))>1)
                {
                    sampleSpace<-table(unlist(traitList))
                    sampleSpace<-sampleSpace^(1-b) #transmission bias
                    currentTraits<-sample(size= length(currentTraits),
                                            x=as.numeric(names(sampleSpace)),replace=TRUE,prob=sampleSpace)
                }	

            ##mutation
            index<-which(runif(length(currentTraits))<mu)
            if (length(index)>0)
                {
                    newTraits<-newTraitCounter+1:length(index)
                    newTraitCounter<-max(newTraits)    
                    currentTraits[index]=newTraits
                }

            ##expression
            traitList[[timeCounter]]<-currentTraits
            timeCounter=timeCounter+1; if(timeCounter>w){timeCounter=1}

            ##store pots
            population<-c(population,currentTraits)
        }

    
    sampledTraits<-sample(population,size=sum(endSeq)) #All pots produced in the second phase

                                        #Analysis:
    traits=unique(c(1:length(startSeq), sampledTraits))

    ## Observed Data ###
    start<-instances(rep(1:length(startSeq),startSeq),traits)
    end<-instances(rep(1:length(endSeq), endSeq),traits)
    startProb<-start/sum(start)
    endProb<-end/sum(end)

    ## Simulation Output ###
    startSim<-instances(unlist(initialTraitList),traits)	
    endSim<-instances(sampledTraits, traits)
    endSimProb<-endSim/sum(endSim)
    
    ##Compute Epsilon#
    starters<-which(startProb>0)
    epsilon = sqrt(sum((endProb[starters]-endSimProb[starters])^2))
    epsilon=epsilon/length(starters)


    if(diagnostic==TRUE)
        {
    return(list(epsilon=epsilon,endSimProb=endSimProb[starters]))
        }
    if(diagnostic==FALSE)
        {
    return(epsilon)
        }
}



            

