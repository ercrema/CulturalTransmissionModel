#################################
#    Equilibrium Version        #
#################################

### PARAMETERS:
## rho ... production rate per household
## eta ... number of transmission events 
## mu  ... innovation rate
## b   ... frequency bias
## w   ... memory (in number of transmission events)
## T   ... number of phases
## warmUp ... number of transmission events for the burn-in phase
## sampleRange ... observed sample size (n) for each phase
## Hhat ... average number of households
## obsData ... observed frequencies of cultural variants
## diagnostic ... if set to TRUE returns predicted frequencies of cultural variants for each phase. 


equilibrium<-function(rho,eta,mu,b,T,warmUp,w,sampleRange,Hhat,obsData,diagnostic=FALSE)
    {

        ##Compute Population Size:
        nu=round(rho*Hhat)

        ##round number of transmission events:
        newTraitCounter=nu

        ##Initialise with all Agents having a different trait
        initialTraits<-1:nu


        ##Outputs#

        traitList<-vector("list",length=w)
        fullSamplePerPhase<-vector("list",length=T)
        timeCounter=2
        if (w==1){timeCounter=1}

        ##Initial Production
        traitList[[1]]=initialTraits
        currentTraits<-initialTraits

### Start of Burn-in Stage ###     


        for (t in 2:warmUp)
            {
                ##transmission
                if (length(unique(unlist(traitList)))>1)
                    {
                        sampleSpace<-table(unlist(traitList))
                        sampleSpace<-sampleSpace^(1-b) #transmission bias
                        currentTraits<-sample(size=nu,x=as.numeric(names(sampleSpace)),replace=TRUE,prob=sampleSpace)
                    }


                ##mutation
                index<-which(runif(nu)<mu)
                if (length(index)>0)
                    {
                        newTraits<-newTraitCounter+1:length(index)
                        newTraitCounter<-max(newTraits)    
                        currentTraits[index]=newTraits
                    }

                ##expression 
                traitList[[timeCounter]]<-currentTraits
                timeCounter=timeCounter+1; if(timeCounter>w){timeCounter=1}
            }

### End of Burn-in Stage ###     

        for (x in 1:T)
            {
                for (y in 1:eta[x])
                    {
                                        #transmission
                        if (length(unique(unlist(traitList)))>1)
                            {
                                sampleSpace<-table(unlist(traitList))
                                sampleSpace<-sampleSpace^(1-b) #transmission bias
                                currentTraits<-sample(size=nu,x=as.numeric(names(sampleSpace)),replace=TRUE,prob=sampleSpace)
                            }

                        index<-which(runif(nu)<mu)
                        if (length(index)>0)
                            {
                                newTraits<-newTraitCounter+1:length(index)
                                newTraitCounter<-max(newTraits)    
                                currentTraits[index]=newTraits
                            }

                                        #expression
                        traitList[[timeCounter]]<-currentTraits
                        timeCounter=timeCounter+1; if(timeCounter>w){timeCounter=1}

                                        #store pots
                        fullSamplePerPhase[[x]]<-c(fullSamplePerPhase[[x]],currentTraits)

                    }
                                        # end of block procedures:
                replace=TRUE
                if (length(fullSamplePerPhase[[x]])>sampleRange[x]) {replace=FALSE}
                fullSamplePerPhase[[x]]<-sample(fullSamplePerPhase[[x]],size=sampleRange[x],replace=replace)
            }


#Post Analysis:
allTraits=unique(unlist(fullSamplePerPhase))
resMatrix<-t(sapply(fullSamplePerPhase,instances,allTraits))
        
#Calculate diversity indices
        if (length(allTraits)>1)
            {
                d.epsilon=sqrt(sum((simpson(resMatrix)-simpson(obsData))^2))
            }
        if (length(allTraits)==1)
            {
                d.epsilon=sqrt(sum((0-simpson(obsData))^2))              
            }
#Reorder resMatrix based on rank at first time-block

        simFirstK=1 #placeholder: k during the first phase
        obsData=reOrder(obsData) #reorder observed data
        freqObs=prop.table(obsData,1) #proportions in the observed data
        obsFirstK=sum(obsData[1,]>0) #observed k during the first phase
        freqObs=freqObs[,1:obsFirstK]  #consider only variants present at the beginning
        freqSim=rep(1,nrow(obsData))
        if (length(allTraits)>1)
            {
                resMatrix=reOrder(resMatrix)
                simFirstK=sum(resMatrix[1,]>0)
                freqSim=prop.table(resMatrix,1)              
                freqSim=freqSim[,1:simFirstK]
            }
        
        if (simFirstK<obsFirstK) #if simulation has smaller initial k, add empty columns
            {
                freqSim=cbind(freqSim,matrix(0,nrow=8,ncol=obsFirstK -simFirstK))
            }

        if (simFirstK>=obsFirstK)  #if simulation has larger or equal initial k, trim the frequency matrix
            {
                freqSim=freqSim[,1:obsFirstK]
            }
        raw.epsilon=sqrt(sum((freqObs-freqSim)^2))

        
        ## Divide by number of control points ##
        raw.epsilon=raw.epsilon/length(freqObs)
        d.epsilon=d.epsilon/nrow(obsData)
        
        
        if (diagnostic==FALSE)
            {
                return(list(d.epsilon=d.epsilon,raw.epsilon=raw.epsilon))
            }

        if (diagnostic==TRUE)
            {
                return(list(raw=freqSim,d.epsilon=d.epsilon,raw.epsilon=raw.epsilon))
            }
    }


