varpop<-function(freqMat,H,rho=2,wy=5,r,mu=0.001,b=0,alpha=1,duration=20,diagnostic=FALSE)
    {
        ##Define number of phases to examine
        T=nrow(freqMat)-1
        ##Define population size of vessels from recovery rate
        ni=apply(freqMat,1,sum)
        Ni=ni*r

        ##Define number of time steps so for given number of potters and pots we match the target population      
        etaSequence=round(((Ni[2:length(H)]*2)/rho - H[2:length(H)] + H[1:c(length(H)-1)])
            / (H[2:length(H)] + H[1:c(length(H)-1)]))
     
        ##Convert Memory parameter to actual number of timesteps
        w=ceiling(sum(etaSequence)/(duration*T)*wy)
        
        ## Define House Sequence (Define number of individuals in the last m steps of the previous phase)
        houseSequence=c(H[1],H)        
        initialPoolSize<-ceiling(seq(houseSequence[1],houseSequence[2],length.out=etaSequence[1])[c(etaSequence[1]-w+1):etaSequence[1]]*rho)	

        ## Generate time-varying sequence of potters 
        Nseq=vector("list",length=T)
        for (x in 1:T)
            {
                Nseq[[x]]=ceiling(seq(houseSequence[x+1],houseSequence[x+2],length.out=etaSequence[x])*rho)
            }

        ## Variants Settings ##
        k_initial=ncol(freqMat)
        colnames(freqMat)=1:k_initial
        index.check=which(freqMat[c(1),]>0) #traits that exist during the first phase
        startSeq=freqMat[1,]
        newTraitCounter=k_initial

        ## Conduct Dirichlet sampling ##
        totalSamplePool=sample(x=1:k_initial,size=sum(initialPoolSize),prob=rdirichlet(1,startSeq+alpha),replace=TRUE)
        totalSamplePool[1:length(index.check)]=index.check #ensure all initial variants exist in the sample pool
        totalSamplePool=sample(totalSamplePool) #shuffle
        traitList<-vector("list",length=w)
        for (x in 1:length(traitList)) #redistribute
            {
                traitList[[x]]=totalSamplePool[1:initialPoolSize[x]]
                totalSamplePool=totalSamplePool[-c(1:initialPoolSize[x])]
            }

        ## Output Storage Declaration ##
        PopulationPerPhase<-vector("list",length=length(etaSequence))
        SamplePerPhase<-vector("list",length=length(etaSequence))
        
        ## Core Simulation Starts Here ##
        timeCounter=1

        for (i in 1:length(etaSequence))
            {
                if(Nseq[[i]][1]<Nseq[[i]][length(Nseq[[i]])]){option="increase"}
                if(Nseq[[i]][1]>Nseq[[i]][length(Nseq[[i]])]){option="decrease"}

                for (t in 1:etaSequence[i])
                    {
                        if(t==1){currentTraits=numeric(length=Nseq[[i]][1])}
                        if (t>1)
                            {
                                        #population decrease:        
                                if (option=="decrease") {currentTraits<-sample(currentTraits,size=Nseq[[i]][t])}
                                        #population increase:        
                                if (option=="increase") {currentTraits<-c(currentTraits,sample(currentTraits,size=Nseq[[i]][t]-length(currentTraits),replace=TRUE))}
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
                        PopulationPerPhase[[i]]<-c(PopulationPerPhase[[i]],currentTraits)

                    }
                SamplePerPhase[[i]]<-sample(PopulationPerPhase[[i]],size=ni[i+1]) ##All pots produced in the second phase

            }
        ## List all traits
        traits=unique(c(1:k_initial, unlist(PopulationPerPhase)))
        ## Count Frequencies per phase
        simFreq=t(sapply(SamplePerPhase,instances,variants=traits))

        simFreq.prop=prop.table(simFreq,1)[,index.check]
        obsFreq.prop=prop.table(freqMat[-1,],1)[,index.check]

        epsilon=sqrt(sum((simFreq.prop-obsFreq.prop)^2))

        ## Divide by number of control points ##
        epsilon=epsilon/length(simFreq.prop)
        
        if (diagnostic==TRUE)
            {
                return(list(simFreq=simFreq.prop,epsilon=epsilon))
            }
        if (diagnostic==FALSE)
            {
                return(epsilon=epsilon)
            }
        
    }
