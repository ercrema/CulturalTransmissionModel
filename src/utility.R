
## Function for Dirichlet sampling ##
rdirichlet<-function (n, alpha) 
{
    l <- length(alpha)
    x <- matrix(rgamma(l * n, alpha), ncol = l, byrow = TRUE)
    sm <- x %*% rep(1, l)
    x/as.vector(sm)
}

## Utility function for counting number of variants ##
instances<-function(x,variants)    
{
    x=c(x,variants)
    res<-table(x)-1
    return(res)    
}

## Utility function for re-ordering frequency matrices ##
reOrder<-function(data)
{
#Remove Empty Columns
    if (any(apply(data,2,sum)==0)){data=data[,-which(apply(data,2,sum)==0)]}
    data = data[,order(data[1,],decreasing=TRUE)]	

    tmp=which(data[1,]==0)
    soFar=min(tmp)-1

    for (x in 2:c(nrow(data)))
        {
            if (any(data[x,tmp]>0))
                {
                    count=length(which(data[x,tmp]>0))
                    data=data[,c(1:soFar,order(data[x,tmp],decreasing=TRUE)+soFar)]
                    soFar= soFar+ count	
                    tmp=which(apply(data[1:x,],2,sum)==0)
                }
        }
    return(data)	
}

## compute simpson's diversity index ##

simpson<-function (x, i = 1, base = exp(1)) 
{
    x <- drop(as.matrix(x))
    if (length(dim(x)) > 1) {
        total <- apply(x, i, sum)
        x <- sweep(x, i, total, "/")
    }
    else {
        x <- x/sum(x)
    }

    x <- x * x
    if (length(dim(x)) > 1) 
        H <- apply(x, 1, sum, na.rm = TRUE)
    else {H <- sum(x, na.rm = TRUE)}

    return(1 - H)
}

