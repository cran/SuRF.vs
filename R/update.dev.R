#' update.dev
#'
#' This function is to derive the deviance distribution based on the permutation method
#' This function is not to be used independently but will be called by the function selectnew()
#' @param data: the variable 'data' within seqcutoff()
#' @param vslist: a vector of selected variables
#' @param C: the number of permutation times
#' @param family: family=stats::gaussian(link="identity"));family=stats::binomial(link="logit");family=list(family="cox");etc.
#' @param weights: In a binomial model, weights: =TRUE: if weighted version is desired; =FALSE, otherwise ; In other models,weights: =vector of weights of the same size as the sample size N: if weighted version is desired;=FALSE, otherwise (other generalized model)
#' @return dev: a vector of deviance after C permutations (length OF this vector is C)
#' @export

update_dev<-function(data,vslist, C,weights,ncores,family){
###calculating weights
    len<-dim(data)[2]
    N=dim(data)[1]

    #resample the data by weights

    if(family$family=="binomial"){
      if(is.null(weights)==TRUE|weights==FALSE){
        #unweighted version
        w=rep(1,N) }else{
          N1=length(which(data[,len]==1))
          N0=length(which(data[,len]==0))

          if(N1>=N0){
            w=ifelse(data[,len]==1,1,(N1/N0))
          }else{
            w=ifelse(data[,len]==1,(N0/N1),1)
          }
        }
    }else{

      if(!is.null(weights)){
        if(length(weights)==1&&weights==FALSE){
          w=rep(1,N)
        }else{
          if(is.numeric(weights)&&length(weights)==N)w=weights else{stop("Weights are not numeric or not in the same length as sample size");}
        }

      }else w=rep(1,N)
    }


    sel=1:(len-1)

    X=data[,sel]
    y=data[,len]
    N=length(y)
    p=dim(X)[2]
    nonzeroind=which(apply(X,2,sum)!=0) # select predictors not of all zero's
    dev=NULL

    xmat2=data.frame(X[,vslist],y=y)
    mod2=stats::glm(y~.,data=xmat2,family=family,weights=w,control=list(maxit=500))

  ##

  #setup parallel backend to use many processors
    if(ncores>1){
        cl <- parallel::makeCluster(ncores) #not to overload your computer
        doParallel::registerDoParallel(cl)
        dev<- foreach::`%dopar%`(foreach::foreach(j=1:C, .combine=cbind),  {

            ind=sample(1:N) #permute index
            XX=X
            yy=y

            if(length(vslist)!=0) XX[,vslist]=X[ind,vslist]

            yy=y[ind]
            ww=w[ind]
            devseq=apply(XX[,setdiff(nonzeroind,vslist)],2,function(x){xmat1=data.frame(XX[,vslist],x,yy=yy);mod1=stats::glm(yy~.,data=xmat1,family=family,weights=ww,control=list(maxit=500));return(summary(mod2)[4]$deviance-summary(mod1)[4]$deviance)})


            dmax=max(devseq)
            dmax


        })

        dev=as.vector(dev)
        parallel::stopCluster(cl)
    }else{
        dev<-rep(0,C)

        for(j in 1:C){
            ind=sample(1:N) #permute index
            XX=X
            yy=y

            if(length(vslist)!=0) XX[,vslist]=X[ind,vslist]

            yy=y[ind]
            ww=w[ind]
            devseq=apply(XX[,setdiff(nonzeroind,vslist)],2,function(x){xmat1=data.frame(XX[,vslist],x,yy=yy);mod1=stats::glm(yy~.,data=xmat1,family=family,weights=ww,control=list(maxit=500));return(summary(mod2)[4]$deviance-summary(mod1)[4]$deviance)})


            dev[j]=max(devseq)
        }
    }
    return(dev)
}




