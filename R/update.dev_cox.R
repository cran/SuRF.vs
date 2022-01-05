#' update.dev_cox
#'
#' For COX proportional model ONLY. This function is to derive the deviance distribution based on the permutation method;This function is not to be used independently but will be called by the function selectnew()
#' @param data: the variable 'data' within seqcutoff()
#' @param vslist: a vector of selected variables
#' @param C: the number of permutation times
#' @param family: family=stats::gaussian(link="identity"));family=stats::binomial(link="logit");family=list(family="cox");etc.
#' @param weights: =TRUE: if weighted version is desired,   =FALSE, otherwise (binomial model); weights: =vector of weights of the same size as the sample size N: if weighted version is desired,   =FALSE, otherwise (other generalized model)
#' @return dev: a vector of deviance after C permutations (length OF this vector is C)
#' @export


update_dev_cox<-function(data,vslist, C,ncores,weights){

    len<-dim(data)[2]
    N=dim(data)[1]





    sel=1:(len-2)

    X=data[,sel]
    y=data[,(len-1):len]
    N=dim(y)[1]
    p=dim(X)[2]
    nonzeroind=which(apply(X,2,sum)!=0) # select predictors not of all zero's
    dev=NULL


    #APPLY WEIGHTS and resample the data by weights
    if(!is.null(weights)){
      if(length(weights)==1&&weights==FALSE){
        w=rep(1,N)
      }else{
        if(is.numeric(weights)&&length(weights)==N)w=weights else{stop("Weights are not numeric or not in the same length as sample size");}
      }

    }else w=rep(1,N)

    yy=Surv(y[,1],y[,2])
    mod2=coxph(yy~NULL,weights=w)


    #setup parallel backend to use many processors

    if(ncores>1){
        cl <- parallel::makeCluster(ncores) #not to overload your computer
        doParallel::registerDoParallel(cl)

        dev<- foreach::`%dopar%`(foreach::foreach(j=1:C, .combine=cbind),  {

            ind=sample(1:N) #permute index
            XX=X
            yy=y

            if(length(vslist)!=0) XX[,vslist]=X[ind,vslist]

            yy=yy[ind,]
            yyy=Surv(yy[,1],yy[,2])

            devseq=apply(XX[,setdiff(nonzeroind,vslist)],2,function(x){xmat1=cbind(XX[,vslist],x);colnames(xmat1)=c(vslist,"x");xmat1=as.matrix(xmat1);mod1=coxph(yyy~xmat1,singular.ok=TRUE,weights=w);return(-2*(mod2$`loglik`[1]-mod1$`loglik`[2]))})


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
            yy=y[ind,]

            if(length(vslist)!=0) XX[,vslist]=X[ind,vslist]


            yyy=Surv(yy[,1],yy[,2])

            devseq=apply(XX[,setdiff(nonzeroind,vslist)],2,function(x){xmat1=cbind(XX[,vslist],x);colnames(xmat1)=c(vslist,"x");xmat1=as.matrix(xmat1);mod1=coxph(yyy~xmat1,singular.ok=TRUE,weights=w);return(-2*(mod2$`loglik`[1]-mod1$`loglik`[2]))})

            dev[j]=max(devseq)
        }
    }
    return(dev)
}
