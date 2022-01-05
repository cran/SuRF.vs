#' SURF
#'
#' This main function is to apply subsampling, ranking forward selection method(SuRF)
#' @param X - count data, need to be converted to proportion
#' @param Xo - other type of predictor variables
#' @param y - response variable, a vecotr for most families. For family="cox", y will should be a matrix of the response variable in column1 and censoring status in column 2.
#' @param fold - number of folds for cross-validation in Lasso
#' @param Alpha - Alpha parameter for elastic net
#' @param prop - proportion of observations left out in subsampling
#' @param weights - use weighted regression: for unbalanced class sizes (bimomial family only) or weighted sample for other families;In a binomial model, weights: =TRUE: if weighted version is desired; =FALSE, otherwise ; In other models,weights: =vector of weights of the same size as the sample size N: if weighted version is desired;=FALSE, otherwise (other generalized model)
#' @param B - number of subsamples to take
#' @param C - number of permutations used to estimate null distribution
#' @param alpha_u - the upper bound of significance level for the permutation test: alpha_u has to be in the range of (0,1). The large of this value, the longer the program will run;
#' @param alpha - the alpha value of interest (alpha >0 and must be <=alpha_u). It can be a single value or a vector.If missing, by default it is 0.05.
#' @param ncores:  whether SuRF should compute in parallel: 1 indicate NOT; anything greater will compute in parallel
#' @param display.progress - whether SuRF should print a message on completion of each
#' @return Bmod: sub-sampling results
#' @return trdata: data frame including both X and y
#' @return ranklist: ranking table
#' @return modpath: variable selection path (along the alpha range)
#' @return selmod: model results at the selected alpha(s)
#' @return family: model family used
#' @keywords variable selection, LASSO, forward selection
#' @import glmnet
#' @import survival
#' @import dplyr
#' @import doParallel
#' @export
#' @examples
#' library(survival)
#' library(glmnet)
#' library(SuRF)
#' N=100;p=200
#' nzc=p/3
#' X=matrix(rnorm(N*p),N,p)
#' beta=rnorm(nzc)
#' fx=x[,seq(nzc)]%*%beta/3
#' hx=exp(fx)
#' ty=rexp(N,hx)
#' tcens=rbinom(n=N,prob=.3,size=1)# censoring indicator (1 or 0)

#' Xo=NULL
#' B=20
#' Alpha=1
#' fold=5
#' ncores=1
#' prop=0.1
#' C=3
#' alpha_u=0.2
#' alpha=seq(0.01,0.1,len=20)
#'
#' #binomial model
#' XX=X[,1:2]
#' f=1+XX%*%c(2,1.5)
#' p=exp(f)/(1+exp(f))
#' y=rbinom(100,1,p)
#' weights=FALSE
#' family=stats::binomial(link="logit")
#' surf_binary=SURF(Xo=X,y=y,X=NULL,fold=5,Alpha=1,prop=0.1,weights=weights,B=10,C=5,ncores=1,display.progress=TRUE,family=family,alpha_u=0.1,alpha=alpha)
#'
#' #linear regression
#' y=1+XX%*%c(0.1,0.2)
#' family=stats::gaussian(link="identity")
#' surf_lm=SURF(Xo=X,y=y,X=NULL,fold=5,Alpha=1,prop=0.1,weights=weights,B=10,C=5,ncores=1,display.progress=TRUE,family=family,alpha_u=0.1,alpha=alpha)
#'
#' #cox proportional model
#' y=cbind(time=ty,status=1-tcens)
#' weights=rep(1,100)
#' rseed=floor(runif(20,1,100))
#' weights[rseed]=2
#' family=list(family="cox")
#' surf_cox=SURF(Xo=X,y=y,X=NULL,fold=5,Alpha=1,prop=0.1,weights=weights,B=10,C=5,ncores=1,display.progress=TRUE,family=family,alpha_u=alpha_u,alpha=alpha)



SURF<-function(Xo,y,X=NULL,fold=10,Alpha=1,prop=0.1,weights=FALSE,B=1000,C=200,ncores=1,display.progress=TRUE,family=stats::binomial(link="logit"),alpha_u=0.1,alpha=0.05){


    if(family$family!="cox"){

        if(family$link!=get(family$family)()$link){
            warning(paste("glmnet uses default link function ",get(family$family)()$link," whereas you have specified ",family$link,". Models used for ranking and variable selection won't match.",sep="\""))
        }
    }
    #pcutoff=1-pval

    if(!is.element("doParallel", utils::installed.packages()[,1])){
        ncores=1
    }



    traindata=dataclean(X.c=X,X.o=Xo,y=y)
    if(display.progress)  print("clean completed")



    tempmod=Subsample_B(B=B,data=traindata,fold=fold,Alpha=Alpha,prop=prop,weights=weights,ncores,family)
    if(display.progress)  print("subsample  completed")

    if(family$family=="cox")ranking=Ranking_cox(data=traindata,model=tempmod) else ranking=Ranking(data=traindata,model=tempmod)

    if(display.progress)    print("ranking completed")
    if(C==0){
        ans<-list("Bmod"=tempmod,"trdata"=traindata,"ranklist"=ranking,"family"=family)
        class(ans)<-"SURF_fit"
        stop("SURF only performs the ranking without the permutation and variable selection step")
        return(ans)
    }else{

    mod=selpath(data=traindata,weights=weights,ranktable=ranking$table,ncores,family=family,C=C,alpha_u=alpha_u)
    if(display.progress)    print("selection path completed")

    ind=length(which(alpha>alpha_u|alpha<=0)) #find invalid alpha values
    if(ind>0)print("warining: alpha value is<=0 or some values are greater than 'alpha_u'; results will be output for valid alpha values only")

    alpha=alpha[which(alpha<=alpha_u&alpha>0)]
    if(length(alpha)==0){print("warining: no supplied alpha value is valid; results will be given at alpha_u");alpha=alpha_u}

    mod_alpha=apply(as.matrix(alpha),1,function(x)selvar_alpha(res=mod,alpha=x))
    for(i in 1:length(alpha)) mod_alpha[[i]]$alpha=alpha[i]
    if(display.progress)    print("select variable at the 'alpha' significance level completed")

    ans<-list("Bmod"=tempmod,"trdata"=traindata,"ranklist"=ranking,"modpath"=mod,"selmod"=mod_alpha,"family"=family)
    class(ans)<-"SURF_fit"
    return(ans)
    }

}
