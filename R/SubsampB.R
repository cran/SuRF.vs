#' Subsample_B
#'
#' This function is to run sub-sampling B times
#' @param B: the number of sub-samplings to run (e.g., B=1000)
#' @param data: the dataframe should be arranged in the way such that columns are X1,X2,X3....,Xp, status.  Where Xi's are variables and status is the outcome(for the logistic regression, the outcome is in terms of 0/1)
#' @param fold: fold used in lasso cross validation to select the tuning parameter
#' @param Type: should use 'class' for classification always
#' @param Alpha: 1 for Lasso,0 for ridgeression
#' @param prop: percentage of samples left out for each subsamping
#' @param weights: In a binomial model, weights: =TRUE: if weighted version is desired; =FALSE, otherwise ; In other models,weights: =vector of weights of the same size as the sample size N: if weighted version is desired;=FALSE, otherwise (other generalized model)
#' @return Class.Err: mis-classification error on the left out ones over B runs. A vector of length B.
#' @return Lambda: tuning parameters selected from B runs. It is a vector of length B
#' @return BETA: It is a matrix used to save the beta coefficients from all B runs
#' #' @export

 Subsample_B<-function(B,data,fold,Alpha,prop,weights,ncores,family){
    if(family$family=="binomial"){
        ###Should add option to use deviance for binomial distributions
        Type="class"
    }else{
        Type="deviance"
    }
    if(ncores>1){
        cl <- parallel::makeCluster(ncores) #not to overload your computer
        doParallel::registerDoParallel(cl)
    }

    ### User should register cluster prior to calling function

    var.list=NULL




    #length(which(substring(colnames(data),1,1)=="X"&(!is.na(as.numeric(substring(colnames(data),2))))))
    #This is the number of predictors


    if(ncores>1){
        RESULTS=foreach::`%dopar%`(foreach::foreach(i = 1:B, .combine = rbind), {
            if(family$family=="cox"){model=Subsample.w_cox(data=data,fold=fold,Alpha=Alpha,prop=prop,weights=weights)}else{
            model=Subsample.w(data=data,fold=fold,Alpha=Alpha,prop=prop,weights=weights,family,Type)}

            beta=model$Beta
            lamb=model$lambda
            #var.list=c(var.list,model$coef[which(abs(as.numeric(model$coef[,2]))>eps),1][-1])
            error=model$Error
            out=c(error,lamb,beta)
            return(out)

        })
    }else{

        if(family$family=="cox"){
            len=dim(data)[2]-2
            RESULTS=matrix(0,B,len+2)
            for(i in 1:B){
                model=Subsample.w_cox(data=data,fold=fold,Alpha=Alpha,prop=prop,weights=weights)
                beta=model$Beta
                lamb=model$lambda

                error=model$Error
                RESULTS[i,]=c(error,lamb,beta)

            }


        }else{
            len=dim(data)[2]-1
            RESULTS=matrix(0,B,len+3)
            for(i in 1:B){
                model=Subsample.w(data=data,fold=fold,Alpha=Alpha,prop=prop,weights=weights,family,Type)

                beta=model$Beta
                lamb=model$lambda
                #var.list=c(var.list,model$coef[which(abs(as.numeric(model$coef[,2]))>eps),1][-1])
                error=model$Error
                RESULTS[i,]=c(error,lamb,beta)

            }

        }
    }
    #  if(in.parallel){
    #      on.exit(parallel::stopCluster(cl))
    #  }
    return(list(Class.Err=sum(RESULTS[,1])/B,Lambda=RESULTS[,2],BETA=RESULTS[,-c(1:2)]))
}
