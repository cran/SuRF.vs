#' Subsample.w_cox
#'
#' This function is to subsample the data and perform LASSO (single time) on the selected samples for cox proportional model
#' @param B: the number of sub-samplings to run (e.g., B=1000)
#' @param data: the dataframe should be arranged in the way such that columns are X1,X2,X3....,Xp, status.  Where Xi's are variables and status is the outcome(for the logistic regression, the outcome is in terms of 0/1)
#' @param fold: fold used in lasso cross validation to select the tuning parameter
#' @param Type: should use 'class' for classification always
#' @param Alpha: 1 for Lasso,0 for ridgeression
#' @param prop: percentage of samples left out for each sub-sampling
#' @param weights: = a vector of weights: if weighted version is desired,   =FALSE, otherwise
#' @return #lambda: the tuning parameter that within 1 sd of the tuning parameter gives the lowest CV error
#' @return coef: a table shows the name of the selected variables by LASSO and its coefficients
#' @return table: there are a equal proportion of samples from each status left out  and we use the model built on the selected subsamples to predict those left out ones. Table contains two columns: column1 is the predicted value and column2 isthe true value of the outcome
#' @return error: misclassification error based on the above table
#' @return Beta: should be a vector of length p+1 and this is the beta coefficients from the LASSO model.

#' @export




Subsample.w_cox<-function(data,fold,Alpha,prop,weights){
    N=dim(data)[1]
    ID=1:N


    len<-dim(data)[2]


    test.ID=sample(ID,N*prop)

    #APPLY WEIGHTS and resample the data by weights
    if(!is.null(weights)){
        if(length(weights)==1&&weights==FALSE){
            w=rep(1,N)
        }else{
        if(is.numeric(weights)&&length(weights)==N)w=weights else{stop("Weights are not numeric or not in the same length as sample size")}
        }

    }else w=rep(1,N)


    train.ID=ID[-test.ID]

    sel=1:(len-2)


    test.X=data[test.ID,sel]
    train.X=data[-test.ID,sel]

    test.Y=data[test.ID,(len-1)]
    train.Y= data[-test.ID,(len-1):len]
    train.Y2=Surv(train.Y[,1], train.Y[,2])
    train.w=w[-test.ID]
    N=length(test.ID)



    P=dim(train.X)[2]

    Beta=numeric(P)
    NonZero_ind=which(apply(train.X,2,sum)!=0)
    Beta[-NonZero_ind]=0
    ### remove predictors which are all zero from data





    model.cv=glmnet::cv.glmnet(as.matrix(train.X[,NonZero_ind]),train.Y2,family="cox",weights =train.w)


    min.s=model.cv$lambda.min
    model.pr=stats::predict(model.cv,as.matrix(test.X[,NonZero_ind]),s="lambda.min",type="link")


    #  Beta[NonZero_ind]=as.vector(stats::coef(model.cv,s=min.s)[-1])
    Beta[NonZero_ind]=as.vector(stats::coef(model.cv,s=min.s))

    #  Beta=c(Beta,stats::coef(model.cv,s=min.s)[1])

    ###Numeric predictor, use MSE
    table=cbind( as.numeric(noquote(model.pr)),test.Y)
    error.w=mean((table[,1]-table[,2])^2)


    COEF=cbind(rownames(stats::coef(model.cv))[-which(as.logical(stats::coef(model.cv)==0))],stats::coef(model.cv)[-which(as.logical(stats::coef(model.cv)==0))])

    return(list(lambda=min.s,coef=COEF,table=table,Error=error.w,Beta=Beta))

}
