#' Subsample.w
#'
#' This function is to subsample the data and perform LASSO (single time) on the selected samples
#' @param B: the number of sub-samplings to run (e.g., B=1000)
#' @param data: the dataframe should be arranged in the way such that columns are X1,X2,X3....,Xp, status.  Where Xi's are variables and status is the outcome(for the logistic regression, the outcome is in terms of 0/1)
#' @param fold: fold used in lasso cross validation to select the tuning parameter
#' @param Type: should use 'class' for classification always
#' @param Alpha: 1 for Lasso,0 for ridgeression
#' @param prop: percentage of samples left out for each subsamping
#' @param weights: =TRUE: if weighted version is desired; =FALSE, otherwise (binomial model);weights: =vector of weights of the same size as the sample size N: if weighted version is desired;=FALSE, otherwise (other generalized model)
#' @return lambda: the tuning parameter that within 1 sd of the tuning parameter gives the lowest CV error
#' @return coef: a table shows the name of the selected variables by LASSO and its coefficients
#' @return table: there are a equal proportion of samples from each status left out  and we use the model built on the selected
#' @return subsamples to predict those left out ones. Table contains two columns: column1 is the predicted value and column2 is the true class
#' @return error: misclassification error based on the above table
#' @return Beta: should be a vector of length p+1 and this is the beta coefficients from the LASSO model; Be aware of that the intercept is placed at the end of this vector

#' @export


Subsample.w<-function(data,fold,Alpha,prop,weights,family,Type){
  N=dim(data)[1]
  ID=1:N
  len<-dim(data)[2]

  #DETERMINE THE WEIGHTS
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
        if(is.numeric(weights)&&length(weights)==N)w=weights else{stop("Weights are not numeric or not in the same length as sample size")}
      }

    }else w=rep(1,N)
  }

  #SELECT THE TRAINING AND TEST SAMPLE
  if(family$family=="binomial"){
###      stratified sampling by classes
      test.ID=NULL
      n.samp1=ceiling(length(ID[data[,len]==1])*prop)
      n.samp0=ceiling(length(ID[data[,len]==0])*prop)
      test.ID=sample(ID[data[,len]==1],n.samp1)

      test.ID=c(test.ID,sample(ID[data[,len]==0],n.samp0))
  }else{
      test.ID=sample(ID,N*prop)

  }

  train.ID=ID[-test.ID]

  sel=1:(len-1)

  test.X=data[test.ID,sel]
  train.X=data[-test.ID,sel]
  test.Y=data[test.ID,len]
  train.Y=data[-test.ID,len]
  test.w=w[test.ID]
  train.w=w[-test.ID]

  N=length(test.ID)
  N1=length(which(test.Y==1))
  N0=length(which(test.Y==0))


  P=dim(train.X)[2]

  Beta=numeric(P)
  NonZero_ind=which(apply(train.X,2,sum)!=0)
  Beta[-NonZero_ind]=0
### remove predictors which are all zero from data

  if(family$family=="binomial"){
      train.Y<-as.factor(train.Y)
      test.Y<-as.factor(test.Y)
  }


  if(is.null(weights)|weights==FALSE){
      model.cv=glmnet::cv.glmnet(as.matrix(train.X[,NonZero_ind]),train.Y,family=family$family,type.measure=Type,nfolds=fold,alpha=Alpha,standardize=T ,parallel=FALSE)
  } else {

    model.cv=glmnet::cv.glmnet(as.matrix(train.X[,NonZero_ind]),train.Y,family=family$family,type.measure=Type,nfolds=fold,alpha=Alpha,standardize=T,weights=train.w,parallel=FALSE)
  }

  min.s=model.cv$lambda.min
  model.pr=stats::predict(model.cv,as.matrix(test.X[,NonZero_ind]),s="lambda.min")


  Beta[NonZero_ind]=as.vector(stats::coef(model.cv,s=min.s)[-1])

  Beta=c(Beta,stats::coef(model.cv,s=min.s)[1])

  if(family$family=="binomial"){
      table=cbind( as.numeric(noquote(model.pr)),c(rep(1,N1),rep(0,N0)))
      colnames(table)=c("Predicted Class","True Class")
      if(N1<=N0){
          error.w=((N0/N1)*table[1,2]+table[2,1])/((N0/N1)*(table[1,2]+table[2,2])+(table[1,1]+table[2,1]))
      }else{
          error.w=((N1/N0)*table[2,1]+table[1,2])/((N1/N0)*(table[1,1]+table[2,1])+(table[1,2]+table[2,2]))
      }
  }else{
###Numeric predictor, use MSE
      table=cbind( as.numeric(noquote(model.pr)),test.Y)
      error.w=mean((table[,1]-table[,2])^2)
  }

  COEF=cbind(rownames(stats::coef(model.cv))[-which(as.logical(stats::coef(model.cv)==0))],stats::coef(model.cv)[-which(as.logical(stats::coef(model.cv)==0))])


  return(list(lambda=min.s,coef=COEF,table=table,Error=error.w,Beta=Beta))

}
