#' dataclean
#'
#' This function is to 1)Scale the count data (count data only) to proportion 2)create a data frame consisting of proportion data, and 3) Keep an variable name list (original variable names and names in terms of X's, e.g.X1,X2,..,etc. )
#' #environmental data (host genome and other information about observations)
#' @param #X.c: data frame that has count data from all levels (only count data will be row scaled)
#' @param X.o: data frame that has other environmental variables (no scaling will be done, those variables will scaled together with proportion data in LASSO step)
#' @param y: a vector representing the outcome (0 or 1 for binomial model)
#' @return data.Xy: a dataframe containing all variables named as X1,X2,...,Xp and the binary outcome (called status)in the last column; this data frame will be used in the other functions for data analysis

#' @export


dataclean<-function(X.c,X.o,y){

  if(length(X.c)!=0)X.c=X.c/apply(X.c,1,sum) # rowscale the count data
  if(length(X.o)==0&length(X.c)!=0)X=data.frame(X.c)
  if(length(X.o)==0&length(X.c)==0)stop("No predictor variables given")
  if(length(X.o)!=0&length(X.c)==0)X=data.frame(X.o)
  if(length(X.o)!=0&length(X.c)!=0)X=data.frame(X.c,X.o)




  Data.Xy=data.frame(X,y) # data frame that has all X's and y (status)

  return(Data.Xy)
}





