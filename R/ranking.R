#' Ranking
#'
#' This function is to rank the variables after B subsampings;It also removes the highly correlated variables from lower level

#' @param data: the data object return from the dataclean function (the last column is the outcome)
#' @param model: model object from sub-sampling procedure
#' @return table: a table shows the ranked variable list  with its frequency (descending order)
#' @return Beta: coefficients flag (1 or 0) indicating if the variable is selected; intercept is not included;
#' @export


Ranking<-function(data,model){


  data.X=data[,-dim(data)[2]]

    binary.beta=(model$BETA!=0)
        #matrix(ifelse(c(model$BETA)==0,0,1),nrow=dim(model$BETA)[1]) #COEFFICIENT MATRIX OF 0 OR 1 FORMAT
  temp=binary.beta[,-dim(binary.beta)[2]] #REMOVE INTERCEPT COLUMN

    freqs<-colSums(temp) #CALCULATE FREQUENCY
    ranking<-order(freqs,decreasing=TRUE)
    table=data.frame("frequency"=freqs[ranking], "var"=as.character(names(data)[ranking]),stringsAsFactors=FALSE)
    if(sum(table$frequency)!=0){

###the following code removes highly correlated variables (redundant) that at lower level
      mat=(stats::cor(as.matrix(data.X[,table[table[,1]>0,2]])))
        remove.seq=NULL
        for(i in seq_len(dim(mat)[1])){
            remove.seq=c(remove.seq,which(mat[i,-c(1:i)]>0.9999999999)+i)
        }
        remove.seq=unique(remove.seq)
        keep.seq=setdiff(seq_len(dim(mat)[1]),remove.seq)
        table=table[table[,1]>0,][keep.seq,]
    }
  #table only show variables that at least appear once over B runs after removing redundant variables
#    table$var_original=names(data)[table$var]
    return(list(table=table,Beta=temp))
}

