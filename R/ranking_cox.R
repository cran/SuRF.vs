#' Ranking_cox
#'
#' This function is to rank the variables after B subsampings for cox proportional model; It also removes the highly correlated variables from lower level

#' @param data: the data object return from the dataclean function(the last column is the outcome)
#' @param model: cox proportion model object from sub-sampling procedure (B times)
#' @return table: a table shows the ranked variable list  with its frequency (descending order)
#' @return Beta: coefficients flag (1 or 0) indicating if the variable is selected; intercept is not included;
#' @export


Ranking_cox<-function(data,model){

    ncol=dim(data)[2]
    data.X=data[,1:(ncol-2)]

    binary.beta=(model$BETA!=0)

    temp=binary.beta
    freqs<-colSums(   binary.beta) #CALCULATE FREQUENCY
    ranking<-order(freqs,decreasing=TRUE)
    table=data.frame("frequency"=freqs[ranking], "var"=as.character(names(data)[ranking]),stringsAsFactors=FALSE)
    if(sum(table$frequency)!=0){

        ###the following code removes highly correlated variables (redundant) that at lower level
        mat=(stats::cor(as.matrix(data.X[,table[table[,1]>0,2]])))
        remove.seq=NULL
        for(i in 1:dim(mat)[1]){
            remove.seq=c(remove.seq,which(mat[i,-c(1:i)]>0.9999999999)+i)
        }
        remove.seq=unique(remove.seq)
        keep.seq=setdiff(1:dim(mat)[1],remove.seq)
        table=table[table[,1]>0,][keep.seq,]
    }

    return(list(table=table,Beta=temp))
}
