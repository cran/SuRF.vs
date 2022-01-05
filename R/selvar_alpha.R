#' selvar_alpha
#'
#' This function is to extract summarize the results at 'alpha' level from 'mod' object to obtain 1)final selected variables 'selvar' #2)pvalue of each selected variable according to the variables in the 'selvar' 3)deviance contributed by each selected variable (given the previous selected variables 4)deviance permutation distribution 5)cutoff value based on (1-alpha)percentile of the deviance permutation distribution;when no variable is selected, only return the last deviance distribution and the cutoff value;this function can be used separately after running selpath(); the alpha value must be >0 and <= alpha_u parameter from SURF()
#' @param res: 'mod' object returned from 'selpath' function
#' @param alpha: alpha level(default alpha=0.05)(a single value up to the value 'alpha_u' sepecified in selpath() function)
#' @return selvar:final selected variable
#' @return pval:pvalue of each selected variable (present if at least 1 var is selected)
#' @return devlist:deviance contributed by each selected variable (given the previous selected variables;present if at least 1 var is selected )
#' @return dist.mat:a list of deviance permutation distributions (including the distribution from the step from which no more variable is added)

#' @importFrom dplyr %>%
#' @export

if(getRversion() >= "2.15.1")  utils::globalVariables(c("a_lower","a_upper"))
### These are actually column names in a data frame.  However, the
### automated check cannot detect this.  This line is to prevent
### automated code-checking from giving an error.  The function is not
### in any way influenced by the existence or otherwise of these
### variables.


selvar_alpha<-function(res,alpha){

				f=res[1]
				len=length(f[[1]]);
				cate=as.data.frame(matrix(unlist(lapply(f[[1]],function(x)x$alpha.range)),ncol=2,byrow=T))
				colnames(cate)=c("a_lower","a_upper")
				cate2=cate[with(cate, order(a_lower, a_upper)),]
				cate3=cate2 %>% distinct(a_lower,a_upper, .keep_all = TRUE)#remove duplicated interval
				cate4=cate3[with(cate3, order(a_lower,-a_upper)),]

				len=dim(cate4)[1]
				selres=vector("list", length = len)

				for(j in 1:len){
					cat=as.vector(unlist(cate4[j,]))
					ind=unlist(lapply(lapply(f[[1]],function(x){arange=x[[1]];return(arange)}),function(x)identical(unlist(x),cat)))
					selres[[j]]=list(arange=cat, vslist=lapply(f[[1]][ind],function(x)x$vslist))
				}

				selres2=lapply(selres,function(x){len=length(x$vslist);return(list(vslist=x$vslist[[len]],alpha.range=x$arange))})
				selres3=lapply(selres2,function(x)if(x[['alpha.range']][1]<=alpha&&x[['alpha.range']][2]>=alpha)return(x))
				selres4=selres3[!sapply(selres3,is.null)]

				len=length(selres4)
				selres5=vector("list", length = len)

				for(k in 1:len){
					a=sum(unlist(lapply(selres4,function(x){len1=length(selres4[[k]][['vslist']]);len2=length(x[['vslist']]);return(identical(selres4[[k]][['vslist']],x[['vslist']][1:len1]))})))
					if(sum(unlist(lapply(selres4,function(x){len1=length(selres4[[k]][['vslist']]);len2=length(x[['vslist']]);return(identical(selres4[[k]][['vslist']],x[['vslist']][1:len1]))})))>1)selres5[[k]]=NULL else selres5[[k]]=selres4[[k]]
				}

				selres6=selres5[!sapply(selres5,is.null)]
				len=length(selres5)
				if(len>1)selres6=selres5[!sapply(lapply(selres5,function(x)x[['vslist']]),is.null)]
				finalres=selres6

var=finalres[[1]][['vslist']]
varlen=length(var)+1
output=vector("list", length = varlen)
varlist=lapply(res$selpoint,function(x)x$vslist)
varlist2=lapply(res$sel.nodes,function(x)x$vslist)


if(!is.null(var)){
for(i in 1:varlen){

if(i<=varlen-1){
output[[i]]$selvar=var[i]
index1=which(lapply(varlist,function(x)paste(x,collapse=" "))==paste(var[1:i],collapse=" "))
output[[i]]$pval=res$selpoint[[index1]]$pval #pvalue
output[[i]]$dev=res$selpoint[[index1]]$dev #deviance value contributed by the corresponding variable
}
if(i==1)index2=which(lapply(varlist2,function(x)is.null(x))==TRUE) else index2=which(lapply(varlist2,function(x)paste(x,collapse=" "))==paste(var[1:(i-1)],collapse=" "))

output[[i]]$cutoff=stats::quantile(res$sel.nodes[[index2]]$dev.dis,1-alpha) #cutoff value based on (1-alpha)percentile of the dev.dist
output[[i]]$dev.dist=res$sel.nodes[[index2]]$dev.dist

}
summary=list(selvar=unlist(lapply(output,function(x)x[["selvar"]])),pval=unlist(lapply(output,function(x)x[["pval"]])),devlist=unlist(lapply(output,function(x)x[["dev"]])),cutoff=unlist(lapply(output,function(x)x[["cutoff"]])),dist.mat=lapply(output,function(x)x[["dev.dist"]]))
}else summary=list(selvar=NULL,dist.mat=lapply(output,function(x)x$dev.dist),cutoff=unlist(lapply(output,function(x)x$cutoff)))

return(summary)
}
