#' selpath
#'
#' #This function is to trace the selection path
#' @param data: the dataframe should be arranged in the way such that columns are X1,X2,X3....,Xp, status.  Where Xi's are variables and status is the outcome(for the logistic regression, the outcome is 0/1)
#' @param ranktable :ranking table from ranking step
#' @param ncores: no of parallel computing cores
#' @param family :generalized model families
#' @param C: the number of permutation times
#' @param alpha_u: the upper significance level
#' @param In a binomial model, weights: =TRUE: if weighted version is desired; =FALSE, otherwise ; In other models,weights: =vector of weights of the same size as the sample size N: if weighted version is desired;=FALSE, otherwise (other generalized model)
#' @return selpoint: a list. it contains each selected variable point,information includes 1)vslist: the variable sect before selecting this variable listed in 'selvar' 2)alpha.range: the variable will be selected within this alpha range 3)pval: pvalue of the variable 4)selvar:selected variable 5)vslist:variable sect after selecting the variable listed in 'selvar'
#' @return sel.nodes: a list. deviance distributions used for selecting the new variable; it includes 1)vslist: the variable sect before the new selection 2)dev.dist: the permutation for selecting the new variable 3)vtlist that has  i)pval: pvalue of the proposed variable ii)selvar: selected variable (the proposed variable is NULL if not selected) iii)dev:deviance contributed by the proposed variable
#' @export


selpath<-function(data,weights,ranktable,ncores,family,C,alpha_u){

alpha_l=0
#alpha_u=0.2
l1=l2=1
#select variables
m=1

sel.nodes=list()
selpoint=list()

#the first selection point starts from a empty vslist and the initial range of alpha (e.g, (0,0.2))
selpoint[[1]]=c(selpoint,list(alpha.range=c(0,alpha_u),selvar=NULL,vslist=NULL)) #individual nodes
vslist=NULL

#repeat the following process until no additional points/leaves generated (all branches are closed)

repeat{
	l.all=0  #the inital sum of leaves starts 0 for every new level(additional variable selection)



 for(i in l1:l2){


		              #only excute it for non-inital point
                                                if(m!=1){
			vslist=selpoint[[i]]$vslist;
			alpha_l=selpoint[[i]]$alpha.range[1];
			alpha_u=selpoint[[i]]$alpha.range[2];
			}

  #check if the ranktable has been the end? no more variables to select from.
  testfinish=ranktable[match(setdiff(ranktable$var,vslist),ranktable$var),]
  if(nrow(testfinish)==0){
                          out=list(vslist=vslist,dev.dist=NULL,vtlist=NULL);
                          sel.nodes[[m]]=out;
                          new=lapply(out$vtlist,function(x){a=unlist(x$selvar);b=unlist(out$vslist);c=c(b,a);return(list(alpha.range=x$alpha.range,pval=x$pval,selvar=a,vslist=c,dev=x$dev))})
                          #selpoint records the info of all leaves
                          selpoint=c(selpoint,new)

                          #l.new mesuare how many leaves generated from a previous point(from a single previous point from upper level)
                          l.new=length(new)

                          #l.all measures new points from all previous points
                          l.all=l.all+l.new

                          #always accumulate the records from the deviance distribution
                          m=m+1
                          l.all=0
                          } else{
                  #out contains the output of the potential new nodes (results from the deviance distribution)
                  out=selectnew(vslist=vslist,ranktable=ranktable,data=data,weights=weights,ncores=ncores,family=family,alpha_l=alpha_l,alpha_u=alpha_u,C=C)

                 #sel.nodes record the deviance distribution results
                 sel.nodes[[m]]=out

                 #new records the new points(new selected leaf)
                 new=lapply(out$vtlist,function(x){a=unlist(x$selvar);b=unlist(out$vslist);c=c(b,a);return(list(alpha.range=x$alpha.range,pval=x$pval,selvar=a,vslist=c,dev=x$dev))})

                 #selpoint records the info of all leaves
                  selpoint=c(selpoint,new)

                 #l.new mesuare how many leaves generated from a previous point(from a single previous point from upper level)
                 l.new=length(new)

                  #l.all measures new points from all previous points
                  l.all=l.all+l.new

                  #always accumulate the records from the deviance distribution
                   m=m+1
                  }
}

if(l.all==0)break  #l.all indicates the total number of points at the current level have no future points, break if 0(0 indicates no more variables)
#loop the points from the additional points
l1=l2+1
l2=l.all+l2
}

return(list(selpoint=selpoint,sel.nodes=sel.nodes))
}



