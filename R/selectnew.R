#' selectnew
#'
#' This function is to add new node (new deviance distribution for adding the new variable)
#'
#' @param data: the dataframe should be arranged in the way such that columns are X1,X2,X3....,Xp, status.  Where Xi's are variables and status is the outcome(for the logistic regression, the outcome is in terms of 0/1)
#' @param ranktable :ranking table from ranking step
#' @param ncores: no of parallel computing cores
#' @param family :generalized model families
#' @param C: the number of permutation times
#' @param In binomial model, weights: =TRUE: if weighted version is desired; =FALSE, otherwise ; In other models,weights: =vector of weights of the same size as the sample size N: if weighted version is desired;=FALSE, otherwise (other generalized model)
#' @param alpha_l is the minimum significance level(>=0)
#' @param alpha_u: the upper significance level
#' @param In a binomial model, weights: =TRUE: if weighted version is desired; =FALSE, otherwise ; In other models,weights: =vector of weights of the same size as the sample size N: if weighted version is desired;=FALSE, otherwise (other generalized model)
#' @return vslist: the updated list of selected variables
#' @return dev.dist: deviance distributions used for selecting the new variable
#' @return vtlist has 1)alpha.range for the newly selected variable,2)selvar: the newly selected variable name,3)pval: pvalue for the newly selected varaible,and 4)dev: the deviance value contributed by the newly selected variable


#' @export


selectnew<-function(vslist,ranktable,data,weights,ncores=1,family,alpha_l,alpha_u,C){
  ecdf_fun <- function(x,perc) stats::ecdf(x)(perc)
len=dim(data)[2]
N=dim(data)[1]



#APPLY WEIGHTS and resample the data by weights

      if(family$family=="binomial"){
        if(is.null(weights)==TRUE|weights==FALSE){
          #unweighted version
          w=rep(1,length(data[,len])) }else{
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
    if(is.numeric(weights)&&length(weights)==N)w=weights else{stop("Weights are not numeric or not in the same length as sample size");}
  }

}else w=rep(1,N)
}
#update the current ranktable by removing the selected variables from the ranktable
c.ranktable=ranktable[match(setdiff(ranktable$var,vslist),ranktable$var),]

#vclist is the current candidate variable list
vclist=c.ranktable[,2]
freq=c.ranktable[,1]

#derive a NEW NULL distribution given the current selected variable list: vslist

if(family$family=="cox")newdev.dist=update_dev_cox(data=data,vslist=vslist,C=C,ncores,weights=weights) else newdev.dist=update_dev(data=data,vslist=vslist,C=C,weights=weights,ncores=ncores,family=family)

#vtlist is a list of proposed models


if(is.null(vslist))vtlist=as.matrix(apply(as.matrix(vclist),1,function(x){(c(vslist,x))})) else{
  vtlist=t(apply(as.matrix(vclist),1,function(x){(c(vslist,x))}))}

colnames(vtlist)=paste("V",1:dim(vtlist)[2],sep="")

if(family$family=="cox") pval_table=data.frame(vtlist,freq,t(apply( vtlist,1,function(y){z=y;

 							        yy=data[,(len-1):len] #response and censored info
								y_status=Surv(yy[,1],yy[,2])#response var in the survival model
                                                                dat=data.frame(y_status,data[,c(z)]);

                                                                colnames(dat)[-1]=c(z)
                                                                vnames=colnames(dat)[-1]
                                                                fmla <- stats::as.formula(paste("y_status ~ ", paste(vnames, collapse= "+")))
                                                                mod=coxph(fmla,data=dat,singular.ok = TRUE,weights = w)
                                                                ANOVA=stats::anova(mod,test="Chisq");
                                                                vardev=ANOVA[dim(vtlist)[2]+1,2];
                                                                squantile=1-ecdf_fun(newdev.dist,vardev);

                                                                temp=cbind(squantile,vardev);
                                                                colnames(temp)=c("pval","deviance");
                                                                return(temp);

}
)
)

) else   pval_table=data.frame(vtlist,freq,t(apply( vtlist,1,function(y){z=y;
                                                                dat=data.frame(data[,c(z)],data[,len]);
                                                                colnames(dat)=c(z,"status");
                                                                mod=stats::glm(status~.,data=dat,family=family,weights=w,control=list(maxit=500));
                                                                ANOVA=stats::anova(mod,test="Chisq");
                                                                vardev=ANOVA[dim(vtlist)[2]+1,2];
                                                                squantile=1-ecdf_fun(newdev.dist,vardev);

                                                                temp=cbind(squantile,vardev);
                                                                colnames(temp)=c("pval","deviance");
                                                                return(temp);

}
)
)

)

colnames(pval_table)[(ncol(pval_table)-1):ncol(pval_table)]=c("pval","deviance")

#order the data by frequency and deviance value, both in descending order;
#this step is to handle the tie senario;
pval_table=pval_table[with(pval_table, order(-freq, -deviance)), ]

#keep rows that have a p-value less than the current alpha_u
#(may we should do: remove any row that contributes a deviance greater than the critical value according to 1-alpha_u )
cutoff=stats::quantile(newdev.dist,1-alpha_u)
#pval_table=pval_table[pval_table$pval<alpha_u,]
pval_table=pval_table[pval_table$deviance>cutoff,]

l.vc=dim(pval_table)[1]

##################################################################################
#if no variable left in the current table, return empty list for 'VTLIST'
#als return the selecti
##################################################################################
if(l.vc==0){sel=vector("list", length=0);
            return(list(vslist=vslist,dev.dist=newdev.dist,vtlist=sel))
            }else{





##################################################################################
#for non-empty list (at least 1 variable left in pval_table currently )
##################################################################################
cummaxdev=0
ind=numeric(l.vc)
ind[1]=1

#entries will be kept if its deviance value is greater than the maximum cumulative deviance value from variables above it
 if(l.vc>1){ for(i in 2:l.vc){
cummaxdev=max(cummaxdev,pval_table$deviance[i-1] )
  if(pval_table$deviance[i]>cummaxdev) ind[i]=1
  }
pval_table=pval_table[which(ind==1),]
}

#if some variable(s) has(have) an  p-value less than the alpha_l;
#only include the variable up to the first variable that has p-value<alpha_l;
ind_alpha_low=which(pval_table$pval<alpha_l)
len_max=length(pval_table$pval)  #check the total number of records in pval_table
if(length(ind_alpha_low))pval_table=pval_table[1:min(ind_alpha_low),]  #




#formatting the data to give the range of significance level
arange=c(alpha_u,pval_table$pval)
#replace the pvalue by the mimimum alpha value if there is a variable that has pvalue<alpha_l
#there will be up to 1 such variable (has a pvalue<alpha_l)
arange[arange<alpha_l]=alpha_l

#pval is a vector of p values of all variables in pval_table
pval=pval_table$pval

#DEV is a vector of deviance values that are contributed by the corresponding variables
DEV=pval_table$deviance



#pvar is a vector of variable names of all variables in pval_table
#pvar=pval_table$vtlist
#print(pval_table)


pvar=t(rev(as.vector(pval_table[colnames(pval_table)[substr(colnames(pval_table),1,1)=="V"]]))[1])


#tentative number of selected new nodes;
size.range=length(arange)-1
#print(size.range)
#sel is a list of selected nodes given that variables that have been selected
sel=vector("list", length = size.range);

#each selected node will contain info: 1)selvar,variable name 2)pvalue of this variable and 3)alpha range,alpha.range 4)deviance contribution from the variable
  for(i in 1:size.range){
      sel[[i]] =list(alpha.range=arange[(i+1):i],selvar=ifelse(is.na(pvar[i]),NA,pvar[i]),pval=pval[i],dev=DEV[i]);

 }


#delete:sel=sel[unlist(lapply(sel,function(x)length(x)))!=0]

#remove the entry if the selected variable is NA
sel=sel[unlist(lapply(sel,function(x)!is.na(x$selvar)))]

#prepare the output of the function
res=list(vslist=vslist,dev.dist=newdev.dist,vtlist=sel)

return(res)
}

}
